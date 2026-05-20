//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ArcLengthSolveObject.h"

#include "Executioner.h"
#include "FEProblemBase.h"
#include "NonlinearSystemBase.h"
#include "AuxiliarySystem.h"
#include "MooseVariableScalar.h"

#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/implicit_system.h"

#include <cmath>

InputParameters
ArcLengthSolveObject::validParams()
{
  InputParameters params = FEProblemSolve::validParams();
  params.addClassDescription(
      "Cylindrical (Riks) arc-length method per Box 4.4 of de Souza Neto et al. "
      "Performs two linear solves per Newton iteration with the same K_T (regular "
      "step + tangential), then solves a scalar quadratic for the iterative load "
      "factor δλ. Predictor at k=1 uses ±Δs/√(δūᵀδū) with sign by secant criterion.");
  params.addRequiredParam<Real>("delta_s",
                                "Arc-length step Δs (incremental displacement L2 norm).");
  params.addParam<unsigned int>("al_max_iter", 30, "Maximum AL Newton iterations per step.");
  params.addParam<Real>("al_rel_tol", 1e-6, "Convergence: ||r|| / ||lambda*q^|| < tol.");
  params.addParam<Real>("al_lambda_abs_tol", 0.0,
                        "Absolute convergence tolerance on |delta_lambda|. "
                        "0 = disabled. Useful for NIPG/DG with stiff k_T.");
  MooseEnum solver_kind("lu gmres", "lu");
  params.addParam<MooseEnum>("al_linear_solver", solver_kind,
                              "Linear solver for the two K_T solves per AL Newton iteration. "
                              "'lu' = direct (MUMPS in parallel); 'gmres' = GMRES + LU PC "
                              "(recommended for non-symmetric K_T, e.g. NIPG DG).");
  params.addRequiredParam<VariableName>(
      "scalar_variable",
      "AuxScalarVariable holding the load parameter λ (NOT a NL Variable).");
  params.addParam<TagName>(
      "arc_length_load_tag", "arc_length_load",
      "Residual tag of the AL-controlled load contributions (used to cache q^).");
  return params;
}

ArcLengthSolveObject::ArcLengthSolveObject(Executioner & ex)
  : FEProblemSolve(ex),
    _delta_s(getParam<Real>("delta_s")),
    _max_iter(getParam<unsigned int>("al_max_iter")),
    _rtol(getParam<Real>("al_rel_tol")),
    _lambda_abs_tol(getParam<Real>("al_lambda_abs_tol")),
    _linear_solver(std::string(getParam<MooseEnum>("al_linear_solver"))),
    _scalar_var_name(getParam<VariableName>("scalar_variable")),
    _al_tag_name(getParam<TagName>("arc_length_load_tag")),
    _al_tag_id(Moose::INVALID_TAG_ID),
    _kappa_var(nullptr),
    _aux_sys(nullptr),
    _lambda_old(0.0),
    _lambda_current(0.0),
    _q_ref_cached(false)
{
}

void
ArcLengthSolveObject::initialSetup()
{
  FEProblemSolve::initialSetup();

  // λ lives in the auxiliary system (not the NL K matrix).
  _aux_sys = &_problem.getAuxiliarySystem();
  _kappa_var = &_aux_sys->getScalarVariable(/*tid=*/0, _scalar_var_name);

  if (!_problem.vectorTagExists(_al_tag_name))
    mooseError("ArcLengthSolveObject: AL load tag '",
               _al_tag_name,
               "' not found. Add `extra_tag_vectors = ",
               _al_tag_name,
               "` under [Problem].");
  _al_tag_id = _problem.getVectorTagID(_al_tag_name);

  _console << "ArcLengthSolveObject ready. λ-var (aux) '" << _scalar_var_name
           << "', AL tag '" << _al_tag_name << "' (id=" << _al_tag_id
           << "), Δs=" << _delta_s << std::endl;
}

bool
ArcLengthSolveObject::solve()
{
  return boxAlgorithm();
}

std::unique_ptr<NumericVector<libMesh::Number>>
ArcLengthSolveObject::makeVec()
{
  // Clone the system's GHOSTED residual vector to get a parallel layout
  // that exactly matches MOOSE's assembly expectations (local dofs + ghost
  // slots from the dof_map send list). Hand-rolling with `init(n_dofs,
  // PARALLEL)` or `solution().clone()` gives a non-ghost-padded vector;
  // MOOSE's residual assembly writes into ghosted indices and segfaults.
  auto & nl = _problem.currentNonlinearSystem();
  auto v = nl.residualGhosted().zero_clone();
  return v;
}

void
ArcLengthSolveObject::writeLambda(Real lambda)
{
  // Mirror in this object so readLambda() can return without indexing the
  // parallel aux solution vector (which crashes on non-owner ranks).
  _lambda_current = lambda;
  // Push to aux system so coupled BCs see the new value.
  _kappa_var->setValues(lambda);
  _kappa_var->insert(_aux_sys->solution());
  _aux_sys->solution().close();
  _aux_sys->system().update();
}

Real
ArcLengthSolveObject::readLambda() const
{
  // Return the locally-cached mirror -- consistent across ranks because
  // writeLambda() sets it uniformly. Indexing the aux parallel solution
  // vector for a SCALAR DOF crashes on non-owner ranks (the DOF lies
  // outside every rank's local range).
  return _lambda_current;
}

void
ArcLengthSolveObject::computeResidual(NumericVector<libMesh::Number> & r)
{
  auto & nl = _problem.currentNonlinearSystem();
  // MOOSE residual computation reads variables from the localized
  // (ghosted) solution view, not the bare parallel `nl.solution()`. Passing
  // the wrong one segfaults during assembly in parallel.
  _problem.computeResidual(*nl.system().current_local_solution.get(), r, nl.number());
}

libMesh::SparseMatrix<libMesh::Number> *
ArcLengthSolveObject::systemMatrix()
{
  auto & nl = _problem.currentNonlinearSystem();
  auto * impl = dynamic_cast<libMesh::ImplicitSystem *>(&nl.system());
  if (!impl || !impl->matrix)
    mooseError("ArcLengthSolveObject: nonlinear system has no implicit matrix.");
  return impl->matrix;
}

void
ArcLengthSolveObject::assembleJacobian()
{
  auto & nl = _problem.currentNonlinearSystem();
  _problem.computeJacobianTag(nl.solution(), *systemMatrix(), nl.systemMatrixTag());
}

bool
ArcLengthSolveObject::linearSolve(const NumericVector<libMesh::Number> & rhs,
                                  NumericVector<libMesh::Number> & dx)
{
  auto * petsc_K = dynamic_cast<libMesh::PetscMatrix<libMesh::Number> *>(systemMatrix());
  if (!petsc_K)
    mooseError("ArcLengthSolveObject: matrix not a PetscMatrix.");
  Mat J = petsc_K->mat();

  KSP ksp;
  PetscErrorCode ierr = KSPCreate(_communicator.get(), &ksp);
  if (ierr) return false;
  ierr = KSPSetOperators(ksp, J, J);
  if (ierr) { (void)KSPDestroy(&ksp); return false; }

  // Pick solver type: LU (direct) or GMRES (with LU as preconditioner).
  // For non-symmetric K_T (NIPG/IIPG DG), GMRES converges robustly; pure
  // LU works too but at higher per-iteration memory cost. The LU PC
  // amortizes the same factorization across the two K_T solves.
  if (_linear_solver == "gmres")
  {
    ierr = KSPSetType(ksp, KSPGMRES);
    if (ierr) { (void)KSPDestroy(&ksp); return false; }
    ierr = KSPSetTolerances(ksp, 1e-10, 1e-50, PETSC_DEFAULT, 1000);
    if (ierr) { (void)KSPDestroy(&ksp); return false; }
  }
  else  // "lu"
  {
    ierr = KSPSetType(ksp, KSPPREONLY);
    if (ierr) { (void)KSPDestroy(&ksp); return false; }
  }
  PC pc;
  ierr = KSPGetPC(ksp, &pc);
  if (ierr) { (void)KSPDestroy(&ksp); return false; }
  ierr = PCSetType(pc, PCLU);
  if (ierr) { (void)KSPDestroy(&ksp); return false; }
  // PCLU's default factor solver is sequential. Pick MUMPS when running in
  // parallel; fall back to PETSc's built-in (single-rank) otherwise.
  if (_communicator.size() > 1)
  {
    ierr = PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
    if (ierr)
    {
      (void)KSPDestroy(&ksp);
      mooseError("ArcLengthSolveObject: parallel linear solve requires MUMPS "
                 "(PCFactorSetMatSolverType MATSOLVERMUMPS failed, ierr=",
                 ierr,
                 "). Rebuild PETSc with --download-mumps, or run in serial.");
    }
  }

  auto * petsc_rhs = dynamic_cast<const libMesh::PetscVector<libMesh::Number> *>(&rhs);
  auto * petsc_dx = dynamic_cast<libMesh::PetscVector<libMesh::Number> *>(&dx);
  if (!petsc_rhs || !petsc_dx)
  {
    (void)KSPDestroy(&ksp);
    mooseError("ArcLengthSolveObject: vectors must be PetscVector.");
  }

  ierr = KSPSolve(ksp, petsc_rhs->vec(), petsc_dx->vec());
  (void)KSPDestroy(&ksp);
  if (ierr) {
    _console << "    linearSolve: KSPSolve failed, ierr=" << ierr << std::endl;
    return false;
  }
  return true;
}

void
ArcLengthSolveObject::cacheReferenceLoad()
{
  auto & nl = _problem.currentNonlinearSystem();
  if (!_q_ref)
    _q_ref = makeVec();

  const Real lambda_save = readLambda();

  auto r0 = makeVec();
  writeLambda(0.0);
  _problem.computeResidual(*nl.system().current_local_solution.get(), *r0, nl.number());

  auto r1 = makeVec();
  writeLambda(1.0);
  _problem.computeResidual(*nl.system().current_local_solution.get(), *r1, nl.number());

  *_q_ref = *r0;
  _q_ref->add(-1.0, *r1);
  _q_ref->close();

  writeLambda(lambda_save);
  _q_ref_cached = true;
}

bool
ArcLengthSolveObject::boxAlgorithm()
{
  auto & nl = _problem.currentNonlinearSystem();
  auto & u = nl.solution();

  // Cache q^ once (proportional loading assumption).
  if (!_q_ref_cached)
    cacheReferenceLoad();
  const Real q_norm = _q_ref->l2_norm();
  if (q_norm <= 0.0)
    mooseError("ArcLengthSolveObject: ||q^|| = 0 — no AL-tagged load found. "
               "Check that the loading BC has `extra_vector_tags = ",
               _al_tag_name,
               "`.");

  // Δu^(0) = 0 for this step.
  if (!_du_inc)
    _du_inc = makeVec();
  _du_inc->zero();

  // λ for this step starts at the last converged λ.
  Real lambda = _lambda_old;
  writeLambda(lambda);

  auto r       = makeVec();
  auto du_star = makeVec();
  auto du_bar  = makeVec();
  auto neg_r   = makeVec();

  for (unsigned int k = 1; k <= _max_iter; ++k)
  {
    // Assemble K_T at current u.
    assembleJacobian();

    // Two linear solves with same K_T.
    computeResidual(*r);
    *neg_r = *r;
    neg_r->scale(-1.0);
    if (!linearSolve(*neg_r, *du_star))
      return false;
    if (!linearSolve(*_q_ref, *du_bar))
      return false;

    Real dlambda = 0.0;
    const Real norm_du_bar_sq = du_bar->dot(*du_bar);
    if (!std::isfinite(norm_du_bar_sq) || norm_du_bar_sq <= 0.0)
    {
      _console << "  AL k=" << k << " ||δū||² invalid (" << norm_du_bar_sq
               << ") — aborting step." << std::endl;
      return false;
    }

    if (k == 1)
    {
      // Predictor: sign by previous step's increment direction.
      Real sign = 1.0;
      if (_du_old)
      {
        const Real s_dot = _du_old->dot(*du_bar);
        if (s_dot < 0.0) sign = -1.0;
      }
      dlambda = sign * _delta_s / std::sqrt(norm_du_bar_sq);
    }
    else
    {
      // Corrector: scalar quadratic in δλ.
      auto temp = makeVec();
      *temp = *_du_inc;
      *temp += *du_star;
      const Real a = norm_du_bar_sq;
      const Real b = 2.0 * temp->dot(*du_bar);
      const Real c = temp->dot(*temp) - _delta_s * _delta_s;
      const Real disc = b * b - 4.0 * a * c;
      if (disc < 0.0)
      {
        _console << "  AL k=" << k << " quadratic disc < 0 — aborting step." << std::endl;
        return false;
      }
      const Real sqd = std::sqrt(disc);
      const Real x1 = (-b + sqd) / (2.0 * a);
      const Real x2 = (-b - sqd) / (2.0 * a);
      // Pick root that minimises ||r|| at the candidate new state. Robust at
      // sharp turning points where Crisfield's max-dot-product criterion
      // would lock onto the wrong branch and oscillate.
      auto try_root = [&](Real xi) -> Real {
        auto du_try = makeVec();
        *du_try = *du_star;
        du_try->add(xi, *du_bar);
        auto u_try = makeVec();
        *u_try = u;
        *u_try += *du_try;
        u_try->close();
        // Temporarily apply (u_try, λ + xi) and measure residual.
        const Real lam_save = readLambda();
        writeLambda(lam_save + xi);
        // Swap solution: copy u_try into nl.solution(), evaluate, restore.
        auto u_save = makeVec();
        *u_save = u;
        u_save->close();
        u = *u_try;
        u.close();
        nl.system().update();
        auto r_try = makeVec();
        computeResidual(*r_try);
        const Real n = r_try->l2_norm();
        u = *u_save;
        u.close();
        nl.system().update();
        writeLambda(lam_save);
        return n;
      };
      const Real n1 = try_root(x1);
      const Real n2 = try_root(x2);
      dlambda = (n1 < n2) ? x1 : x2;
    }

    // Update λ and u.
    lambda += dlambda;
    writeLambda(lambda);

    auto du_k = makeVec();
    *du_k = *du_star;
    du_k->add(dlambda, *du_bar);

    u += *du_k;
    u.close();
    nl.system().update();
    *_du_inc += *du_k;
    _du_inc->close();

    // Convergence check: residual ratio AND/OR |δλ| absolute tolerance.
    computeResidual(*r);
    const Real r_norm = r->l2_norm();
    const Real ratio = r_norm / (q_norm * std::max(std::abs(lambda), 1.0));
    const bool conv_r = (ratio < _rtol);
    const bool conv_l = (_lambda_abs_tol > 0.0 &&
                          std::abs(dlambda) < _lambda_abs_tol);
    _console << "  AL k=" << k << "  λ=" << lambda
             << "  ||r||/||λq^||=" << ratio
             << "  δλ=" << dlambda
             << (conv_r ? "  [r-conv]" : "")
             << (conv_l ? "  [λ-conv]" : "")
             << std::endl;
    if (conv_r || conv_l)
    {
      _lambda_old = lambda;
      if (!_du_old)
        _du_old = makeVec();
      *_du_old = *_du_inc;
      _du_old->close();
      return true;
    }
  }
  _console << "ArcLengthSolveObject: max iterations reached without convergence." << std::endl;
  return false;
}
