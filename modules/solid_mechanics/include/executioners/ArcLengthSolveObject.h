//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FEProblemSolve.h"
#include "MooseTypes.h"

#include <petscsnes.h>

class MooseVariableScalar;
class AuxiliarySystem;

/**
 * Cylindrical (Riks/Crisfield) arc-length method per Box 4.4 of
 * de Souza Neto, Perić, Owen, *Computational Methods for Plasticity*
 * (Wiley 2008), §4.4.2 — non-consistent two-solve scheme.
 *
 *   Per Newton iteration k (within one outer arc-length step Δs):
 *     (a) Assemble K_T at current u^(k-1).
 *     (b) Solve  K_T · δu*  = -r^(k-1)         ← regular Newton step
 *         Solve  K_T · δū   = +q̂                ← tangential to load
 *     (c) Compute δλ^(k):
 *           - if k==1 (predictor): δλ^(1) = ±Δs / √(δūᵀ δū)
 *             sign by sign(Δu_oldᵀ δū), default +1.
 *           - else: scalar quadratic a δλ² + b δλ + c = 0 with
 *               a = δūᵀ δū
 *               b = 2 (Δu^(k-1) + δu*)ᵀ δū
 *               c = (Δu^(k-1) + δu*)ᵀ (Δu^(k-1) + δu*) − Δs²
 *             pick root that maximises (Δu^(k))ᵀ Δu^(k-1).
 *     (d) λ^(k)  = λ^(k-1) + δλ^(k)
 *     (e) δu^(k) = δu* + δλ^(k) · δū
 *     (f) Update u, Δu
 *     (g) Recompute r, check convergence ||r||/||q̂|| < tol.
 *
 * λ is an AuxScalarVariable controlled by this object; it is *not*
 * a DOF of the nonlinear system. The reference load q̂ is cached once
 * on the first solve from the AL-tagged residual contribution at λ=1
 * (proportional loading assumption).
 */
class ArcLengthSolveObject : public FEProblemSolve
{
public:
  static InputParameters validParams();
  ArcLengthSolveObject(Executioner & ex);

  virtual void initialSetup() override;
  virtual bool solve() override;

protected:
  /// Set the aux scalar λ DOF and propagate so BCs read the new value.
  void writeLambda(Real lambda);

  /// Read the current λ from the aux scalar DOF.
  Real readLambda() const;

  /// Compute residual r = f_int(u) − λ·q̂ at the current (u, λ).
  void computeResidual(NumericVector<libMesh::Number> & r);

  /// Assemble the system Jacobian K_T at the current u.
  void assembleJacobian();

  /// Get the underlying libMesh implicit-system matrix.
  libMesh::SparseMatrix<libMesh::Number> * systemMatrix();

  /// Solve K_T · dx = rhs via a fresh KSP (LU). Returns true on success.
  bool linearSolve(const NumericVector<libMesh::Number> & rhs,
                   NumericVector<libMesh::Number> & dx);

  /// Cache the AL reference load q̂ once. Proportional loading: q̂ does not
  /// depend on u or λ, so this is computed at the first solve and reused.
  void cacheReferenceLoad();

  /// One outer arc-length step: predictor + corrector to convergence.
  bool boxAlgorithm();

  /// Allocate a working vector matching the NL system size.
  std::unique_ptr<NumericVector<libMesh::Number>> makeVec();

  /// Configurable parameters
  const Real _delta_s;
  const unsigned int _max_iter;
  const Real _rtol;
  /// Absolute convergence tolerance on |δλ| (per-iteration load-factor
  /// increment). 0 = disabled. Useful for NIPG/DG where the residual ratio
  /// alone can be misleading: a tiny δλ on a stiff problem means we have
  /// effectively converged even if ||r||/||λ q^|| > rtol.
  const Real _lambda_abs_tol;
  /// Linear-solver type for the two K_T·δ = rhs solves inside each AL
  /// Newton iteration. "lu" -> direct LU (MUMPS in parallel); "gmres" ->
  /// GMRES + LU preconditioner (recommended for non-symmetric K_T as in
  /// NIPG/IIPG DG, where the consistency + adjoint terms break symmetry).
  const std::string _linear_solver;
  const VariableName _scalar_var_name;
  const TagName _al_tag_name;
  TagID _al_tag_id;

  /// Aux scalar variable carrying λ (NOT in the NL system).
  MooseVariableScalar * _kappa_var;
  AuxiliarySystem * _aux_sys;

  /// Cumulative load factor from previous converged step.
  Real _lambda_old;

  /// Local mirror of the current λ value (consistent across all MPI ranks).
  /// We push this into the aux system solution so coupled BCs see it; we
  /// also keep it here so readLambda() never has to index the parallel aux
  /// solution vector (which segfaults on non-owner ranks for SCALAR DOFs).
  Real _lambda_current;

  /// q̂ cached after first solve (proportional loading).
  bool _q_ref_cached;

  /// Persistent state across timesteps.
  std::unique_ptr<NumericVector<libMesh::Number>> _q_ref;     // q̂
  std::unique_ptr<NumericVector<libMesh::Number>> _du_old;    // Δu_n (previous converged step)
  std::unique_ptr<NumericVector<libMesh::Number>> _du_inc;    // Δu^(k) within current step
};
