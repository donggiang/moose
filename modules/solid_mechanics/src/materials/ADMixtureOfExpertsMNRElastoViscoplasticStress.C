//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADMixtureOfExpertsMNRElastoViscoplasticStress.h"
#include "MooseException.h"
#include <torch/script.h>
#include <torch/torch.h>
#include <iostream>
#include <memory>
#include <cmath>
#include <algorithm>

registerMooseObject("SolidMechanicsApp", ADMixtureOfExpertsMNRElastoViscoplasticStress);

// ---------------------------------------------------------------------------
// Helper: solve 3×3 linear system A·x = b using Gaussian elimination with
// partial pivoting.  A and b are overwritten.  x is written into x[].
// ---------------------------------------------------------------------------
static void
solve3x3(double A[3][3], double b[3], double x[3])
{
  for (int col = 0; col < 3; ++col)
  {
    // Partial pivot: find row with largest |A[row][col]|
    int pivot = col;
    for (int row = col + 1; row < 3; ++row)
      if (std::abs(A[row][col]) > std::abs(A[pivot][col]))
        pivot = row;

    // Swap rows in A and b
    std::swap(b[col], b[pivot]);
    for (int k = 0; k < 3; ++k)
      std::swap(A[col][k], A[pivot][k]);

    // Eliminate below
    for (int row = col + 1; row < 3; ++row)
    {
      if (std::abs(A[col][col]) < 1e-300)
        continue; // singular or near-singular
      double factor = A[row][col] / A[col][col];
      b[row] -= factor * b[col];
      for (int k = col; k < 3; ++k)
        A[row][k] -= factor * A[col][k];
    }
  }

  // Back substitution
  for (int row = 2; row >= 0; --row)
  {
    x[row] = b[row];
    for (int col = row + 1; col < 3; ++col)
      x[row] -= A[row][col] * x[col];
    x[row] /= (std::abs(A[row][row]) > 1e-300 ? A[row][row] : 1.0);
  }
}

// ===========================================================================
// InputParameters
// ===========================================================================

InputParameters
ADMixtureOfExpertsMNRElastoViscoplasticStress::validParams()
{
  InputParameters params = ADComputeStressBase::validParams();
  params.addClassDescription(
      "Computes stress using the Multivariate Newton Return-Mapping (MNR-RR) "
      "scheme with a Mixture-of-Experts ML constitutive model.  Solves for "
      "[∆λ, ∆ρ_c, ∆ρ_w] simultaneously, giving larger stable time steps than "
      "the univariate (UNR-RR) scheme.");
  params.addRequiredParam<MaterialPropertyName>("mu", "Shear modulus");
  params.addRequiredParam<MaterialPropertyName>("temp", "Temperature");
  params.addRequiredParam<Real>("initial_evm", "Initial equivalent plastic strain");
  params.addRequiredParam<Real>("initial_rhoc", "Initial cell dislocation density [m^-2]");
  params.addRequiredParam<Real>("initial_rhow", "Initial wall dislocation density [m^-2]");
  params.addRequiredParam<MaterialPropertyName>("flux", "Normalised neutron flux [0, 1]");
  params.addParam<Real>("newton_tol_lambda", 5e-7, "Absolute convergence tolerance for the plastic strain increment residual R0");
  params.addParam<Real>("newton_tol_rho", 1e-4, "Relative convergence tolerance for dislocation density residuals R1, R2");
  params.addParam<int>("max_newton_iter", 100, "Maximum number of Newton iterations in the inner loop");
  return params;
}

// ===========================================================================
// Constructor
// ===========================================================================

ADMixtureOfExpertsMNRElastoViscoplasticStress::ADMixtureOfExpertsMNRElastoViscoplasticStress(
    const InputParameters & parameters)
  : ADComputeStressBase(parameters),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getADMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "stress")),
    _mechanical_strain_old(
        getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _mu(getADMaterialProperty<Real>("mu")),
    _temp(getADMaterialProperty<Real>("temp")),
    _rhoc(declareADProperty<Real>("rhoc")),
    _rhow(declareADProperty<Real>("rhow")),
    _rhoc_old(getMaterialPropertyOld<Real>("rhoc")),
    _rhow_old(getMaterialPropertyOld<Real>("rhow")),
    _flux(getADMaterialProperty<Real>("flux")),
    _total_effective_plastic_strain(declareADProperty<Real>("effective_plastic_strain")),
    _total_effective_plastic_strain_old(getMaterialPropertyOld<Real>("effective_plastic_strain")),
    _effective_plastic_strain_rate(declareADProperty<Real>("effective_plastic_strain_rate")),
    _initial_evm(getParam<Real>("initial_evm")),
    _initial_rhoc(getParam<Real>("initial_rhoc")),
    _initial_rhow(getParam<Real>("initial_rhow")),
    _tol_lambda(getParam<Real>("newton_tol_lambda")),
    _tol_rho(getParam<Real>("newton_tol_rho")),
    _max_iter(getParam<int>("max_newton_iter")),
    _incremental_plastic_strain_tensor(
        declareADProperty<RankTwoTensor>("incremental_plastic_strain_tensor")),
    _incremental_plastic_strain_tensor_old(
        getMaterialPropertyOld<RankTwoTensor>("incremental_plastic_strain_tensor")),
    _solve_timer(moose::internal::getPerfGraphRegistry().registerSection(
        "ADMixtureOfExpertsMNRElastoViscoplasticStress::MoE AD solve", 3))
{
  std::string model_path =
      "/Users/huyngd/projects/moose_neml2/modules/solid_mechanics/test/tests/"
      "MoE/model/model2experts.pt";

  if (!std::filesystem::exists(model_path))
  {
    std::cerr << "Model file does not exist: " << model_path << std::endl;
    throw std::runtime_error("MNR: model file not found");
  }

  try
  {
    model = torch::jit::load(model_path);
    std::cout << "MNR: model loaded successfully." << std::endl;
  }
  catch (const c10::Error & e)
  {
    std::cerr << "MNR: error loading model: " << e.what() << std::endl;
    throw;
  }
}

// ===========================================================================
// initQpStatefulProperties
// ===========================================================================

void
ADMixtureOfExpertsMNRElastoViscoplasticStress::initQpStatefulProperties()
{
  _total_effective_plastic_strain[_qp] = _initial_evm;
  _rhoc[_qp] = _initial_rhoc;
  _rhow[_qp] = _initial_rhow;
}

// ===========================================================================
// computeQpStress  —  MNR-RR algorithm (Algorithm 1, Appendix B)
// ===========================================================================

void
ADMixtureOfExpertsMNRElastoViscoplasticStress::computeQpStress()
{
  PerfGuard time_guard(_app.perfGraph(), _solve_timer);

  // ---- 1. Retrieve and bound old dislocation densities ----
  double rhoc_n = static_cast<double>(_rhoc_old[_qp]);
  double rhow_n = static_cast<double>(_rhow_old[_qp]);

  if (rhoc_n <= 0.0)
    rhoc_n = _initial_rhoc;
  if (rhow_n <= 0.0)
    rhow_n = _initial_rhow;
  if (rhoc_n > 8461801410123.313)
    rhoc_n = 8461801410123.313;
  if (rhow_n > 11999567054170.322)
    rhow_n = 11999567054170.322;

  // ---- 2. Trial stress (purely elastic) ----
  ADRankTwoTensor elastic_strain_increment =
      _mechanical_strain[_qp] - _mechanical_strain_old[_qp];
  RankTwoTensor trial_stress = MetaPhysicL::raw_value(
      _stress_old[_qp] + _elasticity_tensor[_qp] * elastic_strain_increment);

  ADRankTwoTensor deviatoric_trial_stress = trial_stress.deviatoric();
  auto norm_dev_sq = deviatoric_trial_stress.doubleContraction(deviatoric_trial_stress);
  auto effective_trial_stress = std::sqrt(1.5 * MetaPhysicL::raw_value(norm_dev_sq));
  const Real tolerance = 1e-14;

  if (effective_trial_stress > (1.0 + tolerance))
  {
    // ---- 3. Constant scalars for this QP / time step ----
    const double G = MetaPhysicL::raw_value(_mu[_qp]);
    const double T = MetaPhysicL::raw_value(_temp[_qp]);
    const double flux = MetaPhysicL::raw_value(_flux[_qp]);
    const double eps_p_n = MetaPhysicL::raw_value(_total_effective_plastic_strain_old[_qp]);
    const double sigma_trial = effective_trial_stress;
    const double dt = (_dt > 0.0) ? _dt : 1.0;

    // ---- 4. MNR-RR Newton loop for ∆q = [∆λ, ∆ρ_c, ∆ρ_w] ----
    double delta_lambda = 0.0;
    double delta_rhoc = 0.0;
    double delta_rhow = 0.0;

    // Normalization denominators set from the first model evaluation.
    // Tolerances are relative for ρ (their scale is ~1e12-1e13) and
    // absolute for λ (scale ~1e-8 to 1e-6).
    double norm_rhoc = 1.0;
    double norm_rhow = 1.0;
    const double tol_lambda = _tol_lambda;
    const double tol_rho    = _tol_rho;
    const int max_iter = _max_iter;

    double ep_rate = 0.0, cell_rate = 0.0, wall_rate = 0.0;

    for (int k = 0; k < max_iter; ++k)
    {
      // Current effective stress and dislocation densities
      double sigma_eff = sigma_trial - 3.0 * G * delta_lambda;
      double rhoc_curr = rhoc_n + delta_rhoc;
      double rhow_curr = rhow_n + delta_rhow;

      // Clamp to model training bounds
      sigma_eff = std::max(sigma_eff, 0.11428862639259);
      rhoc_curr = std::max(std::min(rhoc_curr, 8461801410123.313), 1e8);
      rhow_curr = std::max(std::min(rhow_curr, 11999567054170.322), 1e8);

      // Evaluate MoE model (returns [1,3]: [ep_rate, cell_rate, wall_rate])
      at::Tensor inputs =
          torch::tensor({{sigma_eff, T, eps_p_n, rhoc_curr, rhow_curr, flux}}, torch::kFloat);
      at::Tensor rates = eval(inputs);

      ep_rate = rates.index({0, 0}).item<double>();
      cell_rate = rates.index({0, 1}).item<double>();
      wall_rate = rates.index({0, 2}).item<double>();

      // Set normalization from first evaluation
      if (k == 0)
      {
        norm_rhoc = std::max(std::abs(cell_rate * dt), 1.0);
        norm_rhow = std::max(std::abs(wall_rate * dt), 1.0);
      }

      // Residuals: R_i = ∆q_i − rate_i · ∆t
      double R0 = delta_lambda - ep_rate * dt;
      double R1 = delta_rhoc - cell_rate * dt;
      double R2 = delta_rhow - wall_rate * dt;

      // Convergence check (skip on k=0 to always take at least one step)
      if (k > 0 && std::abs(R0) < tol_lambda &&
          std::abs(R1) / norm_rhoc < tol_rho &&
          std::abs(R2) / norm_rhow < tol_rho)
        break;

      if (k == max_iter - 1)
        throw MooseException("MNR Newton loop failed to converge after " +
                             std::to_string(max_iter) +
                             " iterations at QP " + std::to_string(_qp) +
                             " (dt=" + std::to_string(_dt) + ").");

      // ---- Jacobian:  J = I − D·∆t ----
      //
      // D = [[-3G·∂ε̇/∂σ,   ∂ε̇/∂ρ_c,    ∂ε̇/∂ρ_w  ],
      //      [-3G·∂ρ̇c/∂σ,  ∂ρ̇c/∂ρ_c,   ∂ρ̇c/∂ρ_w ],
      //      [-3G·∂ρ̇w/∂σ,  ∂ρ̇w/∂ρ_c,   ∂ρ̇w/∂ρ_w ]]
      //
      // Physical Jacobian of outputs w.r.t. inputs is [3,6]; we need
      // columns 0 (σ), 3 (ρ_c), 4 (ρ_w).

      at::Tensor J_full = eval_jacobian_full(inputs); // [3, 6]

      double dEP_dSig  = J_full.index({0, 0}).item<double>();
      double dEP_dRhoc = J_full.index({0, 3}).item<double>();
      double dEP_dRhow = J_full.index({0, 4}).item<double>();
      double dRC_dSig  = J_full.index({1, 0}).item<double>();
      double dRC_dRhoc = J_full.index({1, 3}).item<double>();
      double dRC_dRhow = J_full.index({1, 4}).item<double>();
      double dRW_dSig  = J_full.index({2, 0}).item<double>();
      double dRW_dRhoc = J_full.index({2, 3}).item<double>();
      double dRW_dRhow = J_full.index({2, 4}).item<double>();

      // Full consistent Jacobian: ε^vm_{n+1} = ε^vm_n + ∆λ → ∂ε^vm/∂∆λ = 1
      // The -dX_dEps*dt terms in column 0 account for this dependence.
      double J[3][3] = {
          {1.0 + 3.0 * G * dEP_dSig * dt, -dEP_dRhoc * dt, -dEP_dRhow * dt},
          {3.0 * G * dRC_dSig * dt, 1.0 - dRC_dRhoc * dt, -dRC_dRhow * dt},
          {3.0 * G * dRW_dSig * dt, -dRW_dRhoc * dt, 1.0 - dRW_dRhow * dt}};


      // RHS = -R
      double b[3] = {-R0, -R1, -R2};
      double d[3] = {0.0, 0.0, 0.0};

      solve3x3(J, b, d);

      delta_lambda += d[0];
      delta_rhoc += d[1];
      delta_rhow += d[2];
    }

    // ---- 5. Update state after Newton convergence ----
    _total_effective_plastic_strain[_qp] =
        _total_effective_plastic_strain_old[_qp] + delta_lambda;
    _rhoc[_qp] = rhoc_n + delta_rhoc;
    _rhow[_qp] = rhow_n + delta_rhow;
    _effective_plastic_strain_rate[_qp] = delta_lambda / dt;

    // Plastic strain tensor (radial direction in deviatoric space)
    _incremental_plastic_strain_tensor[_qp] =
        deviatoric_trial_stress *
        (1.5 * delta_lambda / effective_trial_stress);

    // Correct elastic strain increment
    elastic_strain_increment -= _incremental_plastic_strain_tensor[_qp];
  }
  else
  {
    // Purely elastic step — carry state forward
    _total_effective_plastic_strain[_qp] = _total_effective_plastic_strain_old[_qp];
    _rhoc[_qp] = rhoc_n;
    _rhow[_qp] = rhow_n;
    _effective_plastic_strain_rate[_qp] = 0.0;
    _incremental_plastic_strain_tensor[_qp] = ADRankTwoTensor();
  }

  _stress[_qp] = _stress_old[_qp] + _elasticity_tensor[_qp] * elastic_strain_increment;
}

// ===========================================================================
// eval — apply transforms, run inference, apply inverse transforms
// Returns [1, 3]: [ep_rate, cell_rate, wall_rate]
// ===========================================================================

torch::Tensor
ADMixtureOfExpertsMNRElastoViscoplasticStress::eval(const torch::Tensor & tensor)
{
  auto varx1 = minmax_transform(tensor.index({0, 0}), -1.0, 1.0,
                                0.11428862639259419753390289997697,
                                299.94931150862146296276478096842766);
  auto varx2 = minmax_transform(tensor.index({0, 1}), -1.0, 1.0,
                                600.04429302445396388066001236438751,
                                1099.98408980407225499220658093690872);
  auto varx3 = logtransform_transform(tensor.index({0, 2}), -1.0, 1.0,
                                      -39.15241605429763183110480895265937,
                                      -3.91206356603903770974284270778298, 1e-20);
  auto varx4 = logtransform_transform(tensor.index({0, 3}), -1.0, 1.0,
                                      12.56521587815730534032354626106098,
                                      29.64417981061014728538793860934675, 0.0);
  auto varx5 = logtransform_transform(tensor.index({0, 4}), -1.0, 1.0,
                                      29.12196678204971789227784029208124,
                                      30.11574121515777946456182689871639, 0.0);
  auto varx6 = minmax_transform(tensor.index({0, 5}), -1.0, 1.0, 0.0, 1.0);

  at::Tensor transformed =
      torch::tensor({{varx1.item<double>(), varx2.item<double>(), varx3.item<double>(),
                      varx4.item<double>(), varx5.item<double>(), varx6.item<double>()}},
                    torch::kFloat);

  at::Tensor output = model.forward({transformed}).toTensor();

  auto out1 = logtransform_inverse_transform(output.index({0, 0}), 0.0, 1.0,
                                             -30.54946341410551013950680498965085,
                                             7.72739242062932873977842973545194, 0.0);
  auto out2 = symlog_inverse_transform(output.index({0, 1}),
                                       -36.71520922998040958873389172367752,
                                       38.45450753607088500984900747425854);
  auto out3 = symlog_inverse_transform(output.index({0, 2}),
                                       -40.36444506932617315442257677204907,
                                       35.04333232507948281408971524797380);

  return torch::tensor({{out1.item<double>(), out2.item<double>(), out3.item<double>()}},
                       torch::kFloat);
}

// ===========================================================================
// eval_jacobian_full — full physical Jacobian [3, 6] via PyTorch autograd
// J[i,j] = d(physical_output_i) / d(physical_input_j)
// ===========================================================================

torch::Tensor
ADMixtureOfExpertsMNRElastoViscoplasticStress::eval_jacobian_full(const torch::Tensor & tensor)
{
  // ---- Forward transform derivatives (diagonal of d(norm_in)/d(phys_in)) ----
  auto jax1 = minmax_forward_derivative(tensor.index({0, 0}), -1.0, 1.0,
                                        0.11428862639259419753390289997697,
                                        299.94931150862146296276478096842766);
  auto jax2 = minmax_forward_derivative(tensor.index({0, 1}), -1.0, 1.0,
                                        600.04429302445396388066001236438751,
                                        1099.98408980407225499220658093690872);
  auto jax3 = logtransform_forward_derivative(tensor.index({0, 2}), -1.0, 1.0,
                                              -39.15241605429763183110480895265937,
                                              -3.91206356603903770974284270778298, 1e-20);
  auto jax4 = logtransform_forward_derivative(tensor.index({0, 3}), -1.0, 1.0,
                                              12.56521587815730534032354626106098,
                                              29.64417981061014728538793860934675, 0.0);
  auto jax5 = logtransform_forward_derivative(tensor.index({0, 4}), -1.0, 1.0,
                                              29.12196678204971789227784029208124,
                                              30.11574121515777946456182689871639, 0.0);
  auto jax6 = minmax_forward_derivative(tensor.index({0, 5}), -1.0, 1.0, 0.0, 1.0);

  // Jin[j] = d(norm_in_j) / d(phys_in_j)  — shape [6]
  at::Tensor Jin =
      torch::tensor({jax1.item<double>(), jax2.item<double>(), jax3.item<double>(),
                     jax4.item<double>(), jax5.item<double>(), jax6.item<double>()},
                    torch::kFloat);

  // ---- Normalised inputs ----
  auto varx1 = minmax_transform(tensor.index({0, 0}), -1.0, 1.0,
                                0.11428862639259419753390289997697,
                                299.94931150862146296276478096842766);
  auto varx2 = minmax_transform(tensor.index({0, 1}), -1.0, 1.0,
                                600.04429302445396388066001236438751,
                                1099.98408980407225499220658093690872);
  auto varx3 = logtransform_transform(tensor.index({0, 2}), -1.0, 1.0,
                                      -39.15241605429763183110480895265937,
                                      -3.91206356603903770974284270778298, 1e-20);
  auto varx4 = logtransform_transform(tensor.index({0, 3}), -1.0, 1.0,
                                      12.56521587815730534032354626106098,
                                      29.64417981061014728538793860934675, 0.0);
  auto varx5 = logtransform_transform(tensor.index({0, 4}), -1.0, 1.0,
                                      29.12196678204971789227784029208124,
                                      30.11574121515777946456182689871639, 0.0);
  auto varx6 = minmax_transform(tensor.index({0, 5}), -1.0, 1.0, 0.0, 1.0);

  // Repeat 3 times so we can compute the full Jacobian in one backward pass
  at::Tensor p =
      torch::tensor({{varx1.item<double>(), varx2.item<double>(), varx3.item<double>(),
                      varx4.item<double>(), varx5.item<double>(), varx6.item<double>()}},
                    torch::kFloat)
          .repeat({3, 1}); // [3, 6]
  p.requires_grad_(true);

  // ---- Model Jacobian in normalised space: Jmod[i,j] = d(norm_out_i)/d(norm_in_j) ----
  torch::Tensor y = model.forward({p}).toTensor(); // [3, 3]
  torch::Tensor eye_matrix = torch::eye(3, torch::kFloat);
  y.backward(eye_matrix);
  torch::Tensor Jmod = p.grad(); // [3, 6]

  // ---- Inverse transform derivatives (diagonal of d(phys_out)/d(norm_out)) ----
  // Evaluate at the *first* row's normalised output (all 3 rows are identical input)
  auto phi = y[0].unsqueeze(0); // [1, 3]

  auto out1_deriv = logtransform_inverse_derivative(phi.index({0, 0}), 0.0, 1.0,
                                                    -30.54946341410551013950680498965085,
                                                    7.72739242062932873977842973545194, 0.0);
  auto out2_deriv = symlog_inverse_derivative(phi.index({0, 1}),
                                             -36.71520922998040958873389172367752,
                                             38.45450753607088500984900747425854);
  auto out3_deriv = symlog_inverse_derivative(phi.index({0, 2}),
                                             -40.36444506932617315442257677204907,
                                             35.04333232507948281408971524797380);

  // Jout[i] = d(phys_out_i) / d(norm_out_i)  — shape [3]
  at::Tensor Jout =
      torch::tensor(
          {out1_deriv.item<double>(), out2_deriv.item<double>(), out3_deriv.item<double>()},
          torch::kFloat);

  // ---- Chain rule: J[i,j] = Jout[i] * Jmod[i,j] * Jin[j] ----
  // Jout.unsqueeze(1) → [3,1]; Jmod*Jin → [3,6] (broadcast); result → [3,6]
  torch::Tensor J = Jout.unsqueeze(1) * (Jmod * Jin);

  return J; // [3, 6]
}

// ===========================================================================
// Transform functions (identical to UNR class)
// ===========================================================================

torch::Tensor
ADMixtureOfExpertsMNRElastoViscoplasticStress::minmax_transform(
    const torch::Tensor & x, double lb, double ub, double xmin, double xmax)
{
  return (ub - lb) * (x - xmin) / (xmax - xmin) + lb;
}

torch::Tensor
ADMixtureOfExpertsMNRElastoViscoplasticStress::logtransform_transform(
    const torch::Tensor & x, double lb, double ub, double logxmin, double logxmax, double eps)
{
  return (ub - lb) * (torch::log(x + eps) - logxmin) / (logxmax - logxmin) + lb;
}

torch::Tensor
ADMixtureOfExpertsMNRElastoViscoplasticStress::logtransform_inverse_transform(
    const torch::Tensor & z, double lb, double ub, double logxmin, double logxmax, double eps)
{
  return torch::exp((logxmax - logxmin) * (z - lb) / (ub - lb) + logxmin) - eps;
}

torch::Tensor
ADMixtureOfExpertsMNRElastoViscoplasticStress::symlog_inverse_transform(
    const torch::Tensor & z_scaled, double zmin, double zmax)
{
  auto z = (z_scaled + 1.0) * (zmax - zmin) / 2.0 + zmin;
  return torch::sign(z) * (torch::exp(torch::abs(z)) - 1.0);
}

torch::Tensor
ADMixtureOfExpertsMNRElastoViscoplasticStress::minmax_forward_derivative(
    const torch::Tensor & /*x*/, double lb, double ub, double xmin, double xmax)
{
  return torch::tensor((ub - lb) / (xmax - xmin));
}

torch::Tensor
ADMixtureOfExpertsMNRElastoViscoplasticStress::logtransform_forward_derivative(
    const torch::Tensor & x, double lb, double ub, double logxmin, double logxmax, double eps)
{
  return (ub - lb) / (logxmax - logxmin) / (x + eps);
}

torch::Tensor
ADMixtureOfExpertsMNRElastoViscoplasticStress::logtransform_inverse_derivative(
    const torch::Tensor & z, double lb, double ub, double logxmin, double logxmax, double eps)
{
  auto u = torch::exp((logxmax - logxmin) * (z - lb) / (ub - lb) + logxmin);
  return u * (logxmax - logxmin) / (ub - lb);
}

torch::Tensor
ADMixtureOfExpertsMNRElastoViscoplasticStress::symlog_inverse_derivative(
    const torch::Tensor & z_scaled, double zmin, double zmax)
{
  auto z = (z_scaled + 1.0) * (zmax - zmin) / 2.0 + zmin;
  return torch::exp(torch::abs(z)) * (zmax - zmin) / 2.0;
}
