//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "ADComputeStressBase.h"
#include "RankTwoTensor.h"
#include "PerfGraph.h"
#include <torch/script.h>
#include <torch/torch.h>

/**
 * Elasto-viscoplastic stress update using the Multivariate Newton Return-Mapping
 * (MNR-RR) algorithm from the MoE manuscript (Algorithm 1, Appendix B).
 *
 * Unlike the univariate (UNR-RR) scheme that solves only for the plastic strain
 * increment ∆λ and updates ρ_c / ρ_w explicitly afterwards, MNR-RR solves for
 * ∆q = [∆λ, ∆ρ_c, ∆ρ_w] simultaneously.  During the inner Newton iterations
 * the MoE model is evaluated at the current (updated) dislocation densities,
 * which allows much larger time steps.
 *
 * The 3×3 Newton Jacobian is (Eq. 42 of the manuscript):
 *
 *   J = I − D·∆t,   D = [[-3G ∂ε̇/∂σ,  ∂ε̇/∂ρ_c,  ∂ε̇/∂ρ_w],
 *                         [-3G ∂ρ̇_c/∂σ, ∂ρ̇_c/∂ρ_c, ∂ρ̇_c/∂ρ_w],
 *                         [-3G ∂ρ̇_w/∂σ, ∂ρ̇_w/∂ρ_c, ∂ρ̇_w/∂ρ_w]]
 *
 * Derivatives are obtained from the full 3×6 Jacobian of the MoE model via
 * PyTorch autograd, then scaled by the forward/inverse transform derivatives.
 */
class ADMixtureOfExpertsMNRElastoViscoplasticStress : public ADComputeStressBase
{
public:
  static InputParameters validParams();

  ADMixtureOfExpertsMNRElastoViscoplasticStress(const InputParameters & parameters);

  // ---- Input normalisation (forward transforms) ----
  torch::Tensor minmax_transform(const torch::Tensor & x,
                                 double lb, double ub, double xmin, double xmax);
  torch::Tensor logtransform_transform(const torch::Tensor & x,
                                       double lb, double ub,
                                       double logxmin, double logxmax, double eps);

  // ---- Output de-normalisation (inverse transforms) ----
  torch::Tensor logtransform_inverse_transform(const torch::Tensor & z,
                                               double lb, double ub,
                                               double logxmin, double logxmax, double eps);
  torch::Tensor symlog_inverse_transform(const torch::Tensor & z_scaled,
                                         double zmin, double zmax);

  // ---- Derivatives of transforms (for Jacobian) ----
  torch::Tensor minmax_forward_derivative(const torch::Tensor & x,
                                          double lb, double ub, double xmin, double xmax);
  torch::Tensor logtransform_forward_derivative(const torch::Tensor & x,
                                                double lb, double ub,
                                                double logxmin, double logxmax, double eps);
  torch::Tensor logtransform_inverse_derivative(const torch::Tensor & z,
                                                double lb, double ub,
                                                double logxmin, double logxmax, double eps);
  torch::Tensor symlog_inverse_derivative(const torch::Tensor & z_scaled,
                                          double zmin, double zmax);

  /**
   * Evaluate the MoE model: apply forward transforms, run inference, apply
   * inverse transforms.  Returns [1, 3]: [ep_rate, cell_rate, wall_rate].
   * Input tensor is [1, 6]: [sigma_eff, T, ep, rhoc, rhow, flux].
   */
  torch::Tensor eval(const torch::Tensor & tensor);

  /**
   * Compute the full physical Jacobian dOutput/dInput, shape [3, 6].
   * Row i = derivative of output i w.r.t. all 6 physical inputs.
   * Used by the MNR Newton loop to extract columns [0, 3, 4] (σ, ρ_c, ρ_w).
   */
  torch::Tensor eval_jacobian_full(const torch::Tensor & tensor);

protected:
  torch::jit::script::Module model;

  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress() override;

  const std::string _elasticity_tensor_name;
  const ADMaterialProperty<RankFourTensor> & _elasticity_tensor;
  const MaterialProperty<RankTwoTensor> & _stress_old;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;

  const ADMaterialProperty<Real> & _mu;
  const ADMaterialProperty<Real> & _temp;

  ADMaterialProperty<Real> & _rhoc;
  ADMaterialProperty<Real> & _rhow;
  const MaterialProperty<Real> & _rhoc_old;
  const MaterialProperty<Real> & _rhow_old;

  const ADMaterialProperty<Real> & _flux;

  ADMaterialProperty<Real> & _total_effective_plastic_strain;
  const MaterialProperty<Real> & _total_effective_plastic_strain_old;
  ADMaterialProperty<Real> & _effective_plastic_strain_rate;

  const Real & _initial_evm;
  const Real & _initial_rhoc;
  const Real & _initial_rhow;

  const Real & _tol_lambda;
  const Real & _tol_rho;
  const int & _max_iter;

  ADMaterialProperty<RankTwoTensor> & _incremental_plastic_strain_tensor;
  const MaterialProperty<RankTwoTensor> & _incremental_plastic_strain_tensor_old;

  PerfID _solve_timer;
};
