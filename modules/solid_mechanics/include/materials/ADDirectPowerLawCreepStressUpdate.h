//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "StressUpdateBase.h"

/**
 * Integrates isotropic power-law creep by solving directly for the six stress components.
 */
class ADDirectPowerLawCreepStressUpdate : public ADStressUpdateBase
{
public:
  static InputParameters validParams();

  ADDirectPowerLawCreepStressUpdate(const InputParameters & parameters);

  void initQpStatefulProperties() override;
  void propagateQpStatefulProperties() override;

  void updateState(ADRankTwoTensor & strain_increment,
                   ADRankTwoTensor & inelastic_strain_increment,
                   const ADRankTwoTensor & rotation_increment,
                   ADRankTwoTensor & stress_new,
                   const RankTwoTensor & stress_old,
                   const ADRankFourTensor & elasticity_tensor,
                   const RankTwoTensor & elastic_strain_old,
                   bool compute_full_tangent_operator = false,
                   RankFourTensor & tangent_operator = _identityTensor) override;

  bool requiresIsotropicTensor() override { return true; }
  bool isIsotropic() override { return true; }

protected:
  /// Compute the viscoplastic strain rate at a supplied stress.
  ADRankTwoTensor computeFlowRate(const ADRankTwoTensor & stress) const;

  /// Evaluate the backward-Euler stress residual and its derivative.
  void computeResidualAndJacobian(const ADRankTwoTensor & stress,
                                  const ADRankTwoTensor & trial_stress,
                                  const ADRankFourTensor & elasticity_tensor,
                                  ADRankTwoTensor & residual,
                                  ADRankFourTensor & jacobian) const;

  /// Convert the fourth-order stress derivative to the MOOSE six-component ordering.
  void computeVoigtJacobian(const ADRankFourTensor & jacobian,
                            DenseMatrix<ADReal> & voigt_jacobian) const;

  const ADVariableValue * const _temperature;
  const Real _coefficient;
  const Real _n_exponent;
  const Real _m_exponent;
  const Real _activation_energy;
  const Real _gas_constant;
  const Real _start_time;
  const Real _absolute_tolerance;
  const Real _relative_tolerance;
  const unsigned int _max_iterations;
  const bool _verbose;

  ADMaterialProperty<RankTwoTensor> & _creep_strain;
  const MaterialProperty<RankTwoTensor> & _creep_strain_old;
  MaterialProperty<Real> & _local_iterations;
};
