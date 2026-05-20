//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADDGCohesiveZoneOutput.h"

#include "metaphysicl/raw_type.h"

#include <cmath>
using std::sqrt;
using std::pow;

registerMooseObject("SolidMechanicsApp", ADDGCohesiveZoneOutput);

InputParameters
ADDGCohesiveZoneOutput::validParams()
{
  InputParameters params = InterfaceMaterial::validParams();
  params.addClassDescription(
      "Diagnostic output for ADDGCohesiveZone: declares normal/shear jump, "
      "damage, normal/shear traction as boundary MaterialProperty<Real> for "
      "visualization via MaterialRealAux. Mirrors the kernel's TSL math.");
  params.addRequiredCoupledVar("displacements", "Vector of displacement variables.");
  params.addRequiredParam<Real>("penalty_stiffness", "Cohesive penalty stiffness K.");
  params.addRequiredParam<Real>("normal_strength", "Peak normal traction T.");
  params.addRequiredParam<Real>("shear_strength", "Peak shear traction S.");
  params.addRequiredParam<Real>("GI_c", "Mode-I fracture energy.");
  params.addRequiredParam<Real>("GII_c", "Mode-II fracture energy.");
  params.addRequiredParam<Real>("eta_BK", "Benzeggagh-Kenane mixed-mode exponent.");
  params.addParam<Real>(
      "residual_stiffness_ratio",
      0.0,
      "Residual-stiffness floor kappa; must match the ADDGCohesiveZone "
      "kernel's value for the diagnostic tractions to be consistent.");
  return params;
}

ADDGCohesiveZoneOutput::ADDGCohesiveZoneOutput(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _ndisp(coupledComponents("displacements")),
    _K(getParam<Real>("penalty_stiffness")),
    _T(getParam<Real>("normal_strength")),
    _S(getParam<Real>("shear_strength")),
    _GIc(getParam<Real>("GI_c")),
    _GIIc(getParam<Real>("GII_c")),
    _eta_BK(getParam<Real>("eta_BK")),
    _kappa(getParam<Real>("residual_stiffness_ratio")),
    _disp(_ndisp),
    _disp_neighbor(_ndisp),
    _normal_jump(declareProperty<Real>("normal_jump")),
    _shear_jump(declareProperty<Real>("shear_jump")),
    _damage(declareProperty<Real>("czm_damage")),
    _normal_traction(declareProperty<Real>("normal_traction")),
    _shear_traction(declareProperty<Real>("shear_traction")),
    _active_zone(declareProperty<Real>("czm_active_zone"))
{
  for (unsigned int j = 0; j < _ndisp; ++j)
  {
    _disp[j] = &adCoupledValue("displacements", j);
    _disp_neighbor[j] = &adCoupledNeighborValue("displacements", j);
  }
}

void
ADDGCohesiveZoneOutput::computeQpProperties()
{
  // jump = u_element - u_neighbor; n = element outward normal
  ADRealVectorValue jump;
  for (unsigned int j = 0; j < _ndisp; ++j)
    jump(j) = (*_disp[j])[_qp] - (*_disp_neighbor[j])[_qp];

  const RealVectorValue & n = _normals[_qp];
  const ADReal jump_dot_n = jump * n;
  const ADReal delta_n_open = -jump_dot_n;

  const ADRealVectorValue jump_n_vec = jump_dot_n * n;
  const ADRealVectorValue jump_s_vec = jump - jump_n_vec;
  const ADReal delta_s = sqrt(jump_s_vec * jump_s_vec + 1.0e-30);

  const ADReal delta_n_pos =
      (MetaPhysicL::raw_value(delta_n_open) > 0.0) ? delta_n_open : ADReal(0.0);
  const ADReal delta_eff =
      sqrt(delta_n_pos * delta_n_pos + delta_s * delta_s + 1.0e-30);

  const Real delta_0_n = _T / _K;
  const Real delta_0_s = _S / _K;

  ADReal delta_0_eff;
  ADReal delta_c_eff;
  if (MetaPhysicL::raw_value(delta_n_pos) > 1.0e-30)
  {
    const ADReal beta2 = (delta_s * delta_s) / (delta_n_pos * delta_n_pos);
    const ADReal num = (delta_0_n * delta_0_n) * (delta_0_s * delta_0_s) * (1.0 + beta2);
    const ADReal den = (delta_0_s * delta_0_s) + beta2 * (delta_0_n * delta_0_n);
    delta_0_eff = sqrt(num / den);
    const ADReal B = beta2 / (1.0 + beta2);
    const ADReal Gc_mix = _GIc + (_GIIc - _GIc) * pow(B, _eta_BK);
    delta_c_eff = 2.0 * Gc_mix / (_K * delta_0_eff + 1.0e-30);
  }
  else
  {
    delta_0_eff = ADReal(delta_0_s);
    delta_c_eff = 2.0 * _GIIc / (_K * delta_0_s + 1.0e-30);
  }

  ADReal damage_ad;
  if (MetaPhysicL::raw_value(delta_eff) <= MetaPhysicL::raw_value(delta_0_eff))
    damage_ad = ADReal(0.0);
  else if (MetaPhysicL::raw_value(delta_eff) < MetaPhysicL::raw_value(delta_c_eff))
    damage_ad = (delta_c_eff * (delta_eff - delta_0_eff)) /
                (delta_eff * (delta_c_eff - delta_0_eff));
  else
    damage_ad = ADReal(1.0);

  // Effective degradation factor with residual-stiffness floor (matches
  // ADDGCohesiveZone). At full damage the cohesive tangent saturates at
  // kappa*K instead of dropping to zero -- keeps Newton/AL well-posed.
  const ADReal factor = (1.0 - damage_ad) * (1.0 - _kappa) + _kappa;

  ADReal T_n_ad;
  if (MetaPhysicL::raw_value(jump_dot_n) < 0.0)
    T_n_ad = factor * _K * delta_n_open;
  else
    T_n_ad = -_K * jump_dot_n;  // = K * delta_n_open (still <= 0 in closing)

  const ADReal T_s_ad = factor * _K * delta_s;

  // Real-cast for aux output
  const Real damage_real = MetaPhysicL::raw_value(damage_ad);

  _normal_jump[_qp]      = MetaPhysicL::raw_value(delta_n_open);
  _shear_jump[_qp]       = MetaPhysicL::raw_value(delta_s);
  _damage[_qp]           = damage_real;
  _normal_traction[_qp]  = MetaPhysicL::raw_value(T_n_ad);
  _shear_traction[_qp]   = MetaPhysicL::raw_value(T_s_ad);
  // Active process zone: 1 where the TSL is actively softening, 0 otherwise.
  _active_zone[_qp]      = (damage_real > 0.0 && damage_real < 1.0) ? 1.0 : 0.0;
}
