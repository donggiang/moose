//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADDGCohesiveZone.h"

#include "metaphysicl/raw_type.h"

#include <cmath>
using std::sqrt;
using std::pow;

registerMooseObject("SolidMechanicsApp", ADDGCohesiveZone);

InputParameters
ADDGCohesiveZone::validParams()
{
  InputParameters params = ADInterfaceKernel::validParams();
  params.addClassDescription(
      "DG cohesive zone (Neumann-type star traction) on an interior face. "
      "Bilinear traction-separation with Benzeggagh-Kenane mixed-mode. "
      "Removes the need for BreakMeshByBlock + CG CohesiveZone block.");
  params.addRequiredParam<unsigned int>("component",
                                        "Displacement component (0=x, 1=y, 2=z).");
  params.addRequiredCoupledVar("displacements", "Vector of displacement variables.");
  params.addRequiredParam<Real>("penalty_stiffness", "Cohesive penalty stiffness K [F/L^3].");
  params.addRequiredParam<Real>("normal_strength", "Peak normal traction T [F/L^2].");
  params.addRequiredParam<Real>("shear_strength", "Peak shear traction S [F/L^2].");
  params.addRequiredParam<Real>("GI_c", "Mode-I fracture energy [F/L].");
  params.addRequiredParam<Real>("GII_c", "Mode-II fracture energy [F/L].");
  params.addRequiredParam<Real>("eta_BK", "Benzeggagh-Kenane mixed-mode exponent.");
  params.addParam<Real>(
      "viscosity",
      0.0,
      "Viscous regularization coefficient [F.s/L^3]. T_visc = -viscosity * "
      "d/dt [[u]] added in parallel to the static TSL. Helps Newton "
      "traverse the post-peak snap. 0 = rate-independent law.");
  params.addParam<Real>(
      "residual_stiffness_ratio",
      0.0,
      "Residual-stiffness floor kappa in [0, 1). Effective degradation "
      "factor becomes (1-d)*(1-kappa) + kappa, so the cohesive tangent "
      "never drops below kappa*K. Keeps Newton well-posed past the snap. "
      "0 = classical fully-degraded TSL.");
  return params;
}

ADDGCohesiveZone::ADDGCohesiveZone(const InputParameters & parameters)
  : ADInterfaceKernel(parameters),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _K(getParam<Real>("penalty_stiffness")),
    _T(getParam<Real>("normal_strength")),
    _S(getParam<Real>("shear_strength")),
    _GIc(getParam<Real>("GI_c")),
    _GIIc(getParam<Real>("GII_c")),
    _eta_BK(getParam<Real>("eta_BK")),
    _viscosity(getParam<Real>("viscosity")),
    _kappa(getParam<Real>("residual_stiffness_ratio")),
    _disp(_ndisp),
    _disp_neighbor(_ndisp),
    _disp_old(_ndisp),
    _disp_neighbor_old(_ndisp)
{
  if (_component >= _ndisp)
    paramError("component",
               "component=",
               _component,
               " out of bounds for displacements size ",
               _ndisp);

  for (unsigned int j = 0; j < _ndisp; ++j)
  {
    _disp[j] = &adCoupledValue("displacements", j);
    _disp_neighbor[j] = &adCoupledNeighborValue("displacements", j);
    _disp_old[j] = &coupledValueOld("displacements", j);
    _disp_neighbor_old[j] = &coupledNeighborValueOld("displacements", j);
  }
}

ADReal
ADDGCohesiveZone::computeQpResidual(Moose::DGResidualType type)
{
  // Jump vector [[u]] = u_element - u_neighbor   (MOOSE convention).
  ADRealVectorValue jump;
  for (unsigned int j = 0; j < _ndisp; ++j)
    jump(j) = (*_disp[j])[_qp] - (*_disp_neighbor[j])[_qp];

  // Element-side outward normal n. n points from element into neighbor;
  // so jump.n > 0 means element moves TOWARD neighbor (closing) and
  // delta_n_open = -(jump.n) > 0 in opening.
  const RealVectorValue & n = _normals[_qp];
  const ADReal jump_dot_n = jump * n;
  const ADReal delta_n_open = -jump_dot_n;

  // Decompose jump into normal and shear vector parts.
  const ADRealVectorValue jump_n_vec = jump_dot_n * n;
  const ADRealVectorValue jump_s_vec = jump - jump_n_vec;
  const ADReal delta_s = sqrt(jump_s_vec * jump_s_vec + 1.0e-30);

  // Effective separation drives damage. Only OPENING contributes to it;
  // closing (delta_n_open < 0) is elastic compression -- no damage.
  const ADReal delta_n_pos =
      (MetaPhysicL::raw_value(delta_n_open) > 0.0) ? delta_n_open : ADReal(0.0);
  const ADReal delta_eff =
      sqrt(delta_n_pos * delta_n_pos + delta_s * delta_s + 1.0e-30);

  // Single-mode onset and critical jumps.
  const Real delta_0_n = _T / _K;
  const Real delta_0_s = _S / _K;

  // Mixed-mode onset and critical separations (Turon 2006, Eqs. 24-27).
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
    // Pure mode II (no opening): use mode-II onset/critical directly.
    delta_0_eff = ADReal(delta_0_s);
    delta_c_eff = 2.0 * _GIIc / (_K * delta_0_s + 1.0e-30);
  }

  // Damage from the bilinear softening law (monotonic; no history tracking).
  ADReal damage;
  if (MetaPhysicL::raw_value(delta_eff) <= MetaPhysicL::raw_value(delta_0_eff))
    damage = ADReal(0.0);
  else if (MetaPhysicL::raw_value(delta_eff) < MetaPhysicL::raw_value(delta_c_eff))
    damage = (delta_c_eff * (delta_eff - delta_0_eff)) /
             (delta_eff * (delta_c_eff - delta_0_eff));
  else
    damage = ADReal(1.0);

  // Effective degradation factor with optional residual-stiffness floor:
  //   factor(d) = (1-d)*(1-kappa) + kappa     [classic TSL when kappa=0]
  // Keeps the cohesive tangent at kappa*K when d=1 instead of dropping to 0,
  // so the snap-back tangent stays well-conditioned for Newton/AL.
  const ADReal factor = (1.0 - damage) * (1.0 - _kappa) + _kappa;

  // Cohesive traction on element side (Newton's 3rd law: opposite on neighbor):
  //   For opening: T_coh_e = -factor * K (jump_n_vec + jump_s_vec)
  //   For closing: T_coh_e = -K jump_n_vec - factor * K jump_s_vec
  // (Normal compression remains elastic; shear is damaged regardless.)
  ADRealVectorValue T_coh_element;
  if (MetaPhysicL::raw_value(jump_dot_n) < 0.0)
    T_coh_element = -factor * _K * (jump_n_vec + jump_s_vec);
  else
    T_coh_element = -_K * jump_n_vec - factor * _K * jump_s_vec;

  // Optional viscous (Kelvin-Voigt) damper in parallel with the static TSL.
  // T_visc = -viscosity * d/dt [[u]]; uses backward-Euler rate against the
  // previous converged step. The static cohesive AD Jacobian is kept exact
  // (jump_old is non-AD).
  if (_viscosity > 0.0 && _dt > 0.0)
  {
    ADRealVectorValue jump_rate;
    for (unsigned int j = 0; j < _ndisp; ++j)
    {
      const Real jump_old_j = (*_disp_old[j])[_qp] - (*_disp_neighbor_old[j])[_qp];
      jump_rate(j) = (jump(j) - jump_old_j) / _dt;
    }
    T_coh_element -= _viscosity * jump_rate;
  }

  // Bilinear-form contribution per side:
  //   element side:  -T_coh_element[c] * test
  //   neighbor side: +T_coh_element[c] * test_neighbor   (Newton's 3rd law)
  if (type == Moose::Element)
    return -T_coh_element(_component) * _test[_i][_qp];
  else
    return T_coh_element(_component) * _test_neighbor[_i][_qp];
}
