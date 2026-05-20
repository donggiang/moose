//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceMaterial.h"

/**
 * Diagnostic-output InterfaceMaterial that mirrors the constitutive law of
 * ADDGCohesiveZone. Runs on the bonded cohesive face, computes the TSL state
 * (normal/shear jump, damage, normal/shear traction) from the displacement
 * jump, and declares them as MaterialProperty<Real> on the interface.
 *
 * Pair with `(AD)MaterialRealAux` (boundary = bonded_iface, check_boundary_
 * restricted = false) to project each value into a `MONOMIAL CONSTANT`
 * boundary-restricted AuxVariable for ParaView visualization. No effect on
 * the residual.
 */
class ADDGCohesiveZoneOutput : public InterfaceMaterial
{
public:
  static InputParameters validParams();
  ADDGCohesiveZoneOutput(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Number of displacement components
  const unsigned int _ndisp;

  /// Cohesive parameters (must match those of ADDGCohesiveZone)
  const Real _K;
  const Real _T;
  const Real _S;
  const Real _GIc;
  const Real _GIIc;
  const Real _eta_BK;
  /// Residual stiffness floor (same as ADDGCohesiveZone kappa). 0 = pure
  /// bilinear TSL; non-zero keeps the tractions reported here at
  /// (1-d)(1-kappa) + kappa of pre-peak instead of dropping to 0 at d=1.
  const Real _kappa;

  /// Coupled displacements on element and neighbor sides (face QPs)
  std::vector<const ADVariableValue *> _disp;
  std::vector<const ADVariableValue *> _disp_neighbor;

  /// Diagnostic output material properties (Real, not AD: aux output only).
  MaterialProperty<Real> & _normal_jump;
  MaterialProperty<Real> & _shear_jump;
  MaterialProperty<Real> & _damage;
  MaterialProperty<Real> & _normal_traction;
  MaterialProperty<Real> & _shear_traction;

  /// 1 where the cohesive zone is actively softening (0 < damage < 1),
  /// 0 elsewhere (intact ahead of the tip OR fully cracked behind the tip).
  /// Threshold on this in ParaView to isolate the process zone width.
  MaterialProperty<Real> & _active_zone;
};
