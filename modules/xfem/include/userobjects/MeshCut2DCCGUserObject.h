//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MeshCut2DUserObjectBase.h"

/**
 * MeshCut2DCCGUserObject:
 * (1) reads in a mesh describing the crack surface,
 * (2) uses the mesh to do initial cutting of 2D elements, and
 * (3) grows the mesh by a growth rate determined by the C integral.
 */
class MeshCut2DCCGUserObject : public MeshCut2DUserObjectBase
{
public:
  static InputParameters validParams();

  MeshCut2DCCGUserObject(const InputParameters & parameters);

  virtual void initialize() override;
  virtual bool isCutterMeshChanged() const override;

  Real getTotalCrackLength() const;

protected:
  virtual void findActiveBoundaryGrowth() override;

private:
  /// Fracture-integral KI and KII values.
  const std::vector<Real> & _ki_vpp;
  const std::vector<Real> & _kii_vpp;

  /// C-integral values.
  const std::vector<Real> & _c_vpp;

  /// Accumulated crack extension length.
  Real _total_crack_length;

  /// Paris-law coefficients for da = A C^m dt.
  const Real _paris_coeff;
  const Real _paris_exponent;
};
