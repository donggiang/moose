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
 * (1) reads in a mesh describing the crack surface
 * (2) uses the mesh to do initial cutting of 2D elements, and
 * (3) grows the mesh by a growth rate determined by the C integral
 */

class MeshCut2DCCGUserObject : public MeshCut2DUserObjectBase
{
public:
  static InputParameters validParams();

  MeshCut2DCCGUserObject(const InputParameters & parameters);

  virtual void initialize() override;
  virtual bool isCutterMeshChanged() const override;

  Real get_total_crack_length() const;

protected:
  virtual void findActiveBoundaryGrowth() override;

private:
  /// amount to grow crack by for each xfem update step
  const Real & _growth_increment;

  /// Pointer fracture integral ki if available
  const std::vector<Real> * const _ki_vpp;
  /// Pointer fracture integral kii if available
  const std::vector<Real> * const _kii_vpp;
  /// Fracture integral C*
  const std::vector<Real> & _c_vpp;
  Real _total_crack_length;
  /// critical k value for crack growth
  const Real _paris_coeff;
  /// Maximum stress criterion threshold for crack growth.
  const Real _paris_exponent;
  Real ci_old;
};
