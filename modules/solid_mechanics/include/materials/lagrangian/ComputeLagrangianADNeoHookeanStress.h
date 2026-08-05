//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "ComputeLagrangianStressPK1.h"

/// Compressible Neo-Hookean PK1 stress and tangent differentiated from a local energy potential
class ComputeLagrangianADNeoHookeanStress : public ComputeLagrangianStressPK1
{
public:
  static InputParameters validParams();
  ComputeLagrangianADNeoHookeanStress(const InputParameters & parameters);

protected:
  void initialSetup() override;
  void computeQpProperties() override;
  void computeQpPK1Stress() override;

  const unsigned int _ndisp;
  const std::vector<const VariableValue *> _disp;
  const std::vector<const VariableGradient *> _grad_disp;
  const Moose::CoordinateSystemType _coord_system;
  const bool _stabilize_strain;
  const MaterialProperty<Real> & _lambda;
  const MaterialProperty<Real> & _mu;
};
