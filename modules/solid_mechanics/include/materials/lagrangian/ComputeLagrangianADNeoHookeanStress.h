//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "ComputeLagrangianStressPK2.h"

/// Compressible Neo-Hookean hyperelasticity differentiated from a local strain-energy potential
class ComputeLagrangianADNeoHookeanStress : public ComputeLagrangianStressPK2
{
public:
  static InputParameters validParams();
  ComputeLagrangianADNeoHookeanStress(const InputParameters & parameters);

protected:
  void computeQpPK2Stress() override;

  const MaterialProperty<Real> & _lambda;
  const MaterialProperty<Real> & _mu;
};
