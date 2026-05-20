//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IntegratedBC.h"

class Function;

/**
 * Weak Dirichlet BC for scalar diffusion in DG following the Abedi/Rabedi
 * star-value formulation (Note 8 eq. 113 / eq. 115) with
 *   T*       = T_bar   (prescribed)
 *   q* . n   = q . n   (no replacement -- use the computed flux directly)
 *
 * Substituting these star values into eq. (113) and using the f(w) = -kappa
 * grad w test multiplier convention (the same convention as DGDiffusion):
 *   B^e_u(v, u) = -INT v (q . n) dS  +  eps INT (kappa grad v . n)(T_bar - u) dS
 *
 * With q = -kappa grad u and bringing the residual into MOOSE's standard
 * sign convention:
 *   r = - kappa grad u . n * test
 *       + eps * (u - T_bar) * kappa grad test . n
 *
 * NO sigma/h * (u - T_bar) penalty term: that contribution is optional in
 * the Abedi derivation (added in DGFunctionDiffusionDirichletBC for SIPG
 * coercivity but not part of the star-value construction itself). This
 * class is the pure Abedi Dirichlet face contribution.
 */
class DGAbediDirichletBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  DGAbediDirichletBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const Function & _func;
  const Real _epsilon;
  const MaterialProperty<Real> & _diff;
};
