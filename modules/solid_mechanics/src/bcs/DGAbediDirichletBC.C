//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DGAbediDirichletBC.h"

#include "Function.h"
#include "MooseVariableFE.h"

registerMooseObject("SolidMechanicsApp", DGAbediDirichletBC);

InputParameters
DGAbediDirichletBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addClassDescription(
      "Weak Dirichlet BC for scalar diffusion in DG following Abedi/Rabedi "
      "star-value formulation (T* = T_bar, q* = q; no sigma/h penalty).");
  params.addRequiredParam<FunctionName>(
      "function", "Prescribed Dirichlet value T_bar(x, t) on the boundary.");
  params.addRequiredParam<Real>(
      "epsilon",
      "DG variant: -1 SIPG, +1 NIPG, 0 IIPG (matches DGDiffusion convention).");
  params.addParam<MaterialPropertyName>(
      "diff", 1, "Diffusion coefficient material property (kappa).");
  return params;
}

DGAbediDirichletBC::DGAbediDirichletBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _func(getFunction("function")),
    _epsilon(getParam<Real>("epsilon")),
    _diff(getMaterialProperty<Real>("diff"))
{
}

Real
DGAbediDirichletBC::computeQpResidual()
{
  const Real fn = _func.value(_t, _q_point[_qp]);

  Real r = 0.0;
  // Consistency: -INT v (q* . n) dS with q* = q = -kappa grad u, so
  //   -INT v (-kappa grad u . n) = +INT v kappa grad u . n
  // MOOSE residual-sign convention (matching DGDiffusion / the existing
  // DGFunctionDiffusionDirichletBC): subtract kappa grad u . n * test.
  r -= _diff[_qp] * _grad_u[_qp] * _normals[_qp] * _test[_i][_qp];

  // Adjoint-consistency (last term of eq. 113) with T* = T_bar:
  //   eps INT (kappa grad v . n) (T_bar - u) dS
  // In MOOSE residual-sign convention: + eps (u - T_bar) kappa grad v . n.
  r += _epsilon * (_u[_qp] - fn) * _diff[_qp] * _grad_test[_i][_qp] * _normals[_qp];

  // NOTE: no sigma/h * (u - T_bar) * v penalty. That term is optional in
  // Abedi's derivation; this class is the pure star-value Dirichlet form.

  return r;
}

Real
DGAbediDirichletBC::computeQpJacobian()
{
  Real r = 0.0;
  r -= _diff[_qp] * _grad_phi[_j][_qp] * _normals[_qp] * _test[_i][_qp];
  r += _epsilon * _phi[_j][_qp] * _diff[_qp] * _grad_test[_i][_qp] * _normals[_qp];
  return r;
}
