//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ALCoupledScaledNeumannBC.h"

#include "Function.h"

registerMooseObject("SolidMechanicsApp", ALCoupledScaledNeumannBC);

InputParameters
ALCoupledScaledNeumannBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addClassDescription(
      "Neumann boundary condition whose magnitude is scaled by a coupled scalar "
      "variable λ. Companion to ArcLengthScalarKernel for input-only arc-length "
      "(Path C). Residual contribution: -λ · function · test.");
  params.addRequiredParam<FunctionName>(
      "function", "Reference external load f_ext(x, t).");
  params.addRequiredCoupledVar(
      "arc_length_scalar", "Coupled scalar variable carrying the load parameter λ.");
  return params;
}

ALCoupledScaledNeumannBC::ALCoupledScaledNeumannBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _func(getFunction("function")),
    _lambda(coupledScalarValue("arc_length_scalar")),
    _lambda_var(coupledScalar("arc_length_scalar"))
{
}

Real
ALCoupledScaledNeumannBC::computeQpResidual()
{
  return -_lambda[0] * _test[_i][_qp] * _func.value(_t, _q_point[_qp]);
}

Real
ALCoupledScaledNeumannBC::computeQpOffDiagJacobianScalar(unsigned int svar)
{
  if (svar == _lambda_var)
    return -_test[_i][_qp] * _func.value(_t, _q_point[_qp]);
  return 0.0;
}
