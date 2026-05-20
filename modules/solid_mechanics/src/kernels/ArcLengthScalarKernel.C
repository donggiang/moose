//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ArcLengthScalarKernel.h"

registerMooseObject("SolidMechanicsApp", ArcLengthScalarKernel);

InputParameters
ArcLengthScalarKernel::validParams()
{
  InputParameters params = KernelScalarBase::validParams();
  params.addClassDescription(
      "Cylindrical (Riks) arc-length constraint kernel. Couples one displacement "
      "component to the load-parameter scalar λ via the constraint "
      "∫|Δu|² dΩ − Δs² = 0. One instance is required per displacement component; "
      "exactly ONE instance per λ must set `primary_constraint = true` so the "
      "−Δs² constant is included.");
  params.renameCoupledVar(
      "scalar_variable", "lambda", "Load-parameter scalar variable (the AL kappa).");
  params.addRequiredParam<Real>("delta_s",
                                "Prescribed arc-length step Δs (per outer step).");
  params.addRequiredParam<PostprocessorName>(
      "volume_pp",
      "Name of a Postprocessor giving the domain volume V (used to spread the "
      "constant −Δs²/V across the per-QP integrand so that ∫(−Δs²/V) dΩ = −Δs²).");
  params.addRequiredParam<bool>(
      "primary_constraint",
      "Set to true on EXACTLY ONE kernel instance per λ — that one carries the "
      "−Δs² term. Other instances (one per additional displacement component) "
      "should set this to false; they only contribute |Δu_c|².");
  return params;
}

ArcLengthScalarKernel::ArcLengthScalarKernel(const InputParameters & parameters)
  : KernelScalarBase(parameters),
    _delta_s(getParam<Real>("delta_s")),
    _volume(getPostprocessorValue("volume_pp")),
    _primary_constraint(getParam<bool>("primary_constraint")),
    _u_old_qp(_var.slnOld())
{
}

Real
ArcLengthScalarKernel::computeScalarQpResidual()
{
  const Real du = _u[_qp] - _u_old_qp[_qp];
  Real r = du * du;
  if (_primary_constraint && _volume > 0.0)
    r -= (_delta_s * _delta_s) / _volume;
  return r;
}

Real
ArcLengthScalarKernel::computeScalarQpOffDiagJacobian(const unsigned int jvar_num)
{
  if (jvar_num != _var.number())
    return 0.0;
  const Real du = _u[_qp] - _u_old_qp[_qp];
  // d/du_j of (u-u_old)^2 = 2(u-u_old) * phi_j(qp)
  return 2.0 * du * _phi[_j][_qp];
}
