//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADDGAbediSidesetReaction.h"

#include "Function.h"
#include "FEProblemBase.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

#include "metaphysicl/raw_type.h"

registerMooseObject("SolidMechanicsApp", ADDGAbediSidesetReaction);

InputParameters
ADDGAbediSidesetReaction::validParams()
{
  InputParameters params = SideIntegralPostprocessor::validParams();
  params.addClassDescription(
      "Reaction force on a DG Dirichlet sideset extracted from the Abedi "
      "star-value bilinear form (consistency + adjoint), parameterized by "
      "a virtual displacement v(x) = direction * f(x).");
  params.addRequiredParam<RealVectorValue>(
      "direction", "Unit vector picking out the reaction component (e.g., '0 1 0').");
  params.addRequiredParam<Real>("epsilon",
                                "DG variant; must match the BC epsilon. "
                                "-1 SIPG, +1 NIPG, 0 IIPG.");
  params.addRequiredCoupledVar("displacements", "Vector of displacement variables.");
  params.addParam<std::vector<FunctionName>>(
      "function",
      "Prescribed displacement u_bar (one Function per component). Defaults to zero.");
  params.addParam<FunctionName>(
      "test_lifting",
      "Optional scalar lifting f(x) for the virtual displacement. "
      "If omitted, f = 1 (constant) and the adjoint term vanishes "
      "(result matches ADSidesetReaction up to sign).");
  params.addParam<MaterialPropertyName>("stress_tensor", "stress", "Stress material property name.");
  params.addParam<MaterialPropertyName>(
      "elasticity_tensor", "elasticity_tensor", "Elasticity tensor material property name.");
  params.addParam<std::string>("base_name", "Base name for material properties.");
  return params;
}

ADDGAbediSidesetReaction::ADDGAbediSidesetReaction(const InputParameters & parameters)
  : SideIntegralPostprocessor(parameters),
    _ndisp(coupledComponents("displacements")),
    _direction(getParam<RealVectorValue>("direction")),
    _epsilon(getParam<Real>("epsilon")),
    _stress(getADMaterialProperty<RankTwoTensor>(
        (isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "") +
        getParam<MaterialPropertyName>("stress_tensor"))),
    _elasticity_tensor(getADMaterialProperty<RankFourTensor>(
        (isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "") +
        getParam<MaterialPropertyName>("elasticity_tensor"))),
    _disp(_ndisp),
    _u_bar(_ndisp, nullptr),
    _test_lifting(isParamValid("test_lifting")
                      ? &_fe_problem.getFunction(getParam<FunctionName>("test_lifting"), _tid)
                      : nullptr)
{
  for (unsigned int j = 0; j < _ndisp; ++j)
    _disp[j] = &adCoupledValue("displacements", j);

  if (isParamValid("function"))
  {
    const auto fn_names = getParam<std::vector<FunctionName>>("function");
    if (fn_names.size() != _ndisp)
      paramError("function",
                 "Must provide exactly ",
                 _ndisp,
                 " function names (one per displacement component); got ",
                 fn_names.size());
    for (unsigned int j = 0; j < _ndisp; ++j)
      _u_bar[j] = &_fe_problem.getFunction(fn_names[j], _tid);
  }
}

Real
ADDGAbediSidesetReaction::computeQpIntegral()
{
  // Deviation from imposed displacement: u_h - u_bar (defaults u_bar=0).
  ADRealVectorValue u_dev;
  for (unsigned int j = 0; j < _ndisp; ++j)
  {
    const Real u_bar_j = _u_bar[j] ? _u_bar[j]->value(_t, _q_point[_qp]) : 0.0;
    u_dev(j) = (*_disp[j])[_qp] - u_bar_j;
  }

  // Virtual displacement v(x) = direction * f(x).
  const Real f_val = _test_lifting ? _test_lifting->value(_t, _q_point[_qp]) : 1.0;
  const RealGradient grad_f =
      _test_lifting ? _test_lifting->gradient(_t, _q_point[_qp]) : RealGradient(0, 0, 0);

  const RealVectorValue v_field = _direction * f_val;

  // grad_v(i,p) = direction(i) * (grad_f)_p ; strain_v = sym(grad_v).
  RankTwoTensor grad_v_real;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int p = 0; p < LIBMESH_DIM; ++p)
      grad_v_real(i, p) = _direction(i) * grad_f(p);

  const RankTwoTensor strain_v_real = 0.5 * (grad_v_real + grad_v_real.transpose());
  // Promote to AD scalar product so we can multiply by AD elasticity tensor.
  const ADRankTwoTensor strain_v(strain_v_real);
  const ADRankTwoTensor stress_v = _elasticity_tensor[_qp] * strain_v;
  const ADRealVectorValue traction_v = stress_v * _normals[_qp];

  // Consistency: -(sigma(u_h) . n) . v_field
  const ADRealVectorValue sigma_n = _stress[_qp] * _normals[_qp];
  ADReal r = 0.0;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    r -= sigma_n(i) * v_field(i);

  // Adjoint: epsilon * (sigma(v) . n) . (u_h - u_bar)
  for (unsigned int j = 0; j < _ndisp; ++j)
    r += _epsilon * traction_v(j) * u_dev(j);

  return MetaPhysicL::raw_value(r);
}
