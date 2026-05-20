//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADDGElasticityDirichletBC.h"
#include "Function.h"
#include "FEProblemBase.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

#include "libmesh/utility.h"

registerMooseObject("SolidMechanicsApp", ADDGElasticityDirichletBC);

InputParameters
ADDGElasticityDirichletBC::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addClassDescription(
      "Weak Dirichlet BC for vector linear elasticity in DG formulation "
      "(SIPG/NIPG/IIPG/OBB), Liu-Wheeler-Dawson 2009.");
  params.addRequiredParam<unsigned int>(
      "component", "Displacement component this BC operates on: 0=x, 1=y, 2=z.");
  params.addRequiredCoupledVar(
      "displacements", "Vector of displacement variables (1 to 3 components).");
  params.addRequiredParam<Real>("epsilon",
                                "DG variant: -1 SIPG, +1 NIPG, 0 IIPG.");
  params.addRequiredParam<Real>(
      "eta", "Dimensionless penalty parameter (Liu's delta_p).");
  params.addRequiredParam<std::vector<FunctionName>>(
      "function",
      "Functions providing the prescribed displacement, one per component "
      "(e.g. 'u_x_bar u_y_bar u_z_bar' in 3D).");
  params.addParam<MaterialPropertyName>(
      "penalty_modulus", "shear_modulus",
      "Material property providing the penalty stiffness scale (Liu's G).");
  params.addParam<std::string>("base_name", "Base name for material properties.");
  return params;
}

ADDGElasticityDirichletBC::ADDGElasticityDirichletBC(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _eta(getParam<Real>("eta")),
    _epsilon(getParam<Real>("epsilon")),
    _stress(getADMaterialProperty<RankTwoTensor>(
        (isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "") + "stress")),
    _elasticity_tensor(getADMaterialProperty<RankFourTensor>(
        (isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "") +
        "elasticity_tensor")),
    _penalty_modulus(getADMaterialProperty<Real>("penalty_modulus")),
    _disp(_ndisp)
{
  if (_component >= _ndisp)
    paramError("component",
               "component=",
               _component,
               " out of bounds for ",
               _ndisp,
               "-component displacement vector.");

  const auto fn_names = getParam<std::vector<FunctionName>>("function");
  if (fn_names.size() != _ndisp)
    paramError("function",
               "Must provide exactly ",
               _ndisp,
               " function names (one per displacement component); got ",
               fn_names.size());

  _u_bar.resize(_ndisp);
  for (unsigned int j = 0; j < _ndisp; ++j)
  {
    _u_bar[j] = &_fe_problem.getFunction(fn_names[j], _tid);
    _disp[j] = &adCoupledValue("displacements", j);
  }
}

ADReal
ADDGElasticityDirichletBC::computeQpResidual()
{
  // Boundary penalty per Wheeler-Liu-Dawson 2009 eq. (24):
  //   delta_p = eta * G / h_F                        (with eta multiplier)
  // The Dirichlet face residual carries this penalty as
  //   + delta_p * (u - u_bar)_component * test.

  // Face characteristic size (matches DGFunctionDiffusionDirichletBC convention)
  const int order = std::max(libMesh::Order(1), _var.order());
  const Real h_F =
      _current_elem_volume / _current_side_volume * 1.0 / Utility::pow<2>(order);

  // Vector deviation (u - u_bar) at QP across all components
  ADRealVectorValue u_dev;
  for (unsigned int j = 0; j < _ndisp; ++j)
  {
    const Real u_bar_j = _u_bar[j]->value(_t, _q_point[_qp]);
    u_dev(j) = (*_disp[j])[_qp] - u_bar_j;
  }

  // BC penalty coefficient: eta * G / h_F (Wheeler 2009 eq. 24)
  const ADReal penalty_coef = _eta * _penalty_modulus[_qp] / h_F;

  // 1) Consistency: -(sigma . n)_component * test
  ADReal r = -(_stress[_qp] * _normals[_qp])(_component) * _test[_i][_qp];

  // 2) Adjoint: + epsilon * (sigma(v) . n) . (u - u_bar)
  // sigma(v) = C : eps(v); v = e_component * test_scalar
  const ADRankTwoTensor grad_v = buildGradV(_grad_test[_i][_qp]);
  const ADRankTwoTensor strain_v = 0.5 * (grad_v + grad_v.transpose());
  const ADRankTwoTensor stress_v = _elasticity_tensor[_qp] * strain_v;
  const ADRealVectorValue traction_v = stress_v * _normals[_qp];
  r += _epsilon * traction_v * u_dev;

  // 3) Penalty (Wheeler 2009 eq. 24): + eta*(G / h_F) * (u - u_bar)_component * test
  r += penalty_coef * u_dev(_component) * _test[_i][_qp];

  return r;
}
