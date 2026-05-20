//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADDGElasticity.h"
#include "MooseVariableFE.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

#include "libmesh/utility.h"

registerMooseObject("SolidMechanicsApp", ADDGElasticity);

InputParameters
ADDGElasticity::validParams()
{
  InputParameters params = ADDGKernel::validParams();
  params.addClassDescription(
      "DG kernel for vector linear elasticity (SIPG / NIPG / IIPG / OBB) "
      "following Liu, Wheeler, Dawson (2009). Penalty is eta * G / h_F.");
  params.addRequiredParam<unsigned int>(
      "component",
      "Displacement component this kernel operates on: 0=x, 1=y, 2=z.");
  params.addRequiredCoupledVar(
      "displacements", "Vector of displacement variables (1 to 3 components).");
  params.addRequiredParam<Real>(
      "epsilon",
      "DG variant: -1 SIPG (symmetric), +1 NIPG (non-symmetric), 0 IIPG (incomplete).");
  params.addRequiredParam<Real>(
      "eta",
      "Dimensionless penalty parameter (Liu 2009 delta_p), O(p^2). "
      "Penalty coefficient = eta * G / h_F.");
  params.addParam<MaterialPropertyName>(
      "penalty_modulus",
      "shear_modulus",
      "Material property providing the penalty's stiffness scale (Liu's G). "
      "Default: 'shear_modulus'.");
  params.addParam<std::string>("base_name", "Base name for material properties.");
  return params;
}

ADDGElasticity::ADDGElasticity(const InputParameters & parameters)
  : ADDGKernel(parameters),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _eta(getParam<Real>("eta")),
    _epsilon(getParam<Real>("epsilon")),
    _stress(getADMaterialProperty<RankTwoTensor>(
        (isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "") + "stress")),
    _stress_neighbor(getNeighborADMaterialProperty<RankTwoTensor>(
        (isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "") + "stress")),
    _elasticity_tensor(getADMaterialProperty<RankFourTensor>(
        (isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "") +
        "elasticity_tensor")),
    _elasticity_tensor_neighbor(getNeighborADMaterialProperty<RankFourTensor>(
        (isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "") +
        "elasticity_tensor")),
    _penalty_modulus(getADMaterialProperty<Real>("penalty_modulus")),
    _penalty_modulus_neighbor(getNeighborADMaterialProperty<Real>("penalty_modulus")),
    _disp(_ndisp),
    _disp_neighbor(_ndisp)
{
  if (_component >= _ndisp)
    paramError("component",
               "Component index ",
               _component,
               " is out of bounds for ",
               _ndisp,
               "-component displacement vector.");

  for (unsigned int j = 0; j < _ndisp; ++j)
  {
    _disp[j] = &adCoupledValue("displacements", j);
    _disp_neighbor[j] = &adCoupledNeighborValue("displacements", j);
  }
}

ADReal
ADDGElasticity::computeQpResidual(Moose::DGResidualType type)
{
  // Face characteristic size (matches DGDiffusion convention).
  const int order = std::max(libMesh::Order(1), _var.order());
  const Real h_F =
      _current_elem_volume / _current_side_volume * 1.0 / Utility::pow<2>(order);

  // Average normal traction at QP from stress materials (consistency contribution).
  // {sigma . n} = 0.5 (sigma_+ + sigma_-) . n_s
  const ADRealVectorValue traction_avg = 0.5 * (_stress[_qp] * _normals[_qp] +
                                                _stress_neighbor[_qp] * _normals[_qp]);

  // Vector jump in displacement at QP: [u]_j = u_+|_j - u_-|_j
  ADRealVectorValue u_jump;
  for (unsigned int j = 0; j < _ndisp; ++j)
    u_jump(j) = (*_disp[j])[_qp] - (*_disp_neighbor[j])[_qp];

  // Average penalty modulus (Liu's G).
  const ADReal G_avg = 0.5 * (_penalty_modulus[_qp] + _penalty_modulus_neighbor[_qp]);
  const ADReal penalty_coef = _eta * G_avg / h_F;

  ADReal r = 0.0;

  switch (type)
  {
    case Moose::Element:
    {
      // Test on element side. Consistency: contribution to v_+
      // -<sigma . n>_component * v_test (component = _component)
      r -= traction_avg(_component) * _test[_i][_qp];

      // Adjoint: + epsilon * 0.5 * (sigma(v_+) . n) . [u]
      // sigma(v_+) = C_elem : eps(v_+); v_+ = e_component * test_scalar
      const ADRankTwoTensor grad_v = buildGradV(_grad_test[_i][_qp]);
      const ADRankTwoTensor strain_v = 0.5 * (grad_v + grad_v.transpose());
      const ADRankTwoTensor stress_v = _elasticity_tensor[_qp] * strain_v;
      const ADRealVectorValue traction_v = stress_v * _normals[_qp];
      r += _epsilon * 0.5 * traction_v * u_jump;

      // Penalty: + (eta * G / h_F) * [u]_component * v_test
      r += penalty_coef * u_jump(_component) * _test[_i][_qp];
      break;
    }

    case Moose::Neighbor:
    {
      // Test on neighbor side. Consistency: contribution to v_-
      // +<sigma . n>_component * v_test_neighbor
      r += traction_avg(_component) * _test_neighbor[_i][_qp];

      // Adjoint: + epsilon * 0.5 * (sigma(v_-) . n) . [u]
      // sigma(v_-) = C_neighbor : eps(v_-); v_- = e_component * test_scalar_neighbor
      const ADRankTwoTensor grad_v = buildGradV(_grad_test_neighbor[_i][_qp]);
      const ADRankTwoTensor strain_v = 0.5 * (grad_v + grad_v.transpose());
      const ADRankTwoTensor stress_v = _elasticity_tensor_neighbor[_qp] * strain_v;
      const ADRealVectorValue traction_v = stress_v * _normals[_qp];
      r += _epsilon * 0.5 * traction_v * u_jump;

      // Penalty: -(eta * G / h_F) * [u]_component * v_test_neighbor
      r -= penalty_coef * u_jump(_component) * _test_neighbor[_i][_qp];
      break;
    }
  }

  return r;
}
