//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "ComputeLagrangianADNeoHookeanStress.h"

#include "metaphysicl/dualnumberarray.h"

#include <array>

registerMooseObject("SolidMechanicsApp", ComputeLagrangianADNeoHookeanStress);

namespace
{
using LocalDerivative = MetaPhysicL::NumberArray<RankTwoTensor::N2, Real>;
using LocalADReal = DualNumber<Real, LocalDerivative, true>;
using LocalHessianDerivative = MetaPhysicL::NumberArray<RankTwoTensor::N2, LocalADReal>;
using LocalHessianADReal = DualNumber<LocalADReal, LocalHessianDerivative, true>;
}

InputParameters
ComputeLagrangianADNeoHookeanStress::validParams()
{
  InputParameters params = ComputeLagrangianStressPK1::validParams();
  params.addClassDescription(
      "Compute a compressible Neo-Hookean PK1 stress and tangent by locally differentiating the "
      "strain-energy potential with respect to the deformation measure derived from the "
      "displacement gradient.");
  params.addRequiredCoupledVar("displacements", "Displacement variables for the problem.");
  params.addParam<bool>("stabilize_strain", false, "Use the F-bar deformation gradient.");
  params.addParam<MaterialPropertyName>(
      "lambda", "lambda", "First Lame parameter material property.");
  params.addParam<MaterialPropertyName>("mu", "mu", "Shear modulus material property.");
  return params;
}

ComputeLagrangianADNeoHookeanStress::ComputeLagrangianADNeoHookeanStress(
    const InputParameters & parameters)
  : ComputeLagrangianStressPK1(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp(coupledValues("displacements")),
    _grad_disp(coupledGradients("displacements")),
    _coord_system(getBlockCoordSystem()),
    _stabilize_strain(getParam<bool>("stabilize_strain")),
    _lambda(getMaterialProperty<Real>(getParam<MaterialPropertyName>("lambda"))),
    _mu(getMaterialProperty<Real>(getParam<MaterialPropertyName>("mu")))
{
  if (_ndisp != _mesh.dimension())
    paramError("displacements",
               "The number of displacement variables must match the mesh dimension.");
}

void
ComputeLagrangianADNeoHookeanStress::initialSetup()
{
  ComputeLagrangianStressPK1::initialSetup();
  if (!_large_kinematics)
    paramError("large_kinematics", "This material requires large kinematics to be enabled.");
}

void
ComputeLagrangianADNeoHookeanStress::computeQpProperties()
{
  ComputeLagrangianStressPK1::computeQpProperties();

  const bool need_jacobian = _fe_problem.currentlyComputingJacobian() ||
                             _fe_problem.currentlyComputingResidualAndJacobian();
  if (!need_jacobian)
    return;

  if (isPropertyActive(_d_nl_fbar.id()))
    _d_nl_fbar[_qp] =
        (_pk1_jacobian_bypass_fbar[_qp] * _d_F_stab_d_F_avg[_qp])
            .singleProductJ(_F_ust[_qp]) /
        _F_ust[_qp].det();
}

void
ComputeLagrangianADNeoHookeanStress::computeQpPK1Stress()
{
  // Seed the deformation measure passed to the constitutive model. With F-bar stabilization the
  // non-AD kernel supplies the complete derivative of F-bar with respect to the nodal variables.
  std::array<LocalHessianADReal, RankTwoTensor::N2> F;
  for (const auto i : make_range(RankTwoTensor::N))
    for (const auto j : make_range(RankTwoTensor::N))
    {
      const auto component = i * RankTwoTensor::N + j;
      LocalHessianADReal Fij(0.0);
      if (_stabilize_strain)
        Fij.value().value() = _F[_qp](i, j);
      else
      {
        Fij.value().value() = i == j ? 1.0 : 0.0;
        if (i < _ndisp)
          Fij.value().value() += (*_grad_disp[i])[_qp](j);
        if (_coord_system == Moose::COORD_RZ && i == 2 && j == 2)
          Fij.value().value() += (*_disp[0])[_qp] / _q_point[_qp](0);
      }
      Fij.value().derivatives()[component] = 1.0;
      Fij.derivatives()[component].value() = 1.0;
      F[component] = Fij;
    }

  std::array<LocalHessianADReal, RankTwoTensor::N2> C;
  for (const auto i : make_range(RankTwoTensor::N))
    for (const auto j : make_range(RankTwoTensor::N))
    {
      const auto component = i * RankTwoTensor::N + j;
      C[component] = 0.0;
      for (const auto k : make_range(RankTwoTensor::N))
        C[component] += F[k * RankTwoTensor::N + i] * F[k * RankTwoTensor::N + j];
    }

  const auto J = F[0] * (F[4] * F[8] - F[5] * F[7]) - F[1] * (F[3] * F[8] - F[5] * F[6]) +
                 F[2] * (F[3] * F[7] - F[4] * F[6]);
  const auto log_J = log(J);
  const auto trace_C = C[0] + C[4] + C[8];
  const auto energy =
      0.5 * _lambda[_qp] * log_J * log_J - _mu[_qp] * log_J + 0.5 * _mu[_qp] * (trace_C - 3.0);

  for (const auto i : make_range(RankTwoTensor::N))
    for (const auto j : make_range(RankTwoTensor::N))
    {
      const auto stress_component = i * RankTwoTensor::N + j;
      _pk1_stress[_qp](i, j) = energy.derivatives()[stress_component].value();
      for (const auto k : make_range(RankTwoTensor::N))
        for (const auto l : make_range(RankTwoTensor::N))
          _pk1_jacobian[_qp](i, j, k, l) =
              energy.derivatives()[stress_component].derivatives()[k * RankTwoTensor::N + l];
    }
}
