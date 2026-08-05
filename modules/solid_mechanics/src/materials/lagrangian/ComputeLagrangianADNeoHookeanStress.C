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
  InputParameters params = ComputeLagrangianStressPK2::validParams();
  params.addClassDescription(
      "Compute a compressible Neo-Hookean stress and material tangent by locally differentiating "
      "the strain-energy potential.");
  params.addParam<MaterialPropertyName>(
      "lambda", "lambda", "First Lame parameter material property.");
  params.addParam<MaterialPropertyName>("mu", "mu", "Shear modulus material property.");
  return params;
}

ComputeLagrangianADNeoHookeanStress::ComputeLagrangianADNeoHookeanStress(
    const InputParameters & parameters)
  : ComputeLagrangianStressPK2(parameters),
    _lambda(getMaterialProperty<Real>(getParam<MaterialPropertyName>("lambda"))),
    _mu(getMaterialProperty<Real>(getParam<MaterialPropertyName>("mu")))
{
  if (!_large_kinematics)
    paramError("large_kinematics", "This material requires large kinematics to be enabled.");
}

void
ComputeLagrangianADNeoHookeanStress::computeQpPK2Stress()
{
  // Seed both AD levels so one energy evaluation provides its gradient and Hessian with respect
  // to the nine Green-Lagrange strain components.
  std::array<LocalHessianADReal, RankTwoTensor::N2> E;
  for (const auto i : make_range(RankTwoTensor::N))
    for (const auto j : make_range(RankTwoTensor::N))
    {
      const auto component = i * RankTwoTensor::N + j;
      LocalHessianADReal Eij(0.0);
      Eij.value().value() = _E[_qp](i, j);
      Eij.value().derivatives()[component] = 1.0;
      Eij.derivatives()[component].value() = 1.0;
      E[component] = Eij;
    }

  std::array<LocalHessianADReal, RankTwoTensor::N2> C;
  for (const auto i : make_range(RankTwoTensor::N))
    for (const auto j : make_range(RankTwoTensor::N))
    {
      const auto component = i * RankTwoTensor::N + j;
      C[component] = 2.0 * E[component] + (i == j ? 1.0 : 0.0);
    }

  const auto det_C = C[0] * (C[4] * C[8] - C[5] * C[7]) -
                     C[1] * (C[3] * C[8] - C[5] * C[6]) +
                     C[2] * (C[3] * C[7] - C[4] * C[6]);
  const auto log_J = 0.5 * log(det_C);
  const auto trace_C = C[0] + C[4] + C[8];
  const auto energy = 0.5 * _lambda[_qp] * log_J * log_J - _mu[_qp] * log_J +
                      0.5 * _mu[_qp] * (trace_C - 3.0);

  for (const auto i : make_range(RankTwoTensor::N))
    for (const auto j : make_range(RankTwoTensor::N))
    {
      const auto stress_component = i * RankTwoTensor::N + j;
      _S[_qp](i, j) = energy.derivatives()[stress_component].value();
      for (const auto k : make_range(RankTwoTensor::N))
        for (const auto l : make_range(RankTwoTensor::N))
          _C[_qp](i, j, k, l) =
              energy.derivatives()[stress_component].derivatives()[k * RankTwoTensor::N + l];
    }
}
