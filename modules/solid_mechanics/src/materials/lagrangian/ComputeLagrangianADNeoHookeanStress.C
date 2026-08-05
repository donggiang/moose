//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "ComputeLagrangianADNeoHookeanStress.h"

#include "ADReal.h"

registerMooseObject("SolidMechanicsApp", ComputeLagrangianADNeoHookeanStress);

InputParameters
ComputeLagrangianADNeoHookeanStress::validParams()
{
  InputParameters params = ComputeLagrangianStressPK2::validParams();
  params.addClassDescription(
      "Compute a compressible Neo-Hookean stress using the standard Lagrangian implementation "
      "and its material tangent using local automatic differentiation.");
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
  const auto I = RankTwoTensor::Identity();
  const RankTwoTensor Cinv = (2.0 * _E[_qp] + I).inverse();
  _S[_qp] = (_lambda[_qp] * log(_F[_qp].det()) - _mu[_qp]) * Cinv + _mu[_qp] * I;

  // The Green-Lagrange strain and PK2 stress above remain ordinary MOOSE tensors. Seed only
  // the constitutive input used to calculate d(PK2)/dE with local forward AD.
  ADRankTwoTensor E;
  for (const auto i : make_range(RankTwoTensor::N))
    for (const auto j : make_range(RankTwoTensor::N))
    {
      const auto component = i * RankTwoTensor::N + j;
      ADReal Eij{};
      Eij.value() = _E[_qp](i, j);
      Moose::derivInsert(Eij.derivatives(), component, 1.0);
      E(i, j) = Eij;
    }

  const auto C = 2.0 * E + ADRankTwoTensor::Identity();
  const auto Cinv_ad = C.inverse();
  const auto log_J = 0.5 * log(C.det());
  const ADRankTwoTensor S =
      (_lambda[_qp] * log_J - _mu[_qp]) * Cinv_ad +
      _mu[_qp] * ADRankTwoTensor::Identity();

  for (const auto i : make_range(RankTwoTensor::N))
    for (const auto j : make_range(RankTwoTensor::N))
      for (const auto k : make_range(RankTwoTensor::N))
        for (const auto l : make_range(RankTwoTensor::N))
          _C[_qp](i, j, k, l) = S(i, j).derivatives()[k * RankTwoTensor::N + l];
}
