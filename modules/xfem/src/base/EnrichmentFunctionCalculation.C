//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EnrichmentFunctionCalculation.h"

#include <cmath>

void
EnrichmentFunctionCalculation::addValidParams(InputParameters & params)
{
  params.addRangeCheckedParam<Real>(
      "creep_exponent",
      1.0,
      "creep_exponent > 0.0",
      "Power-law creep exponent n used in the crack-tip enrichment radial order r^(1/(n+1)). "
      "The default n=1 recovers the elastic sqrt(r) enrichment.");
}

EnrichmentFunctionCalculation::EnrichmentFunctionCalculation(
    const CrackFrontDefinition * crack_front_definition, const Real creep_exponent)
  : _crack_front_definition(*crack_front_definition),
    _creep_exponent(creep_exponent),
    _lambda(1.0 / (_creep_exponent + 1.0))
{
}

unsigned int
EnrichmentFunctionCalculation::crackTipEnrichementFunctionAtPoint(const Point & point,
                                                                  std::vector<Real> & B)
{
  unsigned int crack_front_point_index =
      _crack_front_definition.calculateRThetaToCrackFront(point, _r, _theta);

  if (MooseUtils::absoluteFuzzyEqual(_r, 0.0))
    mooseError("EnrichmentFunctionCalculation: the distance between a point and the crack "
               "tip/front is zero.");

  if (MooseUtils::absoluteFuzzyEqual(_creep_exponent, 1.0))
  {
    Real st = std::sin(_theta);
    Real st2 = std::sin(_theta / 2.0);
    Real ct2 = std::cos(_theta / 2.0);
    Real sr = std::sqrt(_r);

    B[0] = sr * st2;
    B[1] = sr * ct2;
    B[2] = sr * st2 * st;
    B[3] = sr * ct2 * st;

    return crack_front_point_index;
  }

  const Real st = std::sin(_theta);
  const Real st2 = std::sin(_theta / 2.0);
  const Real ct2 = std::cos(_theta / 2.0);
  const Real rl = std::pow(_r, _lambda);

  B[0] = rl * st2;
  B[1] = rl * ct2;
  B[2] = rl * st2 * st;
  B[3] = rl * ct2 * st;

  return crack_front_point_index;
}

unsigned int
EnrichmentFunctionCalculation::crackTipEnrichementFunctionDerivativeAtPoint(
    const Point & point, std::vector<RealVectorValue> & dB)
{
  unsigned int crack_front_point_index =
      _crack_front_definition.calculateRThetaToCrackFront(point, _r, _theta);

  if (MooseUtils::absoluteFuzzyEqual(_r, 0.0))
    mooseError("EnrichmentFunctionCalculation: the distance between a point and the crack "
               "tip/front is zero.");

  if (MooseUtils::absoluteFuzzyEqual(_creep_exponent, 1.0))
  {
    Real st = std::sin(_theta);
    Real ct = std::cos(_theta);
    Real st2 = std::sin(_theta / 2.0);
    Real ct2 = std::cos(_theta / 2.0);
    Real st15 = std::sin(1.5 * _theta);
    Real ct15 = std::cos(1.5 * _theta);
    Real sr = std::sqrt(_r);

    dB[0](0) = -0.5 / sr * st2;
    dB[0](1) = 0.5 / sr * ct2;
    dB[0](2) = 0.0;
    dB[1](0) = 0.5 / sr * ct2;
    dB[1](1) = 0.5 / sr * st2;
    dB[1](2) = 0.0;
    dB[2](0) = -0.5 / sr * st15 * st;
    dB[2](1) = 0.5 / sr * (st2 + st15 * ct);
    dB[2](2) = 0.0;
    dB[3](0) = -0.5 / sr * ct15 * st;
    dB[3](1) = 0.5 / sr * (ct2 + ct15 * ct);
    dB[3](2) = 0.0;

    return crack_front_point_index;
  }

  const Real st = std::sin(_theta);
  const Real ct = std::cos(_theta);
  const Real st2 = std::sin(_theta / 2.0);
  const Real ct2 = std::cos(_theta / 2.0);
  const Real rlm1 = std::pow(_r, _lambda - 1.0);

  const Real F0 = st2;
  const Real F1 = ct2;
  const Real F2 = st2 * st;
  const Real F3 = ct2 * st;

  const Real dF0 = 0.5 * ct2;
  const Real dF1 = -0.5 * st2;
  const Real dF2 = 0.5 * ct2 * st + st2 * ct;
  const Real dF3 = -0.5 * st2 * st + ct2 * ct;

  dB[0](0) = rlm1 * (_lambda * F0 * ct - dF0 * st);
  dB[0](1) = rlm1 * (_lambda * F0 * st + dF0 * ct);
  dB[0](2) = 0.0;
  dB[1](0) = rlm1 * (_lambda * F1 * ct - dF1 * st);
  dB[1](1) = rlm1 * (_lambda * F1 * st + dF1 * ct);
  dB[1](2) = 0.0;
  dB[2](0) = rlm1 * (_lambda * F2 * ct - dF2 * st);
  dB[2](1) = rlm1 * (_lambda * F2 * st + dF2 * ct);
  dB[2](2) = 0.0;
  dB[3](0) = rlm1 * (_lambda * F3 * ct - dF3 * st);
  dB[3](1) = rlm1 * (_lambda * F3 * st + dF3 * ct);
  dB[3](2) = 0.0;

  return crack_front_point_index;
}

void
EnrichmentFunctionCalculation::rotateFromCrackFrontCoordsToGlobal(const RealVectorValue & vector,
                                                                  RealVectorValue & rotated_vector,
                                                                  const unsigned int point_index)
{
  rotated_vector = _crack_front_definition.rotateFromCrackFrontCoordsToGlobal(vector, point_index);
}
