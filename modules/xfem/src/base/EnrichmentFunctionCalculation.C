//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EnrichmentFunctionCalculation.h"

EnrichmentFunctionCalculation::EnrichmentFunctionCalculation(
    const CrackFrontDefinition * crack_front_definition)
  : _crack_front_definition(*crack_front_definition)
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

  // Power-law creep exponent n from Norton law: edot ~ sigma^n
  const Real n = 1.0; // replace with your actual member/property
  const Real lambda = 1.0 / (n + 1.0);

  const Real st  = std::sin(_theta);
  const Real st2 = std::sin(0.5 * _theta);
  const Real ct2 = std::cos(0.5 * _theta);

  const Real rl = std::pow(_r, lambda);

  // Same angular functions as LEFM, but radial order changed from sqrt(r) to r^lambda
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

  // Power-law creep exponent n from Norton law: edot ~ sigma^n
  const Real n = 1.0; // replace with your actual member/property
  const Real lambda = 1.0 / (n + 1.0);

  const Real st  = std::sin(_theta);
  const Real ct  = std::cos(_theta);
  const Real st2 = std::sin(0.5 * _theta);
  const Real ct2 = std::cos(0.5 * _theta);

  // Angular functions
  const Real F0 = st2;
  const Real F1 = ct2;
  const Real F2 = st2 * st;
  const Real F3 = ct2 * st;

  // Angular derivatives dF/dtheta
  const Real dF0 =  0.5 * ct2;
  const Real dF1 = -0.5 * st2;
  const Real dF2 =  0.5 * ct2 * st + st2 * ct;
  const Real dF3 = -0.5 * st2 * st + ct2 * ct;

  // Common radial factor: d/dr of r^lambda gives lambda r^(lambda-1)
  const Real rlm1 = std::pow(_r, lambda - 1.0);

  // General formula:
  // dB/dx = r^(lambda-1) [ lambda F cos(theta) - F'(theta) sin(theta) ]
  // dB/dy = r^(lambda-1) [ lambda F sin(theta) + F'(theta) cos(theta) ]

  // B0 = r^lambda * sin(theta/2)
  dB[0](0)= rlm1 * (lambda * F0 * ct - dF0 * st);
  dB[0](1) = rlm1 * (lambda * F0 * st + dF0 * ct);
  dB[0](2) = 0.0;

  // B1 = r^lambda * cos(theta/2)
  dB[1](0)= rlm1 * (lambda * F1 * ct - dF1 * st);
  dB[1](1) = rlm1 * (lambda * F1 * st + dF1 * ct);
  dB[1](2)= 0.0;

  // B2 = r^lambda * sin(theta/2) * sin(theta)
  dB[2](0)= rlm1 * (lambda * F2 * ct - dF2 * st);
  dB[2](1) = rlm1 * (lambda * F2 * st + dF2 * ct);
  dB[2](2) = 0.0;

  // B3 = r^lambda * cos(theta/2) * sin(theta)
  dB[3](0) = rlm1 * (lambda * F3 * ct - dF3 * st);
  dB[3](1)= rlm1 * (lambda * F3 * st + dF3 * ct);
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
