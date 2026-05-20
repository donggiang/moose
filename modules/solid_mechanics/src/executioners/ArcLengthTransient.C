//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ArcLengthTransient.h"
#include "FixedPointSolve.h"

registerMooseObject("SolidMechanicsApp", ArcLengthTransient);

InputParameters
ArcLengthTransient::validParams()
{
  InputParameters params = Transient::validParams();
  params += ArcLengthSolveObject::validParams();
  params.addClassDescription(
      "Transient executioner whose inner solve is the cylindrical arc-length "
      "method per Box 4.4 of de Souza Neto, Perić, Owen (2008). Each timestep "
      "performs ONE outer arc-length step.");
  return params;
}

ArcLengthTransient::ArcLengthTransient(const InputParameters & params)
  : Transient(params), _al_solve(*this)
{
  // Inject the AL solve as the inner solve of the FixedPointSolve, replacing
  // the default FEProblemSolve that Transient set up.
  _fixed_point_solve->setInnerSolve(_al_solve);
}

void
ArcLengthTransient::init()
{
  Transient::init();
  _al_solve.initialSetup();
}
