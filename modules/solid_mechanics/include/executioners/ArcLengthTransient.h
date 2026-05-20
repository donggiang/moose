//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Transient.h"
#include "ArcLengthSolveObject.h"

/**
 * Transient executioner whose inner Newton solve is replaced by the Box 4.4
 * cylindrical arc-length algorithm (see `ArcLengthSolveObject`).
 *
 * Usage:
 *   [Executioner]
 *     type = ArcLengthTransient
 *     scalar_variable = lambda
 *     delta_s = 0.05
 *     arc_length_load_tag = arc_length_load
 *     ...
 *   []
 *
 * Each timestep performs ONE outer arc-length step (Δs traversed).
 * Inside the step, two linear solves per Newton iteration. Standard
 * Transient time-step controls (dt, end_time) gate how many outer
 * steps are taken — typically `dt = 1` and `end_time = N` for N steps.
 */
class ArcLengthTransient : public Transient
{
public:
  static InputParameters validParams();
  ArcLengthTransient(const InputParameters & params);

  virtual void init() override;

protected:
  /// Owned arc-length solve object. Injected as the inner solve of the
  /// Transient's fixed-point solve at init().
  ArcLengthSolveObject _al_solve;
};
