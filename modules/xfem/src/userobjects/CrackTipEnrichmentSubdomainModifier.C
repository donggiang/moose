//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CrackTipEnrichmentSubdomainModifier.h"

#include "CrackFrontDefinition.h"
#include "MooseTypes.h"

#include "libmesh/elem.h"
#include "libmesh/elem_range.h"
#include "libmesh/mesh.h"
#include "libmesh/node.h"
#include "libmesh/point.h"
#include "XFEM.h"
#include <unordered_set>
#include <vector>

registerMooseObject("XFEMApp", CrackTipEnrichmentSubdomainModifier);

InputParameters
CrackTipEnrichmentSubdomainModifier::validParams()
{
  InputParameters params = ElementSubdomainModifier::validParams();
  params.addClassDescription("Assign a subdomain to elements surrounding the crack tip using a "
                             "radius-based search.");

  params.addRequiredParam<UserObjectName>(
      "crack_front_definition",
      "Name of the CrackFrontDefinition user object that provides the crack front points.");
  params.addRequiredParam<Real>("cut_off_subdomain_id",
                                "Subdomain ID contains crack tip locations");

  params.addRequiredParam<Real>(
      "cut_off_radius",
      "Radius around each crack front point used to assign elements to the enriched subdomain.");

  ExecFlagEnum & exec = params.set<ExecFlagEnum>("execute_on");
  exec.addAvailableFlags(EXEC_XFEM_SUBDOMAIN_MODIFIER);
  params.set<ExecFlagEnum>("execute_on") = {EXEC_XFEM_SUBDOMAIN_MODIFIER};
  params.set<int>("execution_order_group") = 0;
  return params;
}

CrackTipEnrichmentSubdomainModifier::CrackTipEnrichmentSubdomainModifier(
    const InputParameters & parameters)
  : ElementSubdomainModifier(parameters),
    _cut_off_subdomain_id(getParam<Real>("cut_off_subdomain_id")),
    _cut_off_radius(getParam<Real>("cut_off_radius")),
    _crack_front_definition(getUserObject<CrackFrontDefinition>("crack_front_definition"))
{
  auto * fe_problem = dynamic_cast<FEProblemBase *>(&_subproblem);
  if (!fe_problem)
    mooseError("Problem casting _subproblem to FEProblemBase in CutMeshElementSubdomainModifier");

  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (!_xfem)
    mooseError(name(), " should be used together with XFEM.");
}

void
CrackTipEnrichmentSubdomainModifier::initialize()
{
  // Let the base class do its bookkeeping
  ElementSubdomainModifier::initialize();
  // Make sure we start fresh for this XFEM update
  _original_subdomains.clear();
}
SubdomainID
CrackTipEnrichmentSubdomainModifier::computeSubdomainID()
{
  if (_cut_off_radius <= 0.0)
    return Moose::INVALID_BLOCK_ID;

  const std::size_t num_crack_points = _crack_front_definition.getNumCrackFrontPoints();
  if (num_crack_points == 0)
    return Moose::INVALID_BLOCK_ID;

  // collect crack-front points (often just 1 in 2D)
  std::vector<Point> crack_points;
  crack_points.reserve(num_crack_points);
  for (std::size_t i = 0; i < num_crack_points; ++i)
  {
    const Point * point = _crack_front_definition.getCrackFrontPoint(i);
    if (point)
      crack_points.push_back(*point);
  }

  if (crack_points.empty())
    return Moose::INVALID_BLOCK_ID;

  const Real radius_sq = _cut_off_radius * _cut_off_radius;
  const SubdomainID current_id = _current_elem->subdomain_id();
  const SubdomainID enriched_id = _cut_off_subdomain_id; // e.g. 2
  const SubdomainID base_id = 1;                         // grey region

  // --- geometry check: is this element within the tip radius? ---
  bool is_within_radius = false;

  for (unsigned int n = 0; n < _current_elem->n_nodes() && !is_within_radius; ++n)
  {
    const auto * node = _current_elem->node_ptr(n);
    if (!node)
      continue;

    const Point & node_point = *node;

    for (const auto & crack_point : crack_points)
    {
      const RealVectorValue diff = node_point - crack_point;
      if (diff.norm_sq() <= radius_sq)
      {
        is_within_radius = true;
        break;
      }
    }
  }

  // --- assign subdomain based on geometry ---

  if (is_within_radius)
  {
    // we want this element in the enriched block
    if (current_id != enriched_id)
      return enriched_id; // change to red
    else
      return Moose::INVALID_BLOCK_ID; // already red, no change
  }
  else
  {
    // we want this element in the base block
    if (current_id != base_id)
      return base_id; // change back to grey
    else
      return Moose::INVALID_BLOCK_ID; // already grey, no change
  }
}
