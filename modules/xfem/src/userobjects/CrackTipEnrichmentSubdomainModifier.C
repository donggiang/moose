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
  // params.addRequiredParam<SubdomainName>(
  //     "cut_off_subdomain",
  //     "Subdomain to assign to elements whose nodes lie within the crack-tip radius.");
  params.addRequiredParam<UserObjectName>(
      "crack_front_definition",
      "Name of the CrackFrontDefinition user object that provides the crack front points.");
  params.addRequiredParam<Real>(
      "cut_off_subdomain_id",
      "Radius around each crack front point used to assign elements to the enriched subdomain.");

  params.addRequiredParam<Real>(
      "cut_off_radius",
      "Radius around each crack front point used to assign elements to the enriched subdomain.");

  // ExecFlagEnum & exec = params.set<ExecFlagEnum>("execute_on");
  // exec.addAvailableFlags(EXEC_XFEM_MARK);
  // params.set<ExecFlagEnum>("execute_on") = {EXEC_XFEM_MARK};
  // params.set<int>("execution_order_group") = 100;
  ExecFlagEnum & exec = params.set<ExecFlagEnum>("execute_on");
  exec.addAvailableFlags(EXEC_XFEM_SUBDOMAIN_MODIFIER);
  params.set<ExecFlagEnum>("execute_on") = { EXEC_XFEM_SUBDOMAIN_MODIFIER};
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

  std::cout<<"********CrackTipEnrichmentSubdomainModifier::CrackTipEnrichmentSubdomainModifier : call constructure\n";
    auto * fe_problem = dynamic_cast<FEProblemBase *>(&_subproblem);
  if (!fe_problem)
    mooseError("Problem casting _subproblem to FEProblemBase in CutMeshElementSubdomainModifier");

  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (!_xfem)
    mooseError(name(), " should be used together with XFEM.");
}



// void
// CrackTipEnrichmentSubdomainModifier::initialize()
// {
//   std::cout<<"********* CrackTipEnrichmentSubdomainModifier::initialize()\n";
//   ElementSubdomainModifier::initialize();
//   _enriched_elements.clear();

//   if (_cut_off_radius <= 0.0)
//     return;

//   const std::size_t num_crack_points = _crack_front_definition.getNumCrackFrontPoints();
//   if (num_crack_points == 0)
//     return;

//   std::vector<Point> crack_points;
//   crack_points.reserve(num_crack_points);
//   for (std::size_t i = 0; i < num_crack_points; ++i)
//   {
//     const Point * point = _crack_front_definition.getCrackFrontPoint(i);
//     if (point)
//       crack_points.push_back(*point);
//   }

//   if (crack_points.empty())
//     return;
//   std::cout<< "********* CrackTipEnrichmentSubdomainModifier::initialize(), crack tip location: " << crack_points[0] << "\n";
//   const Real radius_sq = _cut_off_radius * _cut_off_radius;

//   ConstElemRange & elem_range = *_mesh.getActiveLocalElementRange();
//   for (const auto & elem : elem_range)
//   {
//     if (elem->processor_id() != processor_id())
//       continue;

//     bool is_within_radius = false;
//     for (unsigned int n = 0; n < elem->n_nodes(); ++n)
//     {
//       const auto * node = elem->node_ptr(n);
//       if (!node)
//         continue;

//       const Point & node_point = *node;
//       for (const auto & crack_point : crack_points)
//       {
//         const RealVectorValue diff = node_point - crack_point;
//         if (diff.norm_sq() <= radius_sq)
//         {
//           is_within_radius = true;
//           break;
//         }
//       }

//       if (is_within_radius)
//         break;
//     }

//     if (is_within_radius)
//       _enriched_elements.insert(elem->id());
//   }

//   _xfem->printElementConnectivity();
// }

void
CrackTipEnrichmentSubdomainModifier::initialize()
{
  // Let the base class do its bookkeeping
  ElementSubdomainModifier::initialize();
  // // Let the base class do its bookkeeping
  // ElementSubdomainModifier::initialize();

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

     // std::cout<<"********CrackTipEnrichmentSubdomainModifier::CrackTipEnrichmentSubdomainModifier : call computeSubdomainID\n";
  const Real        radius_sq   = _cut_off_radius * _cut_off_radius;
  const SubdomainID current_id  = _current_elem->subdomain_id();
  const SubdomainID enriched_id = _cut_off_subdomain_id;   // e.g. 2
  const SubdomainID base_id     = 1;                       // grey region

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

        //         // std::cout << "  node " << n << " within radius: "
        // //           << "r^2=" << radius_sq << ", dist^2=" << dist_sq << "\n";
        //     std::cout<< "****** CrackTipEnrichmentSubdomainModifier::computeSubdomainID, detect new subdomain, crack tip location: " << crack_points[0] << "\n";
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
      return enriched_id;             // change to red
    else
      return Moose::INVALID_BLOCK_ID; // already red, no change
  }
  else
  {
    // we want this element in the base block
    if (current_id != base_id)
      return base_id;                 // change back to grey
    else
      return Moose::INVALID_BLOCK_ID; // already grey, no change
  }
}


// SubdomainID
// CrackTipEnrichmentSubdomainModifier::computeSubdomainID()
// {
//   if (_cut_off_radius <= 0.0)
//     return Moose::INVALID_BLOCK_ID;

//   const std::size_t num_crack_points = _crack_front_definition.getNumCrackFrontPoints();
//   if (num_crack_points == 0)
//     return Moose::INVALID_BLOCK_ID;

//   std::vector<Point> crack_points;
//   crack_points.reserve(num_crack_points);
//   for (std::size_t i = 0; i < num_crack_points; ++i)
//   {
//     const Point * point = _crack_front_definition.getCrackFrontPoint(i);
//     if (point)
//       crack_points.push_back(*point);
//   }

//   if (crack_points.empty())
//     return Moose::INVALID_BLOCK_ID;

//   const Real radius_sq = _cut_off_radius * _cut_off_radius;
//   const SubdomainID current_id = _current_elem->subdomain_id();
//   const dof_id_type elem_id = _current_elem->id();

//   bool is_within_radius = false;

//   for (unsigned int n = 0; n < _current_elem->n_nodes() && !is_within_radius; ++n)
//   {
//     const auto * node = _current_elem->node_ptr(n);
//     if (!node)
//       continue;

//     const Point & node_point = *node;

//     for (const auto & crack_point : crack_points)
//     {
//       const RealVectorValue diff = node_point - crack_point;
//       if (diff.norm_sq() <= radius_sq)
//       {
//         is_within_radius = true;
//         break;
//       }
//     }
//   }

//   if (is_within_radius)
//   {
//     // element should be in enriched block (_cut_off_subdomain_id)
//     // remember original id the first time we touch it
//     auto it = _original_subdomains.find(elem_id);
//     if (it == _original_subdomains.end())
//       _original_subdomains[elem_id] = current_id;

//     if (current_id != _cut_off_subdomain_id)
//       return _cut_off_subdomain_id;          // change block
//     else
//       return Moose::INVALID_BLOCK_ID;        // no change
//   }
//   else
//   {
//     // element should NOT be enriched:
//     // if we previously overrode it, restore the original subdomain
//     auto it = _original_subdomains.find(elem_id);
//     if (it != _original_subdomains.end())
//     {
//       SubdomainID original_id = it->second;
//       _original_subdomains.erase(it);

//       if (current_id != original_id)
//         return original_id;                  // restore original block
//     }

//     // otherwise do nothing
//     return Moose::INVALID_BLOCK_ID;
//   }
// }


// SubdomainID
// CrackTipEnrichmentSubdomainModifier::computeSubdomainID()
// {
//   if (_cut_off_radius <= 0.0)
//     return Moose::INVALID_BLOCK_ID;

//   const std::size_t num_crack_points = _crack_front_definition.getNumCrackFrontPoints();
//   if (num_crack_points == 0)
//     return Moose::INVALID_BLOCK_ID;

//   std::vector<Point> crack_points;
//   crack_points.reserve(num_crack_points);
//   for (std::size_t i = 0; i < num_crack_points; ++i)
//   {
//     const Point * point = _crack_front_definition.getCrackFrontPoint(i);
//     if (point)
//       crack_points.push_back(*point);
//   }

//   if (crack_points.empty())
//     return Moose::INVALID_BLOCK_ID;

//   const Real radius_sq = _cut_off_radius * _cut_off_radius;
//   const SubdomainID current_id = _current_elem->subdomain_id();
//   const dof_id_type elem_id = _current_elem->id();

//   // std::cout << "****************CrackTipEnrichmentSubdomainModifier::computeSubdomainID, element id: "
//   //           << elem_id << "*************\n";

//   bool is_within_radius = false;

//   // ---- distance check based on nodes ----
//   for (unsigned int n = 0; n < _current_elem->n_nodes() && !is_within_radius; ++n)
//   {
//     const auto * node = _current_elem->node_ptr(n);
//     if (!node)
//       continue;

//     const Point & node_point = *node;

//     for (const auto & crack_point : crack_points)
//     {
//       const RealVectorValue diff = node_point - crack_point;
//       const Real dist_sq = diff.norm_sq();

//       if (dist_sq <= radius_sq)
//       {
//         // std::cout << "  node " << n << " within radius: "
//         //           << "r^2=" << radius_sq << ", dist^2=" << dist_sq << "\n";
//             std::cout<< "****** CrackTipEnrichmentSubdomainModifier::computeSubdomainID, detect new subdomain, crack tip location: " << crack_points[0] << "\n";
//         is_within_radius = true;
//         break;
//       }
//     }
//   }

//   // ---- enforce subdomain 1/2 based purely on geometry ----

//   if (is_within_radius)
//   {
//     // // element should be enriched → block 2
//     // if (current_id != _cut_off_subdomain_id)
//     //   std::cout << "  set elem " << elem_id << " to subdomain " << _cut_off_subdomain_id << "\n";

//     // change only if needed
//     return current_id == _cut_off_subdomain_id ? Moose::INVALID_BLOCK_ID : _cut_off_subdomain_id;
//   }
//   else
//   {
//     // element should be standard → block 1
//     const SubdomainID standard_id = 1;

//     // if (current_id != standard_id)
//     //   std::cout << "  restore elem " << elem_id << " to subdomain " << standard_id << "\n";

//     // change only if needed
//     return current_id == standard_id ? Moose::INVALID_BLOCK_ID : standard_id;
//   }
// }
