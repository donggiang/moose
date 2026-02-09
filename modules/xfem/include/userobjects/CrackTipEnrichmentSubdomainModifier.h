//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementSubdomainModifier.h"
#include "ElementSubdomainModifierBase.h"
#include <unordered_map>
#include <unordered_set>

class CrackFrontDefinition;
class XFEM;
/**
 * Assigns a subdomain to elements surrounding the crack tip. Elements with
 * nodes that fall within a specified radius of the crack-front points are
 * assigned to the target subdomain, while elements that move outside of the
 * radius are restored to their original blocks.
 */
class CrackTipEnrichmentSubdomainModifier : public ElementSubdomainModifier
{
public:
  static InputParameters validParams();

  CrackTipEnrichmentSubdomainModifier(const InputParameters & parameters);

  void initialize() override;
  //void finalize() override;
protected:
  virtual SubdomainID computeSubdomainID() override;

private:
  /// Target subdomain for the enriched region
  const SubdomainID _cut_off_subdomain_id;

  /// Radius used to assign elements to the enriched subdomain
  const Real _cut_off_radius;

  /// Crack front definition providing the crack-tip locations
  const CrackFrontDefinition & _crack_front_definition;

  /// Track the original subdomain id for elements moved to the enriched subdomain
  std::unordered_map<dof_id_type, SubdomainID> _original_subdomains;

  /// Cached set of element ids that belong to the enriched subdomain
  std::unordered_set<dof_id_type> _enriched_elements;
  /// Shared pointer to XFEM used to query cut elements
  std::shared_ptr<XFEM> _xfem;
};