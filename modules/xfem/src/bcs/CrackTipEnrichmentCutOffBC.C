//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CrackTipEnrichmentCutOffBC.h"

registerMooseObject("XFEMApp", CrackTipEnrichmentCutOffBC);

InputParameters
CrackTipEnrichmentCutOffBC::validParams()
{
  InputParameters p = DirichletBC::validParams();
  p.addRequiredParam<Real>("cut_off_radius",
                           "The cut off radius of crack tip enrichment functions");
  p.set<bool>("use_displaced_mesh") = false;
  p.addRequiredParam<UserObjectName>("crack_front_definition",
                                     "The CrackFrontDefinition user object name");
  p.set<bool>("preset") = true;
  return p;
}

CrackTipEnrichmentCutOffBC::CrackTipEnrichmentCutOffBC(const InputParameters & parameters)
  : DirichletBC(parameters),
    _cut_off_radius(getParam<Real>("cut_off_radius")),
    _crack_front_definition(getUserObject<CrackFrontDefinition>("crack_front_definition"))
{
    std::cout<<"***********CrackTipEnrichmentCutOffBC::CrackTipEnrichmentCutOffBC: call constructor*****\n";
}

bool
CrackTipEnrichmentCutOffBC::shouldApply() const
{
  // if(_cut_off_radius==0)
  // {
  // //     std::cout<<"******************************************\n";
  // // std::cout<<"***********CrackTipEnrichmentCutOffBC: node: " << (*_current_node) << "\n";
  //   return true;
  // }
  // Real r, theta;
  // _crack_front_definition.calculateRThetaToCrackFront((*_current_node), r, theta);

  // if (r > _cut_off_radius && _cut_off_radius>0)
  // {
  // // std::cout<<"******************************************\n";
  // // std::cout<<"***********CrackTipEnrichmentCutOffBC: node: " << (*_current_node) << "\n";
    return true;
  // }
  // else
  //   return false;
}
