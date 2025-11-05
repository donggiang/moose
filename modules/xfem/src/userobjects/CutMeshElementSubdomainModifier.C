//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CutMeshElementSubdomainModifier.h"

#include "FEProblemBase.h"
#include "MooseTypes.h"
#include "XFEM.h"

registerMooseObject("XFEMApp", CutMeshElementSubdomainModifier);

InputParameters
CutMeshElementSubdomainModifier::validParams()
{
  InputParameters params = ElementSubdomainModifier::validParams();
  params.addClassDescription(
      "Change element subdomain when the element is cut by the XFEM cutter mesh or located at a "
      "crack tip.");
  params.addRequiredParam<SubdomainID>(
      "cut_subdomain_id", "Subdomain ID assigned to elements cut by the cutter mesh.");
  return params;
}

CutMeshElementSubdomainModifier::CutMeshElementSubdomainModifier(const InputParameters & parameters)
  : ElementSubdomainModifier(parameters),
    _cut_subdomain_id(getParam<SubdomainID>("cut_subdomain_id"))
{
  auto * fe_problem = dynamic_cast<FEProblemBase *>(&_subproblem);
  if (!fe_problem)
    mooseError(
        "Problem casting _subproblem to FEProblemBase in CutMeshElementSubdomainModifier");

  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (!_xfem)
    mooseError(name(), " should be used together with XFEM.");
}

SubdomainID
CutMeshElementSubdomainModifier::computeSubdomainID()
{
  const auto elem_uid = _current_elem->unique_id();

  SubdomainID original_subdomain;
  {
    libMesh::Threads::spin_mutex::scoped_lock lock(_original_subdomain_mutex);
    auto insertion =
        _original_subdomain_ids.try_emplace(elem_uid, _current_elem->subdomain_id());
    original_subdomain = insertion.first->second;
  }

  if (_xfem->isElemCut(_current_elem) || _xfem->isElemAtCrackTip(_current_elem))
    return _cut_subdomain_id;

  return original_subdomain;
}