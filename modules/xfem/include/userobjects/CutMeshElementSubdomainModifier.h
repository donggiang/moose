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
	#include "libmesh/libmesh_common.h"
	#include "libmesh/threads.h"
	#include <memory>
	#include <unordered_map>
	class FEProblemBase;
	class XFEM;
	/**
	 * A subdomain modifier that moves elements cut by the XFEM cutter mesh (or located at
	 * crack tips) to a specified subdomain. Elements that are not cut remain on their
	 * original subdomain.
	 */
	class CutMeshElementSubdomainModifier : public ElementSubdomainModifier
	{
	public:
	  static InputParameters validParams();
	  CutMeshElementSubdomainModifier(const InputParameters & parameters);
	protected:
	  virtual SubdomainID computeSubdomainID() override;
	private:

	  std::shared_ptr<XFEM> _xfem;
	};
