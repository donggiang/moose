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
	  /// Target subdomain for cut or crack-tip elements
	  SubdomainID _cut_subdomain_id;
	  /// XFEM controller used to query cut and crack-tip status
	  std::shared_ptr<XFEM> _xfem;
	  // /// Original subdomain id for each element (keyed by unique id)
	  // std::unordered_map<libMesh::unique_id_type, SubdomainID> _original_subdomain_ids;
	  // /// Mutex to protect access to the original subdomain map in threaded execution
	  // mutable libMesh::Threads::spin_mutex _original_subdomain_mutex;
};