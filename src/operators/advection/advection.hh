/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   ATS

   Base interface for a general-purpose advection operator.
   ------------------------------------------------------------------------- */

#ifndef OPERATOR_ADVECTION_ADVECTION_HH_
#define OPERATOR_ADVECTION_ADVECTION_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Mesh.hh"
#include "CompositeVector.hh"
#include "BoundaryFunction.hh"

namespace Amanzi {
namespace Operators {

class Advection {
 public:
  Advection(Teuchos::ParameterList& advect_plist, const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : advect_plist_(advect_plist), mesh_(mesh)
  {}
  virtual ~Advection() = default;

  Teuchos::RCP<const CompositeVector> flux() const { return flux_; }
  virtual void set_flux(const Teuchos::RCP<const CompositeVector>& flux);

  unsigned int num_dofs() const { return num_dofs_; }
  virtual void set_num_dofs(unsigned int num_dofs);

  Teuchos::RCP<CompositeVector> field() { return field_; }

  virtual void Apply(const Teuchos::RCP<Functions::BoundaryFunction>& bc_flux,
                     bool include_bc_fluxes = true) = 0;

 protected:
  unsigned int num_dofs_;
  Teuchos::RCP<const CompositeVector> flux_;
  Teuchos::RCP<CompositeVector> field_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::ParameterList advect_plist_;
};

} // namespace Operators
} // namespace Amanzi

#endif
