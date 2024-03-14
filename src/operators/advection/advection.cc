/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   ATS

   Base class for advection.
   ------------------------------------------------------------------------- */

#include "CompositeVectorSpace.hh"

#include "advection.hh"

namespace Amanzi {
namespace Operators {

void
Advection::set_flux(const Teuchos::RCP<const CompositeVector>& flux)
{
  // check that flux includes FACES and has one dof
  flux_ = flux;
}

void
Advection::set_num_dofs(unsigned int num_dofs)
{
  if (field_ == Teuchos::null || num_dofs_ != num_dofs) {
    num_dofs_ = num_dofs;
    std::vector<int> ndofs_tmp(2, num_dofs_);

    std::vector<std::string> names(2);
    names[0] = "cell";
    names[1] = "face";

    std::vector<AmanziMesh::Entity_kind> locations(2);
    locations[0] = AmanziMesh::Entity_kind::CELL;
    locations[1] = AmanziMesh::Entity_kind::FACE;

    Teuchos::RCP<CompositeVectorSpace> space = Teuchos::rcp(new CompositeVectorSpace());
    space->SetMesh(mesh_)->SetGhosted()->SetComponents(names, locations, ndofs_tmp);
    field_ = Teuchos::rcp(new CompositeVector(*space));
  }
}

} // namespace Operators
} // namespace Amanzi
