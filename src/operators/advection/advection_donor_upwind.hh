/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   ATS

   Donor upwind advection.
   ------------------------------------------------------------------------- */

#ifndef OPERATOR_ADVECTION_ADVECTION_DONOR_UPWIND_HH_
#define OPERATOR_ADVECTION_ADVECTION_DONOR_UPWIND_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_IntVector.h"

#include "Mesh.hh"
#include "CompositeVector.hh"

#include "advection.hh"

namespace Amanzi {
namespace Operators {

class AdvectionDonorUpwind : public Advection {
 public:
  AdvectionDonorUpwind(Teuchos::ParameterList& advect_plist,
                       const Teuchos::RCP<const AmanziMesh::Mesh> mesh);

  virtual void set_flux(const Teuchos::RCP<const CompositeVector>& flux);
  virtual void Apply(const Teuchos::RCP<Functions::BoundaryFunction>& bc_flux,
                     bool include_bc_fluxes = true);

 private:
  void IdentifyUpwindCells_();

  Teuchos::RCP<Epetra_IntVector> upwind_cell_;
  Teuchos::RCP<Epetra_IntVector> downwind_cell_;
};

} // namespace Operators
} // namespace Amanzi

#endif
