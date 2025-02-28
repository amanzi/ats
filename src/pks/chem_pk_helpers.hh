/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

//! A set of helper functions for doing common things in PKs.
#pragma once

#include "Teuchos_TimeMonitor.hpp"
#include "Epetra_MultiVector.hpp"

namespace Amanzi {

namespace AmanziChemistry {
class Chemistry_PK;
}

// -----------------------------------------------------------------------------
// Helper functions for working with Amanzi's Chemistry PK
// -----------------------------------------------------------------------------
void
convertConcentrationToAmanzi(const Epetra_MultiVector& mol_den,
                             int num_aqueous,
                             const Epetra_MultiVector& tcc_ats,
                             Epetra_MultiVector& tcc_amanzi);

void
convertConcentrationToATS(const Epetra_MultiVector& mol_den,
                          int num_aqueous,
                          const Epetra_MultiVector& tcc_ats,
                          Epetra_MultiVector& tcc_amanzi);

bool
advanceChemistry(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
                 double t_old,
                 double t_new,
                 bool reinit,
                 const Epetra_MultiVector& mol_dens,
                 Teuchos::RCP<Epetra_MultiVector> tcc,
                 Teuchos::Time& timer);

} // namespace Amanzi
