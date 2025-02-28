/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

#include "chem_pk_helpers.hh"
#include "Chemistry_PK.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Helper functions for working with Amanzi's Chemistry PK
// -----------------------------------------------------------------------------
void
convertConcentrationToAmanzi(const Epetra_MultiVector& mol_dens,
                             int num_aqueous,
                             const Epetra_MultiVector& tcc_ats,
                             Epetra_MultiVector& tcc_amanzi)
{
  // convert from mole fraction [mol C / mol H20] to [mol C / L]
  for (int k = 0; k != num_aqueous; ++k) {
    for (int c = 0; c != tcc_ats.MyLength(); ++c) {
      // 1.e-3 converts L to m^3
      tcc_amanzi[k][c] = tcc_ats[k][c] * mol_dens[0][c] * 1.e-3;
    }
  }
}

void
convertConcentrationToATS(const Epetra_MultiVector& mol_dens,
                          int num_aqueous,
                          const Epetra_MultiVector& tcc_amanzi,
                          Epetra_MultiVector& tcc_ats)
{
  // convert from [mol C / L] to mol fraction [mol C / mol H20]
  for (int k = 0; k != num_aqueous; ++k) {
    for (int c = 0; c != tcc_amanzi.MyLength(); ++c) {
      tcc_ats[k][c] = tcc_amanzi[k][c] / (mol_dens[0][c] * 1.e-3);
    }
  }
}


bool
advanceChemistry(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
                 double t_old,
                 double t_new,
                 bool reinit,
                 const Epetra_MultiVector& mol_dens,
                 Teuchos::RCP<Epetra_MultiVector> tcc,
                 Teuchos::Time& timer)
{
  bool fail = false;
  int num_aqueous = chem_pk->num_aqueous_components();
  convertConcentrationToAmanzi(mol_dens, num_aqueous, *tcc, *tcc);
  chem_pk->set_aqueous_components(tcc);

  {
    auto monitor = Teuchos::rcp(new Teuchos::TimeMonitor(timer));
    fail = chem_pk->AdvanceStep(t_old, t_new, reinit);
  }
  if (fail) return fail;

  *tcc = *chem_pk->aqueous_components();
  convertConcentrationToATS(mol_dens, num_aqueous, *tcc, *tcc);
  return fail;
}


} //  namespace Amanzi
