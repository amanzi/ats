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
                             const Epetra_MultiVector& tcc_ats,
                             Epetra_MultiVector& tcc_amanzi)
{
  // convert from mole fraction [mol C / mol H20] to [mol C / L]
  tcc_amanzi.Multiply(1.e-3, mol_dens, tcc_ats, 0.);
}

void
convertConcentrationToATS(const Epetra_MultiVector& mol_dens,
                          const Epetra_MultiVector& tcc_amanzi,
                          Epetra_MultiVector& tcc_ats)
{
  // convert from [mol C / L] to mol fraction [mol C / mol H20]
  tcc_ats.ReciprocalMultiply(1.e3, mol_dens, tcc_amanzi, 0.);
}


bool
advanceChemistry(AmanziChemistry::Chemistry_PK& chem_pk,
                 double t_old,
                 double t_new,
                 bool reinit,
                 const Epetra_MultiVector& mol_dens,
                 Epetra_MultiVector& tcc,
                 Teuchos::Time& timer)
{
  bool fail = false;
  convertConcentrationToAmanzi(mol_dens, tcc, tcc);
  {
    auto monitor = Teuchos::rcp(new Teuchos::TimeMonitor(timer));
    fail = chem_pk.AdvanceStep(t_old, t_new, reinit);
  }
  if (fail) return fail;

  convertConcentrationToATS(mol_dens, tcc, tcc);
  return fail;
}


} //  namespace Amanzi
