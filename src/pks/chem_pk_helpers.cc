/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

#include "Epetra_MultiVector.h"

#include "Key.hh"
#include "State.hh"
#include "PK_Helpers.hh"
#include "chem_pk_helpers.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Helper functions for working with Amanzi's Chemistry PK
// -----------------------------------------------------------------------------
void
convertConcentrationToMolFrac(State& S,
        const KeyTag& tcc,
        const KeyTag& mol_frac,
        const KeyTag& mol_dens,
        const std::string& passwd)
{
  const Epetra_MultiVector& tcc_c =
    *S.Get<CompositeVector>(tcc.first, tcc.second).ViewComponent("cell", false);
  Epetra_MultiVector& mol_frac_c =
    *S.GetW<CompositeVector>(mol_frac.first, mol_frac.second, passwd).ViewComponent("cell", false);

  S.GetEvaluator(mol_dens.first, mol_dens.second).Update(S, passwd);
  const Epetra_MultiVector& mol_dens_c =
    *S.Get<CompositeVector>(mol_dens.first, mol_dens.second).ViewComponent("cell", false);

  // convert from mole fraction [mol C / mol H20] to [mol C / L]
  int ierr = mol_frac_c.ReciprocalMultiply(1.e3, mol_dens_c, tcc_c, 0.);
  AMANZI_ASSERT(!ierr);
  changedEvaluatorPrimary(mol_frac.first, mol_frac.second, S);
}

void
convertMolFracToConcentration(State& S,
        const KeyTag& mol_frac,
        const KeyTag& tcc,
        const KeyTag& mol_dens,
        const std::string& passwd)
{
  const Epetra_MultiVector& mol_frac_c =
    *S.Get<CompositeVector>(mol_frac.first, mol_frac.second).ViewComponent("cell", false);
  Epetra_MultiVector& tcc_c =
    *S.GetW<CompositeVector>(tcc.first, tcc.second, passwd).ViewComponent("cell", false);

  S.GetEvaluator(mol_dens.first, mol_dens.second).Update(S, passwd);
  const Epetra_MultiVector& mol_dens_c =
    *S.Get<CompositeVector>(mol_dens.first, mol_dens.second).ViewComponent("cell", false);

  int ierr = tcc_c.Multiply(1.e-3, mol_dens_c, mol_frac_c, 0.);
  AMANZI_ASSERT(!ierr);
  changedEvaluatorPrimary(tcc.first, tcc.second, S);
}

} //  namespace Amanzi
