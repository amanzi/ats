/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   ATS

   Default base with default implementations of methods for a physical PK.
   ------------------------------------------------------------------------- */

#include "EvaluatorPrimary.hh"
#include "StateDefs.hh"
#include "pk_helpers.hh"
#include "pk_physical_default.hh"

namespace Amanzi {

PK_Physical_Default::PK_Physical_Default(Teuchos::ParameterList& pk_tree,
                                         const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                         const Teuchos::RCP<State>& S,
                                         const Teuchos::RCP<TreeVector>& solution)
  : PK(pk_tree, glist, S, solution),
    PK_Physical(pk_tree, glist, S, solution)
{}


// -----------------------------------------------------------------------------
// Construction of data.
// -----------------------------------------------------------------------------



} // namespace Amanzi
