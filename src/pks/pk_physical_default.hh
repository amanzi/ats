/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A base class with default implementations of methods for a leaf of the PK tree (a conservation equation, or similar).
/*!

`PKPhysicalBase` is a base class providing some functionality for PKs which
are defined on a single mesh, and represent a single process model.  Typically
all leaves of the PK tree will inherit from `PKPhysicalBase`.

.. _pk-physical-default-spec:
.. admonition:: pk-physical-default-spec

    * `"domain name`" ``[string]`` Name from the Mesh_ list on which this PK is defined.

    * `"primary variable key`" ``[string]`` The primary variable
      e.g. `"pressure`", or `"temperature`".  Most PKs supply sane defaults.

    * `"initial condition`" ``[initial-conditions-spec]``  See InitialConditions_.

    INCLUDES:

    - ``[pk-spec]`` This *is a* PK_.
    - ``[debugger-spec]`` Uses a Debugger_

*/
#ifndef ATS_PK_PHYSICAL_BASE_HH_
#define ATS_PK_PHYSICAL_BASE_HH_

#include "Teuchos_ParameterList.hpp"
#include "TreeVector.hh"

#include "Debugger.hh"

#include "EvaluatorPrimary.hh"
#include "PK.hh"
#include "PK_Physical.hh"

namespace Amanzi {

//
// CAN THIS FINALLY GO AWAY?!?!?
//
class PK_Physical_Default : public PK_Physical {
 public:
  PK_Physical_Default(Teuchos::ParameterList& pk_tree,
                      const Teuchos::RCP<Teuchos::ParameterList>& glist,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~PK_Physical_Default() = default;

};

} // namespace Amanzi

#endif
