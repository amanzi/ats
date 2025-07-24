/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

An MPC that coordinates the coupling of flow and transport, or of integrated
flow an integrated transport.

This MPC simply sets a few default evaluators used in coupling (e.g. alias and
time-interpolated evaluators) based upon the knonw temporal discretizations and
coupling strategies of flow and transport.  Nothing done here couldn't also be
done in the input file by the user, but this shifts some of the burden away
from the user.

Specifically, this sets (as appropriate for surface-only, subsurface-only, or
integrated cases):

- Ensures that flow is not subcycled.
- Transport's flux field is given by flow's water_flux field, and is temporally
  piecewise constant throughout the flow interval using the value at flow's
  NEXT tag.

If transport is subcycled, then this additionally sets that:

- porosity is an ALIASED evaluator pointing to flow's NEXT tag (this likely
  will change, but for now changing it will break tests)
- saturation & molar density are temporally interpolated evaluators between
  flow's CURRENT and NEXT tags.

Note this always assumes that Flow is the _first_ PK in the `"PKs order`" list.

`"PK type`" = `"coupled flow and transport`"

.. _pk-coupled-flow-and-transport-spec:
.. admonition:: pk-coupled-flow-and-transport-spec

  INCLUDES:

  - ``[mpc-subcycled-spec]``

*/

#pragma once

#include "mpc_subcycled.hh"

namespace Amanzi {

class MPCFlowTransport : public MPCSubcycled {
 public:
  MPCFlowTransport(Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& soln)
    : PK(pk_tree, global_list, S, soln),
      MPCSubcycled(pk_tree, global_list, S, soln),
      chemistry_(false),
      surface_(false),
      subsurface_(false)
  {}

  void parseParameterList() override;
  void CommitStep(double t_old, double t_new, const Tag& tag) override;

 protected:
  bool chemistry_, surface_, subsurface_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCFlowTransport> reg_;
};

} // namespace Amanzi
