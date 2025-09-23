/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

This modifies the diffusion wave equation for overland flow that includes
freeze-thaw processes.

`"PK type`" = `"overland flow with ice`"

.. _pk-overland-flow-with-ice-spec:
.. admonition:: pk-overland-flow-with-ice-spec

    INCLUDES:

    - ``[overland-pressure-spec]`` See `Overland Flow PK`_.

*/

#ifndef PK_FLOW_ICY_OVERLAND_HH_
#define PK_FLOW_ICY_OVERLAND_HH_

#include "Teuchos_TimeMonitor.hpp"

#include "PK_Factory.hh"
#include "overland_pressure.hh"

namespace Amanzi {
namespace ATS_Physics {

namespace Operators {
class Upwinding;
}

namespace Flow {


class IcyOverlandFlow : public OverlandPressureFlow {
 public:
  IcyOverlandFlow(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& solution)
    : PK(pk_tree, global_list, S, solution), OverlandPressureFlow(pk_tree, global_list, S, solution)
  {}

  // Virtual destructor
  virtual ~IcyOverlandFlow() override {}

 protected:
  // setup methods
  virtual void SetupPhysicalEvaluators_() override;

 private:
  // factory registration
  static RegisteredPKFactory<IcyOverlandFlow> reg_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
