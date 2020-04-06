/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Two-phase overland flow equation.

/*!

This modifies the diffusion wave equation for overland flow that includes
freeze-thaw processes.  This class could completely go away, but it does some
error checking on the input file to make sure freeze-thaw processes are done
correctly.  In the future this service should be done by a preprocessor
generating the input file, and this class would go away.

.. _icy-overland-spec:
.. admonition:: icy-overland-spec

    INCLUDES:

    - ``[overland-pressure-spec]`` See `Overland Flow PK`_.

*/

#ifndef PK_FLOW_ICY_OVERLAND_HH_
#define PK_FLOW_ICY_OVERLAND_HH_

#include "Teuchos_TimeMonitor.hpp"

#include "PK_Factory.hh"
#include "overland_pressure.hh"
#include "icy_height_model.hh"

namespace Amanzi {

namespace Operators { class Upwinding; }

namespace Flow {


class IcyOverlandFlow : public OverlandPressureFlow {

 public:

  IcyOverlandFlow(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& solution) :
    PK(pk_tree, global_list, S, solution),
    OverlandPressureFlow(pk_tree, global_list, S, solution) {}

  // Virtual destructor
  virtual ~IcyOverlandFlow() {}

 protected:
  // setup methods
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

 private:
  // factory registration
  static RegisteredPKFactory<IcyOverlandFlow> reg_;
};

}  // namespace AmanziFlow
}  // namespace AmanziFlow

#endif
