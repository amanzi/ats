/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! The effective main of ATS

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!
ATS input accepts an XML list including a few required elements.

* `"mesh`" ``[mesh-typed-spec-list]`` A list of Mesh_ spec objects, with
  domain names given by the name of the sublist.

* `"regions`" ``[region-spec-list]`` A list of geometric Region_
  specs, with names given by the name of the sublist.

* `"cycle driver`" ``[coordinator-spec]``  See the Coordinator_ spec.

* `"visualization`" ``[visualization-spec-list]`` A list of
  Visualization_ specs, one for each mesh/domain, controlling
  simulation output that is sparse in time but across the entirety
  of the domain.

* `"observations`" ``[observation-spec-list]`` An list of Observation_
  specs, output that is sparse/local in space but across the entirety
  of time.

* `"checkpoint`" ``[checkpoint-spec]`` A Checkpoint_ spec, controlling
  output intended for restart.

* `"PKs`" ``[list]``  List of all PKs to be used in the simulation.

* `"state`" ``[state-spec]`` A State_ spec controlling evaluators and
  other data.

 */
  

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"

struct SimulationDriver
  : public Teuchos::VerboseObject<SimulationDriver>
{
  virtual int Run (const Amanzi::Comm_ptr_type& comm,
                   const Teuchos::RCP<Teuchos::ParameterList>& input_parameter_list);

};
