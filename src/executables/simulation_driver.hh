/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
// Runs the top-level simulation.


/* -*-  mode: c++; indent-tabs-mode: nil -*- */


/*!

ATS's top-level main accepts an XML list including a few required elements.

.. _main-spec:
.. admonition:: main-spec

    * `"mesh`" ``[mesh-typed-spec-list]`` A list of Mesh_ objects.
    * `"regions`" ``[region-typedsublist-spec-list]`` A list of Region_ objects.
    * `"cycle driver`" ``[coordinator-spec]``  See Coordinator_.
    * `"visualization`" ``[visualization-spec-list]`` A list of Visualization_ objects.
    * `"observations`" ``[observation-spec-list]`` An list of Observation_ objects.
    * `"checkpoint`" ``[checkpoint-spec]`` See Checkpoint_.
    * `"PKs`" ``[pk-typed-spec-list]`` A list of PK_ objects.
    * `"state`" ``[state-spec]`` See State_.

 */

#pragma once
#include "Teuchos_ParameterList.hpp"
#include "VerboseObject.hh"

namespace ATS {

struct SimulationDriver {
  int Run(const Teuchos::RCP<const Amanzi::Comm_type>& comm,
          Teuchos::ParameterList&       input_parameter_list);

};

} // namespace
