/*
ATSDriver:
Runs top-level time-loop for standalone ATS simulations
*/

#pragma once

#include "Key.hh"
#include "coordinator.hh"


namespace Amanzi {
class State;
};

namespace ATS {

class ATSDriver : public Coordinator {

public:
  ATSDriver(Teuchos::ParameterList& parameter_list,
                 Teuchos::RCP<Amanzi::State>& S,
                 Amanzi::Comm_ptr_type comm);

  ~ATSDriver() = default;

  // methods
  void cycle_driver();

};

} // namespace ATS
