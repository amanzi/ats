/*
ELM-ATS Coordinator:
Contains functions and data for running ATS as 
an ELM external model. Instantiates and initializes
State variables required by ELM and provides a 
method to advance ATS. 

*/

#pragma once

#include "Teuchos_Time.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MpiComm.h"
#include "AmanziComm.hh"
#include "AmanziTypes.hh"
#include "Key.hh"

#include "VerboseObject.hh"

#include "coordinator.hh"


namespace Amanzi {
class State;
};


namespace ATS {

class ELM_ATSCoordinator : public Coordinator {

public:
  ELM_ATSCoordinator(Teuchos::ParameterList& parameter_list,
              Teuchos::RCP<Amanzi::State>& S,
              Amanzi::Comm_ptr_type comm);
  ~ELM_ATSCoordinator() = default;

  // methods
  virtual void setup();
  virtual void initialize();
  virtual void reinit(double start_time=0.0, bool visout=false);
  virtual bool advance(double dt, bool visout=false, bool chkout=false);
  virtual void finalize();

  // for testing
  double get_end_time() {return t1_;}

};

} // namespace ATS
