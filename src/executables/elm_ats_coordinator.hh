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

  // PK methods
  virtual void setup();
  virtual void initialize();
  virtual bool advance(double dt);

private:
Amanzi::Key domain_sub_;
Amanzi::Key domain_srf_;
Amanzi::Key sub_src_key_;
Amanzi::Key srf_src_key_;
Amanzi::Key pres_key_;
Amanzi::Key satl_key_;
Amanzi::Key por_key_;

};

} // namespace ATS
