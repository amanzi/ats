/*


ELM has a 3D grid
ELM data are passed in either as scalars (surface only variables)
or as 1D column arrays (subsurface, subsurface+surface variables)
The number of cells per column is likely known prior to runtime

This class should do the following:

1. Allow factory registration
2. Constructable by either coordinator or an mpc
3. Setup() - setup required data and structures - mark evals as 
   primary and changed when appropriate
4. Initialize() 


Something needs to have control over the other pks

Need different ELM_ATS pk and ELM_ATS driver
ELM_PK : standard surface_balance PK, but only exchanges info. doesn't do any calculation 
of variables, only derivatives and such

ELM_Driver : handles weirdness of ATS getting called - somewhat similar to an mpc?






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

class ELM_ATSDriver : public Coordinator {

public:
  ELM_ATSDriver(Teuchos::ParameterList& parameter_list,
              Teuchos::RCP<Amanzi::State>& S,
              Amanzi::Comm_ptr_type comm);

  // PK methods
  virtual void setup();
  virtual void initialize();
  virtual bool advance();

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
