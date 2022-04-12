/*
ELM-ATS Driver:
Provides an interface to ATS functionality for ELM 
*/
#pragma once

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "VerboseObject.hh"

#include "elm_ats_coordinator.hh"
#include "State.hh"

namespace ATS {

class ELM_ATSDriver {

public:
  // default constructor and destructor
  ELM_ATSDriver() {};
  ~ELM_ATSDriver() = default;

  // methods
  int setup(MPI_Fint *f_comm, char *input_filename);
  void initialize();
  void advance(double *dt);
  void advance_test();

private:
  std::unique_ptr<ELM_ATSCoordinator> elm_coordinator_;
  Teuchos::RCP<Amanzi::State> S_;
};


// include here temporarily during development
// maybe place into AmanziComm.hh
// or leave as local function
#include "AmanziTypes.hh"
#ifdef TRILINOS_TPETRA_STACK

#ifdef HAVE_MPI
#include "Teuchos_MpiComm.hpp"
#else
#include "Teuchos_SerialComm.hpp"
#endif

#else // Epetra stack

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#endif // trilinos stack

inline
Amanzi::Comm_ptr_type setComm(MPI_Comm comm) {
#ifdef TRILINOS_TPETRA_STACK
#ifdef HAVE_MPI
  return Teuchos::rcp(new Teuchos::MpiComm<int>(comm));
#else
  return Teuchos::rcp(new Teuchos::SerialComm<int>());
#endif
#else
#ifdef HAVE_MPI
  return Teuchos::rcp(new Epetra_MpiComm(comm));
#else
  return Teuchos::rcp(new Epetra_SerialComm());
#endif
#endif
}

} // namespace

