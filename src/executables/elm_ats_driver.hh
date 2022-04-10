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
  int setup(char *input_filename);
  void initialize();
  void advance(double *dt);
  void advance_test();

private:
  std::unique_ptr<ELM_ATSCoordinator> elm_coordinator_;
  Teuchos::RCP<Amanzi::State> S_;
};



//extern "C" {
//  ELM_ATSDriver *ELM_ATSDriver__new () {
//    return new ELM_ATSDriver();
//  }
//  int ELM_ATSDriver__area (ELM_ATSDriver *This) {
//    return This->area();
//  }
//  void ELM_ATSDriver__delete (ELM_ATSDriver *This) {
//    delete This;
//  }
//}



} // namespace

