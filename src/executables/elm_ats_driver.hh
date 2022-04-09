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

struct ELM_ATSDriver {

// PK methods
int setup(char *input_filename);
void initialize();
void advance(double *dt);

void advance_test();

std::unique_ptr<ELM_ATSCoordinator> elm_coordinator_;
Teuchos::RCP<Amanzi::State> S_;

};

} // namespace

