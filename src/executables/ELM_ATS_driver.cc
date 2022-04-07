
#include <iostream>
#include <unistd.h>
#include <sys/resource.h>
#include "errors.hh"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "AmanziComm.hh"
#include "AmanziTypes.hh"

#include "InputAnalysis.hh"

#include "Units.hh"
#include "CompositeVector.hh"
#include "TimeStepManager.hh"
#include "Visualization.hh"
#include "VisualizationDomainSet.hh"
#include "IO.hh"
#include "Checkpoint.hh"
#include "UnstructuredObservations.hh"
#include "State.hh"
#include "PK.hh"
#include "TreeVector.hh"
#include "PK_Factory.hh"


#include "pk_helpers.hh"
#include "ELM_ATS_driver.hh"

namespace ATS {


  ELM_ATSDriver::ELM_ATSDriver(Teuchos::ParameterList& parameter_list,
                         Teuchos::RCP<Amanzi::State>& S,
                         Amanzi::Comm_ptr_type comm ) :
    Coordinator(parameter_list, S, comm) 

    {
      domain_sub_ = parameter_list_->get<std::string>("domain name", "domain");
      domain_srf_ = Amanzi::Keys::readDomainHint(*parameter_list_, domain_sub_, "subsurface", "surface");

      // surface and subsurface source_sink keys
      // these are the main coupling terms
      sub_src_key_ = Amanzi::Keys::readKey(*parameter_list_, domain_sub_, "source", "source_sink");
      srf_src_key_ = Amanzi::Keys::readKey(*parameter_list_, domain_srf_, "source", "source_sink");
      
      // soil state
      pres_key_ = Amanzi::Keys::readKey(*parameter_list_, domain_sub_, "pressure", "pressure");
      satl_key_ = Amanzi::Keys::readKey(*parameter_list_, domain_sub_, "saturation_liquid", "saturation_liquid");

      // soil properties - does this change?
      por_key_ = Amanzi::Keys::readKey(*parameter_list_, domain_sub_, "porosity", "porosity");
      
    };

  // PK methods
void ELM_ATSDriver::setup() {

  // common constants
  S_->Require<double>("atmospheric_pressure",
                      Amanzi::Tags::DEFAULT, "coordinator");
  S_->Require<Amanzi::AmanziGeometry::Point>("gravity",
          Amanzi::Tags::DEFAULT, "coordinator");

  // needed other times
  S_->Require<double>("time", Amanzi::Tags::CURRENT, "time");
  S_->Require<double>("time", Amanzi::Tags::NEXT, "time");

  // set initial tags to t0_
  S_->set_time(Amanzi::Tags::CURRENT, t0_);
  S_->set_time(Amanzi::Tags::NEXT, t0_);

  const auto& mesh_subsurf = S_->GetMesh(domain_sub_);
  const auto& mesh_surf = S_->GetMesh(domain_srf_);

  // require primary variables
  // -- subsurface water source
  S_->Require<Amanzi::CompositeVector,Amanzi::CompositeVectorSpace>(sub_src_key_, Amanzi::Tags::NEXT,  sub_src_key_)
    .SetMesh(mesh_subsurf)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  RequireEvaluatorPrimary(sub_src_key_, Amanzi::Tags::NEXT, *S_);
  // -- surface water source-sink
  S_->Require<Amanzi::CompositeVector,Amanzi::CompositeVectorSpace>(srf_src_key_, Amanzi::Tags::NEXT,  srf_src_key_)
    .SetMesh(mesh_surf)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  RequireEvaluatorPrimary(srf_src_key_, Amanzi::Tags::NEXT, *S_);

  // require soil state
  // -- subsurface pressure
  S_->RequireEvaluator(pres_key_, Amanzi::Tags::NEXT);
  S_->Require<Amanzi::CompositeVector,Amanzi::CompositeVectorSpace>(pres_key_, Amanzi::Tags::NEXT)
    .SetMesh(mesh_subsurf)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  // -- liquid saturation
  S_->RequireEvaluator(satl_key_, Amanzi::Tags::NEXT);
  S_->Require<Amanzi::CompositeVector,Amanzi::CompositeVectorSpace>(satl_key_, Amanzi::Tags::NEXT)
    .SetMesh(mesh_subsurf)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);

  // require soil properties
  // -- porosity
  S_->RequireEvaluator(por_key_, Amanzi::Tags::NEXT);
  S_->Require<Amanzi::CompositeVector,Amanzi::CompositeVectorSpace>(por_key_, Amanzi::Tags::NEXT)
    .SetMesh(mesh_subsurf)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);

  // order matters here -- PKs set the leaves, then observations can use those
  // if provided, and setup finally deals with all secondaries and allocates memory
  pk_->Setup();
  for (auto& obs : observations_) obs->Setup(S_.ptr());
  S_->Setup();

}

void ELM_ATSDriver::initialize() {

  // this is probably wrong
  // I think I may need to initialize the evaluators around the same time the pks call Initialize()
  // after State->InitializeFields() and before State does its final checks,  maybe
  // which means fully replacing Coordinator::initialize()
  // but I'm not sure, so I'll try it first 
  Coordinator::initialize(); //?

  S_->GetW<Amanzi::CompositeVector>(sub_src_key_, Amanzi::Tags::NEXT, sub_src_key_).PutScalar(0.);
  S_->GetRecordW(sub_src_key_, Amanzi::Tags::NEXT, sub_src_key_).set_initialized();
  S_->GetW<Amanzi::CompositeVector>(srf_src_key_, Amanzi::Tags::NEXT, srf_src_key_).PutScalar(0.);
  S_->GetRecordW(srf_src_key_, Amanzi::Tags::NEXT, srf_src_key_).set_initialized();

}

bool ELM_ATSDriver::advance() {

  // !!need to get from ELM unless TAGS are controlled by ELM, but use this for now!!
  double dt = S_->Get<double>("dt", Amanzi::Tags::DEFAULT);
  double t_old = S_->get_time(Amanzi::Tags::CURRENT);
  double t_new = S_->get_time(Amanzi::Tags::NEXT);

  // Get incoming state from ELM
  Epetra_MultiVector& srf_water_src = *S_->GetW<Amanzi::CompositeVector>(srf_src_key_, Amanzi::Tags::NEXT, srf_src_key_)
    .ViewComponent("cell", false);
  Epetra_MultiVector& sub_water_src = *S_->GetW<Amanzi::CompositeVector>(sub_src_key_, Amanzi::Tags::NEXT, sub_src_key_)
    .ViewComponent("cell", false);
  // getELMSources() - !!need to create!!
  ChangedEvaluatorPrimary(srf_src_key_, Amanzi::Tags::NEXT, *S_);
  ChangedEvaluatorPrimary(sub_src_key_, Amanzi::Tags::NEXT, *S_);

  // get other ELM outputs: - there will be more
  //S_->GetEvaluator(pres_key_, Amanzi::Tags::NEXT).Update(*S_, name_);
  //Epetra_MultiVector& pres = *S_->GetW<Amanzi::CompositeVector>(pres_key_, Amanzi::Tags::NEXT).ViewComponent("cell", false);
  //S_->GetEvaluator(poro_key_, Amanzi::Tags::NEXT).Update(*S_, name_);
  //const Epetra_MultiVector& poro = *S_->Get<Amanzi::CompositeVector>(poro_key_, Amanzi::Tags::NEXT).ViewComponent("cell", false);
  //S_->GetEvaluator(sl_key_, Amanzi::Tags::NEXT).Update(*S_, name_);
  //const Epetra_MultiVector& sl = *S_->Get<Amanzi::CompositeVector>(sl_key_, Amanzi::Tags::NEXT).ViewComponent("cell", false);


  auto fail = Coordinator::advance(); // advance pks

  if (fail) {
    Errors::Message msg("ELM_ATSDriver: Coordinator advance failed.");
    Exceptions::amanzi_throw(msg);
  }

  // update ATS->ELM data
  S_->GetEvaluator(pres_key_, Amanzi::Tags::NEXT).Update(*S_, pres_key_);
  const Epetra_MultiVector& pres = *S_->Get<Amanzi::CompositeVector>(pres_key_, Amanzi::Tags::NEXT).ViewComponent("cell", false);
  S_->GetEvaluator(por_key_, Amanzi::Tags::NEXT).Update(*S_, por_key_);
  const Epetra_MultiVector& poro = *S_->Get<Amanzi::CompositeVector>(por_key_, Amanzi::Tags::NEXT).ViewComponent("cell", false);
  S_->GetEvaluator(satl_key_, Amanzi::Tags::NEXT).Update(*S_, satl_key_);
  const Epetra_MultiVector& satl = *S_->Get<Amanzi::CompositeVector>(satl_key_, Amanzi::Tags::NEXT).ViewComponent("cell", false);

  // // getATSData() - !!need to create!!

  return fail;


  

  }


} // namespace ATS
