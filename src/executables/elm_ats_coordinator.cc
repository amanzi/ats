
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
#include "elm_ats_coordinator.hh"

namespace ATS {


  ELM_ATSCoordinator::ELM_ATSCoordinator(Teuchos::ParameterList& parameter_list,
                         Teuchos::RCP<Amanzi::State>& S,
                         Amanzi::Comm_ptr_type comm ) :
    Coordinator(parameter_list, S, comm) 

    {
      domain_sub_ = parameter_list_->get<std::string>("domain name", "domain");
      domain_srf_ = Amanzi::Keys::readDomainHint(*parameter_list_, domain_sub_, "subsurface", "surface");

      // surface and subsurface source_sink keys
      // these are the main coupling terms
      sub_src_key_ = Amanzi::Keys::readKey(*parameter_list_, domain_sub_, "subsurface source", "source_sink");
      srf_src_key_ = Amanzi::Keys::readKey(*parameter_list_, domain_srf_, "surface source", "source_sink");
      
      // soil state
      pres_key_ = Amanzi::Keys::readKey(*parameter_list_, domain_sub_, "pressure", "pressure");
      satl_key_ = Amanzi::Keys::readKey(*parameter_list_, domain_sub_, "saturation_liquid", "saturation_liquid");

      // soil properties - does this change?
      por_key_ = Amanzi::Keys::readKey(*parameter_list_, domain_sub_, "porosity", "porosity");
      
    };

  // PK methods
void ELM_ATSCoordinator::setup() {

  Teuchos::TimeMonitor monitor(*setup_timer_);

  // common constants
  S_->Require<double>("atmospheric_pressure",
                      Amanzi::Tags::DEFAULT, "coordinator");
  S_->Require<Amanzi::AmanziGeometry::Point>("gravity",
          Amanzi::Tags::DEFAULT, "coordinator");

  // needed other times
  S_->require_time(Amanzi::Tags::CURRENT);
  S_->require_time(Amanzi::Tags::NEXT);

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

  // order matters here -- PKs set the leaves, then observations can use those
  // if provided, and setup finally deals with all secondaries and allocates memory
  pk_->Setup();
  for (auto& obs : observations_) obs->Setup(S_.ptr());
  S_->Setup();

}

void ELM_ATSCoordinator::initialize() {

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

  // get the intial timestep
  if (!restart_) {
    double dt = get_dt(false);
    S_->Assign<double>("dt", Amanzi::Tags::DEFAULT, "dt", dt);
    std::cout << "dt init:: " << dt << std::endl;
  }


  // visualization at IC
  visualize();
  checkpoint();

  // Make sure times are set up correctly
  AMANZI_ASSERT(std::abs(S_->get_time(Amanzi::Tags::NEXT)
                         - S_->get_time(Amanzi::Tags::CURRENT)) < 1.e-4);
}

bool ELM_ATSCoordinator::advance(double dt_dummy) {

  // use dt from ATS for now
  double dt = S_->Get<double>("dt", Amanzi::Tags::DEFAULT);
  S_->Assign<double>("dt", Amanzi::Tags::DEFAULT, "dt", dt);
  S_->advance_time(Amanzi::Tags::NEXT, dt);

  std::cout << "dt advance:: " << dt << std::endl;

  // !!need to get from ELM unless TAGS are controlled by ELM, but use this for now!!
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
    Errors::Message msg("ELM_ATSCoordinator: Coordinator advance failed.");
    Exceptions::amanzi_throw(msg);
  }

  S_->set_time(Amanzi::Tags::CURRENT, S_->get_time(Amanzi::Tags::NEXT));
  S_->advance_cycle();

  // make observations, vis, and checkpoints
  for (const auto& obs : observations_) obs->MakeObservations(S_.ptr());
  visualize();
  checkpoint(); // checkpoint with the new dt

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
