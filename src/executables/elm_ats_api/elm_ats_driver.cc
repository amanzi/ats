#include <iostream>
#include <unistd.h>
#include <sys/resource.h>
#include <filesystem>
#include "errors.hh"
#include "dbc.hh"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "VerboseObject.hh"

// registration files
#include "ats_registration_files.hh"

#include "AmanziComm.hh"
#include "CompositeVector.hh"
#include "IO.hh"
#include "UnstructuredObservations.hh"
#include "PK_Helpers.hh"
#include "elm_ats_driver.hh"

namespace ATS {

ELM_ATSDriver*
createELM_ATSDriver(MPI_Fint *f_comm, const char *infile, int npfts) {
  // -- create communicator & get process rank
  //auto comm = getDefaultComm();
  auto c_comm = MPI_Comm_f2c(*f_comm);
  auto comm = getComm(c_comm);
  auto teuchos_comm = Teuchos::rcp(new Teuchos::MpiComm<int>(c_comm));
  auto rank = comm->MyPID();

  // convert input file to std::string for easier handling
  // infile must be null-terminated
  std::string input_filename(infile);

  // check validity of input file name
  if (input_filename.empty()) {
    if (rank == 0)
      std::cerr << "ERROR: no input file provided" << std::endl;
  } else if (!std::filesystem::exists(input_filename)) {
    if (rank == 0)
      std::cerr << "ERROR: input file \"" << input_filename << "\" does not exist." << std::endl;
  }

  // -- parse input file
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(input_filename);

  auto wallclock_timer = Teuchos::TimeMonitor::getNewCounter("wallclock duration");
  return new ELM_ATSDriver(plist, wallclock_timer, teuchos_comm, comm, npfts);
}


ELM_ATSDriver::ELM_ATSDriver(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                             const Teuchos::RCP<Teuchos::Time>& wallclock_timer,
                             const Teuchos::RCP<const Teuchos::Comm<int>>& teuchos_comm,
                             const Amanzi::Comm_ptr_type& comm,
                             int npfts)
  : Coordinator(plist, wallclock_timer, teuchos_comm, comm),
    npfts_(npfts),
    ncolumns_(-1),
    ncells_per_col_(-1)
{
  // -- set default verbosity level to no output
  // -- TODO make the verbosity level an input argument
  VerboseObject::global_default_level = Teuchos::VERB_NONE;

  // domain names
  domain_subsurf_ = Keys::readDomain(*plist_, "domain", "domain");
  domain_surf_ = Keys::readDomainHint(*plist_, domain_subsurf_, "subsurface", "surface");

   // meshes
  mesh_subsurf_ = S_->GetMesh(domain_subsurf_);
  mesh_surf_ = S_->GetMesh(domain_surf_);

  // potential sources
  gross_water_source_key_ = Keys::readKey(*plist_, domain_surf_, "gross water source", "gross_water_source");
  key_map_[ELM::VarID::GROSS_SURFACE_WATER_SOURCE] = gross_water_source_key_;
  pot_evap_key_ = Keys::readKey(*plist_, domain_surf_, "potential evaporation mps", "potential_evaporation_mps");
  key_map_[ELM::VarID::POTENTIAL_EVAPORATION] = pot_evap_key_;
  pot_trans_key_ = Keys::readKey(*plist_, domain_surf_, "potential transpiration mps", "potential_transpiration_mps");
  key_map_[ELM::VarID::POTENTIAL_TRANSPIRATION] = pot_trans_key_;

  // water state
  surf_wc_key_ = Keys::readKey(*plist_, domain_surf_, "water content", "water_content");
  key_map_[ELM::VarID::SURFACE_WATER_CONTENT] = surf_wc_key_;
  wc_key_ = Keys::readKey(*plist_, domain_subsurf_, "water content", "water_content");
  key_map_[ELM::VarID::WATER_CONTENT] = wc_key_;

  // actual water fluxes
  evap_key_ = Keys::readKey(*plist_, domain_surf_, "evaporation", "evaporation");
  key_map_[ELM::VarID::EVAPORATION] = evap_key_;
  col_trans_key_ = Keys::readKey(*plist_, domain_surf_, "total transpiration", "total_transpiration");
  key_map_[ELM::VarID::TRANSPIRATION] = col_trans_key_;
  col_baseflow_key_ = Keys::readKey(*plist_, domain_surf_, "column total baseflow", "baseflow");
  key_map_[ELM::VarID::BASEFLOW] = col_baseflow_key_;
  col_runoff_key_ = Keys::readKey(*plist_, domain_surf_, "runoff generation", "runoff_generation");
  key_map_[ELM::VarID::RUNOFF] = col_runoff_key_;

  // keys for fields used to convert ELM units to ATS units
  surf_mol_dens_key_ = Keys::readKey(*plist_, domain_surf_, "surface molar density", "molar_density_liquid");
  surf_mass_dens_key_ = Keys::readKey(*plist_, domain_surf_, "surface mass density", "mass_density_liquid");
  subsurf_mol_dens_key_ = Keys::readKey(*plist_, domain_subsurf_, "molar density", "molar_density_liquid");
  subsurf_mass_dens_key_ = Keys::readKey(*plist_, domain_subsurf_, "mass density", "mass_density_liquid");

  // cell vol keys
  surf_cv_key_ = Keys::getKey(domain_surf_, "cell_volume");
  cv_key_ = Keys::getKey(domain_subsurf_, "cell_volume");

  // -- mesh info
  ncolumns_ = mesh_surf_->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  auto col_zero = mesh_subsurf_->columns.getCells(0);
  ncells_per_col_ = col_zero.size();
}


void
ELM_ATSDriver::setup()
{
  // potential fluxes (ELM -> ATS)
  ATS::requireEvaluatorAtNext(gross_water_source_key_, Amanzi::Tags::NEXT, *S_, gross_water_source_key_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(pot_evap_key_, Amanzi::Tags::NEXT, *S_, pot_evap_key_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(pot_trans_key_, Amanzi::Tags::NEXT, *S_, pot_trans_key_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  ATS::requireEvaluatorAtNext(surf_cv_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(cv_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  ATS::requireEvaluatorAtNext(surf_mol_dens_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(surf_mass_dens_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(subsurf_mol_dens_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(subsurf_mass_dens_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  ATS::requireEvaluatorAtNext(evap_key_, Amanzi::Tags::NEXT, *S_)
   .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(col_trans_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(col_baseflow_key_, Amanzi::Tags::NEXT, *S_, col_baseflow_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(col_runoff_key_, Amanzi::Tags::NEXT, *S_, col_runoff_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", AmanziMesh::CELL, 1);

  Coordinator::setup();

  ATS::requireEvaluatorAtNext(pres_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
}


//
// dz & depth are currently ignored -- presumed identical between ATS & ELM
//
ELM_ATSDriver::MeshInfo ELM_ATSDriver::getMeshInfo()
{
  ELM_ATSDriver::MeshInfo info;
  info.ncols_local = ncolumns_;
  info.ncols_global = mesh_surf_->getMap(AmanziMesh::Entity_kind::CELL, false).NumGlobalElements();
  info.nlevgrnd = ncells_per_col_;
  return info;
}


void ELM_ATSDriver::initialize()
{
  // Assign start time to ATS
  t0_ = S_->get_time();

  // initialize ATS data, commit initial conditions
  Coordinator::initialize();

  // set as zero and tag as initialized
  initZero_(gross_water_source_key_);
  initZero_(pot_trans_key_);
  initZero_(pot_evap_key_);

  // no current evaluators for these, treat as primary
  initZero_(col_baseflow_key_);
  initZero_(col_runoff_key_);

  // visualization at IC -- TODO remove this or place behind flag
  visualize();
  checkpoint();
}



void ELM_ATSDriver::advance(double dt, bool do_chkp, bool do_vis)
{
  S_->Assign<double>("dt", Amanzi::Tags::DEFAULT, "dt", dt);
  S_->advance_time(Amanzi::Tags::NEXT, dt);

  // solve model for a single timestep
  bool fail = Coordinator::advance();
  if (fail) {
    Errors::Message msg("ELM_ATSDriver: advance(dt) failed.  Make ATS subcycle for proper ELM use.");
    Exceptions::amanzi_throw(msg);
  }

  // make observations, vis, and checkpoints
  for (const auto& obs : observations_) obs->MakeObservations(S_.ptr());

  // vis/checkpoint if EITHER ATS or ELM request it
  if (do_vis && !visualize()) visualize(true);
  if (do_chkp && !checkpoint()) checkpoint(true);
}


// simulates external timeloop with dt coming from calling model
void ELM_ATSDriver::advanceTest()
{
  while (S_->get_time() < t1_) {
    // use dt from ATS for testing
    double dt = Coordinator::get_dt(false);
    // call main method
    advance(dt, false, false);
  }
}

void ELM_ATSDriver::finalize()
{
  WriteStateStatistics(*S_, *vo_);
  report_memory();
  Teuchos::TimeMonitor::summarize(*vo_->os());
  Coordinator::finalize();
}


void ELM_ATSDriver::initZero_(const Key& key)
{
  auto& vec = S_->GetW<CompositeVector>(key, Amanzi::Tags::NEXT, key);
  vec.PutScalar(0.);
  S_->GetRecordW(key, Amanzi::Tags::NEXT, key).set_initialized();
}




void ELM_ATSDriver::copyToSurf_(double const * const in, const Key& key)
{
  Key owner = S_->GetRecord(key, Amanzi::Tags::NEXT).owner();

  // surf maps directly into columns
  auto& vec = *S_->GetW<CompositeVector>(key, Amanzi::Tags::NEXT, owner)
    .ViewComponent("cell", false);
  AMANZI_ASSERT(vec.MyLength() == ncolumns_);

  for (int i=0; i!=ncolumns_; ++i) vec[0][i] = in[i];

  changedEvaluatorPrimary(key, Amanzi::Tags::NEXT, *S_);
}


//
// ELM data is defined as var(col,lev), meaning that, due to Fortran ordering,
// the column is fasted varying, not the grid cell.
//
void ELM_ATSDriver::copyToSub_(double const * const in, const Key& key)
{
  Key owner = S_->GetRecord(key, Amanzi::Tags::NEXT).owner();
  auto& vec = *S_->GetW<CompositeVector>(key, Amanzi::Tags::NEXT, owner)
    .ViewComponent("cell", false);

  for (int i=0; i!=ncolumns_; ++i) {
    const auto& cells = mesh_subsurf_->columns.getCells(i);
    for (int j=0; j!=ncells_per_col_; ++j) {
      vec[0][cells[j]] = in[j * ncolumns_ + i];
    }
  }

  changedEvaluatorPrimary(key, Amanzi::Tags::NEXT, *S_);
}


void ELM_ATSDriver::copyFromSurf_(double * const out, const Key& key) const
{
  S_->GetEvaluator(key, Amanzi::Tags::NEXT).Update(*S_, "ELM");
  const auto& vec = *S_->Get<CompositeVector>(key, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);

  for (int i=0; i!=ncolumns_; ++i) out[i] = vec[0][i];
}


//
// ELM data is defined as var(col,lev), meaning that, due to Fortran ordering,
// the column is fasted varying, not the grid cell.
//
void ELM_ATSDriver::copyFromSub_(double * const out, const Key& key) const
{
  S_->GetEvaluator(key, Amanzi::Tags::NEXT).Update(*S_, "ELM");
  const auto& vec = *S_->Get<CompositeVector>(key, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);

  for (int i=0; i!=ncolumns_; ++i) {
    const auto cells = mesh_subsurf_->columns.getCells(i);
    for (int j=0; j!=ncells_per_col_; ++j) {
      out[j * ncolumns_ + i] = vec[0][cells[j]];
    }
  }
}


void ELM_ATSDriver::setScalar(const ELM::VarID& scalar_id, double in)
{
  Key scalar_key = key_map_[scalar_id];
  S_->GetW<double>(scalar_key, Tags::DEFAULT, scalar_key) = in;
}

double ELM_ATSDriver::getScalar(const ELM::VarID& scalar_id)
{
  Key scalar_key = key_map_[scalar_id];
  return S_->Get<double>(scalar_key, Tags::DEFAULT);
}


void ELM_ATSDriver::getField(const ELM::VarID& var_id, double * const out)
{
  switch(var_id) {
    case ELM::VarID::SURFACE_WATER_CONTENT: getField_<ELM::VarID::SURFACE_WATER_CONTENT>(out); break;
    case ELM::VarID::WATER_CONTENT: getField_<ELM::VarID::WATER_CONTENT>(out); break;
    case ELM::VarID::EVAPORATION: getField_<ELM::VarID::EVAPORATION>(out); break;
    case ELM::VarID::TRANSPIRATION: getField_<ELM::VarID::TRANSPIRATION>(out); break;
    case ELM::VarID::RUNOFF: getField_<ELM::VarID::RUNOFF>(out); break;
    case ELM::VarID::BASEFLOW: getField_<ELM::VarID::BASEFLOW>(out); break;
    default: throw std::runtime_error("Ungettable variable");
  }
}



void ELM_ATSDriver::setField(const ELM::VarID& var_id, double const * const out)
{
  switch(var_id) {
    case ELM::VarID::BASE_POROSITY: setField_<ELM::VarID::BASE_POROSITY>(out); break;
    case ELM::VarID::HYDRAULIC_CONDUCTIVITY: setField_<ELM::VarID::HYDRAULIC_CONDUCTIVITY>(out); break;
    case ELM::VarID::CLAPP_HORNBERGER_B: setField_<ELM::VarID::CLAPP_HORNBERGER_B>(out); break;
    case ELM::VarID::CLAPP_HORNBERGER_PSI_SAT: setField_<ELM::VarID::CLAPP_HORNBERGER_PSI_SAT>(out); break;
    case ELM::VarID::RESIDUAL_SATURATION: setField_<ELM::VarID::RESIDUAL_SATURATION>(out); break;
    case ELM::VarID::EFFECTIVE_POROSITY: setField_<ELM::VarID::EFFECTIVE_POROSITY>(out); break;
    case ELM::VarID::ROOT_FRACTION: setField_<ELM::VarID::ROOT_FRACTION>(out); break;
    case ELM::VarID::SURFACE_WATER_CONTENT: setField_<ELM::VarID::SURFACE_WATER_CONTENT>(out); break;
    case ELM::VarID::WATER_CONTENT: setField_<ELM::VarID::WATER_CONTENT>(out); break;
    case ELM::VarID::GROSS_SURFACE_WATER_SOURCE: setField_<ELM::VarID::GROSS_SURFACE_WATER_SOURCE>(out); break;
    case ELM::VarID::POTENTIAL_EVAPORATION: setField_<ELM::VarID::POTENTIAL_EVAPORATION>(out); break;
    case ELM::VarID::POTENTIAL_TRANSPIRATION: setField_<ELM::VarID::POTENTIAL_TRANSPIRATION>(out); break;
    default: throw std::runtime_error("Unsettable variable");
  }
}


double const *
ELM_ATSDriver::getFieldPtr(const ELM::VarID& var_id)
{
  Key var_key = key_map_[var_id];
  return &(*S_->Get<CompositeVector>(var_key, Tags::NEXT).ViewComponent("cell", false))[0][0];
}

double *
ELM_ATSDriver::getFieldPtrW(const ELM::VarID& var_id)
{
  Key var_key = key_map_[var_id];
  Key owner = S_->GetRecord(var_key, Amanzi::Tags::NEXT).owner();
  return &(*S_->GetW<CompositeVector>(var_key, Tags::NEXT, owner).ViewComponent("cell", false))[0][0];
}



template<ELM::VarID var_id>
void ELM_ATSDriver::getField_(double * const out)
{
  Key field_key = key_map_[var_id];
  if (Keys::getDomain(field_key) == domain_subsurf_) {
    copyFromSub_(out, field_key);
  } else {
    copyFromSurf_(out, field_key);
  }
}


template<ELM::VarID var_id>
void ELM_ATSDriver::setField_(double const * const out)
{
  Key field_key = key_map_[var_id];
  if (Keys::getDomain(field_key) == domain_subsurf_) {
    copyToSub_(out, field_key);
  } else {
    copyToSurf_(out, field_key);
  }
}




} // namespace ATS
