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
createELM_ATSDriver(MPI_Fint *f_comm, const char *infile, const char *logfile, int npfts) {
  // -- create communicator & get process rank
  //auto comm = getDefaultComm();
  auto c_comm = MPI_Comm_f2c(*f_comm);
  auto comm = getComm(c_comm);
  auto teuchos_comm = Teuchos::rcp(new Teuchos::MpiComm<int>(c_comm));
  auto rank = comm->MyPID();

  // convert input file to std::string for easier handling
  // infile must be null-terminated
  std::string input_filename(infile);
  std::string logfile_filename(logfile);

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
  return new ELM_ATSDriver(plist, wallclock_timer, teuchos_comm, comm, logfile_filename, npfts);
}


ELM_ATSDriver::ELM_ATSDriver(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                             const Teuchos::RCP<Teuchos::Time>& wallclock_timer,
                             const Teuchos::RCP<const Teuchos::Comm<int>>& teuchos_comm,
                             const Amanzi::Comm_ptr_type& comm,
                             const std::string& logfile,
                             int npfts)
  : Coordinator(plist, wallclock_timer, teuchos_comm, comm),
    npfts_(npfts),
    ncolumns_(-1),
    ncells_per_col_(-1)
{
  // -- set default verbosity level to no output
  // -- TODO make the verbosity level an input argument
  VerboseObject::global_default_level = Teuchos::VERB_HIGH;
  VerboseObject::global_logfile = logfile;
}


void
ELM_ATSDriver::parseParameterList()
{
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning parseParameterList stage..." << std::endl
               << std::flush;
  }
  // parse my parameter list
  // domain names
  domain_subsurf_ = Keys::readDomain(*plist_, "domain", "domain");
  domain_surf_ = Keys::readDomainHint(*plist_, domain_subsurf_, "subsurface", "surface");

   // meshes
  mesh_subsurf_ = S_->GetMesh(domain_subsurf_);
  mesh_surf_ = S_->GetMesh(domain_surf_);

  key_map_[ELM::VarID::TIME] = {"time", Tags::NEXT};

  // parameters
  poro_key_ = Keys::readKey(*plist_, domain_subsurf_, "porosity", "porosity");
  key_map_[ELM::VarID::EFFECTIVE_POROSITY] = { poro_key_, Tags::NEXT };

  // potential sources
  gross_water_source_key_ = Keys::readKey(*plist_, domain_surf_, "gross water source", "gross_water_source");
  key_map_[ELM::VarID::GROSS_SURFACE_WATER_SOURCE] = { gross_water_source_key_, Tags::NEXT };
  pot_evap_key_ = Keys::readKey(*plist_, domain_surf_, "potential evaporation mps", "potential_evaporation_mps");
  key_map_[ELM::VarID::POTENTIAL_EVAPORATION] = { pot_evap_key_, Tags::NEXT };
  pot_trans_key_ = Keys::readKey(*plist_, domain_surf_, "potential transpiration mps", "potential_transpiration_mps");
  key_map_[ELM::VarID::POTENTIAL_TRANSPIRATION] = { pot_trans_key_, Tags::NEXT };

  // water state
  pres_key_ = Keys::readKey(*plist_, domain_subsurf_, "pressure", "pressure");
  key_map_[ELM::VarID::PRESSURE] = { pres_key_, Tags::NEXT };
  sat_key_ = Keys::readKey(*plist_, domain_subsurf_, "saturation liquid", "saturation_liquid");
  key_map_[ELM::VarID::SATURATION_LIQUID] = { sat_key_, Tags::NEXT };
  pd_key_ = Keys::readKey(*plist_, domain_surf_, "ponded depth", "ponded_depth");
  key_map_[ELM::VarID::PONDED_DEPTH] = { pd_key_, Tags::NEXT };

  // actual water fluxes
  evap_key_ = Keys::readKey(*plist_, domain_surf_, "evaporation", "evaporation");
  key_map_[ELM::VarID::EVAPORATION] = { evap_key_, Tags::NEXT };
  col_trans_key_ = Keys::readKey(*plist_, domain_surf_, "total transpiration", "total_transpiration");
  key_map_[ELM::VarID::TRANSPIRATION] = { col_trans_key_, Tags::NEXT };
  col_baseflow_key_ = Keys::readKey(*plist_, domain_surf_, "column total baseflow", "baseflow");
  key_map_[ELM::VarID::BASEFLOW] = { col_baseflow_key_, Tags::NEXT };
  col_runoff_key_ = Keys::readKey(*plist_, domain_surf_, "runoff generation", "runoff_generation");
  key_map_[ELM::VarID::RUNOFF] = { col_runoff_key_, Tags::NEXT };

  // // keys for fields used to convert ELM units to ATS units
  // surf_mol_dens_key_ = Keys::readKey(*plist_, domain_surf_, "surface molar density", "molar_density_liquid");
  // surf_mass_dens_key_ = Keys::readKey(*plist_, domain_surf_, "surface mass density", "mass_density_liquid");
  // subsurf_mol_dens_key_ = Keys::readKey(*plist_, domain_subsurf_, "molar density", "molar_density_liquid");
  // subsurf_mass_dens_key_ = Keys::readKey(*plist_, domain_subsurf_, "mass density", "mass_density_liquid");

  // cell vol keys
  // surf_cv_key_ = Keys::getKey(domain_surf_, "cell_volume");
  // cv_key_ = Keys::getKey(domain_subsurf_, "cell_volume");

  // -- mesh info
  ncolumns_ = mesh_surf_->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  auto col_zero = mesh_subsurf_->columns.getCells(0);
  ncells_per_col_ = col_zero.size();

  // require my primaries
  // parameters set by ELM
  ATS::requireEvaluatorAtNext(poro_key_, Amanzi::Tags::NEXT, *S_, poro_key_)
    .SetMesh(mesh_subsurf_)->SetComponent("cell", AmanziMesh::CELL, 1);

  // potential fluxes (ELM -> ATS)
  ATS::requireEvaluatorAtNext(gross_water_source_key_, Amanzi::Tags::NEXT, *S_, gross_water_source_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(pot_evap_key_, Amanzi::Tags::NEXT, *S_, pot_evap_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(pot_trans_key_, Amanzi::Tags::NEXT, *S_, pot_trans_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", AmanziMesh::CELL, 1);

  Coordinator::parseParameterList();
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("2a: parseParameterList");
  }
}


void
ELM_ATSDriver::setup()
{
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning setup stage..." << std::endl
               << std::flush;
  }

  // and now require our output variables
  // ATS::requireEvaluatorAtNext(surf_cv_key_, Amanzi::Tags::NEXT, *S_)
  //   .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  // ATS::requireEvaluatorAtNext(cv_key_, Amanzi::Tags::NEXT, *S_)
  //   .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // ATS::requireEvaluatorAtNext(surf_mol_dens_key_, Amanzi::Tags::NEXT, *S_)
  //   .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  // ATS::requireEvaluatorAtNext(surf_mass_dens_key_, Amanzi::Tags::NEXT, *S_)
  //   .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  // ATS::requireEvaluatorAtNext(subsurf_mol_dens_key_, Amanzi::Tags::NEXT, *S_)
  //   .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  // ATS::requireEvaluatorAtNext(subsurf_mass_dens_key_, Amanzi::Tags::NEXT, *S_)
  //   .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // output variables
  ATS::requireEvaluatorAtNext(pd_key_, Amanzi::Tags::NEXT, *S_)
   .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(sat_key_, Amanzi::Tags::NEXT, *S_)
   .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(pres_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  ATS::requireEvaluatorAtNext(evap_key_, Amanzi::Tags::NEXT, *S_)
   .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(col_trans_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(col_baseflow_key_, Amanzi::Tags::NEXT, *S_, col_baseflow_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", AmanziMesh::CELL, 1);
  ATS::requireEvaluatorAtNext(col_runoff_key_, Amanzi::Tags::NEXT, *S_, col_runoff_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", AmanziMesh::CELL, 1);


  // This must be last always -- it allocates memory calling State::setup, so
  // all other setup must be done.
  Coordinator::setup();
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("2b: setup");
  }
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

  // compute dzs on one column only -- presumed terrain following!
  info.dzs.resize(ncells_per_col_);
  const auto& col_cells = mesh_subsurf_->columns.getCells(0);
  double surface_area = mesh_surf_->getCellVolume(0);
  for (int i = 0; i != ncells_per_col_; ++i) {
    info.dzs[i] = mesh_subsurf_->getCellVolume(i) / surface_area;
  }

  // surface area
  info.areas.resize(ncolumns_);
  for (int i = 0; i != ncolumns_; ++i) {
    info.areas[i] = mesh_surf_->getCellVolume(i);
  }
  return info;
}


void ELM_ATSDriver::initialize()
{
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning initialize stage..." << std::endl
               << std::flush;
  }

  // some parameters have been initialized in ELM
  S_->GetRecordW(poro_key_, Tags::NEXT, poro_key_).set_initialized();

  // set as zero and tag as initialized
  initValue_(gross_water_source_key_);
  initValue_(pot_trans_key_);
  initValue_(pot_evap_key_);

  // no current evaluators for these, treat as primary
  initValue_(col_baseflow_key_);
  initValue_(col_runoff_key_);

  // initialize ATS data, commit initial conditions
  Coordinator::initialize();

  // visualization at IC -- TODO remove this or place behind flag
  visualize();
  checkpoint();

  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("3: initialize");
  }
}



void ELM_ATSDriver::advance(double dt, bool do_chkp, bool do_vis)
{
  {
    Teuchos::TimeMonitor timer(*timers_.at("4: solve"));

    AMANZI_ASSERT(std::abs(S_->get_time(Amanzi::Tags::NEXT) - S_->get_time(Amanzi::Tags::CURRENT)) <
                  1.e-4);

    if (vo_->os_OK(Teuchos::VERB_LOW)) {
      *vo_->os()
        << "================================================================================"
        << std::endl
        << std::endl
        << "Cycle = " << S_->get_cycle() << ",  Time [days] = " << std::setprecision(16)
        << S_->get_time() / (60 * 60 * 24) << ",  dt [days] = " << std::setprecision(16)
        << dt / (60 * 60 * 24) << std::endl
        << "--------------------------------------------------------------------------------"
        << std::endl;
    }
    S_->Assign<double>("dt", Amanzi::Tags::DEFAULT, "dt", dt);
    S_->advance_time(Amanzi::Tags::NEXT, dt);

    // solve model for a single timestep
    bool fail = Coordinator::advance();
    if (fail) {
      Errors::Message msg("ELM_ATSDriver: advance(dt) failed.  Make ATS subcycle for proper ELM use.");
      Exceptions::amanzi_throw(msg);
    }

    S_->set_time(Amanzi::Tags::CURRENT, S_->get_time(Amanzi::Tags::NEXT));
    S_->advance_cycle();

    // vis/checkpoint if EITHER ATS or ELM request it
    if (do_vis && !visualize()) visualize(true);
    if (do_chkp && !checkpoint()) checkpoint(true);
    observe();
  }
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("4: solve");
  }
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
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning finalize stage..." << std::endl
               << std::flush;
  }
  Coordinator::finalize();
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("5: finalize");
  }
}


void ELM_ATSDriver::initValue_(const Key& key, double value)
{
  auto& vec = S_->GetW<CompositeVector>(key, Tags::NEXT, key);
  vec.PutScalar(value);
  S_->GetRecordW(key, Tags::NEXT, key).set_initialized();
}


void ELM_ATSDriver::setScalar(const ELM::VarID& scalar_id, double in)
{
  KeyTag scalar_key = key_map_.at(scalar_id);
  S_->GetW<double>(scalar_key.first, scalar_key.second, scalar_key.first) = in;
}

double ELM_ATSDriver::getScalar(const ELM::VarID& scalar_id)
{
  KeyTag scalar_key = key_map_.at(scalar_id);
  return S_->Get<double>(scalar_key.first, scalar_key.second);
}


double const *
ELM_ATSDriver::getFieldPtr(const ELM::VarID& var_id)
{
  Amanzi::KeyTag var_key = key_map_.at(var_id);
  return &(*S_->Get<CompositeVector>(var_key.first, var_key.second).ViewComponent("cell", false))[0][0];
}

double *
ELM_ATSDriver::getFieldPtrW(const ELM::VarID& var_id)
{
  Amanzi::KeyTag var_key = key_map_.at(var_id);
  Amanzi::changedEvaluatorPrimary(var_key.first, var_key.second, *S_);

  Amanzi::Key owner = S_->GetRecord(var_key.first, var_key.second).owner();
  return &(*S_->GetW<CompositeVector>(var_key.first, var_key.second, owner).ViewComponent("cell", false))[0][0];
}


} // namespace ATS
