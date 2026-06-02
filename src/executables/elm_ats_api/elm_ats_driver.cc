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
#include "time_advancer.hh"
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

  // Set global logfile before constructing the driver so that the Driver base
  // class constructor (which creates VerboseObjects) already sees the correct
  // logfile path rather than falling back to stdout.
  VerboseObject::global_default_level = Teuchos::VERB_LOW;
  VerboseObject::global_logfile = logfile_filename;

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
  : Driver(plist, wallclock_timer, teuchos_comm, comm),
    npfts_(npfts),
    ncolumns(-1),
    ncells_per_col_(-1),
    elm_plist_(Teuchos::sublist(Teuchos::sublist(plist, "cycle driver"), "ELM driver")),
    elm_cycle_(0)
{

  if (plist_->isSublist("visualization")) {
    auto vis_list = Teuchos::sublist(plist_, "visualization");
    for (auto& entry : *vis_list) {
      std::string domain_name = entry.first;
      if (S_->HasMesh(domain_name)) {
        auto mesh_p = S_->GetMesh(domain_name);
        auto sublist_p = Teuchos::sublist(vis_list, domain_name);
        if (!sublist_p->isParameter("file name base")) {
          if (domain_name.empty() || domain_name == "domain") {
            sublist_p->set<std::string>("file name base", std::string("ats_vis"));
          } else {
            sublist_p->set<std::string>("file name base", std::string("ats_vis_") + domain_name);
          }
        }
        if (S_->HasMesh(domain_name + "_3d") && sublist_p->get<bool>("visualize on 3D mesh", true))
          mesh_p = S_->GetMesh(domain_name + "_3d");
        auto vis = Teuchos::rcp(new Amanzi::Visualization(*sublist_p));
        vis->set_name(domain_name);
        vis->set_mesh(mesh_p);
        vis->CreateFiles(false);
        visualization_.push_back(vis);
      }
    }
  }
}


bool
ELM_ATSDriver::visualize_(bool force)
{
  bool dump = force;
  double time = S_->get_time();
  int cycle = elm_cycle_;

  if (!dump) {
    for (const auto& vis : visualization_) dump |= vis->DumpRequested(cycle, time);
  }
  if (dump) pk_->CalculateDiagnostics(Tags::NEXT);
  for (const auto& vis : visualization_) {
    if (force || vis->DumpRequested(cycle, time)) WriteVis(*vis, *S_);
  }
  return dump;
}



Key
ELM_ATSDriver::setupIntegratedFlux_(const Key& flux_key, ELM::VarID varid)
{
  Key integrated_key = flux_key + "_step_integrated";

  // inject EvaluatorTimeAccumulated into state/evaluators
  auto& ep = S_->GetEvaluatorList(integrated_key);
  ep.set("evaluator type", "time accumulated");
  ep.set("tag", Amanzi::Tags::NEXT.get());
  ep.set("accumulated key", flux_key);
  ep.set("accumulated tag", Amanzi::Tags::NEXT.get());
  ep.set("accumulation type", "integral");

  // append to "time accumulators" in cycle driver plist
  auto& cd_list = plist_->sublist("cycle driver");
  auto acc_list = cd_list.get<Teuchos::Array<std::string>>(
    "time accumulators", Teuchos::Array<std::string>());
  acc_list.push_back(integrated_key + "@");
  cd_list.set("time accumulators", acc_list);

  // update key_map_ to point at the integrated quantity
  key_map_[varid] = { integrated_key, Amanzi::Tags::NEXT };

  return integrated_key;
}


void
ELM_ATSDriver::parseParameterList()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning parseParameterList stage..." << std::endl
               << std::flush;
  }
  // parse my parameter list
  // domain names
  domain_subsurf_ = Keys::readDomain(*elm_plist_, "domain", "domain");
  domain_surf_ = Keys::readDomainHint(*elm_plist_, domain_subsurf_, "subsurface", "surface");

   // meshes
  mesh_subsurf_ = S_->GetMesh(domain_subsurf_);
  mesh_surf_ = S_->GetMesh(domain_surf_);

  key_map_[ELM::VarID::TIME] = {"time", Tags::NEXT};

  // parameters
  base_poro_key_ = Keys::readKey(*elm_plist_, domain_subsurf_, "base porosity", "base_porosity");
  key_map_[ELM::VarID::BASE_POROSITY] = { base_poro_key_, Tags::NEXT };

  // potential sources
  gross_water_source_key_ = Keys::readKey(*elm_plist_, domain_surf_, "gross water source", "gross_water_source");
  key_map_[ELM::VarID::GROSS_SURFACE_WATER_SOURCE] = { gross_water_source_key_, Tags::NEXT };
  pot_evap_key_ = Keys::readKey(*elm_plist_, domain_surf_, "potential evaporation", "potential_evaporation");
  key_map_[ELM::VarID::POTENTIAL_EVAPORATION] = { pot_evap_key_, Tags::NEXT };
  pot_trans_key_ = Keys::readKey(*elm_plist_, domain_surf_, "potential transpiration mps", "potential_transpiration_mps");
  key_map_[ELM::VarID::POTENTIAL_TRANSPIRATION] = { pot_trans_key_, Tags::NEXT };

  // water state
  pres_key_ = Keys::readKey(*elm_plist_, domain_subsurf_, "pressure", "pressure");
  key_map_[ELM::VarID::PRESSURE] = { pres_key_, Tags::NEXT };
  sat_key_ = Keys::readKey(*elm_plist_, domain_subsurf_, "saturation liquid", "saturation_liquid");
  key_map_[ELM::VarID::SATURATION_LIQUID] = { sat_key_, Tags::NEXT };
  zwt_key_ = Keys::readKey(*elm_plist_, domain_surf_, "water table depth", "water_table_depth");
  key_map_[ELM::VarID::DEPTH_TO_WATER_TABLE] = { zwt_key_, Tags::NEXT };
  pd_key_ = Keys::readKey(*elm_plist_, domain_surf_, "ponded depth", "ponded_depth");
  key_map_[ELM::VarID::PONDED_DEPTH] = { pd_key_, Tags::NEXT };

  surf_wc_key_ = Keys::readKey(*elm_plist_, domain_surf_, "surface water content", "water_content");
  key_map_[ELM::VarID::SURFACE_WATER_CONTENT_OLD] = { surf_wc_key_, Tags::CURRENT };
  wc_key_ = Keys::readKey(*elm_plist_, domain_subsurf_, "water content", "water_content");
  key_map_[ELM::VarID::WATER_CONTENT] = { wc_key_, Tags::NEXT };
  key_map_[ELM::VarID::WATER_CONTENT_OLD] = { wc_key_, Tags::CURRENT };

  lat_key_ = Keys::readKey(*elm_plist_, domain_surf_, "latitude", "latitude");
  lon_key_ = Keys::readKey(*elm_plist_, domain_surf_, "longitude", "longitude");

  // actual water fluxes — each is wrapped in an EvaluatorTimeAccumulated so
  // that the value returned to ELM is integrated over the outer timestep
  evap_key_ = setupIntegratedFlux_(
    Keys::readKey(*elm_plist_, domain_surf_, "evaporation", "evaporation"),
    ELM::VarID::EVAPORATION);
  col_trans_key_ = setupIntegratedFlux_(
    Keys::readKey(*elm_plist_, domain_surf_, "surface transpiration", "transpiration"),
    ELM::VarID::TRANSPIRATION);
  col_baseflow_key_ = setupIntegratedFlux_(
    Keys::readKey(*elm_plist_, domain_surf_, "baseflow generation", "baseflow_mps"),
    ELM::VarID::BASEFLOW);
  col_runoff_key_ = setupIntegratedFlux_(
    Keys::readKey(*elm_plist_, domain_surf_, "runoff generation", "runoff_generation_mps"),
    ELM::VarID::RUNOFF);

  // keys for fields used to convert ELM units to ATS units
  surf_mol_dens_key_ = Keys::readKey(*elm_plist_, domain_surf_, "surface molar density", "molar_density_liquid");
  mol_dens_key_ = Keys::readKey(*elm_plist_, domain_subsurf_, "molar density", "molar_density_liquid");
  // surf_mass_dens_key_ = Keys::readKey(*elm_plist_, domain_surf_, "surface mass density", "mass_density_liquid");
  // subsurf_mass_dens_key_ = Keys::readKey(*elm_plist_, domain_subsurf_, "mass density", "mass_density_liquid");

  // cell vol keys
  surf_cv_key_ = Keys::getKey(domain_surf_, "cell_volume");
  cv_key_ = Keys::getKey(domain_subsurf_, "cell_volume");

  // -- mesh info
  ncolumns = mesh_surf_->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  auto col_zero = mesh_subsurf_->columns.getCells(0);
  ncells_per_col_ = col_zero.size();

  // require my primary variables
  // parameters set by ELM
  requireEvaluatorAtNext(base_poro_key_, Amanzi::Tags::NEXT, *S_, base_poro_key_);

  // potential fluxes (ELM -> ATS)
  requireEvaluatorAtNext(gross_water_source_key_, Amanzi::Tags::NEXT, *S_, gross_water_source_key_);
  requireEvaluatorAtNext(pot_evap_key_, Amanzi::Tags::NEXT, *S_, pot_evap_key_);
  requireEvaluatorAtNext(pot_trans_key_, Amanzi::Tags::NEXT, *S_, pot_trans_key_);

  Driver::parseParameterList();

  // construct TimeAdvancer here so time_advancer_->setup() is called in Driver::setup()
  // before State::Setup() allocates memory.
  // ELM controls checkpointing, so no checkpoint sublist.
  time_advancer_ = Teuchos::rcp(new ATS::TimeAdvancer(
            Teuchos::sublist(plist_, "cycle driver"),
            S_, pk_, tsm_,
            Amanzi::Tags::CURRENT, Amanzi::Tags::NEXT,
            vo_, wallclock_timer_));

  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("2a: parseParameterList");
  }
}


void
ELM_ATSDriver::setup()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning setup stage..." << std::endl
               << std::flush;
  }

  // variables needed for unit changes
  requireEvaluatorAtNext(surf_cv_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireEvaluatorAtNext(cv_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  requireEvaluatorAtNext(surf_mol_dens_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireEvaluatorAtNext(mol_dens_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // ELM --> ATS variables
  // -- parameters set by ELM
  requireEvaluatorAtNext(base_poro_key_, Amanzi::Tags::NEXT, *S_, base_poro_key_)
    .SetMesh(mesh_subsurf_)->SetComponent("cell", AmanziMesh::CELL, 1);

  // -- potential fluxes
  requireEvaluatorAtNext(gross_water_source_key_, Amanzi::Tags::NEXT, *S_, gross_water_source_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", AmanziMesh::CELL, 1);
  requireEvaluatorAtNext(pot_evap_key_, Amanzi::Tags::NEXT, *S_, pot_evap_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", AmanziMesh::CELL, 1);
  requireEvaluatorAtNext(pot_trans_key_, Amanzi::Tags::NEXT, *S_, pot_trans_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", AmanziMesh::CELL, 1);

  // -- water contents at the OLD time -- note, these are already primary and
  // -- owned by overland flow.
  requireEvaluatorAtCurrent(wc_key_, Amanzi::Tags::CURRENT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireEvaluatorAtCurrent(surf_wc_key_, Amanzi::Tags::CURRENT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // save a copy for debugging!
  requireEvaluatorAtCurrent("elm_old_water_content", Amanzi::Tags::NEXT, *S_, "elm_old_water_content")
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireEvaluatorAtCurrent("surface-elm_old_water_content", Amanzi::Tags::NEXT, *S_, "surface-elm_old_water_content")
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // ATS --> ELM variables
  requireEvaluatorAtNext(pd_key_, Amanzi::Tags::NEXT, *S_)
   .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireEvaluatorAtNext(zwt_key_, Amanzi::Tags::NEXT, *S_)
   .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireEvaluatorAtNext(sat_key_, Amanzi::Tags::NEXT, *S_)
   .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireEvaluatorAtNext(pres_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  requireEvaluatorAtNext(evap_key_, Amanzi::Tags::NEXT, *S_)
   .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireEvaluatorAtNext(col_trans_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireEvaluatorAtNext(col_baseflow_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireEvaluatorAtNext(col_runoff_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // mesh info -- optional!
  if (S_->HasEvaluatorList(lat_key_)) {
    requireEvaluatorAtNext(lat_key_, Amanzi::Tags::NEXT, *S_)
      .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
    requireEvaluatorAtNext(lon_key_, Amanzi::Tags::NEXT, *S_)
      .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  }

  // This must be last always -- it allocates memory calling State::setup, so
  // all other setup must be done.
  Driver::setup();
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
  info.ncols_local = ncolumns;
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
  info.areas.resize(ncolumns);
  for (int i = 0; i != ncolumns; ++i) {
    info.areas[i] = mesh_surf_->getCellVolume(i);
  }

  // lat/lon - optional
  info.latitudes.resize(ncolumns, -1.);
  info.longitudes.resize(ncolumns, -1.);

  if (S_->HasEvaluator(lat_key_, Tags::NEXT)) {
    S_->GetEvaluator(lat_key_, Tags::NEXT).Update(*S_, "ELM ATS driver");
    S_->GetEvaluator(lon_key_, Tags::NEXT).Update(*S_, "ELM ATS driver");
    const Epetra_MultiVector& lat = *S_->Get<CompositeVector>(lat_key_, Tags::NEXT)
      .ViewComponent("cell", false);
    const Epetra_MultiVector& lon = *S_->Get<CompositeVector>(lon_key_, Tags::NEXT)
      .ViewComponent("cell", false);

    for (int i = 0; i != ncolumns; ++i) {
      info.latitudes[i] = lat[0][i];
      info.longitudes[i] = lon[0][i];
    }
  }
  return info;
}


void ELM_ATSDriver::initialize()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning initialize stage..." << std::endl
               << std::flush;
  }

  // some parameters have been initialized in ELM
  S_->GetRecordW(base_poro_key_, Tags::NEXT, base_poro_key_).set_initialized();

  // for debugging
  initValue_("elm_old_water_content");
  initValue_("surface-elm_old_water_content");

  // set as zero and tag as initialized
  initValue_(gross_water_source_key_);
  initValue_(pot_trans_key_);
  initValue_(pot_evap_key_);

  // no current evaluators for these, treat as primary
  initValue_(col_baseflow_key_);
  initValue_(col_runoff_key_);

  // initialize ATS data, commit initial conditions; TimeAdvancer::initialize()
  // backs up mesh coordinates and dumps IC vis/obs (no IC checkpoint for ELM).
  Driver::initialize();

  // initialize has to convert IC for water content into ELM units
  // convert ATS water content units [mol] to ELM units [m^3 water per m^3 volume] before returning control to ELM
  {
    auto& wc = *S_->GetW<CompositeVector>(wc_key_, Tags::NEXT,
            S_->GetRecord(wc_key_, Tags::NEXT).owner()).ViewComponent("cell", false);
    const auto& cv = *S_->Get<CompositeVector>(cv_key_, Tags::NEXT).ViewComponent("cell", false);
    const auto& n_liq = *S_->Get<CompositeVector>(mol_dens_key_, Tags::NEXT).ViewComponent("cell", false);
    for (int c = 0; c != wc.MyLength(); ++c) {
      wc[0][c] = wc[0][c] / cv[0][c] / n_liq[0][c];
    }
  }

  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("3: initialize");
  }
}



void ELM_ATSDriver::advance(double dt, bool force_chkp, bool force_vis)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  {
    Teuchos::TimeMonitor timer(*timers_.at("4: solve"));

    AMANZI_ASSERT(std::abs(S_->get_time(Amanzi::Tags::NEXT) - S_->get_time(Amanzi::Tags::CURRENT)) <
                  1.e-4);

    double t_now = S_->get_time(Amanzi::Tags::CURRENT);

    if (vo_->os_OK(Teuchos::VERB_LOW)) {
      *vo_->os()
        << "================================================================================"
        << std::endl
        << std::endl
        << "ELM Cycle = " << elm_cycle_ << ",  Time [days] = " << std::setprecision(16)
        << t_now / (60 * 60 * 24) << ",  dt [days] = " << std::setprecision(16)
        << dt / (60 * 60 * 24) << std::endl
        << "================================================================================"
        << std::endl;
    }

    // convert ELM water content units [m^3 water per m^3 volume] to ATS units [mol] before solving.
    //
    // NOTE: these were just marked as changed in a call to getFieldPtrW so no
    // need to mark them as changed again.
    // -- subsurface wc in volumetric water content --> mols
    {
      auto& wc = *S_->GetW<CompositeVector>(wc_key_, Tags::CURRENT,
              S_->GetRecord(wc_key_, Tags::CURRENT).owner()).ViewComponent("cell", false);
      const auto& cv = *S_->Get<CompositeVector>(cv_key_, Tags::NEXT).ViewComponent("cell", false);
      const auto& n_liq = *S_->Get<CompositeVector>(mol_dens_key_, Tags::NEXT).ViewComponent("cell", false);
      for (int c = 0; c != wc.MyLength(); ++c) {
        wc[0][c] = wc[0][c] * cv[0][c] * n_liq[0][c];
      }
    }

    // -- surface wc in m ponded depth --> mols
    {
      auto& surf_wc = *S_->GetW<CompositeVector>(surf_wc_key_, Tags::CURRENT,
              S_->GetRecord(surf_wc_key_, Tags::CURRENT).owner()).ViewComponent("cell", false);
      const auto& surf_cv = *S_->Get<CompositeVector>(surf_cv_key_, Tags::NEXT).ViewComponent("cell", false);
      const auto& surf_n_liq = *S_->Get<CompositeVector>(surf_mol_dens_key_, Tags::NEXT).ViewComponent("cell", false);
      for (int c = 0; c != surf_wc.MyLength(); ++c) {
        surf_wc[0][c] = surf_wc[0][c] * surf_cv[0][c] * surf_n_liq[0][c];
      }
    }

    // -- save the old values in a new vector so they show up in vis for debugging
    {
      const auto& wc = *S_->Get<CompositeVector>(wc_key_, Tags::CURRENT).ViewComponent("cell", false);
      auto& wc_old = *S_->GetW<CompositeVector>("elm_old_water_content", Tags::NEXT, "elm_old_water_content")
        .ViewComponent("cell", false);
      wc_old = wc;

      const auto& surf_wc = *S_->Get<CompositeVector>(surf_wc_key_, Tags::CURRENT).ViewComponent("cell", false);
      auto& surf_wc_old = *S_->GetW<CompositeVector>("surface-elm_old_water_content", Tags::NEXT, "surface-elm_old_water_content")
        .ViewComponent("cell", false);
      surf_wc_old = surf_wc;
    }

    // advance from t_now to t_now+dt; TimeAdvancer handles subcycling internally
    bool fail = time_advancer_->advance(t_now, t_now + dt);
    if (fail) {
      Errors::Message msg("ELM_ATSDriver: advance(dt) failed.  Make ATS subcycle for proper ELM use.");
      Exceptions::amanzi_throw(msg);
    }

    // potentially visualize
    visualize_();

    // convert ATS water content units [mol] to ELM units [m^3 water per m^3 volume] before returning control to ELM
    {
      auto& wc = *S_->GetW<CompositeVector>(wc_key_, Tags::NEXT,
              S_->GetRecord(wc_key_, Tags::NEXT).owner()).ViewComponent("cell", false);
      const auto& cv = *S_->Get<CompositeVector>(cv_key_, Tags::NEXT).ViewComponent("cell", false);
      const auto& n_liq = *S_->Get<CompositeVector>(mol_dens_key_, Tags::NEXT).ViewComponent("cell", false);
      for (int c = 0; c != wc.MyLength(); ++c) {
        wc[0][c] = wc[0][c] / cv[0][c] / n_liq[0][c];
      }
    }
  }
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("4: solve");
  }

  ++elm_cycle_;
}


// simulates external timeloop with dt coming from calling model
void ELM_ATSDriver::advanceTest()
{
  while (S_->get_time() < t1_) {
    // use dt from ATS for testing
    double dt = pk_->get_dt();
    // call main method
    advance(dt, false, false);
  }
}

void ELM_ATSDriver::finalize()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning finalize stage..." << std::endl
               << std::flush;
  }
  Driver::finalize(false);
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
  if (S_->HasEvaluator(var_key.first, var_key.second)) {
    S_->GetEvaluator(var_key.first, var_key.second).Update(*S_, std::string("elm_ats_driver_on_"+domain_subsurf_));
  }

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
