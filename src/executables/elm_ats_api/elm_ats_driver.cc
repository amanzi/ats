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
#include "pk_helpers.hh"
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

  // keys for fields in mesh info
  lat_key_ = Keys::readKey(*plist_, domain_surf_, "latitude", "latitude");
  lon_key_ = Keys::readKey(*plist_, domain_surf_, "longitude", "longitude");
  elev_key_ = Keys::readKey(*plist_, domain_surf_, "elevation", "elevation");

  // soil parameters/properties
  base_poro_key_ = Keys::readKey(*plist_, domain_subsurf_, "base porosity", "base_porosity");
  poro_key_ = Keys::readKey(*plist_, domain_subsurf_, "porosity", "porosity");
  perm_key_ = Keys::readKey(*plist_, domain_subsurf_, "permeability", "permeability");
  ch_b_key_ = Keys::readKey(*plist_, domain_subsurf_, "Clapp and Hornberger b", "clapp_horn_b");
  ch_smpsat_key_ = Keys::readKey(*plist_, domain_subsurf_, "Clapp and Hornberger soil mafic potential at saturation", "clapp_horn_smpsat");
  ch_sr_key_ = Keys::readKey(*plist_, domain_subsurf_, "Clapp and Hornberger residual saturation", "clapp_horn_sr");

  // potential sources
  root_frac_key_ = Keys::readKey(*plist_, domain_subsurf_, "rooting depth fraction", "rooting_depth_fraction");
  pot_infilt_key_ = Keys::readKey(*plist_, domain_surf_, "potential infiltration mps", "potential_infiltration_mps"); // inputs onto surface (rain, snowmelt)
  pot_evap_key_ = Keys::readKey(*plist_, domain_surf_, "potential evaporation mps", "potential_evaporation_mps");
  pot_trans_key_ = Keys::readKey(*plist_, domain_surf_, "potential transpiration mps", "potential_transpiration_mps");

  // water state
  pd_key_ = Keys::readKey(*plist_, domain_surf_, "ponded depth", "ponded_depth");
  wtd_key_ = Keys::readKey(*plist_, domain_surf_, "water table depth", "water_table_depth");
  pres_key_ = Keys::readKey(*plist_, domain_subsurf_, "pressure", "pressure");
  wc_key_ = Keys::readKey(*plist_, domain_subsurf_, "conserved", "water_content");
  pc_key_ = Keys::readKey(*plist_, domain_subsurf_, "capillary_pressure_gas_liq", "capillary_pressure_gas_liq");
  sat_key_ = Keys::readKey(*plist_, domain_subsurf_, "saturation", "saturation_liquid");

  // water fluxes
  infilt_key_ = Keys::readKey(*plist_, domain_surf_, "surface-subsurface flux", "surface_subsurface_flux");
  evap_key_ = Keys::readKey(*plist_, domain_surf_, "evaporation", "evaporation");
  trans_key_ = Keys::readKey(*plist_, domain_subsurf_, "transpiration", "transpiration");

  // keys for fields used to convert ELM units to ATS units
  surf_mol_dens_key_ = Keys::readKey(*plist_, domain_surf_, "surface molar density", "molar_density_liquid");
  surf_mass_dens_key_ = Keys::readKey(*plist_, domain_surf_, "surface mass density", "mass_density_liquid");
  subsurf_mol_dens_key_ = Keys::readKey(*plist_, domain_subsurf_, "molar density", "molar_density_liquid");
  subsurf_mass_dens_key_ = Keys::readKey(*plist_, domain_subsurf_, "mass density", "mass_density_liquid");

  // need to put into observations or explicitly update if value other than 0.0 is desired
  total_trans_key_ = Keys::readKey(*plist_, domain_surf_, "total transpiration", "total_transpiration");

  // cell vol keys
  surf_cv_key_ = Keys::getKey(domain_surf_, "cell_volume");
  cv_key_ = Keys::getKey(domain_subsurf_, "cell_volume");

  // currently unused keys
  //pres_atm_key_ = Keys::readKey(*plist_, domain_surf_, "atmospheric_pressure", "atmospheric_pressure"); // hardwired as 101325
  //sat_gas_key_ = Keys::readKey(*plist_, domain_subsurf_, "saturation gas", "saturation_gas"); // probably never needed
  //sat_ice_key_ = Keys::readKey(*plist_, domain_subsurf_, "saturation ice", "saturation_ice"); // not until energy

  // -- check that number of surface cells = number of columns
  ncolumns_ = mesh_surf_->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  AMANZI_ASSERT(ncolumns_ == mesh_subsurf_->columns.num_columns_owned);

  // -- get num cells per column - include consistency check later need to know
  //    if coupling zone is the entire subsurface mesh (as currently coded) or
  //    a portion of the total depth specified by # of cells into the
  //    subsurface
  auto col_zero = mesh_subsurf_->columns.getCells(0);
  ncells_per_col_ = col_zero.size();
  for (int col=0; col!=ncolumns_; ++col)
    AMANZI_ASSERT(mesh_subsurf_->columns.getCells(col).size() == ncells_per_col_);
}


void
ELM_ATSDriver::setup()
{
  // potential fluxes (ELM -> ATS)
  requireAtNext(pot_infilt_key_, Amanzi::Tags::NEXT, *S_, pot_infilt_key_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(pot_evap_key_, Amanzi::Tags::NEXT, *S_, pot_evap_key_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(pot_trans_key_, Amanzi::Tags::NEXT, *S_, pot_trans_key_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // subsurface properties
  requireAtNext(base_poro_key_, Amanzi::Tags::NEXT, *S_, base_poro_key_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(perm_key_, Amanzi::Tags::NEXT, *S_, perm_key_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // dynamic subsurface properties
  requireAtNext(root_frac_key_, Amanzi::Tags::NEXT, *S_, root_frac_key_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(poro_key_, Amanzi::Tags::NEXT, *S_) // use base_porosity from elm and ATS model for compressibility
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // Clapp and Hornberger water retention params (ELM -> ATS)
  requireAtNext(ch_b_key_, Amanzi::Tags::NEXT, *S_, ch_b_key_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(ch_smpsat_key_, Amanzi::Tags::NEXT, *S_, ch_smpsat_key_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(ch_sr_key_, Amanzi::Tags::NEXT, *S_, ch_sr_key_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // per-column ATS water state (ATS -> ELM)
  requireAtNext(pd_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(wtd_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // per cell ATS water state
  requireAtNext(pc_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(sat_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // mesh data
  //  requireAtNext(lat_key_, Amanzi::Tags::NEXT, *S_)
  //    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  //  requireAtNext(lon_key_, Amanzi::Tags::NEXT, *S_)
  //    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  //  requireAtNext(elev_key_, Amanzi::Tags::NEXT, *S_)
  //    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(surf_cv_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(cv_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  // may not be necessary? any PK that utilizes this should already have density
  requireAtNext(surf_mol_dens_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(surf_mass_dens_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(subsurf_mol_dens_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(subsurf_mass_dens_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(total_trans_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(wc_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(evap_key_, Amanzi::Tags::NEXT, *S_)
   .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(trans_key_, Amanzi::Tags::NEXT, *S_)
   .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(infilt_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(pres_key_, Amanzi::Tags::NEXT, *S_, "flow")
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  Coordinator::setup();
}



//
// dz & depth are currently ignored -- presumed identical between ATS & ELM
//
void ELM_ATSDriver::get_mesh_info(int& ncols_local,
        int& ncols_global,
        double * const lat,
        double * const lon,
        double * const elev,
        double * const surface_area,
        int * const pft,
        int& nlevgrnd,
        double * const depth)
{
  ncols_local = ncolumns_;
  ncols_global = mesh_surf_->getMap(AmanziMesh::Entity_kind::CELL, false).NumGlobalElements();

  // copyFromSurf_(elev, elev_key_);
  // copyFromSurf_(surface_area, surf_cv_key_);
  // copyFromSurf_(lat, lat_key_);
  // copyFromSurf_(lon, lon_key_);

  // NOTE: figure out how to map from ATS LC types to ELM PFT... --etc
  // hard-coded veg type for now...
  for (int i=0; i!=ncolumns_; ++i) pft[i] = 1;

  nlevgrnd = ncells_per_col_;
  const auto& cells_in_col = mesh_subsurf_->columns.getCells(0);
  double ldepth = 0.;
  double top_face_z = mesh_subsurf_->getFaceCentroid(mesh_subsurf_->columns.getFaces(0)[0])[2];
  for (int i=0; i!=ncells_per_col_; ++i) {
    double bottom_face_z = mesh_subsurf_->getFaceCentroid(mesh_subsurf_->columns.getFaces(0)[i+1])[2];
    double dz = top_face_z - bottom_face_z;
    ldepth += dz/2.;
    depth[i] = ldepth;
    top_face_z = bottom_face_z;
    ldepth += dz/2.;
  }

  // hard-coded Toledo OH for now...
  for (int i=0; i!=ncolumns_; ++i) {
    lat[i] = 41.65;
    lon[i] = -83.54;
  }
}


void ELM_ATSDriver::initialize(double t,
                               double const * const elm_water_content,
                               double const * const elm_pressure)
{
  // Assign start time to ATS
  t0_ = t;

  // initialize ATS data, commit initial conditions
  Coordinator::initialize();

  // initialize pressure field
  ELM_ATSDriver::init_pressure_from_wc_(elm_water_content);

  // set as zero and tag as initialized
  initZero_(root_frac_key_);
  initZero_(pot_infilt_key_);
  initZero_(pot_trans_key_);
  initZero_(pot_evap_key_);
  initZero_(infilt_key_);
  initZero_(trans_key_);
  initZero_(evap_key_);
  initZero_(total_trans_key_);

  // visualization at IC -- TODO remove this or place behind flag
  visualize();
  checkpoint();
}

// use incoming water content to initialize pressure field
void ELM_ATSDriver::init_pressure_from_wc_(double const * const elm_water_content)
{
  // gravity, atmospheric pressure, and liquid water density
  // hardwired for now
  const double g = 9.80665;
  const double p_atm = 101325.0;
  const double rho = 1000.0;

  // evaluators
  S_->GetEvaluator(subsurf_mass_dens_key_, Amanzi::Tags::NEXT).Update(*S_, subsurf_mass_dens_key_);
  const auto& mass_d = *S_->Get<CompositeVector>(subsurf_mass_dens_key_, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);
  S_->GetEvaluator(poro_key_, Amanzi::Tags::NEXT).Update(*S_, poro_key_);
  const auto& por = *S_->Get<CompositeVector>(poro_key_, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);
  S_->GetEvaluator(cv_key_, Amanzi::Tags::NEXT).Update(*S_, cv_key_);
  const auto& volume = *S_->Get<CompositeVector>(cv_key_, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);
  S_->GetEvaluator(surf_cv_key_, Amanzi::Tags::NEXT).Update(*S_, surf_cv_key_);
  const auto& area = *S_->Get<CompositeVector>(surf_cv_key_, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);

  // writable to pressure
  auto& pres = *S_->GetW<CompositeVector>(pres_key_, Amanzi::Tags::NEXT, "flow")
    .ViewComponent("cell", false);

  // WRM model
  auto& wrm_eval = S_->GetEvaluator(sat_key_, Amanzi::Tags::NEXT);
  auto wrm_ptr = dynamic_cast<Amanzi::Flow::WRMEvaluator*>(&wrm_eval);
  AMANZI_ASSERT(wrm_ptr != nullptr);
  auto wrms_ = wrm_ptr->get_WRMs();
  AMANZI_ASSERT(wrms_->second.size() == 1); // only supports one WRM for now
  Teuchos::RCP<Flow::WRM> wrm_ = wrms_->second[0];

  // initialize pressure field from ELM water content
  // per-column hydrostatic pressure in areas of continuous total saturation
  // unsaturated areas are considered to be in contact with atmosphere
  for (int i=0; i!=ncolumns_; ++i) {
    const auto& cells_of_col = mesh_subsurf_->columns.getCells(i);
    int top_sat_idx = -1;
    double sat_depth = 0.0;
    for (int j=0; j!=ncells_per_col_; ++j) {
      // convert ELM water content (kg/m2] to saturation of pore space (0 to 1) [-]
      // VWC  =  elm_wc  *  1/dz    *  1/porosity  *  1/mass density
      // [-]  =  [kg/m2] *  [m^-1]  *  [-]         *  [m3/kg]
      const double dz = volume[0][cells_of_col[j]] / area[0][i];
      const double factor = 1 / (dz * por[0][cells_of_col[j]] * mass_d[0][cells_of_col[j]]);
      const double satl = elm_water_content[j * ncolumns_ + i] * factor;
      if (satl < 1.0) {
        pres[0][cells_of_col[j]] = p_atm - wrm_->capillaryPressure(satl);
        top_sat_idx = -1;
      } else {
        if (top_sat_idx == -1) {
          top_sat_idx = j;
          sat_depth = 0.0;
        }
        sat_depth += dz;
        pres[0][cells_of_col[j]] = p_atm + rho * g * (sat_depth - dz/2);
      }
    }
  }

  // mark pressure as changed and update face values
  changedEvaluatorPrimary(pres_key_, Amanzi::Tags::NEXT, *S_);
  DeriveFaceValuesFromCellValues(S_->GetW<CompositeVector>(pres_key_, Amanzi::Tags::NEXT, "flow"));
  S_->GetRecordW(pres_key_, Amanzi::Tags::NEXT, "flow").set_initialized();

  // update saturation and water content
  S_->GetEvaluator(sat_key_, Amanzi::Tags::NEXT).Update(*S_, sat_key_);
  S_->GetEvaluator(wc_key_, Amanzi::Tags::NEXT).Update(*S_, wc_key_);
}


void ELM_ATSDriver::advance(double dt, bool do_chkp, bool do_vis)
{
  double dt_subcycle = dt;
  double t_end = S_->get_time() + dt_subcycle;

  bool fail{false};
  while (S_->get_time() < t_end && dt_subcycle > 0.0) {
    S_->Assign<double>("dt", Amanzi::Tags::DEFAULT, "dt", dt_subcycle);
    S_->advance_time(Amanzi::Tags::NEXT, dt_subcycle);

    // solve model for a single timestep
    fail = Coordinator::advance();

    if (fail) {
      // reset t_new
      S_->set_time(Amanzi::Tags::NEXT, S_->get_time(Amanzi::Tags::CURRENT));
    } else {
      S_->set_time(Amanzi::Tags::CURRENT, S_->get_time(Amanzi::Tags::NEXT));
      S_->advance_cycle();

      // make observations, vis, and checkpoints
      for (const auto& obs : observations_) obs->MakeObservations(S_.ptr());

      // vis/checkpoint if EITHER ATS or ELM request it
      //if (do_vis && !visualize()) visualize(true);
      //if (do_chkp && !checkpoint()) checkpoint(true);

      visualize(do_vis);
      checkpoint(do_chkp);
    }

    dt_subcycle = Coordinator::get_dt(fail);
  } // while

  if (fail) {
    // write one more vis for help debugging
    S_->advance_cycle(Amanzi::Tags::NEXT);
    visualize(true); // force vis

    // flush observations to make sure they are saved
    for (const auto& obs : observations_) obs->Flush();

    // dump a post_mortem checkpoint file for debugging
    checkpoint_->set_filebasename("post_mortem");
    checkpoint_->Write(*S_, Checkpoint::WriteType::POST_MORTEM);

    Errors::Message msg("ELM_ATSDriver: advance(dt) failed.");
    Exceptions::amanzi_throw(msg);
  }
} // advance()


// simulates external timeloop with dt coming from calling model
void ELM_ATSDriver::advance_test()
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


void ELM_ATSDriver::set_soil_hydrologic_parameters(double const * const base_porosity,
        double const * const hydraulic_conductivity,
        double const * const clapp_horn_b,
        double const * const clapp_horn_smpsat,
        double const * const clapp_horn_sr)
{
  // convert Ksat to perm via rho * g / visc using default rho and visc values.
  copyToSub_(hydraulic_conductivity, perm_key_);
  double factor = 8.9e-4 / (1000 * 9.80665);
  S_->GetW<CompositeVector>(perm_key_, Amanzi::Tags::NEXT, perm_key_).Scale(factor);
  S_->GetRecordW(perm_key_, Amanzi::Tags::NEXT, perm_key_).set_initialized();

  copyToSub_(base_porosity, base_poro_key_);
  copyToSub_(clapp_horn_b, ch_b_key_);
  copyToSub_(clapp_horn_smpsat, ch_smpsat_key_);
  copyToSub_(clapp_horn_sr, ch_sr_key_);

  S_->GetRecordW(base_poro_key_, Amanzi::Tags::NEXT, base_poro_key_).set_initialized();
  S_->GetRecordW(ch_b_key_, Amanzi::Tags::NEXT, ch_b_key_).set_initialized();
  S_->GetRecordW(ch_smpsat_key_, Amanzi::Tags::NEXT, ch_smpsat_key_).set_initialized();
  S_->GetRecordW(ch_sr_key_, Amanzi::Tags::NEXT, ch_sr_key_).set_initialized();
}


void ELM_ATSDriver::set_veg_parameters(double const * const mafic_potential_full_turgor,
        double const * const mafic_potential_wilt_point)
{
  // pass for now! FIXME --etc
}


void ELM_ATSDriver::set_soil_hydrologic_properties(double const * const effective_porosity)
{
  // this isn't well defined -- pass for now --etc
}


void ELM_ATSDriver::set_veg_properties(double const * const rooting_fraction)
{
  copyToSub_(rooting_fraction, root_frac_key_);
}


void
ELM_ATSDriver::set_potential_sources(double const * const elm_surface_input,
                                     double const * const elm_evaporation,
                                     double const * const elm_transpiration)
{
  // ELM's infiltration, evaporation, and transpiration are in units of mm s^-1
  // ATS units are in mol m^-2 s^-1

  // kg / m2 / s  * m3/kg * mol/m3 = mol m^-2 s^-1
  // or
  // mm / s       *  m/mm * mol/m3 = mol m^-2 s^-1
  auto& in = *S_->GetW<CompositeVector>(pot_infilt_key_, Amanzi::Tags::NEXT, pot_infilt_key_)
    .ViewComponent("cell", false);
  auto& ev = *S_->GetW<CompositeVector>(pot_evap_key_, Amanzi::Tags::NEXT, pot_evap_key_)
    .ViewComponent("cell", false);
  auto& tr = *S_->GetW<CompositeVector>(pot_trans_key_, Amanzi::Tags::NEXT, pot_trans_key_)
    .ViewComponent("cell", false);
  AMANZI_ASSERT(in.MyLength() == ncolumns_);
  AMANZI_ASSERT(ev.MyLength() == ncolumns_);
  AMANZI_ASSERT(tr.MyLength() == ncolumns_);

  S_->GetEvaluator(surf_mol_dens_key_, Amanzi::Tags::NEXT).Update(*S_, surf_mol_dens_key_);
  const auto& surf_d = *S_->Get<CompositeVector>(surf_mol_dens_key_, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);

  for (int i=0; i!=ncolumns_; ++i) {
    const double factor = 0.001 * surf_d[0][i];
    in[0][i] = elm_surface_input[i] * factor;
    ev[0][i] = elm_evaporation[i]  * factor;
    tr[0][i] = elm_transpiration[i]  * factor;
  }

  changedEvaluatorPrimary(pot_infilt_key_, Amanzi::Tags::NEXT, *S_);
  changedEvaluatorPrimary(pot_evap_key_, Amanzi::Tags::NEXT, *S_);
  changedEvaluatorPrimary(pot_trans_key_, Amanzi::Tags::NEXT, *S_);
}


void
ELM_ATSDriver::get_waterstate(double * const ponded_depth,
        double * const wt_depth,
        double * const soil_water_potential,
        double * const matric_potential,
        double * const sat_liq,
        double * const sat_ice) // remove?
{
  // convert saturation into [kg/m2] from [-]
  S_->GetEvaluator(sat_key_, Amanzi::Tags::NEXT).Update(*S_, sat_key_);
  const auto& satl = *S_->Get<CompositeVector>(sat_key_, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);

  S_->GetEvaluator(poro_key_, Amanzi::Tags::NEXT).Update(*S_, poro_key_);
  const auto& por = *S_->Get<CompositeVector>(poro_key_, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);

  S_->GetEvaluator(subsurf_mass_dens_key_, Amanzi::Tags::NEXT).Update(*S_, subsurf_mass_dens_key_);
  const auto& dens = *S_->Get<CompositeVector>(subsurf_mass_dens_key_, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);

  // TODO look into ELM effective porosity, ATS ice density, ice saturation
  for (int i=0; i!=ncolumns_; ++i) {
    const auto faces = mesh_subsurf_->columns.getFaces(i);
    const auto cells = mesh_subsurf_->columns.getCells(i);
    for (int j=0; j!=ncells_per_col_; ++j) {
      const double dz = mesh_subsurf_->getFaceCentroid(faces[j])[2] - mesh_subsurf_->getFaceCentroid(faces[j + 1])[2];
      sat_liq[j * ncolumns_ + i] = satl[0][cells[j]] * por[0][cells[j]] * dens[0][cells[j]] * dz;
    }
  }

  // Ponded depth
  // convert ATS m to mm
  S_->GetEvaluator(pd_key_, Amanzi::Tags::NEXT).Update(*S_, "ELM");
  const auto& pd = *S_->Get<CompositeVector>(pd_key_, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);
  for (int i=0; i!=ncolumns_; ++i)
    ponded_depth[i] = pd[0][i] * 1000.0;

  // water table depth
  S_->GetEvaluator(wtd_key_, Amanzi::Tags::NEXT).Update(*S_, "ELM");
  copyFromSurf_(wt_depth, wtd_key_);

  // Soil matric potential
  // convert ATS Pa to -mmH2O
//  int z_index = mesh_subsurf_->space_dimension() - 1;
//  const auto& gravity = S_->Get<AmanziGeometry::Point>("gravity", Amanzi::Tags::DEFAULT);
//  const double g_inv = 1.0 / gravity[z_index]; // should be -9.80665 m s^-2
//
//  S_->GetEvaluator(pres_key_, Amanzi::Tags::NEXT).Update(*S_, pres_key_);
//   const auto& pres = *S_->Get<CompositeVector>(pres_key_, Amanzi::Tags::NEXT)
//     .ViewComponent("cell", false);
//
//  S_->GetEvaluator(pc_key_, Amanzi::Tags::NEXT).Update(*S_, "ELM");
//  const auto& pc = *S_->Get<CompositeVector>(pc_key_, Amanzi::Tags::NEXT)
//    .ViewComponent("cell", false);
//
//  for (int i=0; i!=ncolumns_; ++i) {
//    const auto& cells = mesh_subsurf_->columns.getCells(i);
//    for (int j=0; j!=ncells_per_col_; ++j) {
//      matric_potential[j * ncolumns_ + i] = pc[0][cells[j]] * g_inv;
//      soil_water_potential[j * ncolumns_ + i] = 0.101325 - 1.e-6 * pres[0][cells[j]];
//    }
//  }
}


void
ELM_ATSDriver::get_water_fluxes(double * const surf_subsurf_flx,
        double * const evaporation,
        double * const transpiration,
        double * const net_subsurface_fluxes,
        double * const net_runon)
{
  // Convert and send ATS fluxes to ELM

  // Flux units:          ATS             ELM
  // surface-evaporation  [mol / m2 / s]  [mmH2o/s]
  // surface-infiltration [mol / m2 / s]  [mmH2o/s]
  // transpiration        [mol / m3 / s]  [mmH2o/s]

  // Surface fluxes
  S_->GetEvaluator(infilt_key_, Amanzi::Tags::NEXT).Update(*S_, infilt_key_);
  const auto& infilt = *S_->Get<CompositeVector>(infilt_key_, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);

  S_->GetEvaluator(evap_key_, Amanzi::Tags::NEXT).Update(*S_, evap_key_);
  const auto& evap = *S_->Get<CompositeVector>(evap_key_, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);

  S_->GetEvaluator(surf_mol_dens_key_, Amanzi::Tags::NEXT).Update(*S_, surf_mol_dens_key_);
  const auto& surfdens = *S_->Get<CompositeVector>(surf_mol_dens_key_, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);

  // convert mol/m2/s to mmH2O/s for ELM
    // mol/m2/s * m3/mol  = m/s * mm/m = mm/s
  for (int i=0; i!=ncolumns_; ++i) {
    const double mm_per_mol = 1000.0 / surfdens[0][i];
    surf_subsurf_flx[i] = infilt[0][i] * mm_per_mol;
    evaporation[i] = evap[0][i] * mm_per_mol;
  }

  // Subsurface flux
  S_->GetEvaluator(trans_key_, Amanzi::Tags::NEXT).Update(*S_, "ELM");
  const auto& trans = *S_->Get<CompositeVector>(trans_key_, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);

  S_->GetEvaluator(subsurf_mol_dens_key_, Amanzi::Tags::NEXT).Update(*S_, subsurf_mol_dens_key_);
  const auto& subsurfdens = *S_->Get<CompositeVector>(subsurf_mol_dens_key_, Amanzi::Tags::NEXT)
    .ViewComponent("cell", false);

  // convert mol/m3/s to mmH2O/s by integrating over dz - NO?
  // treat the same as surface fluxes?
  for (int i=0; i!=ncolumns_; ++i) {
    const auto& faces = mesh_subsurf_->columns.getFaces(i);
    const auto& cells = mesh_subsurf_->columns.getCells(i);
    for (int j=0; j!=ncells_per_col_; ++j) {
      double dz = mesh_subsurf_->getFaceCentroid(faces[j])[2] - mesh_subsurf_->getFaceCentroid(faces[j + 1])[2];
      AMANZI_ASSERT(dz > 0.);
      const double factor = dz * 1000.0 / subsurfdens[0][cells[j]];
      transpiration[j * ncolumns_ + i] = trans[0][cells[j]] * factor;
    }
  }
  // unclear how to implement net_subsurface_fluxes and net_runon... --ETC
}


void ELM_ATSDriver::initZero_(const Key& key)
{
  auto& vec = S_->GetW<CompositeVector>(key, Amanzi::Tags::NEXT, key);
  vec.PutScalar(0.);
  S_->GetRecordW(key, Amanzi::Tags::NEXT, key).set_initialized();
}


void ELM_ATSDriver::copyToSurf_(double const * const in, const Key& key, Key owner)
{
  if (owner.empty()) owner = key;
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
void ELM_ATSDriver::copyToSub_(double const * const in, const Key& key, Key owner)
{
  if (owner.empty()) owner = key;
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


} // namespace ATS
