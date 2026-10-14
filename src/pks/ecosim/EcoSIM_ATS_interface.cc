/*--------------------------------------------------------------------------
  ATS

  License: see $ATS_DIR/COPYRIGHT
  Author: Andrew Graus

  This is the main PK for the EcoSIM-ATS interface. This code is written
  following the example of Alquimia with some additional code from the
  SimpleBGC code for walking the columns.

  The idea is to take the basic code used by alquimia and repurpose it so
  that it works on a column by column basis instead of a cell by cell basis

  --------------------------------------------------------------------------*/

#include <algorithm>
#include <iostream>
#include <set>
#include <string>

// TPLs
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "errors.hh"
#include "exceptions.hh"
#include "Mesh.hh"

// include custom evaluators here
// #include "hydraulic_conductivity_evaluator.hh"

#include "PK_Helpers.hh"
#include "EcoSIM_ATS_interface.hh"

namespace Amanzi {
namespace EcoSIM {

EcoSIM::EcoSIM(Teuchos::ParameterList& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& global_list,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& solution):
  PK_Physical_Default(pk_tree, global_list, S, solution),
  PK(pk_tree, global_list, S, solution),
  ncells_per_col_(-1),
  saved_time_(0.0)
  {
    //grab the surface and subsurface domains
    domain_ = Keys::readDomain(*plist_, "domain", "domain");
    domain_surface_ = Keys::readDomainHint(*plist_, domain_, "subsurface", "surface");

    // transport
    mole_fraction_key_ = Keys::readKey(*plist_, domain_, "mole fraction", "mole_fraction");
    //mole_fraction components are accessed by mole_fraction[i][c] where i is the component and c is the cell

    //Flow
    porosity_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
    saturation_liquid_key_ = Keys::readKey(*plist_, domain_, "saturation liquid", "saturation_liquid");
    saturation_gas_key_ = Keys::readKey(*plist_,domain_,"saturation gas", "saturation_gas");
    saturation_ice_key_ = Keys::readKey(*plist_,domain_,"saturation ice", "saturation_ice");
    water_content_key_ = Keys::readKey(*plist_,domain_,"water content","water_content");
    //relative_permeability_key_ = Keys::readKey(*plist_,domain_,"relative permeability","relative_permeability");
    //matric_pressure_key_ = Keys::readKey(*plist_,domain_,"matric pressure","matric_pressure");
    cap_pres_key_ = Keys::readKey(*plist_, domain_, "capillary pressure key", "capillary_pressure_gas_liq");
    //liquid_density_key_ = Keys::readKey(*plist_, domain_, "mass density liquid", "mass_density_liquid");
    liquid_density_key_ = Keys::readKey(*plist_, domain_, "molar density liquid", "molar_density_liquid");
    ice_density_key_ = Keys::readKey(*plist_, domain_, "mass density ice", "mass_density_ice");
    gas_density_key_ = Keys::readKey(*plist_, domain_,"mass density gas", "mass_density_gas");
    gas_density_key_test_ = Keys::readKey(*plist_, domain_, "mass density gas", "mass_density_gas");
    rock_density_key_ = Keys::readKey(*plist_, domain_, "density rock", "density_rock");

    //energy
    T_key_ = Keys::readKey(*plist_, domain_, "temperature", "temperature");
    thermal_conductivity_key_ = Keys::readKey(*plist_, domain_, "thermal conductivity", "thermal_conductivity");

    //Sources
    surface_water_source_key_ = Keys::readKey(*plist_, domain_surface_, "surface water source", "water_source");
    //surface_energy_source_key_ =
    //  Keys::readKey(*plist_, domain_surface_, "surface energy source", "total_energy_source");
    subsurface_water_source_key_ =
      Keys::readKey(*plist_, domain_, "subsurface water source", "water_source");
    //subsurface_energy_source_key_ =
    //  Keys::readKey(*plist_, domain_, "subsurface energy source", "total_energy_source");
    surface_energy_source_ecosim_key_ =
      Keys::readKey(*plist_, domain_surface_, "surface energy source ecosim", "ecosim_source");
    surface_water_source_ecosim_key_ =
      Keys::readKey(*plist_, domain_surface_, "surface water source ecosim", "ecosim_water_source");

    subsurface_energy_source_ecosim_key_ =
      Keys::readKey(*plist_, domain_, "subsurface energy source ecosim", "subsurface_ecosim_source");
    subsurface_water_source_ecosim_key_ =
      Keys::readKey(*plist_, domain_, "subsurface water source ecosim", "subsurface_ecosim_water_source");

    //Other
    cell_volume_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");
    //ecosim_aux_data_key_ = Keys::readKey(*plist_, domain_, "ecosim aux data", "ecosim_aux_data");
    f_wp_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
    f_root_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");

    //Custom Evaluator keys
    hydraulic_conductivity_key_ = Keys::readKey(*plist_, domain_, "hydraulic conductivity", "hydraulic_conductivity");
    //bulk_density_key_ = Keys::readKey(*plist_, domain_, "bulk density", "bulk_density");

    //Surface balance items
    T_surf_key_ = Keys::readKey(*plist_, domain_surface_, "temperature", "temperature");
    sw_key_ =
      Keys::readKey(*plist_, domain_surface_, "incoming shortwave radiation", "incoming_shortwave_radiation");
    sw_incident_key_ =
      Keys::readKey(*plist_, domain_surface_, "incident shortwave radiation", "incident_shortwave_radiation");
    lw_key_ =
      Keys::readKey(*plist_,domain_surface_, "incoming longwave radiation", "incoming_longwave_radiation");
    air_temp_key_ = Keys::readKey(*plist_, domain_surface_, "air temperature", "air_temperature");
    vp_air_key_ = Keys::readKey(*plist_, domain_surface_, "vapor pressure air", "vapor_pressure_air");
    wind_speed_key_ = Keys::readKey(*plist_, domain_surface_, "wind speed", "wind_speed");
    p_rain_key_ = Keys::readKey(*plist_, domain_surface_, "precipitation rain", "precipitation_rain");
    p_snow_key_ = Keys::readKey(*plist_, domain_surface_, "precipitation snow", "precipitation_snow");
    p_total_key_ = Keys::readKey(*plist_, domain_surface_, "precipitation total", "precipitation_total");
    elev_key_ = Keys::readKey(*plist_, domain_surface_, "elevation", "elevation");
    aspect_key_ = Keys::readKey(*plist_, domain_surface_, "aspect", "aspect");
    slope_key_ = Keys::readKey(*plist_, domain_surface_, "slope", "slope_magnitude");
    snow_depth_key_ = Keys::readKey(*plist_, domain_surface_, "snow depth", "snow_depth");
    snow_albedo_key_ = Keys::readKey(*plist_, domain_surface_, "snow_albedo", "snow_albedo");
    canopy_snow_key_ = Keys::readKey(*plist_, domain_surface_, "canopy snow", "canopy_snow");
    //snow_temperature_key_ = Keys::readKey(*plist_, domain_surface_, "snow temperature", "snow_temperature");

    //Plant Phenology Datasets
    lai_key_ = Keys::readKey(*plist_, domain_surface_, "LAI", "LAI");
    sai_key_ = Keys::readKey(*plist_, domain_surface_, "SAI", "SAI");
    v_type_key_ = Keys::readKey(*plist_, domain_surface_, "vegetation type", "vegetation_type");

    //Atmospheric abundance keys
    /*atm_n2_ = plist_->get<double>("atmospheric N2");
    atm_o2_ = plist_->get<double>("atmospheric O2");
    atm_co2_ = plist_->get<double>("atmospheric CO2");
    atm_ch4_ = plist_->get<double>("atmospheric CH4");
    atm_n2o_ = plist_->get<double>("atmospheric N2O");
    atm_h2_ = plist_->get<double>("atmospheric H2");
    atm_nh3_ = plist_->get<double>("atmospheric NH3");*/

    //Starting values and parameters for precribed phenology / albedo

    pressure_at_field_capacity = plist_->get<double>("field capacity [Mpa]");
    pressure_at_wilting_point = plist_->get<double>("wilting point [Mpa]");
    p_bool = plist_->get<bool>("EcoSIM precipitation", false);
    a_bool = plist_->get<bool>("prescribe snow albedo", false);
    pheno_bool = plist_->get<bool>("prescribe phenology", false);
    microbe_bool = plist_->get<bool>("microbe model", false);
    num_pfts = plist_->get<int>("number of pfts", 1);

    // EcoSIM internal state entries (by ATS name) to visualize; none by default
    if (plist_->isParameter("visualize ecosim state")) {
      visualize_internal_state_ =
        plist_->get<Teuchos::Array<std::string>>("visualize ecosim state").toVector();
    }
    internal_state_version_key_ =
      Keys::getKey(domain_surface_, "ecosim_internal_state_layout_version");

    //Parameters for times and time of year
    dt_ = plist_->get<double>("initial time step");
    c_m_ = plist_->get<double>("heat capacity [MJ mol^-1 K^-1]");
    day0_ = plist_->get<int>("starting day of year [0-364]");
    year0_ = plist_->get<int>("starting year");

    curr_day_ = day0_;
    curr_year_ = year0_;

    //This initializes the engine (found in BGCEngine.cc) This is the code that
    //actually points to the driver
    if (!plist_->isParameter("engine")) {
      Errors::Message msg;
      msg << "No 'engine' parameter found in the parameter list for 'BGC'.\n";
      Exceptions::amanzi_throw(msg);
    }
    if (!plist_->isParameter("engine input file")) {
      Errors::Message msg;
      msg << "No 'engine input file' parameter found in the parameter list for 'BGC'.\n";
      Exceptions::amanzi_throw(msg);
    }
    std::string engine_name = plist_->get<std::string>("engine");
    engine_inputfile = plist_->get<std::string>("engine input file");
    eco_engine_ = Teuchos::rcp(new EcoEngine(engine_name, engine_inputfile));
  }


// -- Parse the parameter list and declare evaluators other PKs depend on
void EcoSIM::parseParameterList()
{
  PK_Physical_Default::parseParameterList();

  // These sources are primary variables owned by this PK (see Setup). PKs
  // that use them (e.g. overland flow) are set up before this one, so declare
  // the evaluator type here: parseParameterList() runs for every PK before
  // any Setup().
  for (const auto& key : { surface_energy_source_ecosim_key_, surface_water_source_ecosim_key_,
                           subsurface_energy_source_ecosim_key_, subsurface_water_source_ecosim_key_ }) {
    if (S_->HasEvaluatorList(key) &&
        S_->GetEvaluatorList(key).get<std::string>("evaluator type", "primary variable") !=
          "primary variable") {
      Errors::Message msg;
      msg << "EcoSIM: \"" << key << "\" is computed by the EcoSIM PK; remove its evaluator \""
          << S_->GetEvaluatorList(key).get<std::string>("evaluator type")
          << "\" from the state's evaluators list.";
      Exceptions::amanzi_throw(msg);
    }
    S_->GetEvaluatorList(key).set<std::string>("evaluator type", "primary variable");
  }
}

// -- Destroy ansilary data structures.
EcoSIM::~EcoSIM()
  {
    eco_engine_->FreeState(eco_env_, eco_feedback_);
    eco_engine_->FreeInternalState(eco_istate_);
  }

// -- Setup step
void EcoSIM::Setup() {
  PK_Physical_Default::Setup();
  //Need to do some basic setup of the columns:
  mesh_surf_ = S_->GetMesh(domain_surface_);
  mesh_ = S_->GetMesh(domain_);
  int num_columns_ =
    mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  for (unsigned int column = 0; column != num_columns_; ++column) {
    int f = mesh_surf_->getEntityParent(AmanziMesh::Entity_kind::CELL, column);
    auto col_iter = mesh_->columns.getCells(column);
    std::size_t ncol_cells = col_iter.size();

    double column_area = mesh_->getFaceArea(f);

    if (ncells_per_col_ < 0) {
      ncells_per_col_ = ncol_cells;
    } else {
      AMANZI_ASSERT(ncol_cells == ncells_per_col_);
    }
  }

  //Setting records for variables ONLY used by the EcoSIM PK, this includes, weather forcings,
  // ecosim surface and subsurface forces, and surface variables from EcoSIM that are saved
  // over to ATS (evaporation fluxes, snow related variables)
  if (!S_->HasRecord(snow_depth_key_,tag_next_)) {
        S_->Require<CompositeVector, CompositeVectorSpace>(snow_depth_key_, tag_next_, snow_depth_key_)
          .SetMesh(mesh_surf_)
          ->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  
  // EcoSIM internal state: one field per entry of EcoSIM's layout
  // (ATSStateRegistryMod.F90). PRIVATE entries are carried state, checkpointed
  // and visualized only on request; OUTPUT entries are visualized but not
  // checkpointed. Fields are written only by this PK.
  EcoSizes layout_sizes;
  layout_sizes.ncells_per_col_ = ncells_per_col_;
  layout_sizes.num_components = 1;
  layout_sizes.num_columns = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL,
                                                        AmanziMesh::Parallel_kind::OWNED);
  layout_sizes.num_pfts = num_pfts;
  internal_state_layout_ = eco_engine_->InternalStateLayout(layout_sizes);

  for (const auto& name : visualize_internal_state_) {
    auto entry = std::find_if(internal_state_layout_.begin(), internal_state_layout_.end(),
                              [&](const EcoInternalStateEntry& e) { return e.ats_name == name; });
    if (entry == internal_state_layout_.end()) {
      Errors::Message msg;
      msg << "EcoSIM: \"visualize ecosim state\" entry \"" << name
          << "\" is not an EcoSIM internal state variable. Valid names are:";
      for (const auto& e : internal_state_layout_) msg << " " << e.ats_name;
      Exceptions::amanzi_throw(msg);
    }
  }

  internal_state_keys_.clear();
  for (const auto& entry : internal_state_layout_) {
    Key key = Keys::getKey(domain_surface_, entry.ats_name);
    internal_state_keys_.push_back(key);
    S_->Require<CompositeVector, CompositeVectorSpace>(key, tag_next_, key)
      .SetMesh(mesh_surf_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, entry.num_components);
    auto& record = S_->GetRecordW(key, tag_next_, key);
    if (entry.role == kEcoRoleOutput) {
      record.set_io_checkpoint(false);
    } else {
      bool requested = std::find(visualize_internal_state_.begin(),
                                 visualize_internal_state_.end(),
                                 entry.ats_name) != visualize_internal_state_.end();
      record.set_io_vis(requested);
    }
  }
  S_->Require<int>(internal_state_version_key_, tag_next_, name_);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "EcoSIM internal state (layout version "
               << eco_engine_->InternalStateLayoutVersion() << "):" << std::endl
               << "  ATS field | EcoSIM variable | components | units | role" << std::endl;
    for (int i = 0; i != internal_state_layout_.size(); ++i) {
      const auto& e = internal_state_layout_[i];
      *vo_->os() << "  " << internal_state_keys_[i] << " | " << e.ecosim_name << " | "
                 << e.num_components << " | " << e.units << " | "
                 << (e.role == kEcoRoleOutput ? "output" : "private") << std::endl;
    }
  }

  // Sources EcoSIM returns to ATS. This PK owns them as primary evaluators and
  // marks them changed after every write, so evaluators that use them (e.g.
  // overland flow's molar water source) recompute. Inputs must not declare
  // evaluators for these keys.
  requireEvaluatorPrimary(surface_energy_source_ecosim_key_, tag_next_, *S_, name_)
    .SetMesh(mesh_surf_)
    ->SetGhosted(false)
    ->SetComponent("cell", AmanziMesh::CELL, 1);

  requireEvaluatorPrimary(surface_water_source_ecosim_key_, tag_next_, *S_, name_)
    .SetMesh(mesh_surf_)
    ->SetGhosted(false)
    ->SetComponent("cell", AmanziMesh::CELL, 1);

  requireEvaluatorPrimary(subsurface_energy_source_ecosim_key_, tag_next_, *S_, name_)
    .SetMesh(mesh_)
    ->SetGhosted(false)
    ->SetComponent("cell", AmanziMesh::CELL, 1);

  requireEvaluatorPrimary(subsurface_water_source_ecosim_key_, tag_next_, *S_, name_)
    .SetMesh(mesh_)
    ->SetGhosted(false)
    ->SetComponent("cell", AmanziMesh::CELL, 1);

  /*S_->Require<CompositeVector, CompositeVectorSpace>(snow_temperature_key_ , tag_next_, snow_temperature_key_)
          .SetMesh(mesh_)
          ->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, 1);*/

  S_->RequireEvaluator(snow_albedo_key_, tag_next_);
  S_->Require<CompositeVector, CompositeVectorSpace>(snow_albedo_key_, tag_next_).SetMesh(mesh_surf_)
	        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->RequireEvaluator(sw_key_, tag_next_);
  S_->Require<CompositeVector, CompositeVectorSpace>(sw_key_, tag_next_).SetMesh(mesh_surf_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->RequireEvaluator(sw_incident_key_, tag_next_);
  S_->Require<CompositeVector, CompositeVectorSpace>(sw_incident_key_, tag_next_).SetMesh(mesh_surf_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  
  S_->RequireEvaluator(lai_key_, tag_next_);
  S_->Require<CompositeVector, CompositeVectorSpace>(lai_key_, tag_next_).SetMesh(mesh_surf_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_pfts);

  S_->RequireEvaluator(sai_key_, tag_next_);
  S_->Require<CompositeVector, CompositeVectorSpace>(sai_key_, tag_next_).SetMesh(mesh_surf_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_pfts);

  S_->RequireEvaluator(v_type_key_, tag_next_);
  S_->Require<CompositeVector, CompositeVectorSpace>(v_type_key_, tag_next_).SetMesh(mesh_surf_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_pfts);
  
  if (!S_->HasRecord(canopy_snow_key_, tag_next_)) {
    S_->Require<CompositeVector, CompositeVectorSpace>(canopy_snow_key_, tag_next_, canopy_snow_key_)
      .SetMesh(mesh_surf_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, num_pfts);
  }

  Teuchos::OSTab tab = vo_->getOSTab();

  //EcoSIM can do its own precipitation partitioning so you can put in total precipitation if
  // you want EcoSIM to do it, or snow/rain if the forcing is already split.
  if (p_bool) {
     S_->RequireEvaluator(p_total_key_, tag_next_);
     S_->Require<CompositeVector, CompositeVectorSpace>(p_total_key_, tag_next_).SetMesh(mesh_surf_)
       ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  } else {
     S_->RequireEvaluator(p_snow_key_, tag_next_);
     S_->Require<CompositeVector, CompositeVectorSpace>(p_snow_key_, tag_next_).SetMesh(mesh_surf_)
       ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
     S_->RequireEvaluator(p_rain_key_, tag_next_);
     S_->Require<CompositeVector, CompositeVectorSpace>(p_rain_key_, tag_next_).SetMesh(mesh_surf_)
       ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  //Setup custom evaluators for EcoSIM, found in constitutive relations
  requireEvaluatorAtNext(hydraulic_conductivity_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireEvaluatorAtCurrent(hydraulic_conductivity_key_, tag_current_, *S_, name_);

//Setup variables that were owned by ATS SEB, but now are controlled by EcoSIM
// Can remove SEB from the cycle_driver
  requireEvaluatorAtNext(lw_key_, tag_next_, *S_)
    .SetMesh(mesh_surf_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireEvaluatorAtCurrent(lw_key_, tag_current_, *S_, name_);

  requireEvaluatorAtNext(air_temp_key_, tag_next_, *S_)
    .SetMesh(mesh_surf_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireEvaluatorAtCurrent(air_temp_key_, tag_current_, *S_, name_);

  requireEvaluatorAtNext(vp_air_key_, tag_next_, *S_)
    .SetMesh(mesh_surf_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireEvaluatorAtCurrent(vp_air_key_, tag_current_, *S_, name_);

  requireEvaluatorAtNext(wind_speed_key_, tag_next_, *S_)
    .SetMesh(mesh_surf_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireEvaluatorAtCurrent(wind_speed_key_, tag_current_, *S_, name_);

  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("green") << "Setup of PK was successful"
    << vo_->reset() << std::endl << std::endl;
  }
}

// -- Initialize owned (dependent) variables.
void EcoSIM::Initialize() {
  PK_Physical_Default::Initialize();
  //Need to know the number of components to initialize data structures

  //Transport removal:
  int mole_fraction_num = 1;
  if(microbe_bool){
    const Epetra_MultiVector& mole_fraction= *(S_->GetPtr<CompositeVector>(mole_fraction_key_, Tags::DEFAULT)->ViewComponent("cell"));
    mole_fraction_num = mole_fraction.NumVectors();
    std::cout << "mole fraction number: " << mole_fraction.NumVectors() << std::endl;
  }
  Teuchos::OSTab tab = vo_->getOSTab();

  num_columns_ = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  //Now we call the engine's init state function which allocates the data
  eco_engine_->InitState(eco_env_, eco_feedback_, ncells_per_col_, mole_fraction_num, num_columns_, num_pfts);
  eco_engine_->InitInternalState(eco_istate_, internal_state_layout_, num_columns_);

  int ierr = 0;

  if (S_->HasRecord(ice_density_key_, Tags::DEFAULT)) {
    S_->GetEvaluator(saturation_ice_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(ice_density_key_, Tags::DEFAULT).Update(*S_, name_);
    has_ice = true;
  } else {
    Teuchos::OSTab tab = vo_->getOSTab();
    //*vo_->os() << "Did not find ice key" << std::endl;
    has_ice = false;
  }

  //Check for total precipitation and set the record if it's there

  if (p_bool) {
    S_->GetW<CompositeVector>(p_total_key_, Tags::DEFAULT, "surface-precipitation_total").PutScalar(0.0);
    S_->GetRecordW(p_total_key_, Tags::DEFAULT, "surface-precipitation_total").set_initialized();
  } else {
    S_->GetW<CompositeVector>(p_snow_key_, Tags::DEFAULT, "surface-precipitation_snow").PutScalar(0.0);
    S_->GetRecordW(p_snow_key_, Tags::DEFAULT, "surface-precipitation_snow").set_initialized();
    S_->GetW<CompositeVector>(p_rain_key_, Tags::DEFAULT, "surface-precipitation_rain").PutScalar(0.0);
    S_->GetRecordW(p_rain_key_, Tags::DEFAULT, "surface-precipitation_rain").set_initialized();
  }

  S_->GetW<CompositeVector>(snow_depth_key_, Tags::DEFAULT, "surface-snow_depth").PutScalar(0.0);
  S_->GetRecordW(snow_depth_key_, Tags::DEFAULT, "surface-snow_depth").set_initialized();

  S_->GetW<CompositeVector>(canopy_snow_key_, Tags::DEFAULT, canopy_snow_key_).PutScalar(0.0);
  S_->GetRecordW(canopy_snow_key_, Tags::DEFAULT, canopy_snow_key_).set_initialized();

  for (const auto& key : { surface_energy_source_ecosim_key_, surface_water_source_ecosim_key_,
                           subsurface_energy_source_ecosim_key_, subsurface_water_source_ecosim_key_ }) {
    S_->GetW<CompositeVector>(key, tag_next_, name_).PutScalar(0.0);
    S_->GetRecordW(key, tag_next_, name_).set_initialized();
    changedEvaluatorPrimary(key, tag_next_, *S_);
  }

  //S_->GetW<CompositeVector>(snow_temperature_key_, Tags::DEFAULT, snow_temperature_key_).PutScalar(0.0);
  //S_->GetRecordW(snow_temperature_key_, Tags::DEFAULT, snow_temperature_key_).set_initialized();

  // EcoSIM internal state: carried state is filled from EcoSIM's
  // initialization below and overwritten by the checkpoint on restart;
  // outputs stay zero until the first advance.
  for (const auto& key : internal_state_keys_) {
    S_->GetW<CompositeVector>(key, tag_next_, key).PutScalar(0.0);
    S_->GetRecordW(key, tag_next_, key).set_initialized();
  }
  S_->GetW<int>(internal_state_version_key_, tag_next_, name_) =
    eco_engine_->InternalStateLayoutVersion();
  S_->GetRecordW(internal_state_version_key_, tag_next_, name_).set_initialized();

  //Initialize owned evaluators
  S_->GetW<CompositeVector>(hydraulic_conductivity_key_, Tags::DEFAULT, name_).PutScalar(1.0);
  S_->GetRecordW(hydraulic_conductivity_key_, Tags::DEFAULT, name_).set_initialized();

  //S_->GetW<CompositeVector>(bulk_density_key_, Tags::DEFAULT, "bulk_density").PutScalar(1.0);
  //S_->GetRecordW(bulk_density_key_, Tags::DEFAULT, "bulk_density").set_initialized();

  //S_->GetW<CompositeVector>(matric_pressure_key_, Tags::DEFAULT, "matric_pressure").PutScalar(1.0);
  //S_->GetRecordW(matric_pressure_key_, Tags::DEFAULT, "matric_pressure").set_initialized();

  int num_columns_ = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  //loop over processes instead:
  num_columns_global = mesh_surf_->getMap(AmanziMesh::Entity_kind::CELL, false).NumGlobalElements();
  num_columns_local = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  num_columns_global_ptype = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  //Loop over processes and Initalize EcoSIM on that process
  int numProcesses, p_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  for (int k = 0; k < numProcesses; ++k) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (p_rank==k) {
      num_columns_local = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

      InitializeSingleProcess(p_rank);
    }
  }

  // verbose message
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("green") << "Initialization of PK was successful, T="
        << S_->get_time() << vo_->reset() << std::endl << std::endl;
  }
}

void EcoSIM::CommitStep(double t_old, double t_new, const Tag& tag) {

  // I don't know that we will have much to do here. In SimpleBGC they just copy
  // Data to the pfts, which we won't be doing. In Alquimia they just save the time
  // As below.

  saved_time_ = t_new;

}

bool EcoSIM::AdvanceStep(double t_old, double t_new, bool reinit) {
  double dt = t_new - t_old;
  current_time_ = saved_time_ + dt;


  Teuchos::OSTab out = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << S_->get_time(tag_current_)
               << " t1 = " << S_->get_time(tag_next_) << " h = " << dt << std::endl
               << "Current day: " << curr_day_ << "Current year: " << curr_year_ << std::endl
               << "----------------------------------------------------------------" << std::endl;


  // Ensure dependencies are filled
  // Transport removal

  if(microbe_bool) {
  S_->GetEvaluator(mole_fraction_key_, Tags::DEFAULT).Update(*S_, name_);
  }
  S_->GetEvaluator(porosity_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(saturation_liquid_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(water_content_key_, Tags::DEFAULT).Update(*S_, name_);
  //S_->GetEvaluator(relative_permeability_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(liquid_density_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rock_density_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(T_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(cell_volume_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(f_wp_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(f_root_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(cap_pres_key_, Tags::DEFAULT).Update(*S_, name_);
  //S_->GetEvaluator(subsurface_energy_source_key_, Tags::DEFAULT).Update(*S_, name_);
  //S_->GetEvaluator(subsurface_water_source_key_, Tags::DEFAULT).Update(*S_, name_);
  //S_->GetEvaluator(matric_pressure_key_, Tags::DEFAULT).Update(*S_, name_);


  //Surface data
  S_->GetEvaluator(T_surf_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(sw_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(sw_incident_key_, Tags::DEFAULT).Update(*S_, name_);
  //S_->GetEvaluator(lw_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(air_temp_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(vp_air_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(wind_speed_key_, Tags::DEFAULT).Update(*S_, name_);

  S_->GetEvaluator(elev_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(aspect_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(slope_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(snow_albedo_key_, Tags::DEFAULT).Update(*S_, name_);

  S_->GetEvaluator(lai_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(sai_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(v_type_key_, Tags::DEFAULT).Update(*S_, name_);
  //S_->GetEvaluator(canopy_snow_key_, Tags::DEFAULT).Update(*S_, name_);

  if (p_bool){
    S_->GetEvaluator(p_total_key_, Tags::DEFAULT).Update(*S_, name_);
  } else {
    S_->GetEvaluator(p_rain_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(p_snow_key_, Tags::DEFAULT).Update(*S_, name_);
  }

  //S_->GetEvaluator(surface_energy_source_key_, Tags::DEFAULT).Update(*S_, name_);
  //S_->GetEvaluator(surface_water_source_key_, Tags::DEFAULT).Update(*S_, name_);

  /*if (has_gas) {
    S_->GetEvaluator(saturation_gas_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(gas_density_key_, Tags::DEFAULT).Update(*S_, name_);
  }*/

  if (has_ice) {
    S_->GetEvaluator(saturation_ice_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(ice_density_key_, Tags::DEFAULT).Update(*S_, name_);
  }

  S_->GetEvaluator(T_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(thermal_conductivity_key_, Tags::DEFAULT).Update(*S_, name_);

  //Update owned evaluators
  /*Teuchos::RCP<const CompositeVector> hydra_cond = S_->GetPtr<CompositeVector>(hydraulic_conductivity_key_, Tags::DEFAULT);
  S_->GetEvaluator(hydraulic_conductivity_key_, Tags::DEFAULT).Update(*S_, name_);
  const Epetra_MultiVector& hydraulic_conductivity = *(*S_->Get<CompositeVector>("hydraulic_conductivity", tag_next_)
          .ViewComponent("cell",false))(0);*/

  AmanziMesh::Entity_ID num_columns_ = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  // grab the required fields

  S_->GetEvaluator("porosity", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& porosity = *(*S_->Get<CompositeVector>("porosity", tag_next_)
      .ViewComponent("cell",false))(0);

  S_->GetEvaluator("saturation_liquid", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& liquid_saturation = *(*S_->Get<CompositeVector>("saturation_liquid", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("capillary_pressure_gas_liq", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& capillary_pressure = *(*S_->Get<CompositeVector>("capillary_pressure_gas_liq", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("water_content", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& water_content = *(*S_->Get<CompositeVector>("water_content", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("mass_density_liquid", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& liquid_density = *(*S_->Get<CompositeVector>("mass_density_liquid", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("density_rock", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& rock_density = *(*S_->Get<CompositeVector>("density_rock", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("cell_volume", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& cell_volume = *(*S_->Get<CompositeVector>("cell_volume", tag_next_)
          .ViewComponent("cell",false))(0);

  /*if (has_gas) {
    S_->GetEvaluator("saturation_gas", tag_next_).Update(*S_, name_);
    const Epetra_MultiVector& gas_saturation = *(*S_->Get<CompositeVector>("saturation_gas", tag_next_)
            .ViewComponent("cell",false))(0);
  }*/

  //Atm abundances
  S_->GetEvaluator("surface-incoming_shortwave_radiation", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& sw_rad = *(*S_->Get<CompositeVector>("surface-incoming_shortwave_radiation", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("surface-incoming_longwave_radiation", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& lw_rad = *(*S_->Get<CompositeVector>("surface-incoming_longwave_radiation", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("surface-air_temperature", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& t_air = *(*S_->Get<CompositeVector>("surface-air_temperature", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("surface-vapor_pressure_air", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& p_vap = *(*S_->Get<CompositeVector>("surface-vapor_pressure_air", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("surface-wind_speed", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& v_wind = *(*S_->Get<CompositeVector>("surface-wind_speed", tag_next_)
          .ViewComponent("cell",false))(0);

  //Define before loop to prevent scope issues:
  const Epetra_MultiVector* p_tot = nullptr;
  const Epetra_MultiVector* p_rain = nullptr;
  const Epetra_MultiVector* p_snow = nullptr;

  if(p_bool){
    S_->GetEvaluator("surface-precipitation_total", tag_next_).Update(*S_, name_);
    p_tot = &(*(*S_->Get<CompositeVector>("surface-precipitation_total", tag_next_)
    	.ViewComponent("cell",false))(0));
  } else {
    S_->GetEvaluator("surface-precipitation_rain", tag_next_).Update(*S_, name_);
    p_rain = &(*(*S_->Get<CompositeVector>("surface-precipitation_rain", tag_next_)
          .ViewComponent("cell",false))(0));
    S_->GetEvaluator("surface-precipitation_snow", tag_next_).Update(*S_, name_);
    p_snow = &(*(*S_->Get<CompositeVector>("surface-precipitation_snow", tag_next_)
          .ViewComponent("cell",false))(0));
  }

  S_->GetEvaluator("surface-elevation", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& elevation = *S_->Get<CompositeVector>("surface-elevation", tag_next_)
          .ViewComponent("cell",false);

  S_->GetEvaluator("surface-aspect", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& aspect = *S_->Get<CompositeVector>("surface-aspect", tag_next_)
          .ViewComponent("cell",false);

  S_->GetEvaluator("surface-slope_magnitude", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& slope = *S_->Get<CompositeVector>("surface-slope_magnitude", tag_next_)
          .ViewComponent("cell",false);

  S_->GetEvaluator("surface-LAI", tag_next_).Update(*S_, name_);
  Teuchos::RCP<const Epetra_MultiVector> LAI = S_->Get<CompositeVector>("surface-LAI", tag_next_)
          .ViewComponent("cell",false);

  S_->GetEvaluator("surface-SAI", tag_next_).Update(*S_, name_);
  Teuchos::RCP<const Epetra_MultiVector> SAI = S_->Get<CompositeVector>("surface-SAI", tag_next_)
          .ViewComponent("cell",false);

  S_->GetEvaluator("surface-vegetation_type", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& vegetation_type = *(*S_->Get<CompositeVector>("surface-vegetation_type", tag_next_)
          .ViewComponent("cell",false))(0);
  
  /*S_->GetEvaluator("surface-canopy_snow", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& canopy_snow = *(*S_->Get<CompositeVector>("surface-canopy_snow", tag_next_)
          .ViewComponent("cell",false))(0);*/

  if (has_ice) {
    S_->GetEvaluator("mass_density_ice", tag_next_).Update(*S_, name_);
    const Epetra_MultiVector& ice_density = *(*S_->Get<CompositeVector>("mass_density_ice", tag_next_)
            .ViewComponent("cell",false))(0);

    S_->GetEvaluator("saturation_ice", tag_next_).Update(*S_, name_);
    const Epetra_MultiVector& ice_saturation = *(*S_->Get<CompositeVector>("saturation_ice", tag_next_)
            .ViewComponent("cell",false))(0);
  }

  S_->GetEvaluator("temperature", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& temp = *(*S_->Get<CompositeVector>("temperature", tag_next_)
      .ViewComponent("cell",false))(0);

  S_->GetEvaluator("thermal_conductivity", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& thermal_conductivity = *(*S_->Get<CompositeVector>("thermal_conductivity", tag_next_)
      .ViewComponent("cell",false))(0);

  //loop over processes instead:
  num_columns_global = mesh_surf_->getMap(AmanziMesh::Entity_kind::CELL,false).NumGlobalElements();
  num_columns_local = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  num_columns_global_ptype = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  //Trying to loop over processors now:
  int numProcesses, p_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  for (int k = 0; k < numProcesses; ++k) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (p_rank==k) {
      num_columns_local = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

      AdvanceSingleProcess(dt, p_rank);
    }
  }
  // PLACE TIME AND YEAR ITERATOR HERE
  if (curr_day_ != 364) {
    curr_day_ = curr_day_ + 1;
  } else {
    curr_day_ = 0;
    curr_year_ = curr_year_ + 1;
  }

  return false;
}

// helper function for pushing field to column
void EcoSIM::FieldToColumn_(AmanziMesh::Entity_ID column, const Epetra_Vector& vec,
       Teuchos::Ptr<Epetra_SerialDenseVector> col_vec)
{
  auto col_iter = mesh_->columns.getCells(column);

  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    std::size_t vec_index = col_iter[i];

    (*col_vec)[i] = vec[vec_index];
  }
}

void EcoSIM::FieldToColumn_(AmanziMesh::Entity_ID column, const Teuchos::Ptr<Epetra_SerialDenseVector> vec,
       Teuchos::Ptr<Epetra_SerialDenseVector> col_vec)
{
  auto col_iter = mesh_->columns.getCells(column);

  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    std::size_t vec_index = col_iter[i];

    (*col_vec)[i] = (*vec)[vec_index];
  }
}

//Helper function but for datasets that are multivalued in every cell (concentrations)
void EcoSIM::MatrixFieldToColumn_(AmanziMesh::Entity_ID column, const Epetra_MultiVector& m_arr,
  Teuchos::Ptr<Epetra_SerialDenseMatrix> col_arr)
  {
    int n_comp = m_arr.NumVectors();
    auto col_iter = mesh_->columns.getCells(column);

    for (int j=0; j!=n_comp; ++j){
      for (std::size_t i=0; i!=col_iter.size(); ++i) {
        (*col_arr)(i,j) = m_arr[j][col_iter[i]];
      }
    }
  }

// helper function for pushing column back to field
void EcoSIM::ColumnToField_(AmanziMesh::Entity_ID column, Epetra_Vector& vec,
                               Teuchos::Ptr<Epetra_SerialDenseVector> col_vec)
{
  auto col_iter = mesh_->columns.getCells(column);
  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    vec[col_iter[i]] = (*col_vec)[i];
  }
}

void EcoSIM::ColumnToField_(AmanziMesh::Entity_ID column, Teuchos::Ptr<Epetra_SerialDenseVector> vec,
                               Teuchos::Ptr<Epetra_SerialDenseVector> col_vec)
{
  auto col_iter = mesh_->columns.getCells(column);
  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    (*vec)[col_iter[i]] = (*col_vec)[i];
  }
}

void EcoSIM::MatrixColumnToField_(AmanziMesh::Entity_ID column, Epetra_MultiVector& m_arr,
  Teuchos::Ptr<Epetra_SerialDenseMatrix> col_arr) {

    int n_comp = m_arr.NumVectors();
    auto col_iter = mesh_->columns.getCells(column);

    for (int j=0; j!=n_comp; ++j){
      for (std::size_t i=0; i!=col_iter.size(); ++i) {
        m_arr[j][col_iter[i]] = (*col_arr)(i,j);
      }
    }

  }

// helper function for collecting column dz and depth
void EcoSIM::ColDepthDz_(AmanziMesh::Entity_ID column,
                            Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                            Teuchos::Ptr<Epetra_SerialDenseVector> dz) {
  AmanziMesh::Entity_ID f_above = mesh_surf_->getEntityParent(AmanziMesh::Entity_kind::CELL, column);
  auto col_iter = mesh_->columns.getCells(column);
  ncells_per_col_ = col_iter.size();

  AmanziGeometry::Point surf_centroid = mesh_->getFaceCentroid(f_above);
  AmanziGeometry::Point neg_z(3);
  neg_z.set(0.,0.,-1);

  Teuchos::OSTab tab = vo_->getOSTab();

  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    // depth centroid
    (*depth)[i] = surf_centroid[2] - mesh_->getCellCentroid(col_iter[i])[2];

    const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(col_iter[i]);

    // -- mimics implementation of build_columns() in Mesh
    double mindp = 999.0;
    AmanziMesh::Entity_ID f_below = -1;
    for (std::size_t j=0; j!=faces.size(); ++j) {
      AmanziGeometry::Point normal = mesh_->getFaceNormal(faces[j]);
      if (dirs[j] == -1) normal *= -1;
      normal /= AmanziGeometry::norm(normal);

      double dp = -normal * neg_z;
      if (dp < mindp) {
        mindp = dp;
        f_below = faces[j];
      }
    }

    // -- fill the val
    (*dz)[i] = mesh_->getFaceCentroid(f_above)[2] - mesh_->getFaceCentroid(f_below)[2];
    AMANZI_ASSERT( (*dz)[i] > 0. );
    f_above = f_below;
  }
}

// helper function for collecting dz, depth, and volume for a given column
void EcoSIM::VolDepthDz_(AmanziMesh::Entity_ID column,
                            Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                            Teuchos::Ptr<Epetra_SerialDenseVector> dz,
			    Teuchos::Ptr<Epetra_SerialDenseVector> volume) {
  AmanziMesh::Entity_ID f_above = mesh_surf_->getEntityParent(AmanziMesh::Entity_kind::CELL, column);
  auto col_iter = mesh_->columns.getCells(column);
  ncells_per_col_ = col_iter.size();

  AmanziGeometry::Point surf_centroid = mesh_->getFaceCentroid(f_above);
  AmanziGeometry::Point neg_z(3);
  neg_z.set(0.,0.,-1);

  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    // depth centroid
    (*depth)[i] = surf_centroid[2] - mesh_->getCellCentroid(col_iter[i])[2];

    // dz
    // -- find face_below
    //AmanziMesh::Entity_ID_List faces;
    //std::vector<int> dirs;
    //mesh_->cell_get_faces_and_dirs(col_iter[i], &faces, &dirs);

    const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(col_iter[i]);

    //double vol = mesh_->cell_volume(col_iter[i]);
    (*volume)[i] = mesh_->getCellVolume(col_iter[i]);

    // -- mimics implementation of build_columns() in Mesh
    double mindp = 999.0;
    AmanziMesh::Entity_ID f_below = -1;
    for (std::size_t j=0; j!=faces.size(); ++j) {
      AmanziGeometry::Point normal = mesh_->getFaceNormal(faces[j]);
      if (dirs[j] == -1) normal *= -1;
      normal /= AmanziGeometry::norm(normal);

      double dp = -normal * neg_z;
      if (dp < mindp) {
        mindp = dp;
        f_below = faces[j];
      }
    }

    // -- fill the val
    (*dz)[i] = mesh_->getFaceCentroid(f_above)[2] - mesh_->getFaceCentroid(f_below)[2];
    AMANZI_ASSERT( (*dz)[i] > 0. );
    f_above = f_below;
  }
}

// ATS -> EcoSIM, once at setup: run parameters and flags
void EcoSIM::CopyConfigToEcoSIM_(EcoConfig& config)
{
  config.heat_capacity = c_m_;
  config.field_capacity = pressure_at_field_capacity;
  config.wilting_point = pressure_at_wilting_point;
  config.p_bool = p_bool;
  config.a_bool = a_bool;
  config.pheno_bool = pheno_bool;
  config.microbe_bool = microbe_bool;
  config.pft_file = engine_inputfile.data();
}

// ATS -> EcoSIM, at setup and every advance: the environment, the snow depth
// EcoSIM carries (in the feedback container), and the internal state
void EcoSIM::CopyToEcoSIM_process(int proc_rank,
                                  EcoEnvironment& environment,
                                  EcoFeedback& feedback,
                                  const Tag& water_tag)
{
  //This is the copy function for a loop over a single process instead of a single column
  //Fill state with ATS variables that are going to be changed by EcoSIM
  const Epetra_Vector& porosity = *(*S_->Get<CompositeVector>(porosity_key_, water_tag).ViewComponent("cell", false))(0);

  //Transport removal
  int mole_fraction_num = 1;
  const Epetra_MultiVector* mole_fraction = nullptr;
  if(microbe_bool) {
    mole_fraction= S_->GetPtr<CompositeVector>(mole_fraction_key_, water_tag)->ViewComponent("cell").get();
    mole_fraction_num = mole_fraction->NumVectors();
  }
  
  const Epetra_Vector& liquid_saturation = *(*S_->Get<CompositeVector>(saturation_liquid_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& water_content = *(*S_->Get<CompositeVector>(water_content_key_, water_tag).ViewComponent("cell", false))(0);
  //const Epetra_Vector& relative_permeability = *(*S_->Get<CompositeVector>(relative_permeability_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& liquid_density = *(*S_->Get<CompositeVector>(liquid_density_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& rock_density = *(*S_->Get<CompositeVector>(rock_density_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& cell_volume = *(*S_->Get<CompositeVector>(cell_volume_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& hydraulic_conductivity = *(*S_->Get<CompositeVector>(hydraulic_conductivity_key_, water_tag).ViewComponent("cell", false))(0);
  //const Epetra_Vector& matric_pressure = *(*S_->Get<CompositeVector>(matric_pressure_key_, water_tag).ViewComponent("cell", false))(0);
  //const Epetra_Vector& bulk_density = *(*S_->Get<CompositeVector>(bulk_density_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& rooting_depth_fraction = *(*S_->Get<CompositeVector>(f_root_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& plant_wilting_factor = *(*S_->Get<CompositeVector>(f_wp_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& temp = *(*S_->Get<CompositeVector>(T_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& thermal_conductivity = *(*S_->Get<CompositeVector>(thermal_conductivity_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& capillary_pressure = *(*S_->Get<CompositeVector>(cap_pres_key_, water_tag).ViewComponent("cell", false))(0);

  //const auto& shortwave_radiation = *S_.Get<CompositeVector>(sw_key_, water_tag).ViewComponent("cell", false);
  const Epetra_Vector& temp_surf = *(*S_->Get<CompositeVector>(T_surf_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& shortwave_radiation = *(*S_->Get<CompositeVector>(sw_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& shortwave_incident_radiation = *(*S_->Get<CompositeVector>(sw_incident_key_, water_tag).ViewComponent("cell", false))(0);
  //const Epetra_Vector& longwave_radiation = *(*S_->Get<CompositeVector>(lw_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& air_temperature = *(*S_->Get<CompositeVector>(air_temp_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& vapor_pressure_air = *(*S_->Get<CompositeVector>(vp_air_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& wind_speed= *(*S_->Get<CompositeVector>(wind_speed_key_, water_tag).ViewComponent("cell", false))(0);
  //define these outside of the loop to prevent issues:
  const Epetra_Vector* precipitation = nullptr;
  const Epetra_Vector* precipitation_snow = nullptr;
  
  if(p_bool){
    precipitation = &(*(*S_->Get<CompositeVector>(p_total_key_, water_tag).ViewComponent("cell", false))(0));
  } else {
    precipitation = &(*(*S_->Get<CompositeVector>(p_rain_key_, water_tag).ViewComponent("cell", false))(0));
    precipitation_snow = &(*(*S_->Get<CompositeVector>(p_snow_key_, water_tag).ViewComponent("cell", false))(0));
  }

  const Epetra_Vector& elevation = *(*S_->Get<CompositeVector>(elev_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& aspect = *(*S_->Get<CompositeVector>(aspect_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& slope = *(*S_->Get<CompositeVector>(slope_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& snow_albedo = *(*S_->Get<CompositeVector>(snow_albedo_key_, water_tag).ViewComponent("cell", false))(0);

  Teuchos::RCP<const Epetra_MultiVector> LAI = S_->Get<CompositeVector>(lai_key_, water_tag).ViewComponent("cell", false);
  Teuchos::RCP<const Epetra_MultiVector> SAI = S_->Get<CompositeVector>(sai_key_, water_tag).ViewComponent("cell", false);
  //const Epetra_Vector& vegetation_type = *(*S_->Get<CompositeVector>(v_type_key_, water_tag).ViewComponent("cell", false))(0);
  Teuchos::RCP<const Epetra_MultiVector> vegetation_type = S_->Get<CompositeVector>(v_type_key_, tag_next_).ViewComponent("cell", false);

  // the EcoSIM sources (EcoFeedback) are EcoSIM outputs and are not sent in

  auto& snow_depth = *S_->GetW<CompositeVector>(snow_depth_key_,tag_next_,snow_depth_key_).ViewComponent("cell");
  auto& canopy_snow_in = *S_->GetW<CompositeVector>(canopy_snow_key_,tag_next_,canopy_snow_key_).ViewComponent("cell");
  

  auto col_porosity = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_wc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_relative_permeability = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_mat_p = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_r_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_vol = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_g_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_g_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_temp = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_cond = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_h_cond = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_b_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_depth = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_dz = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_wp = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_rf = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_v_type = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_ss_energy_source = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_ss_water_source = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  //auto col_depth_c = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_cap_pres = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));

  auto col_mole_fraction = Teuchos::rcp(new Epetra_SerialDenseMatrix(ncells_per_col_,mole_fraction_num));

  //Gather columns on this process:
  num_columns_global = mesh_surf_->getMap(AmanziMesh::Entity_kind::CELL,false).NumGlobalElements();
  num_columns_local = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  num_columns_global_ptype = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  int p_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  MPI_Barrier(MPI_COMM_WORLD);

  num_columns_local = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  for (int column=0; column!=num_columns_local; ++column) {
    FieldToColumn_(column,porosity,col_porosity.ptr());
    FieldToColumn_(column,liquid_saturation,col_l_sat.ptr());
    FieldToColumn_(column,water_content,col_wc.ptr());
    //FieldToColumn_(column,relative_permeability,col_relative_permeability.ptr());
    FieldToColumn_(column,liquid_density,col_l_dens.ptr());
    FieldToColumn_(column,rock_density,col_r_dens.ptr());
    FieldToColumn_(column,cell_volume,col_vol.ptr());
    FieldToColumn_(column,hydraulic_conductivity,col_h_cond.ptr());
    //FieldToColumn_(column,bulk_density,col_b_dens.ptr());
    FieldToColumn_(column,plant_wilting_factor,col_wp.ptr());
    FieldToColumn_(column,rooting_depth_fraction,col_rf.ptr());
    //setting matric pressure to capillary pressure here
    FieldToColumn_(column,capillary_pressure,col_mat_p.ptr());
    FieldToColumn_(column,temp, col_temp.ptr());
    FieldToColumn_(column,thermal_conductivity,col_cond.ptr());
    //FieldToColumn_(column,capillary_pressure,col_cap_pres.ptr());
    //FieldToColumn_(column,canopy_snow,col_canopy_snow.ptr());
    //FieldToColumn_(column,vegetation_type,col_v_type.ptr());

    if(microbe_bool) {
      MatrixFieldToColumn_(column, *mole_fraction, col_mole_fraction.ptr());
    }
    
    // This is for computing depth
    //ColDepthDz_(column, col_depth.ptr(), col_dz.ptr());

    //Grabbing the cross sectional area in the z direction
    int f = mesh_surf_->getEntityParent(AmanziMesh::Entity_kind::CELL, column);
    auto col_iter = mesh_->columns.getCells(column);
    std::size_t ncol_cells = col_iter.size();

    double column_area = mesh_->getFaceArea(f);
    //std::cout << "column: " << column << " column_area: " << column_area << std::endl;
    environment.column_area.data[column] = column_area;

    VolDepthDz_(column, col_depth.ptr(), col_dz.ptr(), col_vol.ptr());
    double sum = 0.0;
    /*for (int i = ncells_per_col_ - 1; i >= 0; --i) {
        sum += (*col_dz)[i];
        (*col_depth_c)[i] = sum;
    }*/

    for (int pft = 0; pft != num_pfts; ++pft) {
      const Epetra_Vector& pft_vec = *(*vegetation_type)(pft);
      const Epetra_Vector& lai_vec = *(*LAI)(pft);
      const Epetra_Vector& sai_vec = *(*SAI)(pft);

      environment.plant_functional_type.data[column * ncells_per_col_ + pft] = pft_vec[column];
      environment.LAI.data[column * ncells_per_col_ + pft] = lai_vec[column];
      environment.SAI.data[column * ncells_per_col_ + pft] = sai_vec[column];
    }
    
    
    for (int i=0; i < ncells_per_col_; ++i) {
      environment.liquid_density.data[column * ncells_per_col_ + i] = (*col_l_dens)[i];
      environment.rock_density.data[column * ncells_per_col_ + i] = (*col_r_dens)[i];
      environment.porosity.data[column * ncells_per_col_ + i] = (*col_porosity)[i];
      environment.water_content.data[column * ncells_per_col_ + i] = (*col_wc)[i];
      environment.hydraulic_conductivity.data[column * ncells_per_col_ + i] = (*col_h_cond)[i];
      //environment.bulk_density.data[column * ncells_per_col_ + i] = (*col_b_dens)[i];
      environment.matric_pressure.data[column * ncells_per_col_ + i] = (*col_mat_p)[i];
      environment.temperature.data[column * ncells_per_col_ + i] = (*col_temp)[i];
      //environment.plant_functional_type.data[column * ncells_per_col_ + i] = (*col_v_type)[i];
      environment.plant_wilting_factor.data[column * ncells_per_col_ + i] = (*col_wp)[i];
      environment.rooting_depth_fraction.data[column * ncells_per_col_ + i] = (*col_rf)[i];
      environment.liquid_saturation.data[column * ncells_per_col_ + i] = (*col_l_sat)[i];
      environment.relative_permeability.data[column * ncells_per_col_ + i] = (*col_relative_permeability)[i];
      environment.volume.data[column * ncells_per_col_ + i] = (*col_vol)[i];
      environment.depth.data[column * ncells_per_col_ + i] = (*col_depth)[i];
      //environment.depth_c.data[column * ncells_per_col_ + i] = (*col_depth_c)[i];
      environment.dz.data[column * ncells_per_col_ + i] = (*col_dz)[i];

      /*if (has_gas) {
        environment.gas_saturation.data[column * ncells_per_col_ + i] = (*col_g_sat)[i];
        //environment.gas_density.data[column][i] = (*col_g_dens)[i];
      }*/

      if (has_ice) {
        environment.ice_density.data[column * ncells_per_col_ + i] = (*col_i_dens)[i];
        environment.ice_saturation.data[column * ncells_per_col_ + i] = (*col_i_sat)[i];
      }

    }
    //fill surface variables

    //environment.temperature.data[1] = temp_surf[column];
    feedback.snow_depth.data[column] = snow_depth[0][column];
    for (int pft = 0; pft != num_pfts; ++pft) {
      feedback.canopy_snow.data[column * num_pfts + pft] = canopy_snow_in[pft][column];
    }

    environment.shortwave_radiation.data[column] = shortwave_radiation[column];
    //environment.longwave_radiation.data[column] = longwave_radiation[column];
    environment.air_temperature.data[column] = air_temperature[column];
    environment.vapor_pressure_air.data[column] = vapor_pressure_air[column];
    environment.wind_speed.data[column] = wind_speed[column];
    environment.elevation.data[column] = elevation[column];
    environment.aspect.data[column] = aspect[column];
    environment.slope.data[column] = slope[column];
    environment.snow_albedo.data[column] = snow_albedo[column];
    //environment.vegetation_type.data[column] = vegetation_type[column];

    if(p_bool){
       environment.precipitation.data[column] = (*precipitation)[column];
   } else {
       environment.precipitation.data[column] = (*precipitation)[column];
       environment.precipitation_snow.data[column] = (*precipitation_snow)[column];
   }
    
    if(microbe_bool){ 
      for (int i = 0; i < environment.mole_fraction.cells; i++) {
        for (int k = 0; k < environment.mole_fraction.components; k++) {
          //Assuming column is being grabbed correctly this should look like:
          //int index = column * ncells_per_col_ * num_components + i
          int index = i + (column * ncells_per_col_) + (k * ncells_per_col_ * num_columns_);
        
          environment.mole_fraction.data[index] = (*col_mole_fraction)(i,k);
        }
      }
    }
    
  }

  //Fill the atmospheric abundances
  //NOTE: probably want to add an if statement here to only do this only once
  environment.atm_n2 = atm_n2_;
  environment.atm_o2 = atm_o2_;
  environment.atm_co2 = atm_co2_;
  environment.atm_ch4 = atm_ch4_;
  environment.atm_n2o = atm_n2o_;
  environment.atm_h2 = atm_h2_;
  environment.atm_nh3 = atm_nh3_;

  CopyInternalStateToEcoSIM_(eco_istate_);
}

// EcoSIM -> ATS, every advance: the feedback and the internal state
void EcoSIM::CopyFromEcoSIM_process(const int column,
                                    const EcoFeedback& feedback,
                                    const Tag& water_tag)
{

  //Transport removal
  int mole_fraction_num = 1;
  if (microbe_bool) {
  Epetra_MultiVector& mole_fraction= *(S_->GetPtrW<CompositeVector>(mole_fraction_key_, Tags::DEFAULT, "subsurface transport")->ViewComponent("cell",false));
  int mole_fraction_num = mole_fraction.NumVectors();
  }
  
  auto& porosity = *(*S_->GetW<CompositeVector>(porosity_key_, Tags::DEFAULT, porosity_key_).ViewComponent("cell",false))(0);
  auto& liquid_saturation = *(*S_->GetW<CompositeVector>(saturation_liquid_key_, Tags::DEFAULT, saturation_liquid_key_).ViewComponent("cell",false))(0);
  auto& water_content = *(*S_->GetW<CompositeVector>(water_content_key_, Tags::DEFAULT, water_content_key_).ViewComponent("cell",false))(0);
  //auto& suction_head = *(*S_->GetW<CompositeVector>(suc_key_, Tags::DEFAULT, suc_key_).ViewComponent("cell",false))(0);
  //auto& relative_permeability = *(*S_->GetW<CompositeVector>(relative_permeability_key_, Tags::DEFAULT, relative_permeability_key_).ViewComponent("cell",false))(0);
  auto& liquid_density = *(*S_->GetW<CompositeVector>(liquid_density_key_, Tags::DEFAULT, liquid_density_key_).ViewComponent("cell",false))(0);
  auto& rock_density = *(*S_->GetW<CompositeVector>(rock_density_key_, Tags::DEFAULT, rock_density_key_).ViewComponent("cell",false))(0);
  auto& cell_volume = *(*S_->GetW<CompositeVector>(cell_volume_key_, Tags::DEFAULT, cell_volume_key_).ViewComponent("cell",false))(0);

  auto& surface_energy_source = *(*S_->GetW<CompositeVector>(surface_energy_source_ecosim_key_, tag_next_, name_).ViewComponent("cell", false))(0);
  auto& subsurface_energy_source = *(*S_->GetW<CompositeVector>(subsurface_energy_source_ecosim_key_, tag_next_, name_).ViewComponent("cell", false))(0);

  auto& surface_water_source = *(*S_->GetW<CompositeVector>(surface_water_source_ecosim_key_, tag_next_, name_).ViewComponent("cell", false))(0);
  auto& subsurface_water_source = *(*S_->GetW<CompositeVector>(subsurface_water_source_ecosim_key_, tag_next_, name_).ViewComponent("cell", false))(0);
  auto& temp = *(*S_->GetW<CompositeVector>(T_key_, Tags::DEFAULT, "subsurface energy").ViewComponent("cell",false))(0);
  auto& thermal_conductivity = *(*S_->GetW<CompositeVector>(thermal_conductivity_key_, Tags::DEFAULT, thermal_conductivity_key_).ViewComponent("cell",false))(0);
  //auto& snow_temperature = *(*S_->GetW<CompositeVector>(snow_temperature_key_, Tags::DEFAULT, snow_temperature_key_).ViewComponent("cell", false))(0);
  auto& snow_depth = *S_->GetW<CompositeVector>(snow_depth_key_,tag_next_,snow_depth_key_).ViewComponent("cell");
  auto& canopy_snow_out = *S_->GetW<CompositeVector>(canopy_snow_key_,tag_next_,canopy_snow_key_).ViewComponent("cell");

  auto col_porosity = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_wc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_suc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_relative_permeability = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_r_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_vol = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_g_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_g_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_temp = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_cond = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_h_cond = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_b_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_ss_energy_source = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_ss_water_source = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  //auto col_snow_temperature = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));

  auto col_mole_fraction = Teuchos::rcp(new Epetra_SerialDenseMatrix(mole_fraction_num,ncells_per_col_));

  //Gather columns on this process:
  num_columns_global = mesh_surf_->getMap(AmanziMesh::Entity_kind::CELL, false).NumGlobalElements();
  num_columns_local = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  num_columns_global_ptype = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  //Trying to loop over processors now:
  int p_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  MPI_Barrier(MPI_COMM_WORLD);

  for (int col=0; col!=num_columns_local; ++col) {
    if (std::isnan(feedback.surface_water_source.data[col]) ||
        std::isinf(feedback.surface_water_source.data[col])) {
        std::cout << "Process " << p_rank << " found bad value at column "
                  << col << ": " << feedback.surface_water_source.data[col] << std::endl;
    }
  }

  num_columns_local = mesh_surf_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  double energy_source_tot = feedback.surface_energy_source.data[0];
  double water_source_tot = feedback.surface_water_source.data[0];
  double snow_depth_cell = feedback.snow_depth.data[0];

  for (int col=0; col!=num_columns_local; ++col) {
    
    /*for (int i=0; i < ncells_per_col_; ++i) {
      (*col_ss_water_source)[i] = feedback.subsurface_water_source.data[col * ncells_per_col_ + i];
      (*col_ss_energy_source)[i] = feedback.subsurface_energy_source.data[col * ncells_per_col_ + i];
      (*col_snow_temperature)[i] = state.snow_temperature.data[col * ncells_per_col_ + i];
      (*col_canopy_snow)[i] = state.canopy_snow.data[col * ncells_per_col_ + i];
    }

    ColumnToField_(col, subsurface_water_source, col_ss_water_source.ptr());
    ColumnToField_(col, subsurface_energy_source, col_ss_energy_source.ptr());
    ColumnToField_(col, snow_temperature, col_snow_temperature.ptr());
    ColumnToField_(col, canopy_snow, col_canopy_snow.ptr());*/
    
    surface_energy_source[col] = feedback.surface_energy_source.data[col]/(3600.0);
    surface_water_source[col] = feedback.surface_water_source.data[col]/(3600.0);
    snow_depth[0][col] = feedback.snow_depth.data[col];
    for (int pft = 0; pft != num_pfts; ++pft) {
      canopy_snow_out[pft][col] = feedback.canopy_snow.data[col * num_pfts + pft];
    }
  }
   
  for (int col=0; col!=num_columns_local; ++col) {
    for (int i=0; i < ncells_per_col_; ++i) {
      (*col_ss_water_source)[i] = feedback.subsurface_water_source.data[col * ncells_per_col_ + i];
      //(*col_ss_energy_source)[i] = feedback.subsurface_energy_source.data[col * ncells_per_col_ + i];
      //(*col_snow_temperature)[i] = state.snow_temperature.data[col * ncells_per_col_ + i];
      //(*col_canopy_snow)[i] = state.canopy_snow.data[col * ncells_per_col_ + i];
    }

    ColumnToField_(col, subsurface_water_source, col_ss_water_source.ptr());
    //ColumnToField_(col, subsurface_energy_source, col_ss_energy_source.ptr());
    //ColumnToField_(col, snow_temperature, col_snow_temperature.ptr());
    //ColumnToField_(col, canopy_snow, col_canopy_snow.ptr());


  }

  // EcoSIM gives the subsurface water source as water volume per grid cell
  // (i.e. per ATS cell) per hour [m^3 h^-1]. Subsurface flow multiplies its
  // source by cell volume, so it needs [mol m^-3 s^-1]:
  //   q = S * molar density [mol m^-3] / (3600 [s h^-1] * cell volume [m^3])
  for (int c = 0; c != subsurface_water_source.MyLength(); ++c) {
    subsurface_water_source[c] *= liquid_density[c] / (3600.0 * cell_volume[c]);
  }

  CopyInternalStateFromEcoSIM_(eco_istate_);

  // tell the evaluators that use EcoSIM's sources that they changed
  for (const auto& key : { surface_energy_source_ecosim_key_, surface_water_source_ecosim_key_,
                           subsurface_energy_source_ecosim_key_, subsurface_water_source_ecosim_key_ }) {
    changedEvaluatorPrimary(key, tag_next_, *S_);
  }
}

// Copy the EcoSIM internal state from the ATS fields that hold it. All
// entries are copied; EcoSIM only reads back the PRIVATE ones.
void EcoSIM::CopyInternalStateToEcoSIM_(EcoInternalState& internal_state)
{
  const int n = internal_state.values_per_column;
  int offset = 0;
  for (int i = 0; i != internal_state_layout_.size(); ++i) {
    const auto& field =
      *S_->Get<CompositeVector>(internal_state_keys_[i], tag_next_).ViewComponent("cell", false);
    for (int col = 0; col != num_columns_local; ++col) {
      for (int k = 0; k != internal_state_layout_[i].num_components; ++k) {
        internal_state.values.data[col * n + offset + k] = field[k][col];
      }
    }
    offset += internal_state_layout_[i].num_components;
  }
}

// Copy the EcoSIM internal state into the ATS fields that hold it.
void EcoSIM::CopyInternalStateFromEcoSIM_(const EcoInternalState& internal_state)
{
  const int n = internal_state.values_per_column;
  int offset = 0;
  for (int i = 0; i != internal_state_layout_.size(); ++i) {
    auto& field = *S_->GetW<CompositeVector>(internal_state_keys_[i], tag_next_, internal_state_keys_[i])
                     .ViewComponent("cell", false);
    for (int col = 0; col != num_columns_local; ++col) {
      for (int k = 0; k != internal_state_layout_[i].num_components; ++k) {
        field[k][col] = internal_state.values.data[col * n + offset + k];
      }
    }
    offset += internal_state_layout_[i].num_components;
  }
}

int EcoSIM::InitializeSingleProcess(int proc)
{
  int num_iterations = 1;
  int num_columns = 1;

  num_columns = num_columns_local;

  Teuchos::OSTab tab = vo_->getOSTab();

  // C and Fortran must agree on the exchange containers
  eco_engine_->CheckContainerSizes();

  CopyConfigToEcoSIM_(eco_config_);
  CopyToEcoSIM_process(proc, eco_env_, eco_feedback_, Tags::DEFAULT);

  eco_sizes_.num_columns = num_columns;
  eco_sizes_.ncells_per_col_ = ncells_per_col_;
  eco_sizes_.num_components = 1;
  eco_sizes_.num_pfts = num_pfts;

  eco_engine_->Setup(eco_config_, eco_env_, eco_feedback_, eco_istate_, eco_sizes_,
                     num_iterations, num_columns, ncells_per_col_);
  // EcoSIM's initialized carried state; overwritten by the checkpoint on restart
  CopyInternalStateFromEcoSIM_(eco_istate_);
  
  return 0;
}

int EcoSIM::AdvanceSingleProcess(double dt, int proc)
{
  //Function to run EcoSIMs advance function, ecosim is run as an hourly model
  // Every hour the data is taken from ATS state copied to ecosim
  // The model is run and passed back to ATS.

  int num_iterations = 1;
  int num_columns = 1;

  num_columns = num_columns_local;

  // Time tracking variables
  current_time_ = S_->get_time();        //Current time
  static double last_ecosim_time = 0.0;
  int total_days = static_cast<int>(current_time_ / 86400.0);
  int current_day = (day0_ + total_days) % 365;
  int current_year = year0_ + ((day0_ + total_days)/365);

  eco_env_.current_day = current_day;
  eco_env_.current_year = current_year;

  Teuchos::OSTab tab = vo_->getOSTab();

  if (current_time_ - last_ecosim_time >= 3600.0) {
    *vo_->os() << "Hour completed at total_time: " << current_time_
               << ", Year: " << current_year << ", Day: " << current_day << std::endl;
    *vo_->os() << "Running EcoSIM Advance: " << std::endl;

	// the restored internal state must have been written with this EcoSIM layout
	int saved_version = S_->Get<int>(internal_state_version_key_, tag_next_);
	if (saved_version != eco_engine_->InternalStateLayoutVersion()) {
	  Errors::Message msg;
	  msg << "EcoSIM: internal state was saved with layout version " << saved_version
	      << " but this EcoSIM uses layout version " << eco_engine_->InternalStateLayoutVersion()
	      << "; restarting from this checkpoint is not supported.";
	  Exceptions::amanzi_throw(msg);
	}

	CopyToEcoSIM_process(proc, eco_env_, eco_feedback_, Tags::DEFAULT);

	eco_engine_->Advance(dt, eco_env_, eco_feedback_, eco_istate_, eco_sizes_, num_iterations, num_columns);

    CopyFromEcoSIM_process(proc, eco_feedback_, Tags::DEFAULT);

	last_ecosim_time = current_time_;
  }

  return 0;
}


} // namespace EcoSIM
} // namespace Amanzi
