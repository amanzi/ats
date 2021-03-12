/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <vector>

#include "boost/algorithm/string.hpp"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

#include "BCs.hh"
#include "errors.hh"
#include "Explicit_TI_RK.hh"
#include "FieldEvaluator.hh"
#include "Mesh.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_Diffusion.hh"
#include "PDE_Accumulation.hh"
#include "PK_DomainFunctionFactory.hh"
#include "PK_Utils.hh"

#include "MultiscaleTransportPorosityFactory.hh"
#include "TransportDomainFunction.hh"
#include "TransportBoundaryFunction_Alquimia.hh"
#include "TransportSourceFunction_Alquimia.hh"
#include "TransportDomainFunction_UnitConversion.hh"

#include "transport_ats.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* New constructor compatible with new MPC framework.
****************************************************************** */
Transport_ATS::Transport_ATS(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& glist,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& solution) :
  PK(pk_tree, glist, S, solution),
  PK_PhysicalExplicit<Epetra_Vector>(pk_tree, glist, S, solution)
{
  if (plist_->isParameter("component names")) {
    component_names_ = plist_->get<Teuchos::Array<std::string> >("component names").toVector();
    mol_masses_ = plist_->get<Teuchos::Array<double> >("component molar masses").toVector();
  } else {
    Errors::Message msg("Transport PK: parameter \"component names\" is missing.");
    Exceptions::amanzi_throw(msg);
  }

  subcycling_ = plist_->get<bool>("transport subcycling", false);

  water_source_in_meters_ = plist_->get<bool>("water source in meters", true);

  // initialize io
  Teuchos::RCP<Teuchos::ParameterList> units_list = Teuchos::sublist(glist, "units");
  units_.Init(*units_list);
  vo_ = Teuchos::null;
}


/* ******************************************************************
* Setup for Alquimia.
****************************************************************** */
#ifdef ALQUIMIA_ENABLED
void Transport_ATS::SetupAlquimia(Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk,
                                 Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine)
{
  chem_pk_ = chem_pk;
  chem_engine_ = chem_engine;

  if (chem_engine_ != Teuchos::null) {
    // Retrieve the component names (primary and secondary) from the chemistry
    // engine.
    std::vector<std::string> component_names;
    chem_engine_->GetPrimarySpeciesNames(component_names);
    component_names_ = component_names;
    for (int i = 0; i < chem_engine_->NumAqueousComplexes(); ++i) {
      char secondary_name[128];
      snprintf(secondary_name, 127, "secondary_%d", i);
      component_names_.push_back(secondary_name);
    }
  }
}
#endif

void Transport_ATS::set_states(const Teuchos::RCP<State>& S,
                                  const Teuchos::RCP<State>& S_inter,
                                  const Teuchos::RCP<State>& S_next)
{
    //S_ = S;
  S_inter_ = S_inter;
  S_next_ = S_next;
}


/* ******************************************************************
* Define structure of this PK.
****************************************************************** */
void Transport_ATS::Setup(const Teuchos::Ptr<State>& S)
{
  saturation_key_ = Keys::readKey(*plist_, domain_, "saturation liquid", "saturation_liquid");
  prev_saturation_key_ = Keys::readKey(*plist_, domain_, "previous saturation liquid", "prev_saturation_liquid");
  flux_key_ = Keys::readKey(*plist_, domain_, "mass flux", "mass_flux");
  permeability_key_ = Keys::readKey(*plist_, domain_, "permeability", "permeability");
  tcc_key_ = Keys::readKey(*plist_, domain_, "concentration", "total_component_concentration");
  conserve_qty_key_ = Keys::readKey(*plist_, domain_, "conserved quantity", "total_component_quantity");
  porosity_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
  molar_density_key_ = Keys::readKey(*plist_, domain_, "molar density", "molar_density_liquid");
  tcc_matrix_key_ = Keys::readKey(*plist_, domain_, "tcc matrix", "total_component_concentration_matrix");
  solid_residue_mass_key_ =  Keys::readKey(*plist_, domain_, "solid residue", "solid_residue_mass");
  mass_src_key_ = Keys::readKey(*plist_, domain_, "mass source", "mass_source");
  water_content_key_ = Keys::readKey(*plist_, domain_, "water content", "water_content");
  cv_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");

  water_tolerance_ = plist_->get<double>("water tolerance", 1e-6);
  dissolution_ = plist_->get<bool>("allow dissolution", false);
  max_tcc_ = plist_->get<double>("maximum concentration", 0.9);

  dim = mesh_->space_dimension();

  // cross-coupling of PKs
  Teuchos::RCP<Teuchos::ParameterList> physical_models =
      Teuchos::sublist(plist_, "physical models and assumptions");
  bool abs_perm = physical_models->get<bool>("permeability field is required", false);
  std::string multiscale_model = physical_models->get<std::string>("multiscale model", "single porosity");

  // require my field
  int ncomponents = component_names_.size();
  std::vector<std::vector<std::string> > subfield_names(1);
  subfield_names[0] = component_names_;

  if (component_names_.size() == 0) {
    Errors::Message msg;
    msg << "Transport PK: list of solutes is empty.\n";
    Exceptions::amanzi_throw(msg);
  }

  const auto& tcc_fac = S->RequireField(tcc_key_)->SetMesh(mesh_)->SetGhosted(true);
  name_ = "state"; // note: this is required because the chemistry PK is Amanzi code and uses this.
  S->RequireField(tcc_key_, name_, subfield_names)
      ->SetComponent("cell", AmanziMesh::CELL, ncomponents);

  // we may not use this, but having it in vis makes life so much easier
  S->RequireField(cv_key_);
  S->RequireFieldEvaluator(cv_key_);

  // This vector stores the conserved amount (in mols) of ncomponent
  // transported components, plus two for water.  The first water component is
  // given by the water content (in mols) at the old time plus dt * all fluxes
  // treated explictly.  The second water component is given by the water
  // content at the new time plus dt * all fluxes treated implicitly (notably
  // just DomainCoupling fluxes, which must be able to take all the transported
  // quantity.)
  S->RequireField(conserve_qty_key_, name_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, ncomponents+2);

  // other things that I own
  S->RequireField(prev_saturation_key_, name_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(prev_saturation_key_, name_)->set_io_vis(false);

  S->RequireField(solid_residue_mass_key_,  name_)->SetMesh(mesh_)->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::CELL, ncomponents);

  // require state fields
  if (abs_perm) {
    S->RequireField(permeability_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, dim);
    S->RequireFieldEvaluator(permeability_key_);
  }

  S->RequireField(flux_key_)->SetMesh(mesh_)->SetGhosted(true)
    ->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireFieldEvaluator(flux_key_);

  S->RequireField(saturation_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(saturation_key_);

  S->RequireField(porosity_key_, porosity_key_)->SetMesh(mesh_)->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(porosity_key_);

  S->RequireField(molar_density_key_, molar_density_key_)->SetMesh(mesh_)->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(molar_density_key_);

  if (plist_->sublist("source terms").isSublist("geochemical")){
    S->RequireField(mass_src_key_, mass_src_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(mass_src_key_);
  }

  // require multiscale fields
  multiscale_porosity_ = false;
  if (multiscale_model == "dual porosity") {
    multiscale_porosity_ = true;
    Teuchos::RCP<Teuchos::ParameterList>
        msp_list = Teuchos::sublist(plist_, "multiscale models", true);
    msp_ = CreateMultiscaleTransportPorosityPartition(mesh_, msp_list);

    S->RequireField(tcc_matrix_key_, name_, subfield_names)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, ncomponents);
  }

  // Create verbosity object.
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = plist_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject(name_, vlist));
  db_ = Teuchos::rcp(new Debugger(mesh_, name_, *plist_));
}


/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.
****************************************************************** */
void Transport_ATS::Initialize(const Teuchos::Ptr<State>& S)
{
  // Set initial values for transport variables.
  dt_ = dt_debug_ = t_physics_ = 0.0;
  double time = S->time();
  if (time >= 0.0) t_physics_ = time;

  if (plist_->isSublist("initial condition")) {
    S->GetField(tcc_key_,name_)->Initialize(plist_->sublist("initial condition"));
  }

  internal_tests = 0;
  tests_tolerance = TRANSPORT_CONCENTRATION_OVERSHOOT;

  bc_scaling = 0.0;

  Teuchos::OSTab tab = vo_->getOSTab();
  MyPID = mesh_->get_comm()->MyPID();

  // initialize missed fields
  InitializeFields_(S);

  //create copies
  S->RequireFieldCopy(tcc_key_, "subcycling", name_);
  tcc_tmp = S->GetField(tcc_key_, name_)->GetCopy("subcycling", name_)->GetFieldData();

  S->RequireFieldCopy(saturation_key_, "subcycle_start", name_);
  ws_subcycle_start = S->GetFieldCopyData(saturation_key_, "subcycle_start",name_)
    ->ViewComponent("cell");
  S->RequireFieldCopy(saturation_key_, "subcycle_end", name_);
  ws_subcycle_end = S->GetFieldCopyData(saturation_key_, "subcycle_end", name_)
    ->ViewComponent("cell");
  S->RequireFieldCopy(molar_density_key_, "subcycle_start", name_);
  mol_dens_subcycle_start = S->GetFieldCopyData(molar_density_key_, "subcycle_start",name_)->ViewComponent("cell");
  S->RequireFieldCopy(molar_density_key_, "subcycle_end", name_);
  mol_dens_subcycle_end = S->GetFieldCopyData(molar_density_key_, "subcycle_end", name_)->ViewComponent("cell");

  S->RequireFieldCopy(flux_key_, "next_timestep", name_);
  flux_copy_ = S->GetFieldCopyData(flux_key_,  "next_timestep", name_)->ViewComponent("face", true);
  flux_copy_ -> PutScalar(0.);

  // Check input parameters. Due to limited amount of checks, we can do it earlier.
  Policy(S.ptr());

  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  // extract control parameters
  InitializeAll_();

  // state pre-prosessing
  Teuchos::RCP<const CompositeVector> cv;

  ws_ = S->GetFieldData(saturation_key_)->ViewComponent("cell", false);
  ws_prev_ = S->GetFieldData(prev_saturation_key_)->ViewComponent("cell", false);

  phi_ = S->GetFieldData(porosity_key_)->ViewComponent("cell", false);

  mol_dens_ = S->GetFieldData(molar_density_key_)->ViewComponent("cell", false);
  mol_dens_prev_ = S->GetFieldData(molar_density_key_)->ViewComponent("cell", false);

  tcc = S->GetFieldData(tcc_key_, name_);

  flux_ = S->GetFieldData(flux_key_)->ViewComponent("face", true);
  solid_qty_ = S->GetFieldData(solid_residue_mass_key_, name_)->ViewComponent("cell", false);

  //create vector of conserved quatities
  conserve_qty_ = S->GetFieldData(conserve_qty_key_, name_)->ViewComponent("cell", true);

  // upwind
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  upwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));
  downwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));

  IdentifyUpwindCells();

  // advection block initialization
  current_component_ = -1;

  const Epetra_Map& cmap_owned = mesh_->cell_map(false);

  // reconstruction initialization
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  limiter_ = Teuchos::rcp(new Operators::LimiterCell(mesh_));
  lifting_ = Teuchos::rcp(new Operators::ReconstructionCell(mesh_));

  // mechanical dispersion
  flag_dispersion_ = false;
  if (plist_->isSublist("material properties")) {
    Teuchos::RCP<Teuchos::ParameterList>
        mdm_list = Teuchos::sublist(plist_, "material properties");
    mdm_ = CreateMDMPartition(mesh_, mdm_list, flag_dispersion_);
    if (flag_dispersion_) CalculateAxiSymmetryDirection();
  }

  // create boundary conditions
  if (plist_->isSublist("boundary conditions")) {
    // -- try tracer-type conditions
    PK_DomainFunctionFactory<TransportDomainFunction> factory(mesh_);
    Teuchos::ParameterList& conc_bcs_list = plist_->sublist("boundary conditions").sublist("concentration");

    for (const auto& it : conc_bcs_list) {
      std::string name = it.first;
      if (conc_bcs_list.isSublist(name)) {
        Teuchos::ParameterList& bc_list = conc_bcs_list.sublist(name);
        std::string bc_type = bc_list.get<std::string>("spatial distribution method", "none");

        if (bc_type == "domain coupling") {
          // domain couplings are special -- they always work on all components
          Teuchos::RCP<TransportDomainFunction> bc =
            factory.Create(bc_list, "fields", AmanziMesh::FACE, Kxy);

          for (int i = 0; i < component_names_.size(); i++){
            bc->tcc_names().push_back(component_names_[i]);
            bc->tcc_index().push_back(i);
          }
          bc->set_state(S_);
          bcs_.push_back(bc);

        } else if (bc_type == "subgrid") {
          // subgrid domains take a BC from a single entity of a parent mesh --
          // find the GID of that entity.
          Teuchos::Array<std::string> regions(1, domain_);
          std::size_t last_of = domain_.find_last_of(":");
          AMANZI_ASSERT(last_of != std::string::npos);
          int gid = std::stoi(domain_.substr(last_of+1, domain_.size()));
          bc_list.set("entity_gid_out", gid);

          Teuchos::RCP<TransportDomainFunction> bc =
            factory.Create(bc_list, "boundary concentration", AmanziMesh::FACE, Kxy);

          for (int i = 0; i < component_names_.size(); i++){
            bc->tcc_names().push_back(component_names_[i]);
            bc->tcc_index().push_back(i);
          }
          bc->set_state(S_);
          bcs_.push_back(bc);

        } else {
          Teuchos::RCP<TransportDomainFunction> bc =
            factory.Create(bc_list, "boundary concentration function", AmanziMesh::FACE, Kxy);
          bc->set_state(S_);

          std::vector<std::string> tcc_names = bc_list.get<Teuchos::Array<std::string>>("component names").toVector();
          bc->set_tcc_names(tcc_names);

          // set the component indicies
          for (const auto& n : bc->tcc_names()) {
            bc->tcc_index().push_back(FindComponentNumber(n));
          }
          bcs_.push_back(bc);
        }
      }
    }

#ifdef ALQUIMIA_ENABLED
    // -- try geochemical conditions
    Teuchos::ParameterList& glist = plist_->sublist("boundary conditions").sublist("geochemical");

    for (Teuchos::ParameterList::ConstIterator it = glist.begin(); it != glist.end(); ++it) {
      std::string specname = it->first;
      Teuchos::ParameterList& spec = glist.sublist(specname);

      Teuchos::RCP<TransportBoundaryFunction_Alquimia_Units>
        bc = Teuchos::rcp(new TransportBoundaryFunction_Alquimia_Units(spec, mesh_,
                chem_pk_, chem_engine_));

      bc->set_conversion(1000.0, mol_dens_, true);
      std::vector<int>& tcc_index = bc->tcc_index();
      std::vector<std::string>& tcc_names = bc->tcc_names();

      for (int i = 0; i < tcc_names.size(); i++) {
        tcc_index.push_back(FindComponentNumber(tcc_names[i]));
      }

      bcs_.push_back(bc);
    }
#endif

  } else {
    if (vo_->getVerbLevel() > Teuchos::VERB_NONE) {
      *vo_->os() << vo_->color("yellow") << "No BCs were specified." << vo_->reset() << std::endl;
    }
  }

  // boundary conditions initialization
  time = t_physics_;
  VV_CheckInfluxBC();

  // source term initialization: so far only "concentration" is available.
  if (plist_->isSublist("source terms")) {
    PK_DomainFunctionFactory<TransportDomainFunction> factory(mesh_);
    Teuchos::ParameterList& conc_sources_list = plist_->sublist("source terms").sublist("concentration");

    for (const auto& it : conc_sources_list) {
      std::string name = it.first;
      if (conc_sources_list.isSublist(name)) {
        Teuchos::ParameterList& src_list = conc_sources_list.sublist(name);
        std::string src_type = src_list.get<std::string>("spatial distribution method", "none");

        if (src_type == "domain coupling") {
          // domain couplings are special -- they always work on all components
          Teuchos::RCP<TransportDomainFunction> src =
              factory.Create(src_list, "fields", AmanziMesh::CELL, Kxy);

          for (int i = 0; i < component_names_.size(); i++){
            src->tcc_names().push_back(component_names_[i]);
            src->tcc_index().push_back(i);
          }
          src->set_state(S_);
          srcs_.push_back(src);

        } else {
          Teuchos::RCP<TransportDomainFunction> src =
              factory.Create(src_list, "source function", AmanziMesh::CELL, Kxy);

          std::vector<std::string> tcc_names = src_list.get<Teuchos::Array<std::string>>("component names").toVector();
          src->set_tcc_names(tcc_names);

          // set the component indicies
          for (const auto& n : src->tcc_names()) {
            src->tcc_index().push_back(FindComponentNumber(n));
          }

          src->set_state(S_);
          srcs_.push_back(src);
        }
      }
    }

#ifdef ALQUIMIA_ENABLED
    // -- try geochemical conditions
    Teuchos::ParameterList& geochem_list = plist_->sublist("source terms").sublist("geochemical");

    for (auto it = geochem_list.begin(); it != geochem_list.end(); ++it) {
      std::string specname = it->first;
      Teuchos::ParameterList& spec = geochem_list.sublist(specname);

      Teuchos::RCP<TransportSourceFunction_Alquimia_Units>
          src = Teuchos::rcp(new TransportSourceFunction_Alquimia_Units(spec, mesh_, chem_pk_, chem_engine_));

      auto mass_src = S->GetFieldData(mass_src_key_)->ViewComponent("cell",false);
      src->set_conversion(-1000., mass_src, false);

      for (const auto& n : src->tcc_names()) {
        src->tcc_index().push_back(FindComponentNumber(n));
      }

      srcs_.push_back(src);
    }
#endif
  }

  // Temporarily Transport hosts Henry law.
  PrepareAirWaterPartitioning_();

  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    *vo_->os() << "Number of components: " << tcc->size() << std::endl
               << "cfl=" << cfl_ << " spatial/temporal discretization: "
               << spatial_disc_order << " " << temporal_disc_order << std::endl;
    *vo_->os() << vo_->color("green") << "Initalization of PK is complete."
               << vo_->reset() << std::endl << std::endl;
  }
}


/* ******************************************************************
* Initalized fields left by State and other PKs.
****************************************************************** */
void Transport_ATS::InitializeFields_(const Teuchos::Ptr<State>& S)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // set popular default values when flow PK is off
  if (S->HasField(saturation_key_)) {
    if (S->GetField(saturation_key_)->owner() == name_) {
      if (!S->GetField(saturation_key_, name_)->initialized()) {
        S->GetFieldData(saturation_key_, name_)->PutScalar(1.0);
        S->GetField(saturation_key_, name_)->set_initialized();

        if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
            *vo_->os() << "initilized saturation_liquid to value 1.0" << std::endl;
      }
      InitializeFieldFromField_(prev_saturation_key_, saturation_key_, S, true, true);

    } else {
      if (S->GetField(prev_saturation_key_)->owner() == name_) {
        if (!S->GetField(prev_saturation_key_, name_)->initialized()) {
          InitializeFieldFromField_(prev_saturation_key_, saturation_key_, S, true, true);
          S->GetField(prev_saturation_key_, name_)->set_initialized();

          if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
            *vo_->os() << "initilized prev_saturation_liquid from saturation" << std::endl;
        }
      }
    }
  }
  InitializeFieldFromField_(tcc_matrix_key_, tcc_key_, S, false, false);
  S->GetFieldData(solid_residue_mass_key_, name_)->PutScalar(0.0);
  S->GetField(solid_residue_mass_key_, name_)->set_initialized();
  S->GetField(conserve_qty_key_, name_)->set_initialized();
}


/* ****************************************************************
* Auxiliary initialization technique.
**************************************************************** */
void Transport_ATS::InitializeFieldFromField_(const std::string& field0,
                                                 const std::string& field1,
                                                 const Teuchos::Ptr<State>& S,
                                                 bool call_evaluator,
                                                 bool overwrite)
{
  if (S->HasField(field0)) {
    if (S->GetField(field0)->owner() == name_) {
      if ((!S->GetField(field0, name_)->initialized())||(overwrite)) {
        if (call_evaluator)
            S->GetFieldEvaluator(field1)->HasFieldChanged(S.ptr(), name_);

        const CompositeVector& f1 = *S->GetFieldData(field1);
        CompositeVector& f0 = *S->GetFieldData(field0, name_);

        double vmin0, vmax0, vavg0;
        double vmin1, vmax1, vavg1;

        f0 = f1;

        S->GetField(field0, name_)->set_initialized();
        if ((vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)&&(!overwrite)){
          *vo_->os() << "initiliazed " << field0 << " to " << field1 << std::endl;
        }
      }
    }
  }
}


/* *******************************************************************
* Estimation of the time step based on T.Barth (Lecture Notes
* presented at VKI Lecture Series 1994-05, Theorem 4.2.2.
* Routine must be called every time we update a flow field.
*
* Warning: Barth calculates influx, we calculate outflux. The methods
* are equivalent for divergence-free flows and gurantee EMP. Outflux
* takes into account sinks and sources but preserves only positivity
* of an advected mass.
* ***************************************************************** */
double Transport_ATS::StableTimeStep()
{
  S_next_->GetFieldData(flux_key_)->ScatterMasterToGhosted("face");

  flux_ = S_next_->GetFieldData(flux_key_)->ViewComponent("face", true);

  Teuchos::RCP<Epetra_Map> cell_map = Teuchos::rcp(new Epetra_Map(mesh_->cell_map(false)));
  IdentifyUpwindCells();

  tcc = S_inter_->GetFieldData(tcc_key_, name_);
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");

  // loop over faces and accumulate upwinding fluxes
  std::vector<double> total_outflux(ncells_wghost, 0.0);

  for (int f = 0; f < nfaces_wghost; f++) {
    int c = (*upwind_cell_)[f];
    if (c >= 0) {
      total_outflux[c] += fabs((*flux_)[0][f]);
    }
  }

  Sinks2TotalOutFlux(tcc_prev, total_outflux, 0, num_aqueous - 1);

  // modify estimate for other models
  // if (multiscale_porosity_) {
  //   const Epetra_MultiVector& wcm_prev = *S_next_->GetFieldData("prev_water_content_matrix")->ViewComponent("cell");
  //   const Epetra_MultiVector& wcm = *S_next_->GetFieldData("water_content_matrix")->ViewComponent("cell");

  //   double dtg = S_->final_time() - S_->initial_time();
  //   for (int c = 0; c < ncells_owned; ++c) {
  //     double flux_liquid = (wcm[0][c] - wcm_prev[0][c]) / dtg;
  //     msp_->second[(*msp_->first)[c]]->UpdateStabilityOutflux(flux_liquid, &total_outflux[c]);
  //   }
  // }

  // loop over cells and calculate minimal time step
  double vol, outflux, dt_cell;
  double ws_min_dt, outflux_min_dt;
  vol=0;
  dt_ = dt_cell = TRANSPORT_LARGE_TIME_STEP;
  int cmin_dt = 0;
  for (int c = 0; c < ncells_owned; c++) {
    outflux = total_outflux[c];

    if ((outflux > 0) && ((*ws_prev_)[0][c]>0) && ((*ws_)[0][c]>0) && ((*phi_)[0][c] > 0 )) {
      vol = mesh_->cell_volume(c);
      dt_cell = vol * (*mol_dens_)[0][c] * (*phi_)[0][c] * std::min( (*ws_prev_)[0][c], (*ws_)[0][c] ) / outflux;
    }
    if (dt_cell < dt_) {
      // *vo_->os()<<"Stable step: "<<flux_key_<<" cell "<<c<<" out "<<outflux<<"  dt= "<<dt_cell<<"\n";
      dt_ = dt_cell;
      cmin_dt = c;
      ws_min_dt = std::min( (*ws_prev_)[0][c], (*ws_)[0][c] );
      outflux_min_dt = total_outflux[c];
    }
  }

  if (spatial_disc_order == 2) dt_ /= 2;

  // communicate global time step
  double dt_tmp = dt_;
  const Epetra_Comm& comm = ws_prev_->Comm();
  comm.MinAll(&dt_tmp, &dt_, 1);

  // incorporate developers and CFL constraints
  dt_ = std::min(dt_, dt_debug_);
  dt_ *= cfl_;

  // print optional diagnostics using maximum cell id as the filter
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int cmin_dt_unique = (fabs(dt_tmp * cfl_ - dt_) < 1e-6 * dt_) ? cell_map->GID(cmin_dt) : -2;

    int cmin_dt_tmp = cmin_dt_unique;
    comm.MaxAll(&cmin_dt_tmp, &cmin_dt_unique, 1);
    int min_pid=-1;

    double tmp_package[6];

    if (cell_map->GID(cmin_dt) == cmin_dt_unique) {
      const AmanziGeometry::Point& p = mesh_->cell_centroid(cmin_dt);

      Teuchos::OSTab tab = vo_->getOSTab();
      min_pid = comm.MyPID();
      tmp_package[0] = ws_min_dt;
      tmp_package[1] = outflux_min_dt;
      tmp_package[2] = p[0];
      tmp_package[3] = p[1];
      if (p.dim() == 3) tmp_package[4] = p[2];
      else tmp_package[4] = 0.;
      tmp_package[5] = p.dim();

    }

    int min_pid_tmp = min_pid;
    comm.MaxAll(&min_pid_tmp, &min_pid, 1);

    comm.Broadcast(tmp_package, 6, min_pid);

    Teuchos::OSTab tab = vo_->getOSTab();

    *vo_->os() << "Stable time step "<<dt_<< " is computed at ("<< tmp_package[2]<<", " <<tmp_package[3];
    if (fabs(3 - tmp_package[5]) <1e-10) *vo_->os()<<", "<<tmp_package[4];
    *vo_->os() <<")"<<std::endl;

    *vo_->os() << "Stable time step "<<dt_<< " is limited by saturation/ponded_depth "<<tmp_package[0]<<" and "
	       << "output flux "<<tmp_package[1]<<std::endl;
  }
  return dt_;
}


/* *******************************************************************
* Estimate returns last time step unless it is zero.
******************************************************************* */
double Transport_ATS::get_dt()
{
  if (subcycling_) {
    return 1e+99;
  } else {
    StableTimeStep();
    return dt_;
  }
}


/* *******************************************************************
* MPC will call this function to advance the transport state.
* Efficient subcycling requires to calculate an intermediate state of
* saturation only once, which leads to a leap-frog-type algorithm.
******************************************************************* */
bool Transport_ATS::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool failed = false;
  double dt_MPC = t_new - t_old;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_LOW))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << S_inter_->time()
               << " t1 = " << S_next_->time() << " h = " << dt_MPC << std::endl
               << "----------------------------------------------------------------" << std::endl;

  flux_ = S_next_->GetFieldData(flux_key_)->ViewComponent("face", true);
  *flux_copy_ = *flux_; // copy flux vector from S_next_ to S_;

  ws_ = S_next_->GetFieldData(saturation_key_)->ViewComponent("cell", false);
  mol_dens_ = S_next_->GetFieldData(molar_density_key_)->ViewComponent("cell", false);
  solid_qty_ = S_next_->GetFieldData(solid_residue_mass_key_, name_)->ViewComponent("cell", false);

  // We use original tcc and make a copy of it later if needed.
  tcc = S_inter_->GetFieldData(tcc_key_, name_);
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");
  db_->WriteVector("tcc_old", tcc.ptr());

  // calculate stable time step
  double dt_shift = 0.0, dt_global = dt_MPC;
  double time = S_inter_->intermediate_time();
  if (time >= 0.0) {
    t_physics_ = time;
    dt_shift = time - S_inter_->initial_time();
    dt_global = S_inter_->final_time() - S_inter_->initial_time();
  }

  if (subcycling_)
    StableTimeStep();
  else
    dt_ = dt_MPC;
  double dt_stable = dt_;  // advance routines override dt_

  int interpolate_ws = 0;  // (dt_ < dt_global) ? 1 : 0;

  if ((t_old > S_inter_->initial_time())||(t_new < S_inter_->final_time())) interpolate_ws = 1;

  // start subcycling
  double dt_sum = 0.0;
  double dt_cycle;
  if (interpolate_ws) {
    dt_cycle = std::min(dt_stable, dt_MPC);
    InterpolateCellVector(*ws_prev_, *ws_, dt_shift, dt_global, *ws_subcycle_start);
    InterpolateCellVector(*mol_dens_prev_, *mol_dens_, dt_shift, dt_global, *mol_dens_subcycle_start);
    InterpolateCellVector(*ws_prev_, *ws_, dt_shift + dt_cycle, dt_global, *ws_subcycle_end);
    InterpolateCellVector(*mol_dens_prev_, *mol_dens_, dt_shift + dt_cycle, dt_global, *mol_dens_subcycle_end);
    ws_start = ws_subcycle_start;
    ws_end = ws_subcycle_end;
    mol_dens_start = mol_dens_subcycle_start;
    mol_dens_end = mol_dens_subcycle_end;
  } else {
    dt_cycle = dt_MPC;
    ws_start = ws_prev_;
    ws_end = ws_;
    mol_dens_start = mol_dens_prev_;
    mol_dens_end = mol_dens_;
  }

  db_->WriteVector("sat_old", S_inter_->GetFieldData(prev_saturation_key_).ptr());
  for (int c = 0; c < ncells_owned; c++) {
    double vol_phi_ws_den;
    vol_phi_ws_den = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_prev_)[0][c] * (*mol_dens_prev_)[0][c];
    for (int i=0; i<num_aqueous + num_gaseous; i++) {
      mass_solutes_stepstart_[i] = tcc_prev[i][c] * vol_phi_ws_den;
    }
  }

  int ncycles = 0, swap = 1;
  while (dt_sum < dt_MPC - 1e-6) {
    // update boundary conditions
    time = t_physics_ + dt_cycle / 2;
    for (int i = 0; i < bcs_.size(); i++){
      bcs_[i]->Compute(t_physics_, t_physics_+dt_cycle);
    }

    double dt_try = dt_MPC - dt_sum;
    double tol = 1e-10 * (dt_try + dt_stable);
    bool final_cycle = false;

    if (dt_try >= 2 * dt_stable) {
      dt_cycle = dt_stable;
    } else if (dt_try > dt_stable + tol) {
      dt_cycle = dt_try / 2;
    } else {
      dt_cycle = dt_try;
      final_cycle = true;
    }

    t_physics_ += dt_cycle;
    dt_sum += dt_cycle;

    if (interpolate_ws) {
      if (swap) {  // Initial water saturation is in 'start'.
        ws_start = ws_subcycle_start;
        ws_end = ws_subcycle_end;
        mol_dens_start = mol_dens_subcycle_start;
        mol_dens_end = mol_dens_subcycle_end;

        double dt_int = dt_sum + dt_shift;
        InterpolateCellVector(*ws_prev_, *ws_, dt_int, dt_global, *ws_subcycle_end);
        InterpolateCellVector(*mol_dens_prev_, *mol_dens_, dt_int, dt_global, *mol_dens_subcycle_end);
      } else {  // Initial water saturation is in 'end'.
        ws_start = ws_subcycle_end;
        ws_end = ws_subcycle_start;
        mol_dens_start = mol_dens_subcycle_end;
        mol_dens_end = mol_dens_subcycle_start;

        double dt_int = dt_sum + dt_shift;
        InterpolateCellVector(*ws_prev_, *ws_, dt_int, dt_global, *ws_subcycle_start);
        InterpolateCellVector(*mol_dens_prev_, *mol_dens_, dt_int, dt_global, *mol_dens_subcycle_start);
      }
      swap = 1 - swap;
    }

    if (spatial_disc_order == 1) {  // temporary solution (lipnikov@lanl.gov)
      AdvanceDonorUpwind(dt_cycle);
    } else if (spatial_disc_order == 2 && temporal_disc_order == 1) {
      AdvanceSecondOrderUpwindRK1(dt_cycle);
    } else if (spatial_disc_order == 2 && temporal_disc_order == 2) {
      AdvanceSecondOrderUpwindRK2(dt_cycle);
    }

    // add multiscale model
    if (multiscale_porosity_) {
      double t_int1 = t_old + dt_sum - dt_cycle;
      double t_int2 = t_old + dt_sum;
      AddMultiscalePorosity_(t_old, t_new, t_int1, t_int2);
    }

    if (! final_cycle) {  // rotate concentrations (we need new memory for tcc)
      tcc = Teuchos::RCP<CompositeVector>(new CompositeVector(*tcc_tmp));
    }

    ncycles++;
  }

  dt_ = dt_stable;  // restore the original time step (just in case)

  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", false);

  Advance_Dispersion_Diffusion(t_old, t_new);

  // optional Henry Law for the case of gas diffusion
  if (henry_law_) {
    MakeAirWaterPartitioning_();
  }

  // statistics output
  nsubcycles = ncycles;
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    *vo_->os() << ncycles << " sub-cycles, dt_stable=" << units_.OutputTime(dt_stable)
               << " [sec]  dt_MPC=" << units_.OutputTime(dt_MPC) << " [sec]" << std::endl;

    VV_PrintSoluteExtrema(tcc_next, dt_MPC);
  }
  return failed;
}


void Transport_ATS :: Advance_Dispersion_Diffusion(double t_old, double t_new)
{
  double dt_MPC = t_new - t_old;
  // We define tracer as the species #0 as calculate some statistics.
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", false);
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");
  int num_components = tcc_prev.NumVectors();

  bool flag_diffusion(false);
  for (int i = 0; i < 2; i++) {
    if (diffusion_phase_[i] != Teuchos::null) {
      if (diffusion_phase_[i]->values().size() != 0) flag_diffusion = true;
    }
  }
  if (flag_diffusion) {
    // no molecular diffusion if all tortuosities are zero.
    double tau(0.0);
    for (int i = 0; i < mat_properties_.size(); i++) {
      tau += mat_properties_[i]->tau[0] + mat_properties_[i]->tau[1];
    }
    if (tau == 0.0) flag_diffusion = false;
  }

  if (flag_dispersion_ || flag_diffusion) {
    // default boundary conditions (none inside domain and Neumann on its boundary)
    Teuchos::RCP<Operators::BCs> bc_dummy =
        Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
    auto& bc_model = bc_dummy->bc_model();
    auto& bc_value = bc_dummy->bc_value();
    PopulateBoundaryData(bc_model, bc_value, -1);

    // diffusion operator
    Teuchos::ParameterList& op_list = plist_->sublist("diffusion");
    op_list.set("inverse", plist_->sublist("inverse"));

    Operators::PDE_DiffusionFactory opfactory;
    Teuchos::RCP<Operators::PDE_Diffusion> op1 = opfactory.Create(op_list, mesh_, bc_dummy);
    op1->SetBCs(bc_dummy, bc_dummy);
    Teuchos::RCP<Operators::Operator> op = op1->global_operator();
    Teuchos::RCP<Operators::PDE_Accumulation> op2 =
        Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op));

    const CompositeVectorSpace& cvs = op1->global_operator()->DomainMap();
    CompositeVector sol(cvs), factor(cvs), factor0(cvs), source(cvs), zero(cvs);
    zero.PutScalar(0.0);

    // populate the dispersion operator (if any)
    if (flag_dispersion_) {
      CalculateDispersionTensor_(*flux_, *phi_, *ws_, *mol_dens_);
    }

    int phase, num_itrs(0);
    bool flag_op1(true);
    double md_change, md_old(0.0), md_new, residual(0.0);

    // Disperse and diffuse aqueous components
    for (int i = 0; i < num_aqueous; i++) {
      FindDiffusionValue(component_names_[i], &md_new, &phase);
      md_change = md_new - md_old;
      md_old = md_new;

      if (md_change != 0.0) {
        CalculateDiffusionTensor_(md_change, phase, *phi_, *ws_, *mol_dens_);
        flag_op1 = true;
      }

      // set initial guess
      Epetra_MultiVector& sol_cell = *sol.ViewComponent("cell");
      for (int c = 0; c < ncells_owned; c++) {
        sol_cell[0][c] = tcc_next[i][c];
      }
      if (sol.HasComponent("face")) {
        sol.ViewComponent("face")->PutScalar(0.0);
      }

      if (flag_op1) {
        op->Init();
        Teuchos::RCP<std::vector<WhetStone::Tensor> > Dptr = Teuchos::rcpFromRef(D_);
        op1->Setup(Dptr, Teuchos::null, Teuchos::null);
        op1->UpdateMatrices(Teuchos::null, Teuchos::null);

        // add accumulation term
        Epetra_MultiVector& fac = *factor.ViewComponent("cell");
        for (int c = 0; c < ncells_owned; c++) {
          fac[0][c] = (*phi_)[0][c] * (*ws_)[0][c] * (*mol_dens_)[0][c];
        }
        op2->AddAccumulationDelta(sol, factor, factor, dt_MPC, "cell");
        op1->ApplyBCs(true, true, true);

      } else {
        Epetra_MultiVector& rhs_cell = *op->rhs()->ViewComponent("cell");
        for (int c = 0; c < ncells_owned; c++) {
          double tmp = mesh_->cell_volume(c) * (*ws_)[0][c] * (*phi_)[0][c] * (*mol_dens_)[0][c]/ dt_MPC;
          rhs_cell[0][c] = tcc_next[i][c] * tmp;
        }
      }

      CompositeVector& rhs = *op->rhs();
      int ierr = op->ApplyInverse(rhs, sol);

      if (ierr < 0) {
        Errors::Message msg("TransportExplicit_PK solver failed with message: \"");
        msg << op->returned_code_string() << "\"";
        Exceptions::amanzi_throw(msg);
      }

      residual += op->residual();
      num_itrs += op->num_itrs();

      for (int c = 0; c < ncells_owned; c++) {
        tcc_next[i][c] = sol_cell[0][c];
      }
      if (sol.HasComponent("face")){
        if (tcc_tmp->HasComponent("boundary_face")){
          Epetra_MultiVector& tcc_tmp_bf = *tcc_tmp->ViewComponent("boundary_face",false);
          Epetra_MultiVector& sol_faces = *sol.ViewComponent("face",false);
          const Epetra_Map& vandalay_map = mesh_->exterior_face_map(false);
          const Epetra_Map& face_map = mesh_->face_map(false);
          int nbfaces = tcc_tmp_bf.MyLength();
          for (int bf=0; bf!=nbfaces; ++bf) {
            AmanziMesh::Entity_ID f = face_map.LID(vandalay_map.GID(bf));
            tcc_tmp_bf[i][bf] =  sol_faces[i][f];
          }
        }
      }
    }

    // Diffuse aqueous components. We ignore dispersion
    // tensor (D is reset). Inactive cells (s[c] = 1 and D_[c] = 0)
    // are treated with a hack of the accumulation term.
    D_.clear();
    md_old = 0.0;
    for (int i = num_aqueous; i < num_components; i++) {
      FindDiffusionValue(component_names_[i], &md_new, &phase);
      md_change = md_new - md_old;
      md_old = md_new;

      if (md_change != 0.0 || i == num_aqueous) {
        CalculateDiffusionTensor_(md_change, phase, *phi_, *ws_, *mol_dens_);
      }

      // set initial guess
      Epetra_MultiVector& sol_cell = *sol.ViewComponent("cell");
      for (int c = 0; c < ncells_owned; c++) {
        sol_cell[0][c] = tcc_next[i][c];
      }
      if (sol.HasComponent("face")) {
        sol.ViewComponent("face")->PutScalar(0.0);
      }

      op->Init();
      Teuchos::RCP<std::vector<WhetStone::Tensor> > Dptr = Teuchos::rcpFromRef(D_);
      op1->Setup(Dptr, Teuchos::null, Teuchos::null);
      op1->UpdateMatrices(Teuchos::null, Teuchos::null);

      // add boundary conditions and sources for gaseous components
      PopulateBoundaryData(bc_model, bc_value, i);

      Epetra_MultiVector& rhs_cell = *op->rhs()->ViewComponent("cell");
      ComputeAddSourceTerms(t_new, 1.0, rhs_cell, i, i);
      op1->ApplyBCs(true, true, true);

      // add accumulation term
      Epetra_MultiVector& fac1 = *factor.ViewComponent("cell");
      Epetra_MultiVector& fac0 = *factor0.ViewComponent("cell");

      for (int c = 0; c < ncells_owned; c++) {
        fac1[0][c] = (*phi_)[0][c] * (1.0 - (*ws_)[0][c]) * (*mol_dens_)[0][c];
        fac0[0][c] = (*phi_)[0][c] * (1.0 - (*ws_prev_)[0][c]) * (*mol_dens_prev_)[0][c];
        if ((*ws_)[0][c] == 1.0) fac1[0][c] = 1.0 * (*mol_dens_)[0][c];  // hack so far
      }
      op2->AddAccumulationDelta(sol, factor0, factor, dt_MPC, "cell");

      CompositeVector& rhs = *op->rhs();
      int ierr = op->ApplyInverse(rhs, sol);
      if (ierr < 0) {
        Errors::Message msg("Transport_PK solver failed with message: \"");
        msg << op->returned_code_string() << "\"";
        Exceptions::amanzi_throw(msg);
      }

      residual += op->residual();
      num_itrs += op->num_itrs();

      for (int c = 0; c < ncells_owned; c++) {
        tcc_next[i][c] = sol_cell[0][c];
      }
    }

    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "dispersion solver ||r||=" << residual / num_components
                 << " itrs=" << num_itrs / num_components << std::endl;
    }
  }

}


/* *******************************************************************
* Add multiscale porosity model on sub interval [t_int1, t_int2]:
*   d(VWC_f)/dt -= G_s, d(VWC_m) = G_s
*   G_s = G_w C^* + omega_s (C_f - C_m).
******************************************************************* */
void Transport_ATS::AddMultiscalePorosity_(
    double t_old, double t_new, double t_int1, double t_int2)
{
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell");
  Epetra_MultiVector& tcc_matrix =
     *S_inter_->GetFieldData("total_component_concentration_matrix", name_)->ViewComponent("cell");

  const Epetra_MultiVector& wcf_prev = *S_inter_->GetFieldData("prev_water_content")->ViewComponent("cell");
  const Epetra_MultiVector& wcf = *S_inter_->GetFieldData("water_content")->ViewComponent("cell");

  const Epetra_MultiVector& wcm_prev = *S_inter_->GetFieldData("prev_water_content_matrix")->ViewComponent("cell");
  const Epetra_MultiVector& wcm = *S_inter_->GetFieldData("water_content_matrix")->ViewComponent("cell");

  // multi-node matrix requires more input data
  const Epetra_MultiVector& phi_matrix = *S_inter_->GetFieldData("porosity_matrix")->ViewComponent("cell");

  int nnodes(1);
  Teuchos::RCP<Epetra_MultiVector> tcc_matrix_aux;
  if (S_inter_->HasField("total_component_concentration_matrix_aux")) {
    tcc_matrix_aux = S_inter_->GetFieldData("total_component_concentration_matrix_aux", name_)->ViewComponent("cell");
    nnodes = tcc_matrix_aux->NumVectors() + 1;
  }
  WhetStone::DenseVector tcc_m(nnodes);

  double flux_liquid, flux_solute, wcm0, wcm1, wcf0, wcf1;
  double dtg, dts, t1, t2, tmp0, tmp1, tfp0, tfp1, a, b;
  std::vector<AmanziMesh::Entity_ID> block;

  dtg = t_new - t_old;
  dts = t_int2 - t_int1;
  t1 = t_int1 - t_old;
  t2 = t_int2 - t_old;

  for (int c = 0; c < ncells_owned; ++c) {
    wcm0 = wcm_prev[0][c];
    wcm1 = wcm[0][c];
    flux_liquid = (wcm1 - wcm0) / dtg;

    wcf0 = wcf_prev[0][c];
    wcf1 = wcf[0][c];

    a = t1 / dtg;
    b = t2 / dtg;

    tfp0 = a * wcf1 + (1.0 - a) * wcf0;
    tfp1 = b * wcf1 + (1.0 - b) * wcf0;

    tmp0 = a * wcm1 + (1.0 - a) * wcm0;
    tmp1 = b * wcm1 + (1.0 - b) * wcm0;

    double phim = phi_matrix[0][c];

    for (int i = 0; i < num_aqueous; ++i) {
      tcc_m(0) = tcc_matrix[i][c];
      if (tcc_matrix_aux != Teuchos::null) {
        for (int n = 0; n < nnodes - 1; ++n)
          tcc_m(n + 1) = (*tcc_matrix_aux)[n][c];
      }

      flux_solute = msp_->second[(*msp_->first)[c]]->ComputeSoluteFlux(
          flux_liquid, tcc_next[i][c], tcc_m,
          i, dts, tfp0, tfp1, tmp0, tmp1, phim);

      tcc_matrix[i][c] = tcc_m(0);
      if (tcc_matrix_aux != Teuchos::null) {
        for (int n = 0; n < nnodes - 1; ++n)
          (*tcc_matrix_aux)[n][c] = tcc_m(n + 1);
      }
    }
  }
}


/* *******************************************************************
* Copy the advected tcc field to the state.
******************************************************************* */
void Transport_ATS::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  Teuchos::RCP<CompositeVector> tcc_vec_S = S->GetFieldData(tcc_key_, name_);
  *tcc_vec_S = *tcc_tmp;
  InitializeFieldFromField_(prev_saturation_key_, saturation_key_, S.ptr(), false, true);
  ChangedSolutionPK(S.ptr());

  // Copy to S_ as well
  //
  // WHY?  is this required for subcycling?  this used to be S, which was a
  // bug as it duplicates the above line of code. --etc
  Teuchos::RCP<CompositeVector> tcc_vec_Sold = S_->GetFieldData(tcc_key_, name_);
  *tcc_vec_Sold = *tcc_tmp;
}


/* *******************************************************************
 * A simple first-order transport method
 ****************************************************************** */
void Transport_ATS::AdvanceDonorUpwind(double dt_cycle)
{
  dt_ = dt_cycle;  // overwrite the maximum stable transport step
  mass_solutes_source_.assign(num_aqueous + num_gaseous, 0.0);
  mass_solutes_bc_.assign(num_aqueous + num_gaseous, 0.0);

  // populating next state of concentrations
  tcc->ScatterMasterToGhosted("cell");
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);

  // prepare conservative state in master and slave cells
  double mass_start = 0., tmp1, mass;

  // We advect only aqueous components.
  int num_advect = num_aqueous;
  int num_components = tcc_next.NumVectors();
  conserve_qty_->PutScalar(0.);

  for (int c = 0; c < ncells_owned; c++) {
    double vol_phi_ws_den = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_start)[0][c] * (*mol_dens_start)[0][c];
    (*conserve_qty_)[num_components+1][c] = vol_phi_ws_den;

    for (int i = 0; i < num_advect; i++){
      (*conserve_qty_)[i][c] = tcc_prev[i][c] * vol_phi_ws_den;

      if (dissolution_) {
        if (( (*ws_start)[0][c]  > water_tolerance_) && ((*solid_qty_)[i][c] > 0 )) {  // Dissolve solid residual into liquid
          double add_mass = std::min((*solid_qty_)[i][c], max_tcc_* vol_phi_ws_den - (*conserve_qty_)[i][c]);
          (*solid_qty_)[i][c] -= add_mass;
          (*conserve_qty_)[i][c] += add_mass;
        }
      }

      mass_start += (*conserve_qty_)[i][c];
    }
  }

  db_->WriteCellVector("cons (start)", *conserve_qty_);
  tmp1 = mass_start;
  mesh_->get_comm()->SumAll(&tmp1, &mass_start, 1);

  // advance all components at once
  for (int f = 0; f < nfaces_wghost; f++) {  // loop over master and slave faces
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];
    double u = fabs((*flux_)[0][f]);

    if (c1 >=0 && c1 < ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      for (int i = 0; i < num_advect; i++) {
        double tcc_flux = dt_ * u * tcc_prev[i][c1];
        (*conserve_qty_)[i][c1] -= tcc_flux;
        (*conserve_qty_)[i][c2] += tcc_flux;
      }
      (*conserve_qty_)[num_components+1][c1] -= dt_ * u;
      (*conserve_qty_)[num_components+1][c2] += dt_ * u;
    }
    else if (c1 >=0 && c1 < ncells_owned && (c2 >= ncells_owned || c2 < 0)) {
      for (int i = 0; i < num_advect; i++) {
        double tcc_flux = dt_ * u * tcc_prev[i][c1];
        (*conserve_qty_)[i][c1] -= tcc_flux;
        if (c2 < 0) mass_solutes_bc_[i] -= tcc_flux;
        //AmanziGeometry::Point normal = mesh_->face_normal(f);
      }
      (*conserve_qty_)[num_components+1][c1] -= dt_ * u;

    } else if (c1 >= ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      for (int i = 0; i < num_advect; i++) {
        double tcc_flux = dt_ * u * tcc_prev[i][c1];
        (*conserve_qty_)[i][c2] += tcc_flux;
      }
      (*conserve_qty_)[num_components+1][c2] += dt_ * u;

    } else if (c2 < 0 && c1 >= 0 && c1 < ncells_owned) {
      (*conserve_qty_)[num_components+1][c1] -= dt_ * u;

    } else if (c1 < 0 && c2 >= 0 && c2 < ncells_owned) {
      (*conserve_qty_)[num_components+1][c2] += dt_ * u;
    }
  }

  // loop over exterior boundary sets
  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
      int f = it->first;
      std::vector<double>& values = it->second;
      int c2 = (*downwind_cell_)[f];
      int c1 = (*upwind_cell_)[f];

      double u = fabs((*flux_)[0][f]);
      if (c2 >= 0) {
        for (int i = 0; i < ncomp; i++) {
          int k = tcc_index[i];
          if (k < num_advect) {
            double tcc_flux = dt_ * u * values[i];
            (*conserve_qty_)[k][c2] += tcc_flux;
            mass_solutes_bc_[k] += tcc_flux;
          }
        }
      }
    }
  }
  db_->WriteCellVector("cons (adv)", *conserve_qty_);

  // process external sources
  if (srcs_.size() != 0) {
    double time = t_physics_;
    ComputeAddSourceTerms(time, dt_, *conserve_qty_, 0, num_advect - 1);
  }
  db_->WriteCellVector("cons (src)", *conserve_qty_);

  // recover concentration from new conservative state
  for (int c = 0; c < ncells_owned; c++) {
    double water_new = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_end)[0][c] * (*mol_dens_end)[0][c];
    double water_sink = (*conserve_qty_)[num_components][c]; // water at the new time + outgoing domain coupling source
    double water_total = water_new + water_sink;
    AMANZI_ASSERT(water_total >= water_new);
    (*conserve_qty_)[num_components][c] = water_total;

    // if (std::abs((*conserve_qty_)[num_components+1][c] - water_total) > water_tolerance_
    //     && vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    //   *vo_->os() << "Water balance error (cell " << c << "): " << std::endl
    //              << "  water_old + advected = " << (*conserve_qty_)[num_components+1][c] << std::endl
    //              << "  water_sink = " << water_sink << std::endl
    //              << "  water_new = " << water_new << std::endl;
    // }

    for (int i = 0; i < num_advect; i++) {
      if (water_new > water_tolerance_ && (*conserve_qty_)[i][c] > 0) {
        // there is both water and stuff present at the new time
        // this is stuff at the new time + stuff leaving through the domain coupling, divided by water of both
        tcc_next[i][c] = (*conserve_qty_)[i][c] / water_total;
      } else if (water_sink > water_tolerance_ && (*conserve_qty_)[i][c] > 0) {
        // there is water and stuff leaving through the domain coupling, but it all leaves (none at the new time)
        tcc_next[i][c] = 0.;
      } else {
        // there is no water leaving, and no water at the new time.  Change any stuff into solid
        (*solid_qty_)[i][c] += std::max((*conserve_qty_)[i][c], 0.);
        (*conserve_qty_)[i][c] = 0.;
        tcc_next[i][c] = 0.;
      }
    }
  }
  db_->WriteCellVector("tcc_new", tcc_next);

  double mass_final = 0;
  for (int c = 0; c < ncells_owned; c++) {
    for (int i = 0; i < num_advect; i++){
      mass_final += (*conserve_qty_)[i][c];
    }
  }

  tmp1 = mass_final;
  mesh_->get_comm()->SumAll(&tmp1, &mass_final, 1);

  // update mass balance
  for (int i = 0; i < mass_solutes_exact_.size(); i++) {
    mass_solutes_exact_[i] += mass_solutes_source_[i] * dt_;
  }

  if (internal_tests) {
    VV_CheckGEDproperty(*tcc_tmp->ViewComponent("cell"));
  }
}


/* *******************************************************************
 * We have to advance each component independently due to different
 * reconstructions. We use tcc when only owned data are needed and
 * tcc_next when owned and ghost data. This is a special routine for
 * transient flow and uses first-order time integrator.
 ****************************************************************** */
void Transport_ATS::AdvanceSecondOrderUpwindRK1(double dt_cycle)
{
  dt_ = dt_cycle;  // overwrite the maximum stable transport step
  mass_solutes_source_.assign(num_aqueous + num_gaseous, 0.0);

  // work memory
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);

  // distribute vector of concentrations
  S_inter_->GetFieldData(tcc_key_)->ScatterMasterToGhosted("cell");
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);

  // Epetra_Vector ws_ratio(Copy, *ws_start, 0);
  // for (int c = 0; c < ncells_owned; c++){
  //   double vol_phi_ws_den_end = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_end)[0][c] * (*mol_dens_end)[0][c];
  //   if (vol_phi_ws_den_end > water_tolerance_)  {
  //     double vol_phi_ws_den_start = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_start)[0][c] * (*mol_dens_start)[0][c];
  //     if (vol_phi_ws_den_start > water_tolerance_){
  //       ws_ratio[c] = ( (*ws_start)[0][c] * (*mol_dens_start)[0][c] )
  //                   / ( (*ws_end)[0][c]   * (*mol_dens_end)[0][c]   );
  //     } else {
  //       ws_ratio[c] = 1;
  //     }
  //   }
  //   else  ws_ratio[c]=0.;
  // }

  // We advect only aqueous components.
  int num_advect = num_aqueous;
  int num_components = tcc_next.NumVectors();
  conserve_qty_->PutScalar(0.);

  // prepopulate with initial water for better debugging
  for (int c = 0; c < ncells_owned; c++) {
    double vol_phi_ws_den_start = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_start)[0][c] * (*mol_dens_start)[0][c];
    (*conserve_qty_)[num_components+1][c] = vol_phi_ws_den_start;
  }

  for (int i = 0; i < num_advect; i++) {
    current_component_ = i;  // needed by BJ
    double T = t_physics_;
    Epetra_Vector*& component = tcc_prev(i);
    FunctionalTimeDerivative(T, *component, *(*conserve_qty_)(i));
  }
  db_->WriteCellVector("cons (time_deriv)", *conserve_qty_);

  // calculate the new conc
  for (int c = 0; c < ncells_owned; c++) {
    double water_old = (*conserve_qty_)[num_components+1][c];
    double water_new = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_end)[0][c] * (*mol_dens_end)[0][c];
    double water_sink = (*conserve_qty_)[num_components][c];
    double water_total = water_sink + water_new;
    (*conserve_qty_)[num_components][c] = water_total;

    for (int i=0; i!=num_components; ++i) {
      double cons_qty = (tcc_prev[i][c] + dt_ * (*conserve_qty_)[i][c]) * water_old;
      (*conserve_qty_)[i][c] = cons_qty;
      if (water_new > water_tolerance_ && cons_qty > 0) {
        // there is both water and stuff present at the new time
        // this is stuff at the new time + stuff leaving through the domain coupling, divided by water of both
        tcc_next[i][c] = cons_qty / water_total;
      } else if (water_sink > water_tolerance_ && cons_qty > 0) {
        // there is water and stuff leaving through the domain coupling, but it all leaves (none at the new time)
        tcc_next[i][c] = 0.;
      } else {
        // there is no water leaving, and no water at the new time.  Change any stuff into solid
        (*solid_qty_)[i][c] += std::max(cons_qty, 0.);
        (*conserve_qty_)[i][c] = 0.;
        tcc_next[i][c] = 0.;
      }
    }
  }
  db_->WriteCellVector("tcc_new", tcc_next);

  // update mass balance
  for (int i = 0; i < num_aqueous + num_gaseous; i++) {
    mass_solutes_exact_[i] += mass_solutes_source_[i] * dt_;
  }

  if (internal_tests) {
    VV_CheckGEDproperty(*tcc_tmp->ViewComponent("cell"));
  }
}


/* *******************************************************************
 * We have to advance each component independently due to different
 * reconstructions. This is a special routine for transient flow and
 * uses second-order predictor-corrector time integrator.
 ****************************************************************** */
void Transport_ATS::AdvanceSecondOrderUpwindRK2(double dt_cycle)
{
  dt_ = dt_cycle;  // overwrite the maximum stable transport step
  mass_solutes_source_.assign(num_aqueous + num_gaseous, 0.0);

  // work memory
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  Epetra_Vector f_component(cmap_wghost);//,  f_component2(cmap_wghost);

  // distribute old vector of concentrations
  S_inter_->GetFieldData(tcc_key_)->ScatterMasterToGhosted("cell");
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);

  Epetra_Vector ws_ratio(Copy, *ws_start, 0);
  for (int c = 0; c < ncells_owned; c++){
    if ((*ws_end)[0][c] > 1e-10)  {
      if ((*ws_start)[0][c] > 1e-10){
        ws_ratio[c] = ( (*ws_start)[0][c] * (*mol_dens_start)[0][c] )
                    / ( (*ws_end)[0][c]   * (*mol_dens_end)[0][c]   );
      }else{
        ws_ratio[c] = 1;
      }
    }
    else  ws_ratio[c]=0.;
  }

  // We advect only aqueous components.
  int num_advect = num_aqueous;

  // predictor step
  for (int i = 0; i < num_advect; i++) {
    current_component_ = i;  // needed by BJ

    double T = t_physics_;
    Epetra_Vector*& component = tcc_prev(i);
    FunctionalTimeDerivative(T, *component, f_component);

    for (int c = 0; c < ncells_owned; c++) {
      tcc_next[i][c] = (tcc_prev[i][c] + dt_ * f_component[c]) * ws_ratio[c];
      //if (tcc_next[i][c] < 0) tcc_next[i][c] = 0.;

    }
  }

  tcc_tmp->ScatterMasterToGhosted("cell");

  //if (domain_ == "surface") {
  //*vo_->os()<<"after predictor ToTaL "<<domain_<<" :"<<std::setprecision(10)<<ComputeSolute( tcc_next, 0)<<"\n";
  //}

  // corrector step
  for (int i = 0; i < num_advect; i++) {
    current_component_ = i;  // needed by BJ

    double T = t_physics_;
    Epetra_Vector*& component = tcc_next(i);
    FunctionalTimeDerivative(T, *component, f_component);

    for (int c = 0; c < ncells_owned; c++) {
      double value = (tcc_prev[i][c] + dt_ * f_component[c]) * ws_ratio[c];
      tcc_next[i][c] = (tcc_next[i][c] + value) / 2;
      if (tcc_next[i][c] < 0){
        double vol_phi_ws_den = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_end)[0][c] * (*mol_dens_end)[0][c];
        (*solid_qty_)[i][c] += abs(tcc_next[i][c])*vol_phi_ws_den;
        tcc_next[i][c] = 0.;
      }

    }
  }

  // f_component2.Update(-1, f_component, 1.);
  // double diff_norm;
  // f_component2.NormInf(&diff_norm);
  // *vo_->os()<<domain_<<" difference "<<diff_norm<<"\n";
  //if (domain_ == "surface") {
  //*vo_->os()<<"after corrector ToTaL "<<domain_<<" :"<<std::setprecision(10)<<ComputeSolute( tcc_next, 0)<<"\n";
  //}

  // update mass balance
  for (int i = 0; i < num_aqueous + num_gaseous; i++) {
    mass_solutes_exact_[i] += mass_solutes_source_[i] * dt_ / 2;
  }

  if (internal_tests) {
    VV_CheckGEDproperty(*tcc_tmp->ViewComponent("cell"));
  }

}


/* *******************************************************************
* Advance each component independently due to different field
* reconstructions. This routine uses generic explicit time integrator.
******************************************************************* */
// void Transport_ATS::AdvanceSecondOrderUpwindRKn(double dt_cycle)
// {
//   dt_ = dt_cycle;  // overwrite the maximum stable transport step

//   S_inter_->GetFieldData("total_component_concentration")->ScatterMasterToGhosted("cell");
//   Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
//   Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);

//   // define time integration method
//   auto ti_method = Explicit_TI::forward_euler;
//   if (temporal_disc_order == 2) {
//     ti_method = Explicit_TI::heun_euler;
//   } else if (temporal_disc_order == 3) {
//     ti_method = Explicit_TI::kutta_3rd_order;
//   } else if (temporal_disc_order == 3) {
//     ti_method = Explicit_TI::runge_kutta_4th_order;
//   }

//   // We interpolate ws using dt which becomes local time.
//   double T = 0.0;
//   // We advect only aqueous components.
//   int ncomponents = num_aqueous;

//   for (int i = 0; i < ncomponents; i++) {
//     current_component_ = i;  // it is needed in BJ called inside RK:fun

//     Epetra_Vector*& component_prev = tcc_prev(i);
//     Epetra_Vector*& component_next = tcc_next(i);

//     Explicit_TI::RK<Epetra_Vector> TVD_RK(*this, ti_method, *component_prev);
//     TVD_RK.TimeStep(T, dt_, *component_prev, *component_next);
//   }
// }



/* ******************************************************************
* Computes source and sink terms and adds them to vector tcc.
* Returns mass rate for the tracer.
* The routine treats two cases of tcc with one and all components.
****************************************************************** */
void Transport_ATS::ComputeAddSourceTerms(double tp, double dtp,
                                         Epetra_MultiVector& cons_qty, int n0, int n1)
{
  int num_vectors = cons_qty.NumVectors();
  int nsrcs = srcs_.size();

  for (int m = 0; m < nsrcs; m++) {
    double t0 = tp - dtp;
    srcs_[m]->Compute(t0, tp);

    std::vector<int> tcc_index = srcs_[m]->tcc_index();

    for (auto it = srcs_[m]->begin(); it != srcs_[m]->end(); ++it) {
      int c = it->first;
      std::vector<double>& values = it->second;

      if (c >= ncells_owned) continue;

      if (srcs_[m]->name() == "domain coupling" && n0 == 0) {
        (*conserve_qty_)[num_vectors-2][c] += values[num_vectors-2];
      }

      for (int k = 0; k < tcc_index.size(); ++k) {
        int i = tcc_index[k];
        if (i < n0 || i > n1) continue;

        int imap = i;
        if (num_vectors == 1) imap = 0;

        double value = mesh_->cell_volume(c) * values[k];
        cons_qty[imap][c] += dtp * value;
        mass_solutes_source_[i] += value;
      }
    }
  }
}


void Transport_ATS::Sinks2TotalOutFlux(Epetra_MultiVector& tcc_c,
                                          std::vector<double>& total_outflux, int n0, int n1){

  std::vector<double> sink_add(ncells_wghost, 0.0);
  //Assumption that there is only one sink per component per cell
  double t0 = S_inter_->intermediate_time();
  int num_vectors = tcc_c.NumVectors();
  int nsrcs = srcs_.size();
  Key coupled_flux_key = "surface-surface_subsurface_flux";

  for (int m = 0; m < nsrcs; m++) {
    srcs_[m]->Compute(t0, t0);
    std::vector<int> index = srcs_[m]->tcc_index();

    for (auto it = srcs_[m]->begin(); it != srcs_[m]->end(); ++it) {
      int c = it->first;
      std::vector<double>& values = it->second;

      double val = 0;
      for (int k = 0; k < index.size(); ++k) {
        int i = index[k];
        if (i < n0 || i > n1) continue;

        int imap = i;
        if (num_vectors == 1) imap = 0;

        if ((values[k] < 0) && (tcc_c[imap][c] > 1e-16)) {
          if (srcs_[m]->name() == "domain coupling") {
            //val = std::max(val, fabs(values[k])/tcc_c[imap][c]);
            //val = std::max(val, fabs(values[k]));
            const Epetra_MultiVector& flux_interface_ =
              *S_next_->GetFieldData(coupled_flux_key)->ViewComponent("cell", false);
            val = std::max(val, fabs(flux_interface_[0][c]));
          }
        }
      }
      sink_add[c] = std::max(sink_add[c], val);
    }
  }

  for (int c=0;c<ncells_wghost; c++) total_outflux[c] += sink_add[c];
}



/* *******************************************************************
* Populates operators' boundary data for given component.
* Returns true if at least one face was populated.
******************************************************************* */
bool Transport_ATS::PopulateBoundaryData(
    std::vector<int>& bc_model, std::vector<double>& bc_value, int component)
{
  bool flag = false;

  for (int i = 0; i < bc_model.size(); i++) {
    bc_model[i] = Operators::OPERATOR_BC_NONE;
    bc_value[i] = 0.0;
  }

  AmanziMesh::Entity_ID_List cells;
  for (int f = 0; f < nfaces_wghost; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    if (cells.size() == 1) bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
  }

  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
      int f = it->first;
      std::vector<double>& values = it->second;
      for (int i = 0; i < ncomp; i++) {
        int k = tcc_index[i];
        if (k == component) {
          bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value[f] = values[i];
          flag = true;
        }
      }
    }
  }

  return flag;
}


/* *******************************************************************
* Identify flux direction based on orientation of the face normal
* and sign of the  Darcy velocity.
******************************************************************* */
void Transport_ATS::IdentifyUpwindCells()
{
  for (int f = 0; f < nfaces_wghost; f++) {
    (*upwind_cell_)[f] = -1;  // negative value indicates boundary
    (*downwind_cell_)[f] = -1;
  }

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

    for (int i = 0; i < faces.size(); i++) {
      int f = faces[i];
      double tmp = (*flux_)[0][f] * dirs[i];
      if (tmp > 0.0) {
        (*upwind_cell_)[f] = c;
      } else if (tmp < 0.0) {
        (*downwind_cell_)[f] = c;
      } else if (dirs[i] > 0) {
        (*upwind_cell_)[f] = c;
      } else {
        (*downwind_cell_)[f] = c;
      }
    }
  }
}


void Transport_ATS::ComputeVolumeDarcyFlux(Teuchos::RCP<const Epetra_MultiVector> flux,
                                              Teuchos::RCP<const Epetra_MultiVector> molar_density,
                                              Teuchos::RCP<Epetra_MultiVector>& vol_darcy_flux)
{
  AmanziMesh::Entity_ID_List cells;

  for (int f = 0; f < nfaces_wghost ; f++){
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    double n_liq=0.;
    for (int c=0; c<cells.size();c++) n_liq += (*molar_density)[0][c];
    n_liq /= cells.size();
    if (n_liq > 0) (*vol_darcy_flux)[0][f] = (*flux_)[0][f]/n_liq;
    else (*vol_darcy_flux)[0][f] = 0.;
  }
}


/* *******************************************************************
* Interpolate linearly in time between two values v0 and v1. The time
* is measuared relative to value v0; so that v1 is at time dt. The
* interpolated data are at time dt_int.
******************************************************************* */
void Transport_ATS::InterpolateCellVector(
    const Epetra_MultiVector& v0, const Epetra_MultiVector& v1,
    double dt_int, double dt, Epetra_MultiVector& v_int)
{
  double a = dt_int / dt;
  double b = 1.0 - a;
  v_int.Update(b, v0, a, v1, 0.);
}

}  // namespace Transport
}  // namespace Amanzi

