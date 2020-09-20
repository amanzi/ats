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
                           const Teuchos::RCP<TreeVector>& soln) :
  S_(S),
  soln_(soln)
{
  name_ = Keys::cleanPListName(pk_tree.name());

  // Create miscaleneous lists.
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);

  tp_list_ = Teuchos::sublist(pk_list, name_, true);

  if (tp_list_->isParameter("component names")) {
    component_names_ = tp_list_->get<Teuchos::Array<std::string> >("component names").toVector();
    mol_masses_ = tp_list_->get<Teuchos::Array<double> >("component molar masses").toVector();
  } else {
    Errors::Message msg("Transport PK: parameter component names is missing.");
    Exceptions::amanzi_throw(msg);
  }

  subcycling_ = tp_list_->get<bool>("transport subcycling", false);
  special_source_ = tp_list_->get<bool>("special source", false);
  water_source_in_meters_ = tp_list_->get<bool>("water source in meters", true);

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
  passwd_ = "state";  // owner's password

  domain_name_ = tp_list_->get<std::string>("domain name", "domain");

  saturation_key_ = Keys::readKey(*tp_list_, domain_name_, "saturation liquid", "saturation_liquid");
  prev_saturation_key_ = Keys::readKey(*tp_list_, domain_name_, "previous saturation liquid", "prev_saturation_liquid");
  flux_key_ = Keys::readKey(*tp_list_, domain_name_, "mass flux", "mass_flux");
  permeability_key_ = Keys::readKey(*tp_list_, domain_name_, "permeability", "permeability");
  tcc_key_ = Keys::readKey(*tp_list_, domain_name_, "concentration", "total_component_concentration");
  porosity_key_ = Keys::readKey(*tp_list_, domain_name_, "porosity", "porosity");
  molar_density_key_ = Keys::readKey(*tp_list_, domain_name_, "molar density", "molar_density_liquid");
  tcc_matrix_key_ = Keys::readKey(*tp_list_, domain_name_, "tcc matrix", "total_component_concentration_matrix");
  solid_residue_mass_key_ =  Keys::readKey(*tp_list_, domain_name_, "solid residue", "solid_residue_mass");
  mass_src_key_ = Keys::readKey(*tp_list_, domain_name_, "mass source", "mass_source");
  water_content_key_ = Keys::getKey(domain_, "water_content");

  water_tolerance_ = tp_list_->get<double>("water tolerance", 1e-6);
  dissolution_ = tp_list_->get<bool>("allow dissolution", false);
  max_tcc_ = tp_list_->get<double>("maximum concentration", 0.9);


  mesh_ = S->GetMesh(domain_name_);
  dim = mesh_->space_dimension();

  // cross-coupling of PKs
  Teuchos::RCP<Teuchos::ParameterList> physical_models =
      Teuchos::sublist(tp_list_, "physical models and assumptions");
  bool abs_perm = physical_models->get<bool>("permeability field is required", false);
  std::string multiscale_model = physical_models->get<std::string>("multiscale model", "single porosity");

  // require state fields
  if (abs_perm) {
    S->RequireField(permeability_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, dim);
    S->RequireFieldEvaluator(permeability_key_);
  }

  if (!S->HasField(flux_key_)){
    S->RequireField(flux_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);
    S->RequireFieldEvaluator(flux_key_);
  }

  S->RequireField(saturation_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(saturation_key_);

  // prev_sat does not have an evaluator, this is managed by hand.  not sure why
  S->RequireField(prev_saturation_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(prev_saturation_key_, passwd_)->set_io_vis(false);

  if (!S->HasField(porosity_key_)){
    S->RequireField(porosity_key_, porosity_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(porosity_key_);
  }

  int ncomponents = component_names_.size();
  std::vector<std::vector<std::string> > subfield_names(1);
  subfield_names[0] = component_names_;

  if (!S->HasField(solid_residue_mass_key_)){
    S->RequireField(solid_residue_mass_key_,  passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, ncomponents);
  }

  if (!S->HasField(molar_density_key_)){
    S->RequireField(molar_density_key_, molar_density_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(molar_density_key_);
  }

  if (!S->HasField(mass_src_key_) && tp_list_->sublist("source terms").isSublist("geochemical")){
    S->RequireField(mass_src_key_, mass_src_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(mass_src_key_);
  }

  // require state fields when Transport PK is on
  if (component_names_.size() == 0) {
    Errors::Message msg;
    msg << "Transport PK: list of solutes is empty.\n";
    Exceptions::amanzi_throw(msg);
  }

  S->RequireField(tcc_key_, passwd_, subfield_names)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, ncomponents);

  if (!S->HasField(solid_residue_mass_key_)){
    S->RequireField(solid_residue_mass_key_,  passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, ncomponents);
  }

  // require multiscale fields
  multiscale_porosity_ = false;
  if (multiscale_model == "dual porosity") {
    multiscale_porosity_ = true;
    Teuchos::RCP<Teuchos::ParameterList>
        msp_list = Teuchos::sublist(tp_list_, "multiscale models", true);
    msp_ = CreateMultiscaleTransportPorosityPartition(mesh_, msp_list);

    std::vector<std::vector<std::string> > subfield_names(1);
    subfield_names[0] = component_names_;

    S->RequireField(tcc_matrix_key_, passwd_, subfield_names)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, ncomponents);
  }
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

  if (tp_list_->isSublist("initial conditions")) {
    S->GetField(tcc_key_,passwd_)->Initialize(tp_list_->sublist("initial conditions"));
  }

  internal_tests = 0;
  tests_tolerance = TRANSPORT_CONCENTRATION_OVERSHOOT;

  bc_scaling = 0.0;

  // Create verbosity object.
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = tp_list_->sublist("verbose object");
  vo_ =  Teuchos::rcp(new VerboseObject("TransportPK", vlist));

  Teuchos::OSTab tab = vo_->getOSTab();
  MyPID = mesh_->get_comm()->MyPID();

  // initialize missed fields
  InitializeFields_(S);

  //create copies
  S->RequireFieldCopy(tcc_key_, "subcycling", passwd_);
  tcc_tmp = S->GetField(tcc_key_, passwd_)->GetCopy("subcycling", passwd_)->GetFieldData();

  if (special_source_){
    S->RequireFieldCopy(tcc_key_, "with source", passwd_);
    tcc_w_src = S->GetField(tcc_key_, passwd_)->GetCopy("with source", passwd_)->GetFieldData();
  }

  S->RequireFieldCopy(saturation_key_, "subcycle_start", passwd_);
  ws_subcycle_start = S->GetFieldCopyData(saturation_key_, "subcycle_start",passwd_)
    ->ViewComponent("cell");
  S->RequireFieldCopy(saturation_key_, "subcycle_end", passwd_);
  ws_subcycle_end = S->GetFieldCopyData(saturation_key_, "subcycle_end", passwd_)
    ->ViewComponent("cell");
  S->RequireFieldCopy(molar_density_key_, "subcycle_start", passwd_);
  mol_dens_subcycle_start = S->GetFieldCopyData(molar_density_key_, "subcycle_start",passwd_)->ViewComponent("cell");
  S->RequireFieldCopy(molar_density_key_, "subcycle_end", passwd_);
  mol_dens_subcycle_end = S->GetFieldCopyData(molar_density_key_, "subcycle_end", passwd_)->ViewComponent("cell");

  S->RequireFieldCopy(flux_key_, "next_timestep", passwd_);
  flux_copy_ = S->GetFieldCopyData(flux_key_,  "next_timestep", passwd_)->ViewComponent("face", true);
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

  tcc = S->GetFieldData(tcc_key_, passwd_);

  flux_ = S->GetFieldData(flux_key_)->ViewComponent("face", true);
  solid_qty_ = S->GetFieldData(solid_residue_mass_key_, passwd_)->ViewComponent("cell", false);

  //create vector of conserved quatities
  conserve_qty_ = Teuchos::rcp(new Epetra_MultiVector(*(S->GetFieldData(tcc_key_)->ViewComponent("cell", true))));

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
  if (tp_list_->isSublist("material properties")) {
    Teuchos::RCP<Teuchos::ParameterList>
        mdm_list = Teuchos::sublist(tp_list_, "material properties");
    mdm_ = CreateMDMPartition(mesh_, mdm_list, flag_dispersion_);
    if (flag_dispersion_) CalculateAxiSymmetryDirection();
  }

  // create boundary conditions
  if (tp_list_->isSublist("boundary conditions")) {
    // -- try tracer-type conditions
    PK_DomainFunctionFactory<TransportDomainFunction> factory(mesh_);
    Teuchos::ParameterList& clist = tp_list_->sublist("boundary conditions").sublist("concentration");

    for (Teuchos::ParameterList::ConstIterator it = clist.begin(); it != clist.end(); ++it) {
      std::string name = it->first;
      if (clist.isSublist(name)) {
        Teuchos::ParameterList& bc_list = clist.sublist(name);
        if (name == "coupling") {
          Teuchos::ParameterList::ConstIterator it1 = bc_list.begin();
          std::string specname = it1->first;
          Teuchos::ParameterList& spec = bc_list.sublist(specname);
          Teuchos::RCP<TransportDomainFunction>
            bc = factory.Create(spec, "boundary concentration", AmanziMesh::FACE, Kxy);

          for (int i = 0; i < component_names_.size(); i++){
            bc->tcc_names().push_back(component_names_[i]);
            bc->tcc_index().push_back(i);
          }
          bc->set_state(S_);
          bcs_.push_back(bc);
        } else if (name == "subgrid") {
          Teuchos::ParameterList::ConstIterator it1 = bc_list.begin();
          std::string specname = it1->first;
          Teuchos::ParameterList& spec = bc_list.sublist(specname);
          Teuchos::Array<std::string> regions(1, domain_name_);

          std::size_t last_of = domain_name_.find_last_of("_");
          AMANZI_ASSERT(last_of != std::string::npos);
          int gid = std::stoi(domain_name_.substr(last_of+1, domain_name_.size()));
          spec.set("entity_gid_out", gid);
          Teuchos::RCP<TransportDomainFunction>
            bc = factory.Create(spec, "boundary concentration", AmanziMesh::FACE, Kxy);

          for (int i = 0; i < component_names_.size(); i++){
            bc->tcc_names().push_back(component_names_[i]);
            bc->tcc_index().push_back(i);
          }
          bc->set_state(S_);
          bcs_.push_back(bc);

        } else {
          for (Teuchos::ParameterList::ConstIterator it1 = bc_list.begin(); it1 != bc_list.end(); ++it1) {
            std::string specname = it1->first;
            Teuchos::ParameterList& spec = bc_list.sublist(specname);
            Teuchos::RCP<TransportDomainFunction>
              bc = factory.Create(spec, "boundary concentration", AmanziMesh::FACE, Kxy);

            std::vector<int>& tcc_index = bc->tcc_index();
            std::vector<std::string>& tcc_names = bc->tcc_names();
            bc->set_state(S_);

            tcc_names.push_back(name);
            tcc_index.push_back(FindComponentNumber(name));

            bcs_.push_back(bc);
          }
        }
      }
    }

#ifdef ALQUIMIA_ENABLED
    // -- try geochemical conditions
    Teuchos::ParameterList& glist = tp_list_->sublist("boundary conditions").sublist("geochemical");

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
  // for (int i = 0; i < bcs_.size(); i++) {
  //   bcs_[i]->Compute(time, time);
  // }

  VV_CheckInfluxBC();

  // source term initialization: so far only "concentration" is available.
  if (tp_list_->isSublist("source terms")) {
    PK_DomainFunctionFactory<TransportDomainFunction> factory(mesh_);
    //if (domain_name_ == "domain")  PKUtils_CalculatePermeabilityFactorInWell(S_.ptr(), Kxy);

    Teuchos::ParameterList& clist = tp_list_->sublist("source terms").sublist("concentration");
    for (Teuchos::ParameterList::ConstIterator it = clist.begin(); it != clist.end(); ++it) {
      std::string name = it->first;
      if (clist.isSublist(name)) {
        Teuchos::ParameterList& src_list = clist.sublist(name);
        if (name=="coupling") {
          Teuchos::ParameterList::ConstIterator it1 = src_list.begin();
          std::string specname = it1->first;
          Teuchos::ParameterList& spec = src_list.sublist(specname);
          Teuchos::RCP<TransportDomainFunction> src =
              factory.Create(spec, "sink", AmanziMesh::CELL, Kxy);

          for (int i = 0; i < component_names_.size(); i++){
            src->tcc_names().push_back(component_names_[i]);
            src->tcc_index().push_back(i);
          }
          src->set_state(S_);
          srcs_.push_back(src);

        } else {
          for (Teuchos::ParameterList::ConstIterator it1 = src_list.begin(); it1 != src_list.end(); ++it1) {
            std::string specname = it1->first;
            Teuchos::ParameterList& spec = src_list.sublist(specname);
            Teuchos::RCP<TransportDomainFunction> src =
              factory.Create(spec, "sink", AmanziMesh::CELL, Kxy);
            std::vector<std::string> tcc = spec.get<Teuchos::Array<std::string>>("component names").toVector();
            src->set_tcc_names(tcc);
            // sets the default component name to be the parameterlist's name
            if (src->tcc_names().size() == 0) {
              src->tcc_names().push_back(name);
            }

            // set the component indicies
            for (const auto& n : src->tcc_names()) {
              src->tcc_index().push_back(FindComponentNumber(n));
            }

            src->set_state(S_);
            srcs_.push_back(src);
          }
        }
      }
    }

#ifdef ALQUIMIA_ENABLED
    // -- try geochemical conditions
    Teuchos::ParameterList& glist = tp_list_->sublist("source terms").sublist("geochemical");

    for (auto it = glist.begin(); it != glist.end(); ++it) {
      std::string specname = it->first;
      Teuchos::ParameterList& spec = glist.sublist(specname);

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

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
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
    if (S->GetField(saturation_key_)->owner() == passwd_) {
      if (!S->GetField(saturation_key_, passwd_)->initialized()) {
        S->GetFieldData(saturation_key_, passwd_)->PutScalar(1.0);
        S->GetField(saturation_key_, passwd_)->set_initialized();

        if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
            *vo_->os() << "initilized saturation_liquid to value 1.0" << std::endl;
      }
      InitializeFieldFromField_(prev_saturation_key_, saturation_key_, S, true, true);

    } else {
      if (S->GetField(prev_saturation_key_)->owner() == passwd_) {
        if (!S->GetField(prev_saturation_key_, passwd_)->initialized()) {
          // S->GetFieldData(prev_saturation_key_, passwd_)->PutScalar(1.0);
          // S->GetField(prev_saturation_key_, passwd_)->set_initialized();
          // if (S->HasFieldEvaluator(Keys::getKey(domain_,saturation_key_))){
          //   S->GetFieldEvaluator(Keys::getKey(domain_,saturation_key_))->HasFieldChanged(S.ptr(), "transport");
          // }
          InitializeFieldFromField_(prev_saturation_key_, saturation_key_, S, true, true);
          S->GetField(prev_saturation_key_, passwd_)->set_initialized();

          if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
            *vo_->os() << "initilized prev_saturation_liquid from saturation" << std::endl;
        }
      }
    }
  }

  InitializeFieldFromField_(tcc_matrix_key_, tcc_key_, S, false, false);

  //InitializeFieldFromField_(flux_key_, darcy_flux_key_, S, false, false);

  // if (vol_flux_conversion_){
  //   vol_flux = S->GetFieldData(flux_key_, passwd_)->ViewComponent("face", true);
  //   ComputeVolumeDarcyFlux(S->GetFieldData(darcy_flux_key_)->ViewComponent("face", true),
  //                          S->GetFieldData(molar_density_key_)->ViewComponent("cell", true),
  //                          vol_flux);
  // }

  S->GetFieldData(solid_residue_mass_key_, passwd_)->PutScalar(0.0);
  S->GetField(solid_residue_mass_key_, passwd_)->set_initialized();
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
    if (S->GetField(field0)->owner() == passwd_) {
      if ((!S->GetField(field0, passwd_)->initialized())||(overwrite)) {
        if (call_evaluator)
            S->GetFieldEvaluator(field1)->HasFieldChanged(S.ptr(), passwd_);

        const CompositeVector& f1 = *S->GetFieldData(field1);
        CompositeVector& f0 = *S->GetFieldData(field0, passwd_);

        double vmin0, vmax0, vavg0;
        double vmin1, vmax1, vavg1;

        f0 = f1;

        S->GetField(field0, passwd_)->set_initialized();
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

  tcc = S_inter_->GetFieldData(tcc_key_, passwd_);
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
    //  flux_ = S_next_->GetFieldData(flux_key_)->ViewComponent("face", true);
    // *flux_copy_ = *flux_; // copy flux vector from S_next_ to S_;
    // double norm = 0.;
    // flux_->NormInf(&norm);
    // *vo_->os()<< name()<<" "<<"flux is copied norm:"<<norm<<"\n";
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

  flux_ = S_next_->GetFieldData(flux_key_)->ViewComponent("face", true);
  *flux_copy_ = *flux_; // copy flux vector from S_next_ to S_;

  ws_ = S_next_->GetFieldData(saturation_key_)->ViewComponent("cell", false);
  mol_dens_ = S_next_->GetFieldData(molar_density_key_)->ViewComponent("cell", false);
  solid_qty_ = S_next_->GetFieldData(solid_residue_mass_key_, passwd_)->ViewComponent("cell", false);

  // We use original tcc and make a copy of it later if needed.
  tcc = S_inter_->GetFieldData(tcc_key_, passwd_);
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");

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

  for (int c = 0; c < ncells_owned; c++) {
    double vol_phi_ws_den;
    vol_phi_ws_den = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_prev_)[0][c] * (*mol_dens_prev_)[0][c];
    for (int i=0; i<num_aqueous + num_gaseous; i++){
      mass_solutes_stepstart_[i] = tcc_prev[i][c] * vol_phi_ws_den;
    }
  }

  for (int c = 0; c < ncells_owned; c++) {
    double vol_phi_ws_den;
    vol_phi_ws_den = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_prev_)[0][c] * (*mol_dens_prev_)[0][c];
    for (int i=0; i<num_aqueous + num_gaseous; i++){
      mass_solutes_stepstart_[i] = tcc_prev[i][c] * vol_phi_ws_den;
    }
  }

  int ncycles = 0, swap = 1;
  while (dt_sum < dt_MPC - 1e-6) {
    // update boundary conditions
    time = t_physics_ + dt_cycle / 2;
    for (int i = 0; i < bcs_.size(); i++){
      bcs_[i]->Compute(time, time);
    }

    double dt_try = dt_MPC - dt_sum;
    double tol = 1e-10 * (dt_try + dt_stable);
    bool final_cycle = false;

    // *vo_->os() <<std::setprecision(10)<<"dt_MPC "<<dt_MPC<<" dt_cycle "<<dt_cycle<<" dt_sum "<<dt_sum<<" dt_stable "<<
    //   dt_stable<<" dt_try "<<dt_try<<" "<<dt_try - (dt_stable + tol)<<" tol "<<tol<<"\n";
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
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
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
    Teuchos::ParameterList& op_list =
        tp_list_->sublist("operators").sublist("diffusion operator").sublist("matrix");

    Teuchos::RCP<Operators::BCs> bc_dummy =
        Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

    // default boundary conditions (none inside domain and Neumann on its boundary)
    auto& bc_model = bc_dummy->bc_model();
    auto& bc_value = bc_dummy->bc_value();
    PopulateBoundaryData(bc_model, bc_value, -1);

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

    // Diffuse gaseous components. We ignore dispersion
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
     *S_inter_->GetFieldData("total_component_concentration_matrix", passwd_)->ViewComponent("cell");

  const Epetra_MultiVector& wcf_prev = *S_inter_->GetFieldData("prev_water_content")->ViewComponent("cell");
  const Epetra_MultiVector& wcf = *S_inter_->GetFieldData("water_content")->ViewComponent("cell");

  const Epetra_MultiVector& wcm_prev = *S_inter_->GetFieldData("prev_water_content_matrix")->ViewComponent("cell");
  const Epetra_MultiVector& wcm = *S_inter_->GetFieldData("water_content_matrix")->ViewComponent("cell");

  // multi-node matrix requires more input data
  const Epetra_MultiVector& phi_matrix = *S_inter_->GetFieldData("porosity_matrix")->ViewComponent("cell");

  int nnodes(1);
  Teuchos::RCP<Epetra_MultiVector> tcc_matrix_aux;
  if (S_inter_->HasField("total_component_concentration_matrix_aux")) {
    tcc_matrix_aux = S_inter_->GetFieldData("total_component_concentration_matrix_aux", passwd_)->ViewComponent("cell");
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
  Teuchos::RCP<CompositeVector> tcc;
  tcc = S->GetFieldData(tcc_key_, passwd_);
  *tcc = *tcc_tmp;
  InitializeFieldFromField_(prev_saturation_key_, saturation_key_, S.ptr(), false, true);

  // Copy to S_ as well
  tcc = S->GetFieldData(tcc_key_, passwd_);
  *tcc = *tcc_tmp;

  ChangedSolutionPK(S.ptr());
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
  double vol_phi_ws_den, tcc_flux;
  double mass_start = 0., tmp1, mass;

  // We advect only aqueous components.
  int num_advect = num_aqueous;

  for (int c = 0; c < ncells_owned; c++) {
    vol_phi_ws_den = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_start)[0][c] * (*mol_dens_start)[0][c];
    for (int i = 0; i < num_advect; i++){
      (*conserve_qty_)[i][c] = tcc_prev[i][c] * vol_phi_ws_den;

      if (dissolution_){
        if (( (*ws_start)[0][c]  > water_tolerance_) && ((*solid_qty_)[i][c] > 0 )){  // Dissolve solid residual into liquid
          double add_mass = std::min((*solid_qty_)[i][c], max_tcc_* vol_phi_ws_den - (*conserve_qty_)[i][c]);
          (*solid_qty_)[i][c] -= add_mass;
          (*conserve_qty_)[i][c] += add_mass;
        }
      }

      mass_start += (*conserve_qty_)[i][c];
    }
  }

  tmp1 = mass_start;
  mesh_->get_comm()->SumAll(&tmp1, &mass_start, 1);

  // advance all components at once
  for (int f = 0; f < nfaces_wghost; f++) {  // loop over master and slave faces
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];

    double u = fabs((*flux_)[0][f]);

    if (c1 >=0 && c1 < ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      for (int i = 0; i < num_advect; i++) {
        tcc_flux = dt_ * u * tcc_prev[i][c1];
        (*conserve_qty_)[i][c1] -= tcc_flux;
        (*conserve_qty_)[i][c2] += tcc_flux;
      }

    }
    else if (c1 >=0 && c1 < ncells_owned && (c2 >= ncells_owned || c2 < 0)) {
      for (int i = 0; i < num_advect; i++) {
        tcc_flux = dt_ * u * tcc_prev[i][c1];
        (*conserve_qty_)[i][c1] -= tcc_flux;
        if (c2 < 0) mass_solutes_bc_[i] -= tcc_flux;
        //AmanziGeometry::Point normal = mesh_->face_normal(f);
      }

    } else if (c1 >= ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      for (int i = 0; i < num_advect; i++) {
        tcc_flux = dt_ * u * tcc_prev[i][c1];
        (*conserve_qty_)[i][c2] += tcc_flux;
      }

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
      if (c2 >= 0) {
        double u = fabs((*flux_)[0][f]);
        for (int i = 0; i < ncomp; i++) {
          int k = tcc_index[i];
          if (k < num_advect) {
            tcc_flux = dt_ * u * values[i];
            (*conserve_qty_)[k][c2] += tcc_flux;
            mass_solutes_bc_[k] += tcc_flux;
          }
        }
      }
    }
  }

  // process external sources
  if (srcs_.size() != 0) {
    double time = t_physics_;
    ComputeAddSourceTerms(time, dt_, *conserve_qty_, 0, num_advect - 1);
  }

  // recover concentration from new conservative state
  for (int c = 0; c < ncells_owned; c++) {
    vol_phi_ws_den = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_end)[0][c] * (*mol_dens_end)[0][c];
    for (int i = 0; i < num_advect; i++) {

      if (vol_phi_ws_den > water_tolerance_ && (*conserve_qty_)[i][c] > 0) {
        tcc_next[i][c] = (*conserve_qty_)[i][c] / vol_phi_ws_den;
      }
      else  {
        (*solid_qty_)[i][c] += std::max((*conserve_qty_)[i][c], 0.);
        tcc_next[i][c] = 0.;
      }
    }
  }

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

    // if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH){
    //   tmp1 = mass_solutes_bc_[i];
    //   mesh_->get_comm()->SumAll(&tmp1, &mass_solutes_bc_[i], 1);
    //   *vo_->os() << "*****************\n";
    //   if (domain_name_ == "surface") *vo_->os()<<"Surface mass BC "<<mass_solutes_bc_[i]<<"\n";
    //   else *vo_->os() <<"Subsurface mass BC "<<mass_solutes_bc_[i]<<"\n";
    //   tmp1 = mass_solutes_source_[i];
    //   mesh_->get_comm()->SumAll(&tmp1, &mass_solutes_source_[i], 1);
    //   if (domain_name_ == "surface") *vo_->os()<<"Surface mass_solutes source "<<mass_solutes_source_[i]*dt_<<"\n";
    //   else *vo_->os() << "Subsurface mass_solutes source "<<mass_solutes_source_[i]*dt_<<"\n";
    //   *vo_->os() << "*****************\n";
    // }
  }

  if (internal_tests) {
    VV_CheckGEDproperty(*tcc_tmp->ViewComponent("cell"));
  }

  // if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH){
  //   if (domain_name_ == "surface")  *vo_->os()<<"Surface mass final "<<mass_final<<"\n";
  //   else  *vo_->os()<<"Subsurface mass final "<<mass_final<<"\n";
  // }

  // if (domain_name_ == "surface")  {
  //   Epetra_MultiVector& tcc_w_src_vec = *tcc_w_src->ViewComponent("cell");
  //   *vo_->os()<<"Surface mass final\n";
  //   std::cout<<" conserve_qty "<<" ws "<< "tcc "<< "tcc_w_src\n";
  //   for (int c = 0; c < ncells_owned; c++)
  //     std::cout<<(*conserve_qty_)[0][c]<<" "<<(*ws_end)[0][c]<<" "<<tcc_next[0][c]<<" "<<tcc_w_src_vec[0][c]<<"\n";
  // }

  // if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
  //   *vo_->os()<<"mass error "<<abs(mass_final - (mass_start + mass_solutes_bc_[0] + mass_solutes_source_[0]*dt_) )<<"\n";
  //if (abs(mass_final - (mass_start + mass_solutes_bc_[0] + mass_solutes_source_[0]*dt_) )/mass_final > 1e-6) exit(-1);
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
  Epetra_Vector f_component(cmap_wghost);

  // distribute vector of concentrations
  S_inter_->GetFieldData(tcc_key_)->ScatterMasterToGhosted("cell");
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);


  Epetra_Vector ws_ratio(Copy, *ws_start, 0);
  for (int c = 0; c < ncells_owned; c++){
    double vol_phi_ws_den_end = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_end)[0][c] * (*mol_dens_end)[0][c];
    if (vol_phi_ws_den_end > water_tolerance_)  {
      double vol_phi_ws_den_start = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_start)[0][c] * (*mol_dens_start)[0][c];
      if (vol_phi_ws_den_start > water_tolerance_){
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

  for (int i = 0; i < num_advect; i++) {
    current_component_ = i;  // needed by BJ

    double T = t_physics_;
    Epetra_Vector*& component = tcc_prev(i);
    FunctionalTimeDerivative(T, *component, f_component);

    for (int c = 0; c < ncells_owned; c++) {
      tcc_next[i][c] = (tcc_prev[i][c] + dt_ * f_component[c]) * ws_ratio[c];

      if (tcc_next[i][c] < 0){
        double vol_phi_ws_den = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_end)[0][c] * (*mol_dens_end)[0][c];
        (*solid_qty_)[i][c] += abs(tcc_next[i][c])*vol_phi_ws_den;
        tcc_next[i][c] = 0.;
      }
    }
  }

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

  //if (domain_name_ == "surface") {
  //*vo_->os()<<"after predictor ToTaL "<<domain_name_<<" :"<<std::setprecision(10)<<ComputeSolute( tcc_next, 0)<<"\n";
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
  // *vo_->os()<<domain_name_<<" difference "<<diff_norm<<"\n";
  //if (domain_name_ == "surface") {
  //*vo_->os()<<"after corrector ToTaL "<<domain_name_<<" :"<<std::setprecision(10)<<ComputeSolute( tcc_next, 0)<<"\n";
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

      for (int k = 0; k < tcc_index.size(); ++k) {
        int i = tcc_index[k];
        if (i < n0 || i > n1) continue;

        int imap = i;
        if (num_vectors == 1) imap = 0;

        double value;
        if (srcs_[m]->name() == "domain coupling") {
          value = values[k];
        } else {
          value =  mesh_->cell_volume(c) * values[k];
          // correction for non-SI concentration units
          // if (srcs_[m]->name() == "volume" || srcs_[m]->name() == "weight")
          //     value /= units_.concentration_factor();
        }

        cons_qty[imap][c] += dtp * value;
        mass_solutes_source_[i] += value;
      }
    }
  }
}


/* ******************************************************************
*
*
*
****************************************************************** */
void Transport_ATS::MixingSolutesWthSources(double told, double tnew)

{

  if (!special_source_) return;

  Key mass_source_key =  Keys::readKey(*tp_list_, domain_name_, "water source", "mass_source");
  if (!S_->HasFieldEvaluator(mass_source_key)) return;

  const Epetra_MultiVector& water_source =
    *S_->GetFieldData(mass_source_key)->ViewComponent("cell",false);

  Teuchos::RCP<const Epetra_MultiVector> nliq1_s;

  if (water_source_in_meters_) {
      // External source term is in [m water / s], not in [mols / s], so a
      // density is required.  This density should be upwinded.
    nliq1_s =
      S_->GetFieldData(Keys::getKey(domain_name_,"source_molar_density"))->ViewComponent("cell",false);
  }

  const Epetra_MultiVector& tcc_prev_vec = *tcc->ViewComponent("cell");
  Epetra_MultiVector& tcc_w_src_vec = *tcc_w_src->ViewComponent("cell");
  int num_vectors = tcc_w_src_vec.NumVectors();
  int nsrcs = srcs_.size();
  double dt = (tnew - told);

  for (int m = 0; m < nsrcs; m++) {

    srcs_[m]->Compute(told, tnew);

    std::vector<int> index = srcs_[m]->tcc_index();
    for (auto it = srcs_[m]->begin(); it != srcs_[m]->end(); ++it) {

      if (srcs_[m]->name() == "domain coupling") continue;

      int c = it->first;
      if (c >= ncells_owned) continue;

      std::vector<double>& values = it->second;
      double vol_phi_ws_den = mesh_->cell_volume(c) * (*phi_)[0][c] * (*ws_prev_)[0][c] * (*mol_dens_prev_)[0][c];
      double add_water = 0.;

      if (water_source_in_meters_) {
        double s1 = water_source[0][c] > 0. ? water_source[0][c] * (*nliq1_s)[0][c] : 0.;
        add_water = mesh_->cell_volume(c) * s1 * dt;
      } else {
        double s1 =  water_source[0][c] > 0. ? water_source[0][c] : 0.;
        add_water = mesh_->cell_volume(c) * s1 * dt;
      }

      for (int k = 0; k < index.size(); ++k) {
        int i = index[k];

        int imap = i;
        if (num_vectors == 1) imap = 0;

        double tcc_src_value =  mesh_->cell_volume(c) * values[k] * (tnew - told);
        double consv_qty = tcc_prev_vec[imap][c]*vol_phi_ws_den + tcc_src_value;

        if (vol_phi_ws_den + add_water > water_tolerance_){
          tcc_w_src_vec[imap][c] = consv_qty/(vol_phi_ws_den + add_water);
        }else{
          // if (vo_->getVerbLevel() > Teuchos::VERB_HIGH) {
          //   Teuchos::OSTab tab = vo_->getOSTab();
          //   *vo_->os() << vo_->color("yellow") << "No liquid phase in the sources. No mixing occured." << vo_->reset() << std::endl;
          // }
        }
      }
    }
  }

  // for (int c = 0; c < ncells_owned; c++) {
  //   mass2 += tcc[0][c];
  // }
  // tmp1 = mass2;
  // mesh_->get_comm()->SumAll(&tmp1, &mass2, 1);
  // tmp1 = add_mass;
  // mesh_->get_comm()->SumAll(&tmp1, &add_mass, 1);


}

void Transport_ATS::Sinks2TotalOutFlux(Epetra_MultiVector& tcc,
                                          std::vector<double>& total_outflux, int n0, int n1){

  std::vector<double> sink_add(ncells_wghost, 0.0);
  //Assumption that there is only one sink per component per cell
  double t0 = S_inter_->intermediate_time();
  int num_vectors = tcc.NumVectors();
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

        if ((values[k] < 0)&&(tcc[imap][c]>1e-16)) {
          if (srcs_[m]->name() == "domain coupling") {
            //val = std::max(val, fabs(values[k])/tcc[imap][c]);
            //val = std::max(val, fabs(values[k]));
            const Epetra_MultiVector& flux_interface_ = *S_next_->GetFieldData(coupled_flux_key)->ViewComponent("cell", false);
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
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
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

