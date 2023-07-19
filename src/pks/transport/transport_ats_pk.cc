/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

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
#include "Evaluator.hh"
#include "Mesh.hh"
#include "Mesh_Algorithms.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_Diffusion.hh"
#include "PDE_Accumulation.hh"
#include "PK_DomainFunctionFactory.hh"
#include "PK_Utils.hh"
#include "pk_helpers.hh"

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
                             const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                             const Teuchos::RCP<State>& S,
                             const Teuchos::RCP<TreeVector>& solution)
  : PK(pk_tree, global_plist, S, solution),
    PK_PhysicalExplicit<Epetra_Vector>(pk_tree, global_plist, S, solution)
{
  passwd_ = "state"; // this is what Amanzi uses

  if (plist_->isParameter("component names")) {
    component_names_ = plist_->get<Teuchos::Array<std::string>>("component names").toVector();
    num_components = component_names_.size();
    // otherwise we hopefully get them from chemistry
  }

  if (plist_->isParameter("component molar masses")) {
    mol_masses_ = plist_->get<Teuchos::Array<double>>("component molar masses").toVector();
  } else {
    Errors::Message msg("Transport PK: parameter \"component molar masses\" is missing.");
    Exceptions::amanzi_throw(msg);
  }

  // are we subcycling internally?
  subcycling_ = plist_->get<bool>("transport subcycling", false);
  tag_flux_next_ts_ = Tag{ name() + "_flux_next_ts" }; // what is this for? --ETC

  // initialize io
  units_.Init(global_plist->sublist("units"));

  // keys, dependencies, etc
  saturation_key_ = Keys::readKey(*plist_, domain_, "saturation liquid", "saturation_liquid");
  flux_key_ = Keys::readKey(*plist_, domain_, "water flux", "water_flux");
  permeability_key_ = Keys::readKey(*plist_, domain_, "permeability", "permeability");
  tcc_key_ = Keys::readKey(*plist_, domain_, "concentration", "total_component_concentration");
  conserve_qty_key_ =
    Keys::readKey(*plist_, domain_, "conserved quantity", "total_component_quantity");
  porosity_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
  molar_density_key_ = Keys::readKey(*plist_, domain_, "molar density", "molar_density_liquid");
  tcc_matrix_key_ =
    Keys::readKey(*plist_, domain_, "tcc matrix", "total_component_concentration_matrix");
  solid_residue_mass_key_ = Keys::readKey(*plist_, domain_, "solid residue", "solid_residue_mass");
  water_src_key_ = Keys::readKey(*plist_, domain_, "water source", "water_source");
  geochem_src_factor_key_ =
    Keys::readKey(*plist_, domain_, "geochem source factor", "geochem_src_factor");
  water_content_key_ = Keys::readKey(*plist_, domain_, "water content", "water_content");
  cv_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");
  key_ = tcc_key_;

  // other parameters
  water_tolerance_ = plist_->get<double>("water tolerance", 1e-6);
  dissolution_ = plist_->get<bool>("allow dissolution", false);
  max_tcc_ = plist_->get<double>("maximum concentration", 0.9);
  dim = mesh_->getSpaceDimension();

  db_ = Teuchos::rcp(new Debugger(mesh_, name_, *plist_));
}

void
Transport_ATS::set_tags(const Tag& current, const Tag& next)
{
  PK_PhysicalExplicit<Epetra_Vector>::set_tags(current, next);
  if (subcycling_) {
    tag_subcycle_current_ = Tag{ Keys::cleanName(name() + "_inner_subcycling_current") };
    tag_subcycle_next_ = Tag{ Keys::cleanName(name() + "_inner_subcycling_next") };
  } else {
    tag_subcycle_current_ = tag_current_;
    tag_subcycle_next_ = tag_next_;
  }
}


/* ******************************************************************
* Setup for Alquimia.
****************************************************************** */
#ifdef ALQUIMIA_ENABLED
void
Transport_ATS::SetupAlquimia(Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk,
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
    //
    // DOES TCC include secondaries?
    //
    // for (int i = 0; i < chem_engine_->NumAqueousComplexes(); ++i) {
    //   char secondary_name[128];
    //   snprintf(secondary_name, 127, "secondary_%d", i);
    //   component_names_.push_back(secondary_name);
    // }
    num_components = component_names_.size();
  }
}
#endif


/* ******************************************************************
* Define structure of this PK.
****************************************************************** */
void
Transport_ATS::Setup()
{
  // cross-coupling of PKs
  Teuchos::RCP<Teuchos::ParameterList> physical_models =
    Teuchos::sublist(plist_, "physical models and assumptions");

  if (num_components == 0) {
    Errors::Message msg("Transport PK: list of solutes is empty.");
    Exceptions::amanzi_throw(msg);
  }

  requireAtNext(tcc_key_, tag_subcycle_next_, *S_, passwd_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_components);
  S_->GetRecordSetW(tcc_key_).set_subfieldnames(component_names_);
  requireAtCurrent(tcc_key_, tag_subcycle_current_, *S_, passwd_);

  // CellVolume is required here -- it may not be used in this PK, but having
  // it makes vis nicer
  requireAtNext(cv_key_, tag_next_, *S_).SetMesh(mesh_)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // Raw data, no evaluator?
  //std::vector<std::string> primary_names(component_names_.begin(), component_names_.begin() + num_primary);
  auto primary_names = component_names_;
  // S_->Require<CompositeVector,CompositeVectorSpace>(solid_residue_mass_key_, tag_next_, name_)
  requireAtNext(solid_residue_mass_key_, tag_subcycle_next_, *S_, name_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, num_components);
  S_->GetRecordSetW(solid_residue_mass_key_).set_subfieldnames(primary_names);

  // This vector stores the conserved amount (in mols) of ncomponent
  // transported components, plus two for water.  The first water component is
  // given by the water content (in mols) at the old time plus dt * all fluxes
  // treated explictly.  The second water component is given by the water
  // content at the new time plus dt * all fluxes treated implicitly (notably
  // just DomainCoupling fluxes, which must be able to take all the transported
  // quantity.)
  //
  // Note that component_names includes secondaries, but we only need primaries
  primary_names.emplace_back("H2O_old");
  primary_names.emplace_back("H2O_new");
  //S_->Require<CompositeVector,CompositeVectorSpace>(conserve_qty_key_, tag_next_, name_)
  requireAtNext(conserve_qty_key_, tag_subcycle_next_, *S_, name_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, num_components + 2);
  S_->GetRecordSetW(conserve_qty_key_).set_subfieldnames(primary_names);

  // dependencies:
  // -- permeability
  bool abs_perm = physical_models->get<bool>("permeability field is required", false);
  if (abs_perm) {
    requireAtNext(permeability_key_, tag_next_, *S_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, dim);
  }

  // HACK ALERT -- FIXME --ETC
  //
  // This PK is liberally sprinkled with hard-coded Tags::NEXT and
  // Tags::CURRENT, forcing all things provided by FLOW to be provided at that
  // tag and not at tag_current and tag_next as it should be.  This is because
  // we don't have a good way of aliasing everything we need yet.  In
  // particular, aliases needed to be introduced between Setup() on flow and
  // Setup() on transport, and this was not possible when the quantity of
  // interest (porosity)'s evaluator was not required directly (only
  // indirectly) in flow PK.
  //
  // This will need to be fixed in amanzi/amanzi#646 somehow....? --ETC
  // -- water flux
  requireAtNext(flux_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  S_->Require<CompositeVector, CompositeVectorSpace>(flux_key_, tag_flux_next_ts_, name_);

  // -- water saturation
  requireAtNext(saturation_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  // Require a copy of saturation at the old time tag
  requireAtCurrent(saturation_key_, Tags::CURRENT, *S_);
  if (subcycling_) {
    S_->Require<CompositeVector, CompositeVectorSpace>(
      saturation_key_, tag_subcycle_current_, name_);
    S_->Require<CompositeVector, CompositeVectorSpace>(saturation_key_, tag_subcycle_next_, name_);
    // S_->RequireEvaluator(saturation_key_, tag_subcycle_current_); // for the future...
    // S_->RequireEvaluator(saturation_key_, tag_subcycle_next_); // for the future...
  }

  requireAtNext(porosity_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  requireAtNext(molar_density_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  requireAtCurrent(molar_density_key_, Tags::CURRENT, *S_);
  if (subcycling_) {
    S_->Require<CompositeVector, CompositeVectorSpace>(
      molar_density_key_, tag_subcycle_current_, name_);
    S_->Require<CompositeVector, CompositeVectorSpace>(
      molar_density_key_, tag_subcycle_next_, name_);
    // S_->RequireEvaluator(molar_density_key_, tag_subcycle_current_); // for the future...
    // S_->RequireEvaluator(molar_density_key_, tag_subcycle_next_); // for the future...
  }


  has_water_src_key_ = false;
  if (plist_->sublist("source terms").isSublist("geochemical")) {
    requireAtNext(water_src_key_, Tags::NEXT, *S_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    has_water_src_key_ = true;
    water_src_in_meters_ = plist_->get<bool>("water source in meters", false);

    if (water_src_in_meters_) {
      geochem_src_factor_key_ = water_src_key_;
    } else {
      // set the coefficient as water source / water density
      Teuchos::ParameterList& wc_eval = S_->GetEvaluatorList(geochem_src_factor_key_);
      wc_eval.set<std::string>("evaluator type", "reciprocal evaluator");
      std::vector<std::string> dep{ water_src_key_, molar_density_key_ };
      wc_eval.set<Teuchos::Array<std::string>>("dependencies", dep);
      wc_eval.set<std::string>("reciprocal", dep[1]);

      requireAtNext(geochem_src_factor_key_, Tags::NEXT, *S_)
        .SetMesh(mesh_)
        ->SetGhosted(true)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    }
  }

  // this is the not-yet-existing source, and is dead code currently! (What does this mean? --ETC)
  if (plist_->sublist("source terms").isSublist("component concentration source")) {
    requireAtNext(water_src_key_, Tags::NEXT, *S_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    has_water_src_key_ = true;
    water_src_in_meters_ = plist_->get<bool>("water source in meters", false);
  }


  if (plist_->sublist("source terms").isSublist("component mass source")) {
    Teuchos::ParameterList& conc_sources_list =
      plist_->sublist("source terms").sublist("component mass source");

    for (const auto& it : conc_sources_list) {
      std::string name = it.first;
      if (conc_sources_list.isSublist(name)) {
        Teuchos::ParameterList& src_list = conc_sources_list.sublist(name);
        std::string src_type = src_list.get<std::string>("spatial distribution method", "none");
        if ((src_type == "field") && (src_list.isSublist("field"))) {
          solute_src_key_ = src_list.sublist("field").get<std::string>("field key");
          requireAtNext(solute_src_key_, Tags::NEXT, *S_)
            .SetMesh(mesh_)
            ->SetGhosted(true)
            ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_components);
          S_->Require<CompositeVector, CompositeVectorSpace>(
            solute_src_key_, tag_next_, solute_src_key_);
        }
      }
    }
  }

  // // alias to next for subcycled cases -- revisit this in state subcycling
  // // revision --ETC
  // if (tag_next_ != Tags::NEXT) {
  //   aliasVector(*S_, tcc_key_, tag_next_, Tags::NEXT);
  //   aliasVector(*S_, conserve_qty_key_, tag_next_, Tags::NEXT);
  //   aliasVector(*S_, tcc_matrix_key_, tag_next_, Tags::NEXT);
  //   aliasVector(*S_, solid_residue_mass_key_, tag_next_, Tags::NEXT);
  // }
}


/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.
****************************************************************** */
void
Transport_ATS::Initialize()
{
  // Set initial values for transport variables.
  dt_ = dt_debug_ = t_physics_ = 0.0;
  double time = S_->get_time();
  if (time >= 0.0) t_physics_ = time;

  if (plist_->isSublist("initial condition")) {
    S_->GetRecordW(tcc_key_, tag_subcycle_next_, passwd_)
      .Initialize(plist_->sublist("initial condition"));
  }

  internal_tests = 0;
  tests_tolerance = TRANSPORT_CONCENTRATION_OVERSHOOT;

  bc_scaling = 0.0;

  Teuchos::OSTab tab = vo_->getOSTab();
  MyPID = mesh_->getComm()->MyPID();

  // initialize missed fields
  InitializeFields_();

  // make this go away -- local pointers to data
  tcc_tmp = S_->GetPtrW<CompositeVector>(tcc_key_, tag_subcycle_next_, passwd_);
  tcc = S_->GetPtrW<CompositeVector>(tcc_key_, tag_subcycle_current_, passwd_);
  *tcc = *tcc_tmp;

  ws_ = S_->Get<CompositeVector>(saturation_key_, Tags::NEXT).ViewComponent("cell", false);
  ws_prev_ = S_->Get<CompositeVector>(saturation_key_, Tags::CURRENT).ViewComponent("cell", false);

  mol_dens_ = S_->Get<CompositeVector>(molar_density_key_, Tags::NEXT).ViewComponent("cell", false);
  mol_dens_prev_ =
    S_->Get<CompositeVector>(molar_density_key_, Tags::CURRENT).ViewComponent("cell", false);

  if (subcycling_) {
    ws_subcycle_current = S_->GetW<CompositeVector>(saturation_key_, tag_subcycle_current_, name_)
                            .ViewComponent("cell");
    ws_subcycle_next =
      S_->GetW<CompositeVector>(saturation_key_, tag_subcycle_next_, name_).ViewComponent("cell");

    mol_dens_subcycle_current =
      S_->GetW<CompositeVector>(molar_density_key_, tag_subcycle_current_, name_)
        .ViewComponent("cell");
    mol_dens_subcycle_next =
      S_->GetW<CompositeVector>(molar_density_key_, tag_subcycle_next_, name_)
        .ViewComponent("cell");
  }

  flux_ = S_->Get<CompositeVector>(flux_key_, Tags::NEXT).ViewComponent("face", true);
  flux_copy_ =
    S_->GetW<CompositeVector>(flux_key_, tag_flux_next_ts_, name_).ViewComponent("face", true);
  flux_copy_->PutScalar(0.);

  phi_ = S_->Get<CompositeVector>(porosity_key_, Tags::NEXT).ViewComponent("cell", false);
  solid_qty_ = S_->GetW<CompositeVector>(solid_residue_mass_key_, tag_next_, name_)
                 .ViewComponent("cell", false);
  conserve_qty_ =
    S_->GetW<CompositeVector>(conserve_qty_key_, tag_next_, name_).ViewComponent("cell", true);

  // Check input parameters. Due to limited amount of checks, we can do it earlier.
  Policy(tag_next_);

  ncells_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  ncells_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  nfaces_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  nfaces_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
  nnodes_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);

  // extract control parameters
  InitializeAll_();

  // upwind
  const Epetra_Map& fmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);
  upwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));
  downwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));
  IdentifyUpwindCells();

  // advection block initialization
  current_component_ = -1;

  // reconstruction initialization
  limiter_ = Teuchos::rcp(new Operators::LimiterCell(mesh_));
  lifting_ = Teuchos::rcp(new Operators::ReconstructionCellLinear(mesh_));

  // mechanical dispersion
  flag_dispersion_ = false;
  if (plist_->isSublist("material properties")) {
    auto mdm_list = Teuchos::sublist(plist_, "material properties");
    mdm_ = CreateMDMPartition(mesh_, mdm_list, flag_dispersion_);
    if (flag_dispersion_) CalculateAxiSymmetryDirection();
  }

  // create boundary conditions
  if (plist_->isSublist("boundary conditions")) {
    // -- try tracer-type conditions
    PK_DomainFunctionFactory<TransportDomainFunction> factory(mesh_, S_);
    Teuchos::ParameterList& conc_bcs_list =
      plist_->sublist("boundary conditions").sublist("concentration");

    for (const auto& it : conc_bcs_list) {
      std::string name = it.first;
      if (conc_bcs_list.isSublist(name)) {
        Teuchos::ParameterList& bc_list = conc_bcs_list.sublist(name);
        std::string bc_type = bc_list.get<std::string>("spatial distribution method", "none");

        if (bc_type == "domain coupling") {
          // See amanzi ticket #646 -- this should probably be tag_subcycle_current_?
          // domain couplings are special -- they always work on all components
          Teuchos::RCP<TransportDomainFunction> bc =
            factory.Create(bc_list, "fields", AmanziMesh::Entity_kind::FACE, Kxy, tag_current_);

          for (int i = 0; i < num_components; i++) {
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
          int gid = std::stoi(domain_.substr(last_of + 1, domain_.size()));
          bc_list.set("entity_gid_out", gid);

          // See amanzi ticket #646 -- this should probably be tag_subcycle_current_?
          Teuchos::RCP<TransportDomainFunction> bc =
            factory.Create(bc_list, "boundary concentration", AmanziMesh::Entity_kind::FACE, Kxy, tag_current_);

          for (int i = 0; i < num_components; i++) {
            bc->tcc_names().push_back(component_names_[i]);
            bc->tcc_index().push_back(i);
          }
          bc->set_state(S_);
          bcs_.push_back(bc);

        } else {
          // See amanzi ticket #646 -- this should probably be tag_subcycle_current_?
          Teuchos::RCP<TransportDomainFunction> bc = factory.Create(
            bc_list, "boundary concentration function", AmanziMesh::Entity_kind::FACE, Kxy, tag_current_);
          bc->set_state(S_);

          std::vector<std::string> tcc_names =
            bc_list.get<Teuchos::Array<std::string>>("component names").toVector();
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
    Teuchos::ParameterList& geochem_plist =
      plist_->sublist("boundary conditions").sublist("geochemical");

    for (Teuchos::ParameterList::ConstIterator it = geochem_plist.begin();
         it != geochem_plist.end();
         ++it) {
      std::string specname = it->first;
      Teuchos::ParameterList& spec = geochem_plist.sublist(specname);

      Teuchos::RCP<TransportBoundaryFunction_Alquimia_Units> bc = Teuchos::rcp(
        new TransportBoundaryFunction_Alquimia_Units(spec, mesh_, chem_pk_, chem_engine_));

      bc->set_conversion(1000.0, mol_dens_, true);
      //bc->set_name("alquimia bc");
      std::vector<int>& tcc_index = bc->tcc_index();
      std::vector<std::string>& tcc_names = bc->tcc_names();

      for (int i = 0; i < tcc_names.size(); i++) {
        tcc_index.push_back(FindComponentNumber(tcc_names[i]));
      }

      bcs_.push_back(bc);
    }
#endif

  } else {
    if (vo_->os_OK(Teuchos::VERB_NONE)) {
      *vo_->os() << vo_->color("yellow") << "No BCs were specified." << vo_->reset() << std::endl;
    }
  }

  // boundary conditions initialization
  time = t_physics_;
  VV_CheckInfluxBC();

  // source term initialization: so far only "concentration" is available.
  if (plist_->isSublist("source terms")) {
    PK_DomainFunctionFactory<TransportDomainFunction> factory(mesh_, S_);
    Teuchos::ParameterList& conc_sources_list =
      plist_->sublist("source terms").sublist("component mass source");

    for (const auto& it : conc_sources_list) {
      std::string name = it.first;
      if (conc_sources_list.isSublist(name)) {
        Teuchos::ParameterList& src_list = conc_sources_list.sublist(name);
        std::string src_type = src_list.get<std::string>("spatial distribution method", "none");

        if (src_type == "domain coupling") {
          // domain couplings are special -- they always work on all components
          Teuchos::RCP<TransportDomainFunction> src =
            factory.Create(src_list, "fields", AmanziMesh::Entity_kind::CELL, Kxy);

          for (int i = 0; i < num_components; i++) {
            src->tcc_names().push_back(component_names_[i]);
            src->tcc_index().push_back(i);
          }
          src->set_state(S_);
          srcs_.push_back(src);
        } else if (src_type == "field") {
          // domain couplings are special -- they always work on all components
          Teuchos::RCP<TransportDomainFunction> src =
            factory.Create(src_list, "field", AmanziMesh::Entity_kind::CELL, Kxy);

          for (int i = 0; i < component_names_.size(); i++) {
            src->tcc_names().push_back(component_names_[i]);
            src->tcc_index().push_back(i);
          }
          src->set_state(S_);
          srcs_.push_back(src);
        } else {
          // See amanzi ticket #646 -- this should probably be tag_subcycle_current_?
          Teuchos::RCP<TransportDomainFunction> src =
            factory.Create(src_list, "source function", AmanziMesh::Entity_kind::CELL, Kxy, tag_current_);

          std::vector<std::string> tcc_names =
            src_list.get<Teuchos::Array<std::string>>("component names").toVector();
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

      Teuchos::RCP<TransportSourceFunction_Alquimia_Units> src = Teuchos::rcp(
        new TransportSourceFunction_Alquimia_Units(spec, mesh_, chem_pk_, chem_engine_));

      if (S_->HasEvaluator(geochem_src_factor_key_, Tags::NEXT)) {
        S_->GetEvaluator(geochem_src_factor_key_, Tags::NEXT).Update(*S_, name_);
      }

      auto src_factor =
        S_->Get<CompositeVector>(geochem_src_factor_key_, Tags::NEXT).ViewComponent("cell", false);
      src->set_conversion(-1000., src_factor, false);

      for (const auto& n : src->tcc_names()) { src->tcc_index().push_back(FindComponentNumber(n)); }

      srcs_.push_back(src);
    }
#endif
  }
  // Temporarily Transport hosts Henry law.
  PrepareAirWaterPartitioning_();

  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    *vo_->os() << "Number of components: " << tcc->size() << std::endl
               << "cfl=" << cfl_ << " spatial/temporal discretization: " << spatial_disc_order
               << " " << temporal_disc_order << std::endl;
    *vo_->os() << vo_->color("green") << "Initalization of PK is complete." << vo_->reset()
               << std::endl
               << std::endl;
  }

  // ETC BEGIN HACKING
  StableTimeStep();
}


/* ******************************************************************
* Initalized fields left by State and other PKs.
****************************************************************** */
void
Transport_ATS::InitializeFields_()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  InitializeFieldFromField_(tcc_matrix_key_, tag_next_, tcc_key_, tag_next_, false, false);
  S_->GetW<CompositeVector>(solid_residue_mass_key_, tag_next_, name_).PutScalar(0.0);
  S_->GetRecordW(solid_residue_mass_key_, tag_next_, name_).set_initialized();
  S_->GetRecordW(conserve_qty_key_, tag_next_, name_).set_initialized();
}

/* ****************************************************************
* Auxiliary initialization technique.
**************************************************************** */
void
Transport_ATS::InitializeFieldFromField_(const Key& field0,
                                         const Tag& tag0,
                                         const Key& field1,
                                         const Tag& tag1,
                                         bool call_evaluator,
                                         bool overwrite)
{
  if (S_->HasRecord(field0, tag0)) {
    if (S_->GetRecord(field0, tag0).owner() == name_) {
      if ((!S_->GetRecord(field0, tag0).initialized()) || overwrite) {
        if (call_evaluator) S_->GetEvaluator(field1, tag0).Update(*S_, name_);

        const CompositeVector& f1 = S_->Get<CompositeVector>(field1, tag1);
        CompositeVector& f0 = S_->GetW<CompositeVector>(field0, tag0, name_);
        f0 = f1;

        S_->GetRecordW(field0, tag0, name_).set_initialized();
        if (vo_->os_OK(Teuchos::VERB_MEDIUM) && (!overwrite)) {
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
double
Transport_ATS::StableTimeStep()
{
  S_->Get<CompositeVector>(flux_key_, Tags::NEXT).ScatterMasterToGhosted("face");
  flux_ = S_->Get<CompositeVector>(flux_key_, Tags::NEXT).ViewComponent("face", true);

  Teuchos::RCP<Epetra_Map> cell_map = Teuchos::rcp(new Epetra_Map(mesh_->getMap(AmanziMesh::Entity_kind::CELL,false)));
  IdentifyUpwindCells();

  tcc = S_->GetPtrW<CompositeVector>(tcc_key_, tag_current_, passwd_);
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");

  // loop over faces and accumulate upwinding fluxes
  std::vector<double> total_outflux(ncells_wghost, 0.0);

  for (int f = 0; f < nfaces_wghost; f++) {
    int c = (*upwind_cell_)[f];
    if (c >= 0) { total_outflux[c] += fabs((*flux_)[0][f]); }
  }

  Sinks2TotalOutFlux(tcc_prev, total_outflux, 0, num_aqueous - 1);

  // loop over cells and calculate minimal time step
  double vol = 0.;
  double ws_min_dt = 0.;
  double outflux_min_dt = 0.;
  dt_ = TRANSPORT_LARGE_TIME_STEP;
  double dt_cell = TRANSPORT_LARGE_TIME_STEP;
  int cmin_dt = 0;
  for (int c = 0; c < ncells_owned; c++) {
    double outflux = total_outflux[c];

    if ((outflux > 0) && ((*ws_prev_)[0][c] > 0) && ((*ws_)[0][c] > 0) && ((*phi_)[0][c] > 0)) {
      vol = mesh_->getCellVolume(c);
      dt_cell = vol * (*mol_dens_)[0][c] * (*phi_)[0][c] *
                std::min((*ws_prev_)[0][c], (*ws_)[0][c]) / outflux;
    }
    if (dt_cell < dt_) {
      dt_ = dt_cell;
      cmin_dt = c;
      ws_min_dt = std::min((*ws_prev_)[0][c], (*ws_)[0][c]);
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
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    int cmin_dt_unique = (std::abs(dt_tmp * cfl_ - dt_) < 1e-6 * dt_) ? cell_map->GID(cmin_dt) : -2;

    int cmin_dt_tmp = cmin_dt_unique;
    comm.MaxAll(&cmin_dt_tmp, &cmin_dt_unique, 1);
    int min_pid = -1;

    double tmp_package[6];
    if (cell_map->GID(cmin_dt) == cmin_dt_unique) {
      const AmanziGeometry::Point& p = mesh_->getCellCentroid(cmin_dt);

      min_pid = comm.MyPID();
      tmp_package[0] = ws_min_dt;
      tmp_package[1] = outflux_min_dt;
      tmp_package[2] = p[0];
      tmp_package[3] = p[1];
      if (p.dim() == 3)
        tmp_package[4] = p[2];
      else
        tmp_package[4] = 0.;
      tmp_package[5] = p.dim();
    }

    int min_pid_tmp = min_pid;
    comm.MaxAll(&min_pid_tmp, &min_pid, 1);
    comm.Broadcast(tmp_package, 6, min_pid);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Stable time step " << dt_ << " is computed at (" << tmp_package[2] << ", "
               << tmp_package[3];
    if (fabs(3 - tmp_package[5]) < 1e-10) *vo_->os() << ", " << tmp_package[4];
    *vo_->os() << ")" << std::endl;
    *vo_->os() << "Stable time step " << dt_ << " is limited by saturation/ponded_depth "
               << tmp_package[0] << " and "
               << "output flux " << tmp_package[1] << std::endl;
  }
  return dt_;
}


/* *******************************************************************
* Estimate returns last time step unless it is zero.
******************************************************************* */
double
Transport_ATS::get_dt()
{
  if (subcycling_) {
    return std::numeric_limits<double>::max();
  } else {
    //StableTimeStep();
    return dt_;
  }
}


/* *******************************************************************
* MPC will call this function to advance the transport state.
* Efficient subcycling requires to calculate an intermediate state of
* saturation only once, which leads to a leap-frog-type algorithm.
******************************************************************* */
bool
Transport_ATS::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool failed = false;
  double dt_MPC = t_new - t_old;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_LOW))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << S_->get_time(tag_current_)
               << " t1 = " << S_->get_time(tag_next_) << " h = " << dt_MPC << std::endl
               << "----------------------------------------------------------------" << std::endl;

  // NOTE: these "flow" variables are hard-coded as Tag::NEXT assuming that
  // flow is supercycled relative to transport and therefore we must
  // interpolate the flow variables from the "global" CURRENT+NEXT to the
  // subcycled current + next.  This would be fixed by having evaluators that
  // interpolate in time, allowing transport to not have to know how flow is
  // being integrated... FIXME --etc
  S_->GetEvaluator(flux_key_, Tags::NEXT).Update(*S_, name_);

  // why are we re-assigning all of these?  The previous pointers shouldn't have changed... --ETC
  flux_ = S_->Get<CompositeVector>(flux_key_, Tags::NEXT).ViewComponent("face", true);
  // why are we copying this?  This should result in constant flux, no need to copy? --ETC
  *flux_copy_ = *flux_; // copy flux vector from S_next_ to S_;

  S_->GetEvaluator(saturation_key_, Tags::NEXT).Update(*S_, name_);
  ws_ = S_->Get<CompositeVector>(saturation_key_, Tags::NEXT).ViewComponent("cell", false);
  S_->GetEvaluator(saturation_key_, Tags::CURRENT).Update(*S_, name_);
  ws_prev_ = S_->Get<CompositeVector>(saturation_key_, Tags::CURRENT).ViewComponent("cell", false);

  S_->GetEvaluator(molar_density_key_, Tags::NEXT).Update(*S_, name_);
  mol_dens_ = S_->Get<CompositeVector>(molar_density_key_, Tags::NEXT).ViewComponent("cell", false);
  S_->GetEvaluator(molar_density_key_, Tags::CURRENT).Update(*S_, name_);
  mol_dens_prev_ =
    S_->Get<CompositeVector>(molar_density_key_, Tags::CURRENT).ViewComponent("cell", false);

  //if (subcycling_) S_->set_time(tag_subcycle_current_, t_old);

  // this is locally created and has no evaluator -- should get a primary
  // variable FE owned by this PK --ETC
  solid_qty_ = S_->GetW<CompositeVector>(solid_residue_mass_key_, tag_next_, name_)
                 .ViewComponent("cell", false);

#ifdef ALQUIMIA_ENABLED
  if (plist_->sublist("source terms").isSublist("geochemical")) {
    for (auto& src : srcs_) {
      if (src->name() == "alquimia source") {
        // src_factor = water_source / molar_density_liquid, both flow
        // quantities, see note above.
        S_->GetEvaluator(geochem_src_factor_key_, Tags::NEXT).Update(*S_, name_);
        auto src_factor = S_->Get<CompositeVector>(geochem_src_factor_key_, Tags::NEXT)
                            .ViewComponent("cell", false);
        Teuchos::RCP<TransportSourceFunction_Alquimia_Units> src_alq =
          Teuchos::rcp_dynamic_cast<TransportSourceFunction_Alquimia_Units>(src);
        src_alq->set_conversion(-1000, src_factor, false);
      }
    }
  }

  if (plist_->sublist("boundary conditions").isSublist("geochemical")) {
    for (auto& bc : bcs_) {
      if (bc->name() == "alquimia bc") {
        Teuchos::RCP<TransportBoundaryFunction_Alquimia_Units> bc_alq =
          Teuchos::rcp_dynamic_cast<TransportBoundaryFunction_Alquimia_Units>(bc);
        bc_alq->set_conversion(1000.0, mol_dens_, true);
      }
    }
  }
#endif

  // We use original tcc and make a copy of it later if needed.
  tcc = S_->GetPtrW<CompositeVector>(tcc_key_, tag_current_, passwd_);
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");
  db_->WriteVector("tcc_old", tcc.ptr());

  // calculate stable time step
  double dt_shift = 0.0, dt_global = dt_MPC;
  double time = t_old;
  if (time >= 0.0) {
    t_physics_ = time;
    dt_shift = time - S_->get_time(tag_current_);
    dt_global = S_->get_time(tag_next_) - S_->get_time(tag_current_);
    AMANZI_ASSERT(std::abs(dt_global - dt_MPC) < 1.e-4);
  }

  if (subcycling_)
    StableTimeStep();
  else
    dt_ = dt_MPC;
  double dt_stable = dt_; // advance routines override dt_

  int interpolate_ws = 0; // (dt_ < dt_global) ? 1 : 0;
  if ((t_old > S_->get_time(tag_current_)) || (t_new < S_->get_time(tag_next_))) interpolate_ws = 1;

  double dt_sum = 0.0;
  double dt_cycle;
  Tag water_tag_current, water_tag_next;
  if (interpolate_ws) {
    dt_cycle = std::min(dt_stable, dt_MPC);
    InterpolateCellVector(*ws_prev_, *ws_, dt_shift, dt_global, *ws_subcycle_current);
    InterpolateCellVector(
      *mol_dens_prev_, *mol_dens_, dt_shift, dt_global, *mol_dens_subcycle_current);
    InterpolateCellVector(*ws_prev_, *ws_, dt_shift + dt_cycle, dt_global, *ws_subcycle_next);
    InterpolateCellVector(
      *mol_dens_prev_, *mol_dens_, dt_shift + dt_cycle, dt_global, *mol_dens_subcycle_next);
    ws_current = ws_subcycle_current;
    ws_next = ws_subcycle_next;
    mol_dens_current = mol_dens_subcycle_current;
    mol_dens_next = mol_dens_subcycle_next;
    water_tag_current = tag_subcycle_current_;
    water_tag_next = tag_subcycle_next_;
  } else {
    dt_cycle = dt_MPC;
    ws_current = ws_prev_;
    ws_next = ws_;
    mol_dens_current = mol_dens_prev_;
    mol_dens_next = mol_dens_;
    water_tag_current = Tags::CURRENT;
    water_tag_next = Tags::NEXT;
  }

  db_->WriteVector("sat_old",
                   S_->GetPtr<CompositeVector>(saturation_key_, water_tag_current).ptr());
  db_->WriteVector("sat_new", S_->GetPtr<CompositeVector>(saturation_key_, water_tag_next).ptr());
  db_->WriteVector("mol_dens_old",
                   S_->GetPtr<CompositeVector>(molar_density_key_, water_tag_current).ptr());
  db_->WriteVector("mol_dens_new",
                   S_->GetPtr<CompositeVector>(molar_density_key_, water_tag_next).ptr());
  db_->WriteVector("poro", S_->GetPtr<CompositeVector>(porosity_key_, Tags::NEXT).ptr());

  for (int c = 0; c < ncells_owned; c++) {
    double vol_phi_ws_den;
    vol_phi_ws_den =
      mesh_->getCellVolume(c) * (*phi_)[0][c] * (*ws_prev_)[0][c] * (*mol_dens_prev_)[0][c];
    for (int i = 0; i < num_aqueous + num_gaseous; i++) {
      mass_solutes_stepstart_[i] = tcc_prev[i][c] * vol_phi_ws_den;
    }
  }

  int ncycles = 0, swap = 1;
  while (dt_sum < dt_MPC - 1e-6) {
    // update boundary conditions
    time = t_physics_ + dt_cycle / 2;
    for (int i = 0; i < bcs_.size(); i++) { bcs_[i]->Compute(t_physics_, t_physics_ + dt_cycle); }

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
      if (swap) { // Initial water saturation is in 'start'.
        ws_current = ws_subcycle_current;
        ws_next = ws_subcycle_next;
        mol_dens_current = mol_dens_subcycle_current;
        mol_dens_next = mol_dens_subcycle_next;

        double dt_int = dt_sum + dt_shift;
        InterpolateCellVector(*ws_prev_, *ws_, dt_int, dt_global, *ws_subcycle_next);
        InterpolateCellVector(
          *mol_dens_prev_, *mol_dens_, dt_int, dt_global, *mol_dens_subcycle_next);
      } else { // Initial water saturation is in 'end'.
        ws_current = ws_subcycle_next;
        ws_next = ws_subcycle_current;
        mol_dens_current = mol_dens_subcycle_next;
        mol_dens_next = mol_dens_subcycle_current;

        double dt_int = dt_sum + dt_shift;
        InterpolateCellVector(*ws_prev_, *ws_, dt_int, dt_global, *ws_subcycle_current);
        InterpolateCellVector(
          *mol_dens_prev_, *mol_dens_, dt_int, dt_global, *mol_dens_subcycle_current);
      }
      swap = 1 - swap;
    }

    if (spatial_disc_order == 1) { // temporary solution (lipnikov@lanl.gov)
      AdvanceDonorUpwind(dt_cycle);
    } else if (spatial_disc_order == 2 && temporal_disc_order == 1) {
      AdvanceSecondOrderUpwindRK1(dt_cycle);
    } else if (spatial_disc_order == 2 && temporal_disc_order == 2) {
      AdvanceSecondOrderUpwindRK2(dt_cycle);
    }

    if (!final_cycle) { // rotate concentrations (we need new memory for tcc)
      // should not be allocating here, we have tons of memory for tcc --ETC
      tcc = Teuchos::RCP<CompositeVector>(new CompositeVector(*tcc_tmp));
    }

    ncycles++;
  }

  dt_ = dt_stable; // restore the original time step (just in case)

  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", false);
  Advance_Dispersion_Diffusion(t_old, t_new);
  // optional Henry Law for the case of gas diffusion
  if (henry_law_) MakeAirWaterPartitioning_();

  // statistics output
  nsubcycles = ncycles;
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    *vo_->os() << ncycles << " sub-cycles, dt_stable=" << units_.OutputTime(dt_stable)
               << " [sec]  dt_MPC=" << units_.OutputTime(dt_MPC) << " [sec]" << std::endl;

    VV_PrintSoluteExtrema(tcc_next, dt_MPC);
  }

  // ETC BEGIN HACKING
  StableTimeStep();
  return failed;
}


void
Transport_ATS ::Advance_Dispersion_Diffusion(double t_old, double t_new)
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
      Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
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
      Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, op));

    const CompositeVectorSpace& cvs = op1->global_operator()->DomainMap();
    CompositeVector sol(cvs), factor(cvs), factor0(cvs), source(cvs), zero(cvs);
    zero.PutScalar(0.0);

    // populate the dispersion operator (if any)
    if (flag_dispersion_) { CalculateDispersionTensor_(*flux_, *phi_, *ws_, *mol_dens_); }

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
      for (int c = 0; c < ncells_owned; c++) { sol_cell[0][c] = tcc_next[i][c]; }
      if (sol.HasComponent("face")) { sol.ViewComponent("face")->PutScalar(0.0); }

      if (flag_op1) {
        op->Init();
        Teuchos::RCP<std::vector<WhetStone::Tensor>> Dptr = Teuchos::rcpFromRef(D_);
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
          double tmp =
            mesh_->getCellVolume(c) * (*ws_)[0][c] * (*phi_)[0][c] * (*mol_dens_)[0][c] / dt_MPC;
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

      for (int c = 0; c < ncells_owned; c++) { tcc_next[i][c] = sol_cell[0][c]; }
      if (sol.HasComponent("face")) {
        if (tcc_tmp->HasComponent("boundary_face")) {
          Epetra_MultiVector& tcc_tmp_bf = *tcc_tmp->ViewComponent("boundary_face", false);
          Epetra_MultiVector& sol_faces = *sol.ViewComponent("face", false);
          const Epetra_Map& vandalay_map = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false);
          const Epetra_Map& face_map = mesh_->getMap(AmanziMesh::Entity_kind::FACE,false);
          int nbfaces = tcc_tmp_bf.MyLength();
          for (int bf = 0; bf != nbfaces; ++bf) {
            AmanziMesh::Entity_ID f = face_map.LID(vandalay_map.GID(bf));
            tcc_tmp_bf[i][bf] = sol_faces[i][f];
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
      for (int c = 0; c < ncells_owned; c++) { sol_cell[0][c] = tcc_next[i][c]; }
      if (sol.HasComponent("face")) { sol.ViewComponent("face")->PutScalar(0.0); }

      op->Init();
      Teuchos::RCP<std::vector<WhetStone::Tensor>> Dptr = Teuchos::rcpFromRef(D_);
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
        if ((*ws_)[0][c] == 1.0) fac1[0][c] = 1.0 * (*mol_dens_)[0][c]; // hack so far
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

      for (int c = 0; c < ncells_owned; c++) { tcc_next[i][c] = sol_cell[0][c]; }
    }

    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "dispersion solver ||r||=" << residual / num_components
                 << " itrs=" << num_itrs / num_components << std::endl;
    }
  }
}


/* *******************************************************************
* Copy the advected tcc field to the state.
******************************************************************* */
void
Transport_ATS::CommitStep(double t_old, double t_new, const Tag& tag_next)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Commiting state @ " << tag_next << std::endl;

  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;

  assign(tcc_key_, tag_current, tag_next, *S_);
  if (tag_next == Tags::NEXT) {
    assign(saturation_key_, tag_current, tag_next, *S_);
    assign(molar_density_key_, tag_current, tag_next, *S_);
  }
}


/* *******************************************************************
 * A simple first-order transport method
 ****************************************************************** */
void
Transport_ATS::AdvanceDonorUpwind(double dt_cycle)
{
  IdentifyUpwindCells();
  dt_ = dt_cycle; // overwrite the maximum stable transport step
  mass_solutes_source_.assign(num_aqueous + num_gaseous, 0.0);
  mass_solutes_bc_.assign(num_aqueous + num_gaseous, 0.0);

  // populating next state of concentrations
  tcc->ScatterMasterToGhosted("cell");
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);

  // prepare conservative state in master and slave cells
  double mass_current = 0., tmp1, mass;

  int num_components = tcc_next.NumVectors();
  conserve_qty_->PutScalar(0.);

  for (int c = 0; c < ncells_owned; c++) {
    double vol_phi_ws_den =
      mesh_->getCellVolume(c) * (*phi_)[0][c] * (*ws_current)[0][c] * (*mol_dens_current)[0][c];
    (*conserve_qty_)[num_components + 1][c] = vol_phi_ws_den;

    for (int i = 0; i < num_advect; i++) {
      (*conserve_qty_)[i][c] = tcc_prev[i][c] * vol_phi_ws_den;

      if (dissolution_) {
        if (((*ws_current)[0][c] > water_tolerance_) &&
            ((*solid_qty_)[i][c] > 0)) { // Dissolve solid residual into liquid
          double add_mass =
            std::min((*solid_qty_)[i][c], max_tcc_ * vol_phi_ws_den - (*conserve_qty_)[i][c]);
          (*solid_qty_)[i][c] -= add_mass;
          (*conserve_qty_)[i][c] += add_mass;
        }
      }

      mass_current += (*conserve_qty_)[i][c];
    }
  }

  db_->WriteCellVector("cons (start)", *conserve_qty_);
  tmp1 = mass_current;
  mesh_->getComm()->SumAll(&tmp1, &mass_current, 1);

  // advance all components at once
  for (int f = 0; f < nfaces_wghost; f++) { // loop over master and slave faces
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];
    double u = fabs((*flux_)[0][f]);

    if (c1 >= 0 && c1 < ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      for (int i = 0; i < num_advect; i++) {
        double tcc_flux = dt_ * u * tcc_prev[i][c1];
        (*conserve_qty_)[i][c1] -= tcc_flux;
        (*conserve_qty_)[i][c2] += tcc_flux;
      }
      (*conserve_qty_)[num_components + 1][c1] -= dt_ * u;
      (*conserve_qty_)[num_components + 1][c2] += dt_ * u;
    } else if (c1 >= 0 && c1 < ncells_owned && (c2 >= ncells_owned || c2 < 0)) {
      for (int i = 0; i < num_advect; i++) {
        double tcc_flux = dt_ * u * tcc_prev[i][c1];
        (*conserve_qty_)[i][c1] -= tcc_flux;
        if (c2 < 0) mass_solutes_bc_[i] -= tcc_flux;
        //AmanziGeometry::Point normal = mesh_->getFaceNormal(f);
      }
      (*conserve_qty_)[num_components + 1][c1] -= dt_ * u;

    } else if (c1 >= ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      for (int i = 0; i < num_advect; i++) {
        double tcc_flux = dt_ * u * tcc_prev[i][c1];
        (*conserve_qty_)[i][c2] += tcc_flux;
      }
      (*conserve_qty_)[num_components + 1][c2] += dt_ * u;

    } else if (c2 < 0 && c1 >= 0 && c1 < ncells_owned) {
      (*conserve_qty_)[num_components + 1][c1] -= dt_ * u;

    } else if (c1 < 0 && c2 >= 0 && c2 < ncells_owned) {
      (*conserve_qty_)[num_components + 1][c2] += dt_ * u;
    }
  }

  Epetra_MultiVector* tcc_tmp_bf = nullptr;
  if (tcc_tmp->HasComponent("boundary_face")) {
    tcc_tmp_bf = &(*tcc_tmp->ViewComponent("boundary_face", false));
  }

  // loop over exterior boundary sets
  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
      int f = it->first;
      int bf = tcc_tmp_bf ? AmanziMesh::getFaceOnBoundaryBoundaryFace(*mesh_, f) : -1;

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

            if (tcc_tmp_bf) (*tcc_tmp_bf)[i][bf] = values[i];
          }
        }
      }
    }
  }
  db_->WriteCellVector("cons (adv)", *conserve_qty_);

  // process external sources
  if (srcs_.size() != 0) {
    double time = t_physics_;
    ComputeAddSourceTerms(time, dt_, *conserve_qty_, 0, num_aqueous - 1);
  }
  db_->WriteCellVector("cons (src)", *conserve_qty_);

  // recover concentration from new conservative state
  for (int c = 0; c < ncells_owned; c++) {
    double water_new =
      mesh_->getCellVolume(c) * (*phi_)[0][c] * (*ws_next)[0][c] * (*mol_dens_next)[0][c];
    double water_sink =
      (*conserve_qty_)[num_components]
                      [c]; // water at the new time + outgoing domain coupling source
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
  // tcc_next.Print(std::cout);
  VV_PrintSoluteExtrema(tcc_next, dt_);

  double mass_final = 0;
  for (int c = 0; c < ncells_owned; c++) {
    for (int i = 0; i < num_advect; i++) { mass_final += (*conserve_qty_)[i][c]; }
  }

  tmp1 = mass_final;
  mesh_->getComm()->SumAll(&tmp1, &mass_final, 1);

  // update mass balance
  for (int i = 0; i < mass_solutes_exact_.size(); i++) {
    mass_solutes_exact_[i] += mass_solutes_source_[i] * dt_;
  }

  if (internal_tests) { VV_CheckGEDproperty(*tcc_tmp->ViewComponent("cell")); }
}


/* *******************************************************************
 * We have to advance each component independently due to different
 * reconstructions. We use tcc when only owned data are needed and
 * tcc_next when owned and ghost data. This is a special routine for
 * transient flow and uses first-order time integrator.
 ****************************************************************** */
void
Transport_ATS::AdvanceSecondOrderUpwindRK1(double dt_cycle)
{
  dt_ = dt_cycle; // overwrite the maximum stable transport step
  mass_solutes_source_.assign(num_aqueous + num_gaseous, 0.0);

  // work memory
  const Epetra_Map& cmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::CELL,true);

  // distribute vector of concentrations
  tcc->ScatterMasterToGhosted("cell");
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);

  // Epetra_Vector ws_ratio(Copy, *ws_current, 0);
  // for (int c = 0; c < ncells_owned; c++) {
  //   double vol_phi_ws_den_next = mesh_->getCellVolume(c) * (*phi_)[0][c] * (*ws_next)[0][c] * (*mol_dens_next)[0][c];
  //   if (vol_phi_ws_den_next > water_tolerance_)  {
  //     double vol_phi_ws_den_current = mesh_->getCellVolume(c) * (*phi_)[0][c] * (*ws_current)[0][c] * (*mol_dens_current)[0][c];
  //     if (vol_phi_ws_den_current > water_tolerance_) {
  //       ws_ratio[c] = ( (*ws_current)[0][c] * (*mol_dens_current)[0][c] )
  //                   / ( (*ws_next)[0][c]   * (*mol_dens_next)[0][c]   );
  //     } else {
  //       ws_ratio[c] = 1;
  //     }
  //   }
  //   else  ws_ratio[c]=0.;
  // }

  // We advect only aqueous components.
  int num_components = tcc_next.NumVectors();
  conserve_qty_->PutScalar(0.);

  // prepopulate with initial water for better debugging
  for (int c = 0; c < ncells_owned; c++) {
    double vol_phi_ws_den_current =
      mesh_->getCellVolume(c) * (*phi_)[0][c] * (*ws_current)[0][c] * (*mol_dens_current)[0][c];
    (*conserve_qty_)[num_components + 1][c] = vol_phi_ws_den_current;
  }

  for (int i = 0; i < num_advect; i++) {
    current_component_ = i; // needed by BJ
    double T = t_physics_;
    Epetra_Vector*& component = tcc_prev(i);
    FunctionalTimeDerivative(T, *component, *(*conserve_qty_)(i));
  }
  db_->WriteCellVector("cons (time_deriv)", *conserve_qty_);

  // calculate the new conc
  for (int c = 0; c < ncells_owned; c++) {
    double water_old = (*conserve_qty_)[num_components + 1][c];
    double water_new =
      mesh_->getCellVolume(c) * (*phi_)[0][c] * (*ws_next)[0][c] * (*mol_dens_next)[0][c];
    double water_sink = (*conserve_qty_)[num_components][c];
    double water_total = water_sink + water_new;
    (*conserve_qty_)[num_components][c] = water_total;

    for (int i = 0; i != num_components; ++i) {
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

  if (internal_tests) { VV_CheckGEDproperty(*tcc_tmp->ViewComponent("cell")); }
}


/* *******************************************************************
 * We have to advance each component independently due to different
 * reconstructions. This is a special routine for transient flow and
 * uses second-order predictor-corrector time integrator.
 ****************************************************************** */
void
Transport_ATS::AdvanceSecondOrderUpwindRK2(double dt_cycle)
{
  dt_ = dt_cycle; // overwrite the maximum stable transport step
  mass_solutes_source_.assign(num_aqueous + num_gaseous, 0.0);

  // work memory
  const Epetra_Map& cmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::CELL,true);
  Epetra_Vector f_component(cmap_wghost); //,  f_component2(cmap_wghost);

  // distribute old vector of concentrations
  tcc->ScatterMasterToGhosted("cell");
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);

  Epetra_Vector ws_ratio(Copy, *ws_current, 0);
  for (int c = 0; c < ncells_owned; c++) {
    if ((*ws_next)[0][c] > 1e-10) {
      if ((*ws_current)[0][c] > 1e-10) {
        ws_ratio[c] = ((*ws_current)[0][c] * (*mol_dens_current)[0][c]) /
                      ((*ws_next)[0][c] * (*mol_dens_next)[0][c]);
      } else {
        ws_ratio[c] = 1;
      }
    } else
      ws_ratio[c] = 0.;
  }

  // predictor step
  for (int i = 0; i < num_advect; i++) {
    current_component_ = i; // needed by BJ

    double T = t_physics_;
    Epetra_Vector*& component = tcc_prev(i);
    FunctionalTimeDerivative(T, *component, f_component);

    for (int c = 0; c < ncells_owned; c++) {
      tcc_next[i][c] = (tcc_prev[i][c] + dt_ * f_component[c]) * ws_ratio[c];
    }
  }

  tcc_tmp->ScatterMasterToGhosted("cell");

  // corrector step
  for (int i = 0; i < num_advect; i++) {
    current_component_ = i; // needed by BJ

    double T = t_physics_;
    Epetra_Vector*& component = tcc_next(i);
    FunctionalTimeDerivative(T, *component, f_component);

    for (int c = 0; c < ncells_owned; c++) {
      double value = (tcc_prev[i][c] + dt_ * f_component[c]) * ws_ratio[c];
      tcc_next[i][c] = (tcc_next[i][c] + value) / 2;
      if (tcc_next[i][c] < 0) {
        double vol_phi_ws_den =
          mesh_->getCellVolume(c) * (*phi_)[0][c] * (*ws_next)[0][c] * (*mol_dens_next)[0][c];
        (*solid_qty_)[i][c] += abs(tcc_next[i][c]) * vol_phi_ws_den;
        tcc_next[i][c] = 0.;
      }
    }
  }

  // update mass balance
  for (int i = 0; i < num_aqueous + num_gaseous; i++) {
    mass_solutes_exact_[i] += mass_solutes_source_[i] * dt_ / 2;
  }

  if (internal_tests) { VV_CheckGEDproperty(*tcc_tmp->ViewComponent("cell")); }
}


/* ******************************************************************
* Computes source and sink terms and adds them to vector tcc.
* Returns mass rate for the tracer.
* The routine treats two cases of tcc with one and all components.
****************************************************************** */
void
Transport_ATS::ComputeAddSourceTerms(double tp,
                                     double dtp,
                                     Epetra_MultiVector& cons_qty,
                                     int n0,
                                     int n1)
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
        (*conserve_qty_)[num_vectors - 2][c] += values[num_vectors - 2];
      }

      for (int k = 0; k < tcc_index.size(); ++k) {
        int i = tcc_index[k];
        if (i < n0 || i > n1) continue;

        int imap = i;
        if (num_vectors == 1) imap = 0;
        double value = mesh_->getCellVolume(c) * values[k];
        cons_qty[imap][c] += dtp * value;
        mass_solutes_source_[i] += value;
      }
    }
  }
}


void
Transport_ATS::Sinks2TotalOutFlux(Epetra_MultiVector& tcc_c,
                                  std::vector<double>& total_outflux,
                                  int n0,
                                  int n1)
{
  std::vector<double> sink_add(ncells_wghost, 0.0);
  //Assumption that there is only one sink per component per cell
  double t0 = S_->get_time(tag_current_);
  int num_vectors = tcc_c.NumVectors();
  int nsrcs = srcs_.size();

  // YUCK no requires, just random data access and hard-coded names --ETC
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
            const Epetra_MultiVector& flux_interface_ =
              *S_->Get<CompositeVector>(coupled_flux_key, Tags::NEXT).ViewComponent("cell", false);
            val = std::max(val, fabs(flux_interface_[0][c]));
          }
        }
      }
      sink_add[c] = std::max(sink_add[c], val);
    }
  }

  for (int c = 0; c < ncells_wghost; c++) total_outflux[c] += sink_add[c];
}


/* *******************************************************************
* Populates operators' boundary data for given component.
* Returns true if at least one face was populated.
******************************************************************* */
bool
Transport_ATS::PopulateBoundaryData(std::vector<int>& bc_model,
                                    std::vector<double>& bc_value,
                                    int component)
{
  bool flag = false;

  for (int i = 0; i < bc_model.size(); i++) {
    bc_model[i] = Operators::OPERATOR_BC_NONE;
    bc_value[i] = 0.0;
  }

  for (int f = 0; f < nfaces_wghost; f++) {
    auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
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
void
Transport_ATS::IdentifyUpwindCells()
{
  for (int f = 0; f < nfaces_wghost; f++) {
    (*upwind_cell_)[f] = -1; // negative value indicates boundary
    (*downwind_cell_)[f] = -1;
  }

  for (int c = 0; c < ncells_wghost; c++) {
    const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);

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


void
Transport_ATS::ComputeVolumeDarcyFlux(Teuchos::RCP<const Epetra_MultiVector> flux,
                                      Teuchos::RCP<const Epetra_MultiVector> molar_density,
                                      Teuchos::RCP<Epetra_MultiVector>& vol_darcy_flux)
{

  for (int f = 0; f < nfaces_wghost; f++) {
    auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    double n_liq = 0.;
    for (int c = 0; c < cells.size(); c++) n_liq += (*molar_density)[0][c];
    n_liq /= cells.size();
    if (n_liq > 0)
      (*vol_darcy_flux)[0][f] = (*flux_)[0][f] / n_liq;
    else
      (*vol_darcy_flux)[0][f] = 0.;
  }
}


/* *******************************************************************
* Interpolate linearly in time between two values v0 and v1. The time
* is measuared relative to value v0; so that v1 is at time dt. The
* interpolated data are at time dt_int.
******************************************************************* */
void
Transport_ATS::InterpolateCellVector(const Epetra_MultiVector& v0,
                                     const Epetra_MultiVector& v1,
                                     double dt_int,
                                     double dt,
                                     Epetra_MultiVector& v_int)
{
  double a = dt_int / dt;
  double b = 1.0 - a;
  v_int.Update(b, v0, a, v1, 0.);
}

} // namespace Transport
} // namespace Amanzi
