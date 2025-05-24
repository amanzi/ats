/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoont@ornl.gov)
           Phong Le (lepv@ornl.gov)
*/

/*
  Transport PK

*/

#include <algorithm>
#include <vector>

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

#include "errors.hh"
#include "Explicit_TI_RK.hh"
#include "Evaluator.hh"
#include "Mesh.hh"
#include "MeshAlgorithms.hh"
#include "TensorVector.hh"
#include "OperatorDefs.hh"
#include "BCs.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_Diffusion.hh"
#include "PDE_Accumulation.hh"

#include "PK_Utils.hh"
#include "PK_Helpers.hh"

#include "PK_DomainFunctionFactory.hh"
#include "TransportDomainFunction.hh"
#include "TransportBoundaryFunction_Alquimia.hh"
#include "TransportSourceFunction_Alquimia.hh"
#include "TransportDomainFunction_UnitConversion.hh"

#include "PK_DomainFunction.hh"

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
    PK_Physical_Default(pk_tree, global_plist, S, solution),
    dt_stable_(-1.0),
    dt_max_(-1.0),
    has_diffusion_(false),
    has_dispersion_(false)
{
  // initialize io
  units_.Init(global_plist->sublist("units"));
}

void
Transport_ATS::parseParameterList()
{
  if (!plist_->isParameter("primary variable key suffix")) {
    plist_->set<std::string>("primary variable key suffix", "mole_fraction");
  }

  // with subfield names, the header width is often insufficient
  if (!plist_->sublist("verbose object").isParameter("debug cell header width"))
    plist_->sublist("verbose object").set("debug cell header width", 34);
  PK_Physical_Default::parseParameterList();

  // protect user from old naming convention
  if (Keys::getVarName(key_) == "total_component_concentration") {
    Errors::Message msg;
    msg << "Transport_ATS PK \"" << name() << "\": primary variable can no longer be called "
        << "\"total_component_concentration\", but should instead be left blank (to use "
        << "\"molar_fraction\") or provided something else.  Transport units are "
        << "[mol-C mol-H2O^-1], not [mol L^-1], and therefore should not be called concentration.";
    Exceptions::amanzi_throw(msg);
  }

  if (component_names_.size() == 0) {
    // not set by chemistry... must get set by user
    component_names_ = plist_->get<Teuchos::Array<std::string>>("component names").toVector();
    num_components_ = component_names_.size();
  }

  // NOTE: names MUST be aqueous, solid, gaseous
  num_aqueous_ = plist_->get<int>("number of aqueous components", component_names_.size());

  // parameters
  molar_masses_ = readParameterMapByComponent(plist_->sublist("component molar masses [kg / mol C]"), 1.0);
  tcc_max_ = readParameterMapByComponent(plist_->sublist("component maximum concentration [mol C / mol H2O]"), -1.0);

  if (plist_->isSublist("molecular diffusivity [m^2 s^-1]")) {
    has_diffusion_ = true;
    molec_diff_ = readParameterMapByComponent(plist_->sublist("molecular diffusivity [m^2 s^-1]"), 0.);

    tortuosity_ = readParameterMapByPhase(plist_->sublist("tortuosity [-]"), 1.);
  }

  // keys, dependencies, etc
  // -- flux -- only needed at new time, evaluator controlled elsewhere
  flux_key_ = Keys::readKey(*plist_, domain_, "water flux", "water_flux");

  mass_flux_key_ = Keys::readKey(*plist_, domain_, "mass flux", "mass_flux");
  requireEvaluatorAtNext(mass_flux_key_, tag_next_, *S_, name_);

  // -- liquid water content - need at new time, copy at current time
  lwc_key_ = Keys::readKey(*plist_, domain_, "liquid water content", "water_content");
  requireEvaluatorAtCurrent(lwc_key_, tag_current_, *S_, name_);

  water_src_key_ = Keys::readKey(*plist_, domain_, "water source", "water_source");
  cv_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");

  // workspace, no evaluator
  conserve_qty_key_ =
    Keys::readKey(*plist_, domain_, "conserved quantity", "total_component_quantity");
  requireEvaluatorAtNext(conserve_qty_key_, tag_next_, *S_, name_);

  solid_residue_mass_key_ = Keys::readKey(*plist_, domain_, "solid residue", "solid_residue_quantity");

  geochem_src_factor_key_ =
    Keys::readKey(*plist_, domain_, "geochem source factor", "geochem_src_factor");

  if (chem_engine_ != Teuchos::null) {
    // needed by geochemical bcs
    molar_dens_key_ = Keys::readKey(*plist_, domain_, "molar density liquid", "molar_density_liquid");
  }

  // dispersion coefficient tensor
  dispersion_tensor_key_ = Keys::readKey(*plist_, domain_, "dispersion coefficient", "dispersion_coefficient");
  has_dispersion_ = S_->HasEvaluatorList(dispersion_tensor_key_);

  // other parameters
  // -- a small amount of water, used to define when we are going to completely dry out a grid cell
  water_tolerance_ = plist_->get<double>("water tolerance [mol H2O / m^d]", 1e-6);

  // global transport parameters
  cfl_ = plist_->get<double>("cfl", 1.0);
  dt_max_ = plist_->get<double>("maximum timestep [s]", TRANSPORT_LARGE_TIME_STEP);

  adv_spatial_disc_order_ = plist_->get<int>("advection spatial discretization order", 1);
  if (adv_spatial_disc_order_ < 1 || adv_spatial_disc_order_ > 2) {
    Errors::Message msg;
    msg << "Transport_ATS: \"advection spatial discretization order\" must be 1 or 2, not "
        << adv_spatial_disc_order_;
    Exceptions::amanzi_throw(msg);
  }

  temporal_disc_order_ = plist_->get<int>("temporal discretization order", 1);
  if (temporal_disc_order_ < 1 || temporal_disc_order_ > 2) {
    Errors::Message msg;
    msg << "Transport_ATS: \"temporal discretization order\" must be 1 or 2, not " << temporal_disc_order_;
    Exceptions::amanzi_throw(msg);
  }

  db_ = Teuchos::rcp(new Debugger(mesh_, name_, *plist_));
}


Transport_ATS::ParameterMap
Transport_ATS::readParameterMapByComponent(Teuchos::ParameterList& plist,
        double default_val)
{
  Transport_ATS::ParameterMap map;
  for (int i = 0; i != num_components_; ++i) {
    map[component_names_[i]] = plist.get<double>(component_names_[i], default_val);
  }
  return map;
}


Transport_ATS::ParameterMap
Transport_ATS::readParameterMapByPhase(Teuchos::ParameterList& plist,
        double default_val)
{
  Transport_ATS::ParameterMap map;
  map["aqueous"] = plist.get<double>("aqueous", default_val);
  map["solid"] = plist.get<double>("solid", default_val);
  map["gaseous"] = plist.get<double>("gaseous", default_val);
  return map;
}


/* ******************************************************************
* Setup for Alquimia.
****************************************************************** */
#ifdef ALQUIMIA_ENABLED
void
Transport_ATS::setChemEngine(Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk)
{
  chem_pk_ = chem_pk;
  chem_engine_ = chem_pk->getChemEngine();

  if (chem_engine_ != Teuchos::null) {
    // Retrieve the component names (primary and secondary) from the chemistry
    // engine.
    chem_engine_->GetPrimarySpeciesNames(component_names_);
    num_components_ = component_names_.size();
  }
}
#endif


/* ******************************************************************
* Define structure of this PK.
****************************************************************** */
void
Transport_ATS::Setup()
{
  PK_Physical_Default::Setup();
  SetupTransport_();
  SetupPhysicalEvaluators_();
}


/* ******************************************************************
* Define structure of this PK.
****************************************************************** */
void
Transport_ATS::SetupTransport_()
{
  // upwind and downwind vectors
  const Epetra_Map& fmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::FACE, true);
  upwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));
  downwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));

  if (adv_spatial_disc_order_ == 2) {
    // reconstruction initialization
    Teuchos::ParameterList& recon_list = plist_->sublist("reconstruction");

    // check and set defaults
    if (!recon_list.isParameter("limiter extension for transport"))
      recon_list.set<bool>("limiter extension for transport", true);
    if (!recon_list.isParameter("limiter"))
      recon_list.set<std::string>("limiter", "tensorial");

    lifting_ = Teuchos::rcp(new Operators::ReconstructionCellLinear(mesh_));
    lifting_->Init(recon_list);

    limiter_ = Teuchos::rcp(new Operators::LimiterCell(mesh_));
    limiter_->Init(recon_list);
  }

  adv_bcs_ = Teuchos::rcp(
    new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));

  // workspace for diffusion and dispersion solve
  if (has_dispersion_) {
    // note this space has the wrong number of DoFs, but that will be corrected
    // by the evaluator later.  The rest of the info (name, location, and mesh)
    // are needed.
    CompositeVectorSpace disp_space;
    disp_space.SetMesh(mesh_)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    S_->Require<TensorVector, TensorVector_Factory>(dispersion_tensor_key_, tag_next_)
      .set_map(disp_space);
    S_->RequireEvaluator(dispersion_tensor_key_, tag_next_);
  }

  // operator and boundary conditions for diffusion/dispersion solve
  if (has_dispersion_ || has_diffusion_) {
    // default boundary conditions (none inside domain and Neumann on its boundary)
    diff_bcs_ = Teuchos::rcp(
      new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
    PopulateBoundaryData_(-1, *diff_bcs_);

    // diffusion operator
    Operators::PDE_DiffusionFactory opfactory;
    Teuchos::ParameterList& op_list = plist_->sublist("diffusion");
    diff_op_ = opfactory.Create(op_list, mesh_, diff_bcs_);
    diff_global_op_ = diff_op_->global_operator();
    diff_acc_op_ = Teuchos::rcp(
      new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, diff_global_op_));

    // diffusion workspace
    const CompositeVectorSpace& cvs = diff_global_op_->DomainMap();
    diff_sol_ = Teuchos::rcp(new CompositeVector(cvs));
  }


  // NOTE: these to go away! ETC
  //
  // source term setup
  // --------------------------------------------------------------------------------
  if (plist_->isSublist("source terms")) {
    auto sources_list = Teuchos::sublist(plist_, "source terms");

    // sources of mass of C
    if (sources_list->isSublist("component mass source")) {
      auto conc_sources_list = Teuchos::sublist(sources_list, "component mass source");
      PK_DomainFunctionFactory<TransportDomainFunction> factory(mesh_, S_);

      for (const auto& it : *conc_sources_list) {
        std::string name = it.first;
        if (conc_sources_list->isSublist(name)) {
          auto src_list = Teuchos::sublist(conc_sources_list, name);

          std::string src_type = src_list->get<std::string>("spatial distribution method", "none");

          if (src_type == "domain coupling") {
            Teuchos::RCP<TransportDomainFunction> src =
              factory.Create(*src_list, "fields", AmanziMesh::Entity_kind::CELL, Teuchos::null, tag_current_);

            // domain couplings functions is special -- always work on all components
            for (int i = 0; i < num_components_; i++) {
              src->tcc_names().push_back(component_names_[i]);
              src->tcc_index().push_back(i);
            }
            src->set_state(S_);
            srcs_.push_back(src);

          } else if (src_type == "field") {
            Teuchos::RCP<TransportDomainFunction> src =
              factory.Create(*src_list, "field", AmanziMesh::Entity_kind::CELL, Teuchos::null, tag_current_);

            for (int i = 0; i < num_components_; i++) {
              src->tcc_names().push_back(component_names_[i]);
              src->tcc_index().push_back(i);
            }
            src->set_state(S_);
            srcs_.push_back(src);

            Teuchos::ParameterList flist = src_list->sublist("field");
            // if there are more than one field
            if (flist.isParameter("number of fields")) {
              if (flist.isType<int>("number of fields")) {
                int num_fields = flist.get<int>("number of fields");

                if (num_fields < 1) { // ERROR -- invalid number of fields
                  AMANZI_ASSERT(0);
                }

                for (int fid = 1; fid != (num_fields + 1); ++fid) {
                  std::stringstream sublist_name;
                  sublist_name << "field " << fid << " info";
                  Key field_key = Keys::readKey(flist.sublist(sublist_name.str()), domain_, "field");
                  Tag field_tag = Keys::readTag(flist.sublist(sublist_name.str()));
                  requireEvaluatorAtNext(field_key, field_tag, *S_)
                    .SetMesh(mesh_)
                    ->SetGhosted(true)
                    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
                }
              }
            } else { // if only one field
              Key field_key = Keys::readKey(src_list->sublist("field"), domain_, "field");
              Tag field_tag = Keys::readTag(src_list->sublist("field"), "field");
              requireEvaluatorAtNext(field_key, field_tag, *S_)
                .SetMesh(mesh_)
                ->SetGhosted(true)
                ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_components_);
              // NOTE: this code should be moved to the PK_DomainFunctionField, and all other
              // PK_DomainFunction* should be updated to make sure they require their data! --ETC
              // Problem: This requires the pk_helper.hh to be included in PK_DomainFunctionField.hh --PL
            }

          } else { // all others work on a subset of components
            Teuchos::RCP<TransportDomainFunction> src = factory.Create(
              *src_list, "source function", AmanziMesh::Entity_kind::CELL, Teuchos::null, tag_current_);
            src->set_tcc_names(src_list->get<Teuchos::Array<std::string>>("component names").toVector());
            for (const auto& n : src->tcc_names()) {
              src->tcc_index().push_back(FindComponentNumber_(n));
            }

            src->set_state(S_);
            srcs_.push_back(src);
          }
        }
      }
    }

    // sources of water that include C at a known concentration
    if (sources_list->isSublist("geochemical")) {
      // note these are computed at the flow PK's NEXT tag, which assumes all
      // sources are dealt with implicitly (backward Euler).  This could be relaxed --ETC
      requireEvaluatorAtNext(water_src_key_, tag_next_, *S_)
        .SetMesh(mesh_)
        ->SetGhosted(true)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
      has_water_src_key_ = true;

      // this flag is just for convenience -- some flow PKs accept a water
      // source in m/s not in mol/m^d/s.
      water_src_in_meters_ = plist_->get<bool>("water source in meters", false);

      if (water_src_in_meters_) {
        geochem_src_factor_key_ = water_src_key_;
      } else {
        // set the coefficient as water source / water density
        Teuchos::ParameterList& wc_eval = S_->GetEvaluatorList(geochem_src_factor_key_);
        wc_eval.set<std::string>("evaluator type", "reciprocal evaluator");
        std::vector<std::string> dep{ water_src_key_, molar_dens_key_ };
        wc_eval.set<Teuchos::Array<std::string>>("dependencies", dep);
        wc_eval.set<std::string>("reciprocal", dep[1]);
      }
      requireEvaluatorAtNext(geochem_src_factor_key_, tag_next_, *S_)
        .SetMesh(mesh_)
        ->SetGhosted(true)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    }
  }

  //
  // create boundary conditions
  // --------------------------------------------------------------------------------
  if (plist_->isSublist("boundary conditions")) {
    // -- try tracer-type conditions
    PK_DomainFunctionFactory<TransportDomainFunction> factory(mesh_, S_);
    auto bcs_list = Teuchos::sublist(plist_, "boundary conditions");
    auto conc_bcs_list = Teuchos::sublist(bcs_list, "mole fraction");

    for (const auto& it : *conc_bcs_list) {
      std::string name = it.first;
      if (conc_bcs_list->isSublist(name)) {
        Teuchos::ParameterList& bc_list = conc_bcs_list->sublist(name);
        std::string bc_type = bc_list.get<std::string>("spatial distribution method", "none");

        if (bc_type == "domain coupling") {
          // domain couplings are special -- they always work on all components
          Teuchos::RCP<TransportDomainFunction> bc =
            factory.Create(bc_list, "fields", AmanziMesh::Entity_kind::FACE, Teuchos::null, tag_current_);

          for (int i = 0; i < num_components_; i++) {
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
          bc_list.set("entity GID", gid);

          Teuchos::RCP<TransportDomainFunction> bc = factory.Create(
            bc_list, "boundary mole fraction", AmanziMesh::Entity_kind::FACE, Teuchos::null, tag_current_);

          for (int i = 0; i < num_components_; i++) {
            bc->tcc_names().push_back(component_names_[i]);
            bc->tcc_index().push_back(i);
          }
          bc->set_state(S_);
          bcs_.push_back(bc);

        } else {
          Teuchos::RCP<TransportDomainFunction> bc =
            factory.Create(bc_list,
                           "boundary mole fraction function",
                           AmanziMesh::Entity_kind::FACE,
                           Teuchos::null,
                           tag_current_);
          bc->set_state(S_);

          std::vector<std::string> tcc_names =
            bc_list.get<Teuchos::Array<std::string>>("component names").toVector();
          bc->set_tcc_names(tcc_names);

          // set the component indicies
          for (const auto& n : bc->tcc_names()) {
            bc->tcc_index().push_back(FindComponentNumber_(n));
          }
          bcs_.push_back(bc);
        }
      }
    }

#ifdef ALQUIMIA_ENABLED
    // -- try geochemical conditions
    auto geochem_list = Teuchos::sublist(bcs_list, "geochemical");

    for (const auto& it : *geochem_list) {
      std::string specname = it.first;
      Teuchos::ParameterList& spec = geochem_list->sublist(specname);
      Teuchos::RCP<TransportBoundaryFunction_Alquimia_Units> bc = Teuchos::rcp(
        new TransportBoundaryFunction_Alquimia_Units(spec, mesh_, chem_pk_, chem_engine_));

      std::vector<int>& tcc_index = bc->tcc_index();
      std::vector<std::string>& tcc_names = bc->tcc_names();

      for (int i = 0; i < tcc_names.size(); i++) {
        tcc_index.push_back(FindComponentNumber_(tcc_names[i]));
      }

      bcs_.push_back(bc);
    }
#endif

  } else {
    if (vo_->os_OK(Teuchos::VERB_NONE)) {
      *vo_->os() << vo_->color("yellow") << "No BCs were specified." << vo_->reset() << std::endl;
    }
  }
}


void
Transport_ATS::SetupPhysicalEvaluators_()
{
  // primary variable
  S_->Require<CompositeVector, CompositeVectorSpace>(key_, tag_next_, passwd_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_components_);
  S_->GetRecordSetW(key_).set_subfieldnames(component_names_);

  // -- water flux
  requireEvaluatorAtNext(flux_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  // -- water state
  requireEvaluatorAtNext(lwc_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  requireEvaluatorAtCurrent(lwc_key_, tag_current_, *S_, name_);

  if (!molar_dens_key_.empty()) {
    requireEvaluatorAtNext(molar_dens_key_, tag_next_, *S_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    requireEvaluatorAtCurrent(molar_dens_key_, tag_current_, *S_);
  }

  // CellVolume it may not be used in this PK, but having it makes vis nicer
  requireEvaluatorAtNext(cv_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // Need to figure out primary vs secondary -- are both in component names? --ETC
  std::vector<std::string> subfield_names(num_aqueous_);
  for (int i = 0; i != num_aqueous_; ++i) subfield_names[i] = component_names_[i];
  requireEvaluatorAtNext(solid_residue_mass_key_, tag_next_, *S_, name_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, num_components_);
  S_->GetRecordSetW(solid_residue_mass_key_).set_subfieldnames(subfield_names);

  // This vector stores the conserved amount (in mols) of num_components_
  // transported solutes, plus two for water.
  //
  // The first water component is given by dt * all water fluxes for which the
  // resulting component flux is treated implicitly (notably just
  // DomainCoupling fluxes like infiltration, which must be able to take all
  // the transported quantity.)  This is used to invert for the new
  // concentration.
  //
  // The second water component is given by dt * all water fluxes for which the
  // resulting component flux is treated explicitly (advection, most
  // sources/sinks of water).  Then, we can compute:
  //
  //   (lwc_new - lwc_old) - (Q_implicit + Q_explicit) * dt
  //
  // The result is a "water balance error" -- in many problems it will be
  // zero, but in problems where there are sources/sinks of water with no
  // corresponding transport of chemical species, it will be this amount.
  // Examples of this include evaporation, freezing of liquid --> solid water,
  // etc.
  //

  // Note that component_names includes secondaries, but we only need primaries
  subfield_names.emplace_back("H2O_sources_implicit");
  subfield_names.emplace_back("H2O_sources_explicit");
  requireEvaluatorAtNext(conserve_qty_key_, tag_next_, *S_, name_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, num_aqueous_ + 2);
  S_->GetRecordSetW(conserve_qty_key_).set_subfieldnames(subfield_names);

  // -- solute mass flux
  requireEvaluatorAtNext(mass_flux_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->SetComponent("face", AmanziMesh::Entity_kind::FACE, num_aqueous_);  
}


//
// Set initial conditions
//
void
Transport_ATS::Initialize()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  PK_Physical_Default::Initialize();

  // initialize missed fields
  InitializeFields_();

  // Move to Setup() with other sources? --ETC
  //
  // This must be called after S_->setup() since "water_source" data not
  // created before this step. --PL
  if (plist_->isSublist("source terms")) {
    auto sources_list = Teuchos::sublist(plist_, "source terms");
    if (sources_list->isSublist("geochemical")) {
#ifdef ALQUIMIA_ENABLED
      // -- try geochemical conditions
      auto geochem_list = Teuchos::sublist(sources_list, "geochemical");

      for (const auto& it : *geochem_list) {
        std::string specname = it.first;
        Teuchos::ParameterList& spec = geochem_list->sublist(specname);
        Teuchos::RCP<TransportSourceFunction_Alquimia_Units> src = Teuchos::rcp(
          new TransportSourceFunction_Alquimia_Units(spec, mesh_, chem_pk_, chem_engine_));

        if (S_->HasEvaluator(geochem_src_factor_key_, tag_next_)) {
          S_->GetEvaluator(geochem_src_factor_key_, tag_next_).Update(*S_, name_);
        }

        auto src_factor =
          S_->Get<CompositeVector>(geochem_src_factor_key_, tag_next_).ViewComponent("cell", false);
        src->set_conversion(-1000., src_factor, false);

        for (const auto& n : src->tcc_names()) { src->tcc_index().push_back(FindComponentNumber_(n)); }
        srcs_.push_back(src);
      }
#endif
    }
  }

  // can also now setup the joint diffusion/dispersion workspace tensor
  int D_rank = -1;
  int D_dim = mesh_->getSpaceDimension();
  if (has_diffusion_) D_rank = 1; // scalar
  if (has_dispersion_) {
    // dispersion rank is 1 or 2
    D_rank = S_->Require<TensorVector, TensorVector_Factory>(dispersion_tensor_key_, tag_next_).get_rank();
  }
  if (D_rank >= 0) {
    CompositeVectorSpace D_space;
    D_space.SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    D_ = Teuchos::rcp(new TensorVector(D_space, D_dim, D_rank, false));
  }

  // compute the stable dt for the initial timestep
  dt_stable_ = ComputeStableTimeStep_();

  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    *vo_->os() << "Number of components: " << num_components_ << std::endl
               << "  aqueous: " << num_aqueous_ << std::endl << "    ";
    for (int i = 0; i != num_aqueous_; ++i) *vo_->os() << component_names_[i] << ", ";

    *vo_->os() << "cfl=" << cfl_ << " spatial/temporal discretization: " << adv_spatial_disc_order_
               << " " << temporal_disc_order_ << std::endl << std::endl;
  }
}


/* ******************************************************************
* Initalized fields left by State and other PKs.
****************************************************************** */
void
Transport_ATS::InitializeFields_()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  S_->GetW<CompositeVector>(solid_residue_mass_key_, tag_next_, name_).PutScalar(0.0);
  S_->GetRecordW(solid_residue_mass_key_, tag_next_, name_).set_initialized();

  // initialize conserved quantity
  S_->GetEvaluator(lwc_key_, tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& lwc = *S_->Get<CompositeVector>(lwc_key_, tag_next_).ViewComponent("cell", false);
  const Epetra_MultiVector& tcc = *S_->Get<CompositeVector>(key_, tag_next_).ViewComponent("cell", false);
  Epetra_MultiVector& conserve_qty = *S_->GetW<CompositeVector>(conserve_qty_key_, tag_next_, name_).ViewComponent("cell", false);
  for (int i = 0; i != num_aqueous_; ++i) {
    conserve_qty(i)->Multiply(1., *lwc(0), *tcc(i), 0.);
  }
  S_->GetRecordW(conserve_qty_key_, tag_next_, name_).set_initialized();
}


/* *******************************************************************
* Estimation of the timestep based on T.Barth (Lecture Notes
* presented at VKI Lecture Series 1994-05, Theorem 4.2.2.
* Routine must be called every time we update a flow field.
*
* Warning: Barth calculates influx, we calculate outflux. The methods
* are equivalent for divergence-free flows and guarantee EMP. Outflux
* takes into account sinks and sources but preserves only positivity
* of an advected mass.
* ***************************************************************** */
double
Transport_ATS::ComputeStableTimeStep_()
{
  double dt = TRANSPORT_LARGE_TIME_STEP;

  // Get flux at faces for time NEXT
  IdentifyUpwindCells_();

  int ncells_owned = S_->GetMesh(domain_)->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  // flux at next tag
  S_->GetEvaluator(flux_key_, tag_next_).Update(*S_, name_);
  const CompositeVector& flux_cv = S_->Get<CompositeVector>(flux_key_, tag_next_);
  flux_cv.ScatterMasterToGhosted();
  const Epetra_MultiVector& flux = *flux_cv.ViewComponent("face", true);

  // extensive liquid water content at start and end of step
  S_->GetEvaluator(lwc_key_, tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& lwc_old = *S_->Get<CompositeVector>(lwc_key_, tag_current_).ViewComponent("cell");
  const Epetra_MultiVector& lwc_new = *S_->Get<CompositeVector>(lwc_key_, tag_next_).ViewComponent("cell");

  // loop over ALL faces and accumulate outgoing fluxes from each OWNED cell
  std::vector<double> total_outflux(ncells_owned, 0.0);

  for (int f = 0; f < flux.MyLength(); f++) {
    int c = (*upwind_cell_)[f];
    if (c >= 0 && c < ncells_owned) {
      total_outflux[c] += std::abs(flux[0][f]);
    }
  }

  // loop over cells and calculate minimal timestep
  double min_dt_lwc = 0.;
  double min_dt_outflux = 0.;
  int min_dt_cell = 0;

  for (int c = 0; c != ncells_owned; ++c) {
    double outflux = total_outflux[c];
    double min_lwc = std::min(lwc_old[0][c], lwc_new[0][c]);
    double dt_cell = TRANSPORT_LARGE_TIME_STEP;

    if (outflux > 0 && min_lwc > 0) {
      dt_cell = min_lwc / outflux;
    }

    if (dt_cell < dt) {
      dt = dt_cell;
      min_dt_cell = c;
      min_dt_lwc = min_lwc;
      min_dt_outflux = total_outflux[c];
    }
  }

  if (adv_spatial_disc_order_ == 2) dt /= 2;

  // communicate global timestep
  double dt_tmp = dt;
  auto comm = mesh_->getComm();
  comm->MinAll(&dt_tmp, &dt, 1);

  // incorporate developers and CFL constraints
  dt = std::min(dt, dt_max_);
  dt *= cfl_;
  return dt;
}


/* *******************************************************************
* Estimate returns last timestep unless it is zero.
******************************************************************* */
double
Transport_ATS::get_dt()
{
  return dt_stable_;
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
  double dt = t_new - t_old;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_LOW))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << t_old
               << " t1 = " << t_new << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;
  AMANZI_ASSERT(std::abs(S_->get_time(tag_current_) - t_old) < 1.e-4);
  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t_new) < 1.e-4);
  db_->WriteCellInfo(true);

  // check the stable step size again -- flow can now be computed at the
  // correct, new time, and may have changed, resulting in a smaller dt.
  //
  // By failing the step, it will get repeated with the smaller dt, again
  // potentially recomputing a flow field, but should be closer.
  dt_stable_ = ComputeStableTimeStep_();
  if (dt > dt_stable_ + 1.e-4) {
    if (vo_->os_OK(Teuchos::VERB_LOW))
      *vo_->os() << "Failed step: requested dt = " << dt << " > stable dt = " << dt_stable_ << std::endl;
    return true;
  }

  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Water state:" << std::endl;
  std::vector<std::string> vnames{ "lwc_old", "lwc_new" };
  std::vector<Teuchos::Ptr<const CompositeVector>> vecs{
    S_->GetPtr<CompositeVector>(lwc_key_, tag_current_).ptr(),
    S_->GetPtr<CompositeVector>(lwc_key_, tag_next_).ptr(),
  };
  db_->WriteVectors(vnames, vecs);
  db_->WriteVector("mol_ratio_old",
                   S_->GetPtr<CompositeVector>(key_, tag_current_).ptr(),
                   true,
                   S_->GetRecordSet(key_).subfieldnames());

  // update geochemical condition conversions
#ifdef ALQUIMIA_ENABLED
  if (plist_->sublist("source terms").isSublist("geochemical")) {
    for (auto& src : srcs_) {
      if (src->getType() == DomainFunction_kind::ALQUIMIA) {
        // src_factor = water_source / molar_density_liquid, both flow
        // quantities, see note above.
        S_->GetEvaluator(geochem_src_factor_key_, tag_next_).Update(*S_, name_);
        auto src_factor = S_->Get<CompositeVector>(geochem_src_factor_key_, tag_next_)
                            .ViewComponent("cell", false);
        Teuchos::RCP<TransportSourceFunction_Alquimia_Units> src_alq =
          Teuchos::rcp_dynamic_cast<TransportSourceFunction_Alquimia_Units>(src);
        src_alq->set_conversion(-1000, src_factor, false);
      }
    }
  }

  if (plist_->sublist("boundary conditions").isSublist("geochemical")) {
    for (auto& bc : bcs_) {
      if (bc->getType() == DomainFunction_kind::ALQUIMIA) {
        Teuchos::RCP<TransportBoundaryFunction_Alquimia_Units> bc_alq =
          Teuchos::rcp_dynamic_cast<TransportBoundaryFunction_Alquimia_Units>(bc);

        S_->GetEvaluator(molar_dens_key_, tag_next_).Update(*S_, name_);
        auto molar_dens = S_->Get<CompositeVector>(molar_dens_key_, tag_next_)
          .ViewComponent("cell", false);
        bc_alq->set_conversion(1000.0, molar_dens, true);
      }
    }
  }
#endif

  // update boundary conditions
  for (int i = 0; i < bcs_.size(); i++) {
    bcs_[i]->Compute(t_old, t_new);
  }

  // advance advection term
  if (temporal_disc_order_ == 1) {
    AdvanceAdvectionSources_RK1_(t_old, t_new, adv_spatial_disc_order_);
  } else if (temporal_disc_order_ == 2) {
    AdvanceAdvectionSources_RK2_(t_old, t_new, adv_spatial_disc_order_);
  } else {
    AMANZI_ASSERT(false);
  }

  // advance dispersion and diffusion term
  AdvanceDispersionDiffusion_(t_old, t_new);

  // statistics output
  const Epetra_MultiVector& tcc_new = *S_->Get<CompositeVector>(key_, tag_next_)
    .ViewComponent("cell");
  ChangedSolutionPK(tag_next_);
  db_->WriteVector("mol_ratio_new",
                   S_->GetPtr<CompositeVector>(key_, tag_next_).ptr(),
                   true,
                   S_->GetRecordSet(key_).subfieldnames());
  return failed;
}


void
Transport_ATS ::AdvanceDispersionDiffusion_(double t_old, double t_new)
{
  if (!has_diffusion_ && !has_dispersion_) return;
  double dt = t_new - t_old;

  Epetra_MultiVector& tcc_new = *S_->GetW<CompositeVector>(key_, tag_next_, passwd_)
    .ViewComponent("cell", false);

  // needed for diffusion coefficent and for accumulation term
  const Epetra_MultiVector& lwc = *S_->Get<CompositeVector>(lwc_key_, tag_next_)
    .ViewComponent("cell", false);
  const Epetra_MultiVector& cv = *S_->Get<CompositeVector>(cv_key_, tag_next_)
    .ViewComponent("cell", false);

  //
  // first step  -- aqueous dispersion + diffusion
  //
  if (has_dispersion_) {
    //
    // The dispersion tensor is consistent for all aqueous phase components.
    S_->GetEvaluator(dispersion_tensor_key_, tag_next_).Update(*S_, name_);
    (*D_) = S_->Get<TensorVector>(dispersion_tensor_key_, tag_next_);

    // scale by (volumetric) liquid water content
    for (int c = 0; c != lwc.MyLength(); ++c) {
      (*D_)[c] *= (lwc[0][c] / cv[0][c]);
    }
  } else {
    D_->PutScalar(0.);
  }

  // we track only the difference in molecular diffusion from component to
  // component -- if they don't exist or stay the same, the matrices do not
  // change and can be reused.
  double md_old(0.0);

  for (int i = 0; i != num_aqueous_; ++i) {
    // add molecular diffusion to the dispersion tensor
    bool changed_tensor(false);
    if (has_diffusion_) {
      double md_new = molec_diff_[component_names_[i]] * tortuosity_["aqueous"];
      double md_change = md_new - md_old;
      if (std::abs(md_change) > 1.e-12) {
        // shift the tensor diagonal by the lwc * delta diff coef
        for (int c = 0; c != cv.MyLength(); ++c) {
          (*D_)[c] += md_change * lwc[0][c] / cv[0][c];
        }
        md_old = md_new;
        changed_tensor = true;
      }
    }

    //
    // apply the diffusion operator
    //
    // -- set initial guess
    Epetra_MultiVector& diff_sol_cell = *diff_sol_->ViewComponent("cell", false);
    *diff_sol_cell(0) = *tcc_new(i);
    if (diff_sol_->HasComponent("face")) diff_sol_->ViewComponent("face", false)->PutScalar(0.0);

    // -- build the matrices if needed
    if (changed_tensor || i == 0) {
      // update mass, stiffness matrices of diffusion operator
      diff_global_op_->Init();
      diff_op_->SetTensorCoefficient(Teuchos::rcpFromRef(D_->data));
      diff_op_->UpdateMatrices(Teuchos::null, Teuchos::null);

      // add accumulation term
      diff_acc_op_->AddAccumulationTerm(S_->Get<CompositeVector>(lwc_key_, tag_next_),
              dt, "cell", false);
    }

    // whether or not diffusion operator is changed, RHS is different
    Epetra_MultiVector& rhs_cell = *diff_global_op_->rhs()->ViewComponent("cell");
    for (int c = 0; c < rhs_cell.MyLength(); c++) {
      rhs_cell[0][c] = lwc[0][c] / dt * tcc_new[i][c];
    }
    // apply BCs -- must always do
    diff_op_->ApplyBCs(true, true, true);

    // -- apply the inverse
    CompositeVector& rhs = *diff_global_op_->rhs();
    int ierr = diff_global_op_->ApplyInverse(rhs, *diff_sol_);

    if (ierr != 0) {
      Errors::Message msg("TransportExplicit_PK solver failed with message: \"");
      msg << diff_global_op_->returned_code_string() << "\"";
      Exceptions::amanzi_throw(msg);
    }

    // -- copy back the solution
    *tcc_new(i) = *(*diff_sol_->ViewComponent("cell", false))(0);

  }

  // NOTE: here was gas diffusion! --ETC
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

  PK_Physical_Default::CommitStep(t_old, t_new, tag_next);

  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;
  assign(lwc_key_, tag_current, tag_next, *S_);
}


/* *******************************************************************
 * Advance RK1
 ****************************************************************** */
void
Transport_ATS::AdvanceAdvectionSources_RK1_(double t_old,
        double t_new,
        int spatial_order)
{
  double dt = t_new - t_old;

  // distribute vector of concentrations
  // scattering total concentration from master to others
  S_->Get<CompositeVector>(key_, tag_current_).ScatterMasterToGhosted();
  const Epetra_MultiVector& tcc_old = *S_->Get<CompositeVector>(key_, tag_current_)
    .ViewComponent("cell", true);

  // old and new water contents -- note these were updated in StableStep
  const Epetra_MultiVector& lwc_old = *S_->Get<CompositeVector>(lwc_key_, tag_current_)
    .ViewComponent("cell", false);

  // populating conserved quantity (unit: molC)
  // The conserve_qty (M) has `num_components_+2` vectors, the extras for tracking water
  Epetra_MultiVector& conserve_qty =
    *S_->GetW<CompositeVector>(conserve_qty_key_, tag_next_, name_)
    .ViewComponent("cell", false);

  // populating solute mass flux (unit: molC)
  Epetra_MultiVector& mass_flux = 
    *S_->GetW<CompositeVector>(mass_flux_key_, tag_next_, name_)
    .ViewComponent("face", false);  

  // Mass <-- Concentration * water content
  // conserved component quantity [mol C] = (mol C / mol H20) * mol H2O
  for (int i = 0; i != num_aqueous_; ++i) {
    conserve_qty(i)->Multiply(1., *lwc_old(0), *tcc_old(i), 0.);
  }
  conserve_qty(num_aqueous_)->PutScalar(0.);
  conserve_qty(num_aqueous_ + 1)->PutScalar(0.);
  db_->WriteCellVector("qnty (old)", conserve_qty,
                       S_->GetRecordSet(conserve_qty_key_).subfieldnames());

  // advection: M <-- M + dt * div q * C0
  if (spatial_order == 1) {
    // conserve_qty = conserve_qty + div qC
    AddAdvection_FirstOrderUpwind_(t_old, t_new, tcc_old, conserve_qty, mass_flux);
  } else if (spatial_order == 2) {
    // conserve_qty = conserve_qty + div qC
    AddAdvection_SecondOrderUpwind_(t_old, t_new, tcc_old, conserve_qty, mass_flux);
  }
  db_->WriteCellVector("qnty (adv)", conserve_qty,
                       S_->GetRecordSet(conserve_qty_key_).subfieldnames());

  // process external sources: M <-- M + dt * Q
  AddSourceTerms_(t_old, t_new, conserve_qty, 0, num_aqueous_ - 1);
  db_->WriteCellVector("qnty (src)", conserve_qty,
                       S_->GetRecordSet(conserve_qty_key_).subfieldnames());

  // invert for C1: C1 <-- M / WC1, also deals with dissolution/precipitation
  // tcc_new, the new solution
  Epetra_MultiVector& tcc_new = *S_->GetW<CompositeVector>(key_, tag_next_, passwd_)
    .ViewComponent("cell", false);

  // solid quantity (unit: molC) stores extra solute mass
  // solid_qty has `num_components_` vectors only (solute mass)
  Epetra_MultiVector& solid_qty =
    *S_->GetW<CompositeVector>(solid_residue_mass_key_, tag_next_, name_)
    .ViewComponent("cell", false);

  InvertTccNew_(conserve_qty, tcc_new, &solid_qty, true);
  db_->WriteCellVector("mol_ratio (pre-diff)", tcc_new,
                       S_->GetRecordSet(key_).subfieldnames());
}


/* *******************************************************************
 * We have to advance each component independently due to different
 * reconstructions.
 ****************************************************************** */
void
Transport_ATS::AdvanceAdvectionSources_RK2_(double t_old,
        double t_new,
        int spatial_order)
{
  double dt = t_new - t_old;

  // distribute vector of concentrations
  // scattering total concentration from master to others
  S_->Get<CompositeVector>(key_, tag_current_).ScatterMasterToGhosted();
  const Epetra_MultiVector& tcc_old = *S_->Get<CompositeVector>(key_, tag_current_)
    .ViewComponent("cell", true);

  // old and new water contents -- note these were updated in StableStep
  const Epetra_MultiVector& lwc_old = *S_->Get<CompositeVector>(lwc_key_, tag_current_)
    .ViewComponent("cell", false);

  // populating conserved quantity (unit: molC)
  // The conserve_qty (M) has `num_components_+2` vectors, the extras for tracking water
  Epetra_MultiVector& conserve_qty =
    *S_->GetW<CompositeVector>(conserve_qty_key_, tag_next_, name_)
    .ViewComponent("cell", false);

  // populating solute mass flux (unit: molC)
  Epetra_MultiVector& mass_flux = 
    *S_->GetW<CompositeVector>(mass_flux_key_, tag_next_, name_)
    .ViewComponent("face", false);

  // -- M <-- C0 * W0
  // -- mol H2O * (mol C / mol H20) --> mol C, the conserved component quantity
  for (int i = 0; i != num_aqueous_; ++i) {
    conserve_qty(i)->Multiply(1., *lwc_old(0), *tcc_old(i), 0.);
  }
  conserve_qty(num_aqueous_)->PutScalar(0.);
  conserve_qty(num_aqueous_ + 1)->PutScalar(0.);
  db_->WriteCellVector("qnty (start)", conserve_qty,
                       S_->GetRecordSet(conserve_qty_key_).subfieldnames());

  // Predictor Step:
  // -- advection: M <-- M + dt * div q * C0
  if (spatial_order == 1) {
    AddAdvection_FirstOrderUpwind_(t_old, t_new, tcc_old, conserve_qty, mass_flux);
  } else if (spatial_order == 2) {
    AddAdvection_SecondOrderUpwind_(t_old, t_new, tcc_old, conserve_qty, mass_flux);
  }
  db_->WriteCellVector("qnty (pred adv)", conserve_qty,
                       S_->GetRecordSet(conserve_qty_key_).subfieldnames());

  // -- process external sources: M <-- M + dt * Q(t0)
  AddSourceTerms_(t_old, t_new, conserve_qty, 0, num_aqueous_ - 1);
  db_->WriteCellVector("qnty (pred src)", conserve_qty,
                       S_->GetRecordSet(conserve_qty_key_).subfieldnames());

  // -- invert for C': C' <-- M / WC1, note no dissolution/precip
  {
    Epetra_MultiVector& tcc_new = *S_->GetW<CompositeVector>(key_, tag_next_, passwd_)
      .ViewComponent("cell", false);

    InvertTccNew_(conserve_qty, tcc_new, nullptr, false);
    db_->WriteCellVector("mol_ratio (pred)", tcc_new,
                       S_->GetRecordSet(key_).subfieldnames());
  }

  // Corrector Step:
  // -- M <-- M/2 + M0 C0 / 2
  //      <-- (C0 W0 + dt div q C0 + dt Q0) / 2 + W0 C0 / 2
  for (int i = 0; i != num_aqueous_; ++i) {
    conserve_qty(i)->Multiply(0.5, *lwc_old(0), *tcc_old(i), 0.5);
  }
  conserve_qty(num_aqueous_)->PutScalar(0.);
  conserve_qty(num_aqueous_ + 1)->PutScalar(0.);
  db_->WriteCellVector("qnty (corr start)", conserve_qty,
                       S_->GetRecordSet(conserve_qty_key_).subfieldnames());

  // -- advect the predicted C'
  {
    S_->Get<CompositeVector>(key_, tag_next_).ScatterMasterToGhosted();
    const Epetra_MultiVector& tcc_new = *S_->Get<CompositeVector>(key_, tag_next_)
      .ViewComponent("cell", true);

    //   M <-- M + dt/2 * div q * C'
    //     <-- (C0 W0 + dt div q C0 + dt Q0) / 2 + W0 C0 / 2 + dt div q C' / 2
    if (spatial_order == 1) {
      AddAdvection_FirstOrderUpwind_(t_old + dt/2., t_new, tcc_new, conserve_qty, mass_flux);
    } else if (spatial_order == 2) {
      AddAdvection_SecondOrderUpwind_(t_old + dt/2., t_new, tcc_new, conserve_qty, mass_flux);
    }
    db_->WriteCellVector("qnty (corr adv)", conserve_qty,
                         S_->GetRecordSet(conserve_qty_key_).subfieldnames());

  }

  // -- add sources at the new time, predicted C'
  //   M <-- M + dt/2 Q(t1)
  //     <-- (C0 W0 + dt div q C0 + dt Q(t0)) / 2 + W0 C0 / 2 + dt div q C' / 2 + dt Q(t1) / 2
  //     <-- C0 W0 + dt * (div q C0 + div q C') / 2 + dt (Q(t0) + Q(t1)) / 2
  AddSourceTerms_(t_old + dt/2., t_new, conserve_qty, 0, num_aqueous_ - 1);
  db_->WriteCellVector("qnty (corr src)", conserve_qty,
                       S_->GetRecordSet(conserve_qty_key_).subfieldnames());


  // -- Invert to get C1, this time with dissolution/solidification
  {
    // solid quantity (unit: molC) stores extra solute mass
    // solid_qty has `num_components_` vectors only (solute mass)
    Epetra_MultiVector& solid_qty =
      *S_->GetW<CompositeVector>(solid_residue_mass_key_, tag_next_, name_)
      .ViewComponent("cell", false);

    Epetra_MultiVector& tcc_new = *S_->GetW<CompositeVector>(key_, tag_next_, passwd_)
      .ViewComponent("cell", true);

    InvertTccNew_(conserve_qty, tcc_new, &solid_qty, true);
    db_->WriteCellVector("mol_ratio_new", tcc_new,
                         S_->GetRecordSet(key_).subfieldnames());
  }
}



/* *******************************************************************
* Populates operators' boundary data for given component.
* Returns true if at least one face was populated.
******************************************************************* */
bool
Transport_ATS::PopulateBoundaryData_(int component, Operators::BCs& bc)
{
  int nfaces_all =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
  bool flag = false;

  auto& bc_model = bc.bc_model();
  auto& bc_value = bc.bc_value();

  for (int i = 0; i < bc_model.size(); i++) {
    bc_model[i] = Operators::OPERATOR_BC_NONE;
    bc_value[i] = 0.0;
  }

  // set the default BC
  for (int f = 0; f < nfaces_all; f++) {
    auto cells = mesh_->getFaceCells(f);
    if (cells.size() == 1) bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
  }

  if (component >= 0) {
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
  }

  return flag;
}


/* *******************************************************************
* Identify flux direction based on orientation of the face normal
* and sign of the  Darcy velocity.
******************************************************************* */
void
Transport_ATS::IdentifyUpwindCells_()
{
  int ncells_all =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  upwind_cell_->PutValue(-1);
  downwind_cell_->PutValue(-1);

  S_->Get<CompositeVector>(flux_key_, tag_next_).ScatterMasterToGhosted("face");
  const Epetra_MultiVector& flux = *S_->Get<CompositeVector>(flux_key_, tag_next_).ViewComponent("face", true);

  // identify upwind and downwind cell of each face
  for (int c = 0; c != ncells_all; ++c) {
    const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);

    for (int i = 0; i != faces.size(); ++i) {
      int f = faces[i];
      double tmp = flux[0][f] * dirs[i];
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


} // namespace Transport
} // namespace Amanzi
