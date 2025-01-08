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

#include "BCs.hh"
#include "errors.hh"
#include "Explicit_TI_RK.hh"
#include "Evaluator.hh"
#include "Mesh.hh"
#include "MeshAlgorithms.hh"
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
    PK_PhysicalExplicit<Epetra_Vector>(pk_tree, global_plist, S, solution),
    has_water_src_key_(false),
    flow_tag_(Tags::NEXT),
    passwd_("state"),
    internal_tests(0),
    tests_tolerance(TRANSPORT_CONCENTRATION_OVERSHOOT),
    dt_(0.0),
    dt_max_(0.0),
    t_physics_(0.0),
    current_component_(-1),
    flag_dispersion_(false)
{
  // initialize io
  units_.Init(global_plist->sublist("units"));
}

void
Transport_ATS::parseParameterList()
{
  if (plist_->isParameter("component molar masses")) {
    mol_masses_ = plist_->get<Teuchos::Array<double>>("component molar masses").toVector();
  } else {
    Errors::Message msg("Transport PK: parameter \"component molar masses\" is missing.");
    Exceptions::amanzi_throw(msg);
  }

  // primary variable
  tcc_key_ = Keys::readKey(*plist_, domain_, "concentration", "total_component_concentration");
  requireAtNext(tcc_key_, tag_next_, *S_, passwd_);
  requireAtCurrent(tcc_key_, tag_current_, *S_, passwd_);

  // keys, dependencies, etc
  saturation_key_ = Keys::readKey(*plist_, domain_, "saturation liquid", "saturation_liquid");
  flux_key_ = Keys::readKey(*plist_, domain_, "water flux", "water_flux");
  permeability_key_ = Keys::readKey(*plist_, domain_, "permeability", "permeability");
  conserve_qty_key_ =
    Keys::readKey(*plist_, domain_, "conserved quantity", "total_component_quantity");
  porosity_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
  molar_density_key_ = Keys::readKey(*plist_, domain_, "molar density", "molar_density_liquid");
  tcc_matrix_key_ =
    Keys::readKey(*plist_, domain_, "tcc matrix", "total_component_concentration_matrix");
  solid_residue_mass_key_ = Keys::readKey(*plist_, domain_, "solid residue", "solid_residue_mass");
  water_src_key_ = Keys::readKey(*plist_, domain_, "water source", "water_source");
  water_src_tile_key_ = Keys::readKey(*plist_, domain_, "water source tile", "water_source_tile");
  geochem_src_factor_key_ =
    Keys::readKey(*plist_, domain_, "geochem source factor", "geochem_src_factor");
  water_content_key_ = Keys::readKey(*plist_, domain_, "water content", "water_content");
  cv_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");
  key_ = tcc_key_;

  // other parameters
  water_tolerance_ = plist_->get<double>("water tolerance", 1e-6);
  dissolution_ = plist_->get<bool>("allow dissolution", false);
  max_tcc_ = plist_->get<double>("maximum concentration", -1.0);

  // global transport parameters
  cfl_ = plist_->get<double>("cfl", 1.0);
  dt_max_ = plist_->get<double>("maximum timestep", TRANSPORT_LARGE_TIME_STEP);

  spatial_disc_order_ = plist_->get<int>("spatial discretization order", 1);
  if (spatial_disc_order_ < 1 || spatial_disc_order_ > 2) {
    Errors::Message msg;
    msg << "Transport_ATS: \"spatial discretization order\" must be 1 or 2, not " << spatial_disc_order_;
    Exceptions::amanzi_throw(msg);
  }
  temporal_disc_order_ = plist_->get<int>("temporal discretization order", 1);
  if (temporal_disc_order_ < 1 || temporal_disc_order_ > 2) {
    Errors::Message msg;
    msg << "Transport_ATS: \"temporal discretization order\" must be 1 or 2, not " << temporal_disc_order_;
    Exceptions::amanzi_throw(msg);
  }


  if (plist_->isParameter("runtime diagnostics: regions")) {
    runtime_regions_ =
      plist_->get<Teuchos::Array<std::string>>("runtime diagnostics: regions").toVector();
  }

  internal_tests = plist_->get<bool>("enable internal tests", false);
  tests_tolerance =
    plist_->get<double>("internal tests tolerance", TRANSPORT_CONCENTRATION_OVERSHOOT);

  db_ = Teuchos::rcp(new Debugger(mesh_, name_, *plist_));
}

void
Transport_ATS::set_tags(const Tag& current, const Tag& next)
{
  PK_PhysicalExplicit<Epetra_Vector>::set_tags(current, next);
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
  SetupTransport_();
  SetupPhysicalEvaluators_();
}


/* ******************************************************************
* Define structure of this PK.
****************************************************************** */
void
Transport_ATS::SetupTransport_()
{
  // a few last things that ought to be in the constructor, but cannot be there
  // because they depend upon component_names, which may not get set in the
  // constructor (as they get set by the chemical engine).
  if (component_names_.size() == 0) {
    // not set by chemistry... must get set by user
    component_names_ = plist_->get<Teuchos::Array<std::string>>("component names").toVector();
    num_components_ = component_names_.size();
  }
  num_aqueous_ = plist_->get<int>("number of aqueous components", component_names_.size());
  num_advect_ = plist_->get<int>("number of aqueous components advected", num_aqueous_);
  num_gaseous_ = plist_->get<int>("number of gaseous components", 0);

  tcc_max_.resize(num_aqueous_, 1.e10);
  tcc_max_ = plist_->get<Teuchos::Array<double>>("component max concentrations", tcc_max_);

  // statistics of solutes
  if (plist_->isParameter("runtime diagnostics: solute names")) {
    runtime_solutes_ =
      plist_->get<Teuchos::Array<std::string>>("runtime diagnostics: solute names").toVector();
    if (runtime_solutes_.size() == 1 && runtime_solutes_[0] == "all") {
      runtime_solutes_ = component_names_;
    }
  }
  mass_solutes_exact_.assign(num_aqueous_ + num_gaseous_, 0.0);
  mass_solutes_source_.assign(num_aqueous_ + num_gaseous_, 0.0);
  mass_solutes_bc_.assign(num_aqueous_ + num_gaseous_, 0.0);
  mass_solutes_stepstart_.assign(num_aqueous_ + num_gaseous_, 0.0);

  // cross-coupling of PKs
  Teuchos::RCP<Teuchos::ParameterList> physical_models =
    Teuchos::sublist(plist_, "physical models and assumptions");

  if (num_components_ == 0) {
    Errors::Message msg("Transport PK: list of solutes is empty.");
    Exceptions::amanzi_throw(msg);
  }

  // upwind
  const Epetra_Map& fmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::FACE, true);
  upwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));
  downwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));

  // reconstruction initialization
  limiter_ = Teuchos::rcp(new Operators::LimiterCell(mesh_));
  lifting_ = Teuchos::rcp(new Operators::ReconstructionCellLinear(mesh_));

  // dependencies:
  // -- permeability
  bool abs_perm = physical_models->get<bool>("permeability field is required", false);
  if (abs_perm) {
    requireAtNext(permeability_key_, tag_next_, *S_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, mesh_->getSpaceDimension());
  }

  // HACK ALERT -- FIXME --ETC
  // This PK is liberally sprinkled with hard-coded Tags::NEXT and
  // Tags::CURRENT, forcing all things provided by FLOW to be provided at that
  // tag and not at tag_current and tag_next as it should be.  This is because
  // we don't have a good way of aliasing everything we need yet.  In
  // particular, aliases needed to be introduced between Setup() on flow and
  // Setup() on transport, and this was not possible when the quantity of
  // interest (porosity)'s evaluator was not required directly (only
  // indirectly) in flow PK.


  // source term setup
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

          convert_to_field_[name] = src_list->get<bool>("convert to field", false);
          std::string src_type = src_list->get<std::string>("spatial distribution method", "none");

          if (src_type == "domain coupling") {
            Teuchos::RCP<TransportDomainFunction> src =
              factory.Create(*src_list, "fields", AmanziMesh::Entity_kind::CELL, Kxy_, tag_current_);

            // domain couplings functions is special -- always work on all components
            for (int i = 0; i < num_components_; i++) {
              src->tcc_names().push_back(component_names_[i]);
              src->tcc_index().push_back(i);
            }
            src->set_state(S_);
            srcs_.push_back(src);
          } else if (src_type == "field") {
            Teuchos::RCP<TransportDomainFunction> src =
              factory.Create(*src_list, "field", AmanziMesh::Entity_kind::CELL, Kxy_, tag_current_);

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
                  auto field_key = flist.sublist(sublist_name.str()).get<std::string>("field key");
                  auto field_tag = Keys::readTag(flist.sublist(sublist_name.str()), "tag");
                  requireAtNext(field_key, field_tag, *S_)
                    .SetMesh(mesh_)
                    ->SetGhosted(true)
                    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
                }
              }
            } else { // if only one field
              auto field_key = src_list->sublist("field").get<std::string>("field key");
              auto field_tag = Keys::readTag(src_list->sublist("field"), "tag");
              requireAtNext(field_key, field_tag, *S_)
                .SetMesh(mesh_)
                ->SetGhosted(true)
                ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_components_);
              // NOTE: this code should be moved to the PK_DomainFunctionField, and all other
              // PK_DomainFunction* should be updated to make sure they require their data! --ETC
              // Problem: This requires the pk_helper.hh to be included in PK_DomainFunctionField.hh --PL
            }

          } else { // all others work on a subset of components
            Teuchos::RCP<TransportDomainFunction> src = factory.Create(
              *src_list, "source function", AmanziMesh::Entity_kind::CELL, Kxy_, tag_current_);
            src->set_tcc_names(src_list->get<Teuchos::Array<std::string>>("component names").toVector());
            for (const auto& n : src->tcc_names()) {
              src->tcc_index().push_back(FindComponentNumber_(n));
            }

            src->set_state(S_);
            srcs_.push_back(src);
          }

          if (convert_to_field_[name]) {
            name = Keys::cleanName(name);
            if (Keys::getDomain(name)!=domain_){
              name = Keys::getKey(domain_, name);
            }
            requireAtNext(name, Tags::NEXT, *S_, name)
            .SetMesh(mesh_)
            ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, num_components_);
          }
        }
      }
    }

    // sources of water that include C at a known concentration
    if (sources_list->isSublist("geochemical")) {
      // note these are computed at the flow PK's NEXT tag, which assumes all
      // sources are dealt with implicitly (backward Euler).  This could be relaxed --ETC
      requireAtNext(water_src_key_, flow_tag_, *S_)
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
        std::vector<std::string> dep{ water_src_key_, molar_density_key_ };
        wc_eval.set<Teuchos::Array<std::string>>("dependencies", dep);
        wc_eval.set<std::string>("reciprocal", dep[1]);
      }
      requireAtNext(geochem_src_factor_key_, Tags::NEXT, *S_)
        .SetMesh(mesh_)
        ->SetGhosted(true)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    }
  }

  // material properties
  if (plist_->isSublist("material properties")) {
    auto mat_prop_list = Teuchos::sublist(plist_, "material properties");
    mdm_ = CreateMDMPartition(mesh_, mat_prop_list, flag_dispersion_);

    int nblocks = 0;
    for (Teuchos::ParameterList::ConstIterator i = mat_prop_list->begin(); i != mat_prop_list->end(); i++) {
      if (mat_prop_list->isSublist(mat_prop_list->name(i))) nblocks++;
    }

    mat_properties_.resize(nblocks);

    int iblock = 0;
    for (Teuchos::ParameterList::ConstIterator i = mat_prop_list->begin(); i != mat_prop_list->end(); i++) {
      if (mat_prop_list->isSublist(mat_prop_list->name(i))) {
        mat_properties_[iblock] = Teuchos::rcp(new MaterialProperties());

        Teuchos::ParameterList& model_list = mat_prop_list->sublist(mat_prop_list->name(i));

        mat_properties_[iblock]->tau[0] = model_list.get<double>("aqueous tortuosity", 0.0);
        mat_properties_[iblock]->tau[1] = model_list.get<double>("gaseous tortuosity", 0.0);
        mat_properties_[iblock]->regions =
          model_list.get<Teuchos::Array<std::string>>("regions").toVector();
        iblock++;
      }
    }
  }

  // transport diffusion (default is none)
  diffusion_phase_.resize(TRANSPORT_NUMBER_PHASES, Teuchos::null);

  if (plist_->isSublist("molecular diffusion")) {
    auto diff_list = Teuchos::sublist(plist_, "molecular diffusion");
    if (diff_list->isParameter("aqueous names")) {
      diffusion_phase_[0] = Teuchos::rcp(new DiffusionPhase());
      diffusion_phase_[0]->names() =
        diff_list->get<Teuchos::Array<std::string>>("aqueous names").toVector();
      diffusion_phase_[0]->values() =
        diff_list->get<Teuchos::Array<double>>("aqueous values").toVector();
    }

    if (diff_list->isParameter("gaseous names")) {
      diffusion_phase_[1] = Teuchos::rcp(new DiffusionPhase());
      diffusion_phase_[1]->names() =
        diff_list->get<Teuchos::Array<std::string>>("gaseous names").toVector();
      diffusion_phase_[1]->values() =
        diff_list->get<Teuchos::Array<double>>("gaseous values").toVector();
    }
  }

  // create boundary conditions
  if (plist_->isSublist("boundary conditions")) {
    // -- try tracer-type conditions
    PK_DomainFunctionFactory<TransportDomainFunction> factory(mesh_, S_);
    auto bcs_list = Teuchos::sublist(plist_, "boundary conditions");
    auto conc_bcs_list = Teuchos::sublist(bcs_list, "concentration");

    for (const auto& it : *conc_bcs_list) {
      std::string name = it.first;
      if (conc_bcs_list->isSublist(name)) {
        Teuchos::ParameterList& bc_list = conc_bcs_list->sublist(name);
        std::string bc_type = bc_list.get<std::string>("spatial distribution method", "none");

        if (bc_type == "domain coupling") {
          // domain couplings are special -- they always work on all components
          Teuchos::RCP<TransportDomainFunction> bc =
            factory.Create(bc_list, "fields", AmanziMesh::Entity_kind::FACE, Kxy_, tag_current_);

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
          bc_list.set("entity_gid_out", gid);

          Teuchos::RCP<TransportDomainFunction> bc = factory.Create(
            bc_list, "boundary concentration", AmanziMesh::Entity_kind::FACE, Kxy_, tag_current_);

          for (int i = 0; i < num_components_; i++) {
            bc->tcc_names().push_back(component_names_[i]);
            bc->tcc_index().push_back(i);
          }
          bc->set_state(S_);
          bcs_.push_back(bc);

        } else {
          Teuchos::RCP<TransportDomainFunction> bc =
            factory.Create(bc_list,
                           "boundary concentration function",
                           AmanziMesh::Entity_kind::FACE,
                           Kxy_,
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

      bc->set_conversion(1000.0, mol_dens_, true);
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
  // -- water flux
  requireAtNext(flux_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  // -- water saturation
  requireAtNext(saturation_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  requireAtCurrent(saturation_key_, Tags::CURRENT, *S_);

  requireAtNext(porosity_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  requireAtNext(molar_density_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  requireAtCurrent(molar_density_key_, Tags::CURRENT, *S_);

  requireAtNext(tcc_key_, tag_next_, *S_, passwd_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_components_);
  S_->GetRecordSetW(tcc_key_).set_subfieldnames(component_names_);
  requireAtCurrent(tcc_key_, tag_current_, *S_, passwd_);

  // CellVolume it may not be used in this PK, but having it makes vis nicer
  requireAtNext(cv_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // Need to figure out primary vs secondary -- are both in component names? --ETC
  std::vector<std::string> primary_names = component_names_;
  requireAtNext(solid_residue_mass_key_, tag_next_, *S_, name_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, num_components_);
  S_->GetRecordSetW(solid_residue_mass_key_).set_subfieldnames(primary_names);
  // This vector stores the conserved amount (in mols) of num_components_ transported solutes, plus two for water.
  // - The first water component is given by the water content (in mols) at the current (old) time plus dt * all fluxes treated explictly.
  // - The second water component is given by the water content at the new time plus dt * all fluxes treated implicitly
  //   (notably just DomainCoupling fluxes, which must be able to take all the transported quantity.)

  // Note that component_names includes secondaries, but we only need primaries
  primary_names.emplace_back("H2O_old");
  primary_names.emplace_back("H2O_new");
  requireAtNext(conserve_qty_key_, tag_next_, *S_, name_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, num_components_ + 2);
  S_->GetRecordSetW(conserve_qty_key_).set_subfieldnames(primary_names);
}


//
// Set initial conditions
//
void
Transport_ATS::Initialize()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // Set initial values for transport variables.
  if (plist_->isSublist("initial condition")) {
    S_->GetRecordW(tcc_key_, tag_next_, passwd_).Initialize(plist_->sublist("initial condition"));
  }

  // initialize missed fields
  InitializeFields_();

  // make this go away -- local pointers to data are a no-no! --ETC
  tcc_tmp = S_->GetPtrW<CompositeVector>(tcc_key_, tag_next_, passwd_);
  tcc = S_->GetPtrW<CompositeVector>(tcc_key_, tag_current_, passwd_);

  ws_ = S_->Get<CompositeVector>(saturation_key_, Tags::NEXT).ViewComponent("cell", false);
  ws_prev_ = S_->Get<CompositeVector>(saturation_key_, Tags::CURRENT).ViewComponent("cell", false);
  mol_dens_ = S_->Get<CompositeVector>(molar_density_key_, Tags::NEXT).ViewComponent("cell", false);
  mol_dens_prev_ =
    S_->Get<CompositeVector>(molar_density_key_, Tags::CURRENT).ViewComponent("cell", false);

  // upwind
  IdentifyUpwindCells_();

  // mechanical dispersion
  if (flag_dispersion_) CalculateAxiSymmetryDirection_();

  // boundary conditions initialization
  double time = t_physics_;

  const Epetra_MultiVector& flux = *S_->Get<CompositeVector>(flux_key_, Tags::NEXT).ViewComponent("face", true);
  CheckInfluxBC_(flux);

  // Move to Setup() with other sources? --ETC
  // This must be called after S_->setup() since "water_source" data not created before this step. --PL
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

        if (S_->HasEvaluator(geochem_src_factor_key_, Tags::NEXT)) {
          S_->GetEvaluator(geochem_src_factor_key_, Tags::NEXT).Update(*S_, name_);
        }

        auto src_factor =
          S_->Get<CompositeVector>(geochem_src_factor_key_, Tags::NEXT).ViewComponent("cell", false);
        src->set_conversion(-1000., src_factor, false);

        for (const auto& n : src->tcc_names()) { src->tcc_index().push_back(FindComponentNumber_(n)); }

        srcs_.push_back(src);
      }
#endif
    }
  }

  // Temporarily Transport hosts Henry law.
  PrepareAirWaterPartitioning_();

  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    *vo_->os() << "Number of components: " << tcc->size() << std::endl
               << "cfl=" << cfl_ << " spatial/temporal discretization: " << spatial_disc_order_
               << " " << temporal_disc_order_ << std::endl;
    *vo_->os() << vo_->color("green") << "Initalization of PK is complete." << vo_->reset()
               << std::endl
               << std::endl;
  }
  ComputeStableTimeStep_();
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
  S_->GetW<CompositeVector>(conserve_qty_key_, tag_next_, name_).PutScalar(0.0);
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
  // Get flux at faces for time NEXT
  IdentifyUpwindCells_();

  // Total concentration at current tag
  tcc = S_->GetPtrW<CompositeVector>(tcc_key_, tag_current_, passwd_);
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");

  // flux at flow tag
  const Epetra_MultiVector& flux = *S_->Get<CompositeVector>(flux_key_, Tags::NEXT).ViewComponent("face", true);

  int ncells_all =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  int nfaces_all =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  // loop over faces and accumulate upwinding fluxes
  std::vector<double> total_outflux(ncells_all, 0.0);

  for (int f = 0; f < nfaces_all; f++) {
    int c = (*upwind_cell_)[f];
    if (c >= 0) { total_outflux[c] += std::abs(flux[0][f]); }
  }

  ComputeSinks2TotalOutFlux_(tcc_prev, total_outflux, 0, num_aqueous_ - 1);

  // loop over cells and calculate minimal timestep
  double vol = 0.;
  double ws_min_dt = 0.;
  double outflux_min_dt = 0.;
  dt_ = TRANSPORT_LARGE_TIME_STEP;
  double dt_cell = TRANSPORT_LARGE_TIME_STEP;
  int cmin_dt = 0;

  S_->GetEvaluator(porosity_key_, Tags::NEXT).Update(*S_, name_);
  const Epetra_MultiVector& phi = *S_->Get<CompositeVector>(porosity_key_, Tags::NEXT).ViewComponent("cell", false);

  for (int c = 0; c < tcc_prev.MyLength(); c++) {
    double outflux = total_outflux[c];

    if ((outflux > 0) && ((*ws_prev_)[0][c] > 0) && ((*ws_)[0][c] > 0) && (phi[0][c] > 0)) {
      vol = mesh_->getCellVolume(c);
      dt_cell = vol * (*mol_dens_)[0][c] * phi[0][c] *
                std::min((*ws_prev_)[0][c], (*ws_)[0][c]) / outflux;
    }
    if (dt_cell < dt_) {
      dt_ = dt_cell;
      cmin_dt = c;
      ws_min_dt = std::min((*ws_prev_)[0][c], (*ws_)[0][c]);
      outflux_min_dt = total_outflux[c];
    }
  }

  if (spatial_disc_order_ == 2) dt_ /= 2;

  // communicate global timestep
  double dt_tmp = dt_;
  const Epetra_Comm& comm = ws_prev_->Comm();
  comm.MinAll(&dt_tmp, &dt_, 1);

  // incorporate developers and CFL constraints
  dt_ = std::min(dt_, dt_max_);
  dt_ *= cfl_;

  // print optional diagnostics using maximum cell id as the filter
  auto& cell_map = mesh_->getMap(AmanziMesh::Entity_kind::CELL, false);
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    int cmin_dt_unique = (std::abs(dt_tmp * cfl_ - dt_) < 1e-6 * dt_) ? cell_map.GID(cmin_dt) : -2;

    int cmin_dt_tmp = cmin_dt_unique;
    comm.MaxAll(&cmin_dt_tmp, &cmin_dt_unique, 1);
    int min_pid = -1;

    double tmp_package[6];
    if (cell_map.GID(cmin_dt) == cmin_dt_unique) {
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
    *vo_->os() << "Stable timestep " << dt_ << " is computed at (" << tmp_package[2] << ", "
               << tmp_package[3];
    if (std::abs(3 - tmp_package[5]) < 1e-10) *vo_->os() << ", " << tmp_package[4];
    *vo_->os() << ")" << std::endl;
    *vo_->os() << "Stable timestep " << dt_ << " is limited by saturation/ponded_depth "
               << tmp_package[0] << " and "
               << "output flux " << tmp_package[1] << std::endl;
  }
  return dt_;
}


/* *******************************************************************
* Estimate returns last timestep unless it is zero.
******************************************************************* */
double
Transport_ATS::get_dt()
{
  return dt_;
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
  S_->GetEvaluator(saturation_key_, Tags::NEXT).Update(*S_, name_);
  S_->GetEvaluator(saturation_key_, Tags::CURRENT).Update(*S_, name_);
  S_->GetEvaluator(molar_density_key_, Tags::NEXT).Update(*S_, name_);
  S_->GetEvaluator(molar_density_key_, Tags::CURRENT).Update(*S_, name_);

#ifdef ALQUIMIA_ENABLED
  if (plist_->sublist("source terms").isSublist("geochemical")) {
    for (auto& src : srcs_) {
      if (src->getType() == DomainFunction_kind::ALQUIMIA) {
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
      if (bc->getType() == DomainFunction_kind::ALQUIMIA) {
        Teuchos::RCP<TransportBoundaryFunction_Alquimia_Units> bc_alq =
          Teuchos::rcp_dynamic_cast<TransportBoundaryFunction_Alquimia_Units>(bc);
        bc_alq->set_conversion(1000.0, mol_dens_, true);
      }
    }
  }
#endif

  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");
  db_->WriteVector("tcc_old", tcc.ptr());

  // calculate stable timestep
  double dt_shift = 0.0, dt_global = dt_MPC;
  double time = t_old;
  if (time >= 0.0) {
    t_physics_ = time;
    dt_shift = time - S_->get_time(tag_current_);
    dt_global = S_->get_time(tag_next_) - S_->get_time(tag_current_);
    AMANZI_ASSERT(std::abs(dt_global - dt_MPC) < 1.e-4);
  }

  ComputeStableTimeStep_();
  double dt_stable = dt_; // advance routines override dt_
  double dt_sum = 0.0;
  double dt_cycle;
  dt_cycle = std::min(dt_stable, dt_MPC);

  Tag water_tag_current = Tags::CURRENT;
  Tag water_tag_next = Tags::NEXT;

  S_->GetEvaluator(porosity_key_, Tags::NEXT).Update(*S_, name_);
  const Epetra_MultiVector& phi = *S_->Get<CompositeVector>(porosity_key_, Tags::NEXT).ViewComponent("cell", false);

  db_->WriteVector("sat_old",
                   S_->GetPtr<CompositeVector>(saturation_key_, water_tag_current).ptr());
  db_->WriteVector("sat_new", S_->GetPtr<CompositeVector>(saturation_key_, water_tag_next).ptr());
  db_->WriteVector("mol_dens_old",
                   S_->GetPtr<CompositeVector>(molar_density_key_, water_tag_current).ptr());
  db_->WriteVector("mol_dens_new",
                   S_->GetPtr<CompositeVector>(molar_density_key_, water_tag_next).ptr());
  db_->WriteVector("poro", S_->GetPtr<CompositeVector>(porosity_key_, Tags::NEXT).ptr());

  for (int c = 0; c < tcc_prev.MyLength(); c++) {
    double vol_phi_ws_den;
    vol_phi_ws_den =
      mesh_->getCellVolume(c) * phi[0][c] * (*ws_prev_)[0][c] * (*mol_dens_prev_)[0][c];
    for (int i = 0; i < num_aqueous_ + num_gaseous_; i++) {
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

    if (spatial_disc_order_ == 1) { // temporary solution (lipnikov@lanl.gov)
      AdvanceDonorUpwind_(dt_cycle);
    } else if (spatial_disc_order_ == 2 && temporal_disc_order_ == 1) {
      AdvanceSecondOrderUpwindRK1_(dt_cycle);
    } else if (spatial_disc_order_ == 2 && temporal_disc_order_ == 2) {
      AdvanceSecondOrderUpwindRK2_(dt_cycle);
    } else {
      AMANZI_ASSERT(false);
    }

    if (!final_cycle) { // rotate concentrations (we need new memory for tcc)
      // should not be allocating here, we have tons of memory for tcc --ETC
      tcc = Teuchos::RCP<CompositeVector>(new CompositeVector(*tcc_tmp));
    }

    ncycles++;
  }

  dt_ = dt_stable; // restore the original timestep (just in case)

  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", false);
  Advance_Dispersion_Diffusion_(t_old, t_new);
  // optional Henry Law for the case of gas diffusion
  if (henry_law_) MakeAirWaterPartitioning_();

  // statistics output
  nsubcycles = ncycles;
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    *vo_->os() << ncycles << " sub-cycles, dt_stable=" << units_.OutputTime(dt_stable)
               << " [sec]  dt_MPC=" << units_.OutputTime(dt_MPC) << " [sec]" << std::endl;

    PrintSoluteExtrema(tcc_next, dt_MPC);
  }

  ComputeStableTimeStep_();

  ChangedSolutionPK(tag_next_);
  return failed;
}


void
Transport_ATS ::Advance_Dispersion_Diffusion_(double t_old, double t_new)
{
  double dt_MPC = t_new - t_old;
  // We define tracer as the species #0 as calculate some statistics.
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", false);
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");
  int num_components_ = tcc_prev.NumVectors();

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
    Teuchos::RCP<Operators::BCs> bc_dummy = Teuchos::rcp(
      new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
    auto& bc_model = bc_dummy->bc_model();
    auto& bc_value = bc_dummy->bc_value();
    PopulateBoundaryData_(bc_model, bc_value, -1);

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

    S_->GetEvaluator(porosity_key_, Tags::NEXT).Update(*S_, name_);
    const Epetra_MultiVector& phi = *S_->Get<CompositeVector>(porosity_key_, Tags::NEXT).ViewComponent("cell", false);

    // populate the dispersion operator (if any)
    if (flag_dispersion_) {
      const Epetra_MultiVector& flux = *S_->Get<CompositeVector>(flux_key_, Tags::NEXT).ViewComponent("face", true);
      CalculateDispersionTensor_(flux, phi, *ws_, *mol_dens_);
    }

    int phase;  // transport phase -- 0: liquid, 1: gas
    int num_itrs(0);
    bool flag_op1(true);
    double md_change, md_old(0.0), md_new;  // molecular diffusion
    double residual(0.0);

    // Disperse and diffuse aqueous components
    for (int i = 0; i < num_aqueous_; i++) {
      FindDiffusionValue_(component_names_[i], &md_new, &phase);
      md_change = md_new - md_old;
      md_old = md_new;

      if (md_change != 0.0) {
        CalculateDiffusionTensor_(md_change, phase, phi, *ws_, *mol_dens_);
        flag_op1 = true;
      }

      // set initial guess
      Epetra_MultiVector& sol_cell = *sol.ViewComponent("cell");
      for (int c = 0; c < sol_cell.MyLength(); c++) { sol_cell[0][c] = tcc_next[i][c]; }
      if (sol.HasComponent("face")) { sol.ViewComponent("face")->PutScalar(0.0); }

      if (flag_op1) {
        op->Init();
        Teuchos::RCP<std::vector<WhetStone::Tensor>> Dptr = Teuchos::rcpFromRef(D_);
        op1->Setup(Dptr, Teuchos::null, Teuchos::null);
        op1->UpdateMatrices(Teuchos::null, Teuchos::null);

        // add accumulation term
        Epetra_MultiVector& fac = *factor.ViewComponent("cell");
        for (int c = 0; c < fac.MyLength(); c++) {
          fac[0][c] = phi[0][c] * (*ws_)[0][c] * (*mol_dens_)[0][c];
        }
        op2->AddAccumulationDelta(sol, factor, factor, dt_MPC, "cell");
        op1->ApplyBCs(true, true, true);

      } else {
        Epetra_MultiVector& rhs_cell = *op->rhs()->ViewComponent("cell");
        for (int c = 0; c < rhs_cell.MyLength(); c++) {
          double tmp =
            mesh_->getCellVolume(c) * (*ws_)[0][c] * phi[0][c] * (*mol_dens_)[0][c] / dt_MPC;
          rhs_cell[0][c] = tcc_next[i][c] * tmp;
        }
      }

      CompositeVector& rhs = *op->rhs();
      int ierr = op->ApplyInverse(rhs, sol);

      if (ierr != 0) {
        Errors::Message msg("TransportExplicit_PK solver failed with message: \"");
        msg << op->returned_code_string() << "\"";
        Exceptions::amanzi_throw(msg);
      }

      residual += op->residual();
      num_itrs += op->num_itrs();

      for (int c = 0; c < sol_cell.MyLength(); c++) { tcc_next[i][c] = sol_cell[0][c]; }
      if (sol.HasComponent("face")) {
        if (tcc_tmp->HasComponent("boundary_face")) {
          Epetra_MultiVector& tcc_tmp_bf = *tcc_tmp->ViewComponent("boundary_face", false);
          Epetra_MultiVector& sol_faces = *sol.ViewComponent("face", false);
          const Epetra_Map& vandalay_map =
            mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false);
          const Epetra_Map& face_map = mesh_->getMap(AmanziMesh::Entity_kind::FACE, false);
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
    for (int i = num_aqueous_; i < num_components_; i++) {
      FindDiffusionValue_(component_names_[i], &md_new, &phase);
      md_change = md_new - md_old;
      md_old = md_new;

      if (md_change != 0.0 || i == num_aqueous_) {
        CalculateDiffusionTensor_(md_change, phase, phi, *ws_, *mol_dens_);
      }

      // set initial guess
      Epetra_MultiVector& sol_cell = *sol.ViewComponent("cell");
      for (int c = 0; c < sol_cell.MyLength(); c++) { sol_cell[0][c] = tcc_next[i][c]; }
      if (sol.HasComponent("face")) { sol.ViewComponent("face")->PutScalar(0.0); }

      op->Init();
      Teuchos::RCP<std::vector<WhetStone::Tensor>> Dptr = Teuchos::rcpFromRef(D_);
      op1->Setup(Dptr, Teuchos::null, Teuchos::null);
      op1->UpdateMatrices(Teuchos::null, Teuchos::null);

      // add boundary conditions and sources for gaseous components
      PopulateBoundaryData_(bc_model, bc_value, i);

      Epetra_MultiVector& rhs_cell = *op->rhs()->ViewComponent("cell");
      ComputeAddSourceTerms_(t_new, 1.0, rhs_cell, i, i);
      op1->ApplyBCs(true, true, true);

      // add accumulation term
      Epetra_MultiVector& fac1 = *factor.ViewComponent("cell");
      Epetra_MultiVector& fac0 = *factor0.ViewComponent("cell");

      for (int c = 0; c < fac0.MyLength(); c++) {
        fac1[0][c] = phi[0][c] * (1.0 - (*ws_)[0][c]) * (*mol_dens_)[0][c];
        fac0[0][c] = phi[0][c] * (1.0 - (*ws_prev_)[0][c]) * (*mol_dens_prev_)[0][c];
        if ((*ws_)[0][c] == 1.0) fac1[0][c] = 1.0 * (*mol_dens_)[0][c]; // hack so far
      }
      op2->AddAccumulationDelta(sol, factor0, factor, dt_MPC, "cell");

      CompositeVector& rhs = *op->rhs();
      int ierr = op->ApplyInverse(rhs, sol);
      if (ierr != 0) {
        Errors::Message msg("Transport_PK solver failed with message: \"");
        msg << op->returned_code_string() << "\"";
        Exceptions::amanzi_throw(msg);
      }

      residual += op->residual();
      num_itrs += op->num_itrs();

      for (int c = 0; c < tcc_next.MyLength(); c++) { tcc_next[i][c] = sol_cell[0][c]; }
    }

    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "dispersion solver ||r||=" << residual / num_components_
                 << " itrs=" << num_itrs / num_components_ << std::endl;
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
Transport_ATS::AdvanceDonorUpwind_(double dt_cycle)
{
  IdentifyUpwindCells_();
  dt_ = dt_cycle; // overwrite the maximum stable transport step
  mass_solutes_source_.assign(num_aqueous_ + num_gaseous_, 0.0);
  mass_solutes_bc_.assign(num_aqueous_ + num_gaseous_, 0.0);

  // scattering total concentration from master to others
  tcc->ScatterMasterToGhosted("cell");
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);       // tag current
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);   // tag next

  S_->GetEvaluator(porosity_key_, Tags::NEXT).Update(*S_, name_);
  const Epetra_MultiVector& phi = *S_->Get<CompositeVector>(porosity_key_, Tags::NEXT).ViewComponent("cell", false);

  double mass_current = 0.;
  double tmp1, mass;
  int num_components_ = tcc_next.NumVectors();

  // populating conserved quantity (unit: molC)
  // The conserve_qty has `num_components_+2` vectors
  Epetra_MultiVector& conserve_qty =
    *S_->GetW<CompositeVector>(conserve_qty_key_, tag_next_, name_).ViewComponent("cell", false);
  conserve_qty.PutScalar(0.);     // set all values in conserve_qty to 0 -- why?

  // populating solid quantity (unit: molC)
  // solid_qty has `num_components_` vectors only (solute mass)
  Epetra_MultiVector& solid_qty =
    *S_->GetW<CompositeVector>(solid_residue_mass_key_, tag_next_, name_).ViewComponent("cell", false);

  for (int c = 0; c < conserve_qty.MyLength(); c++) {
    double vol_phi_ws_den =
      mesh_->getCellVolume(c) * phi[0][c] * (*ws_prev_)[0][c] * (*mol_dens_prev_)[0][c];
    (conserve_qty)[num_components_ + 1][c] = vol_phi_ws_den;

    for (int i = 0; i < num_advect_; i++) {
      (conserve_qty)[i][c] = tcc_prev[i][c] * vol_phi_ws_den;       // get and store current solute mass

      if (dissolution_) {
        if (((*ws_prev_)[0][c] > water_tolerance_) &&
            (solid_qty[i][c] > 0)) { // Dissolve solid residual into liquid
          double add_mass =
            std::min(solid_qty[i][c], tcc_max_[i] * vol_phi_ws_den - conserve_qty[i][c]);
          solid_qty[i][c] -= add_mass;
          conserve_qty[i][c] += add_mass;
        }
      }

      mass_current += conserve_qty[i][c]; // this includes total current mass of all solutes
    }
  }

  db_->WriteCellVector("cons (start)", conserve_qty);
  tmp1 = mass_current;
  mesh_->getComm()->SumAll(&tmp1, &mass_current, 1);    // sum all mass_current from all processors
  int nfaces_all =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  const Epetra_MultiVector& flux = *S_->Get<CompositeVector>(flux_key_, Tags::NEXT).ViewComponent("face", true);

  // advance all components at once
  for (int f = 0; f < nfaces_all; f++) {
    // flow moves from upwind cell (c1) to downwind cell (c2).
    // If ci < 0 || ci > ncells_owned -> indicates boundary or halo cells (i=1,2)
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];
    double u = std::abs(flux[0][f]);

    if (c1 >= 0 && c1 < conserve_qty.MyLength() && c2 >= 0 && c2 < conserve_qty.MyLength()) {
      // Here c1 & c2 are inside local domain. Update solute fluxes for both cells
      for (int i = 0; i < num_advect_; i++) {
        double tcc_flux = dt_ * u * tcc_prev[i][c1];
        conserve_qty[i][c1] -= tcc_flux;
        conserve_qty[i][c2] += tcc_flux;
      }
      // Update (tag next) water fluxes for both cells
      conserve_qty[num_components_ + 1][c1] -= dt_ * u;
      conserve_qty[num_components_ + 1][c2] += dt_ * u;

    } else if (c1 >= 0 && c1 < conserve_qty.MyLength() && (c2 >= conserve_qty.MyLength() || c2 < 0)) {
      // downind cell c2 is boundary or belong to another domain owned by other processors
      // Update solute flux for c1
      for (int i = 0; i < num_advect_; i++) {
        double tcc_flux = dt_ * u * tcc_prev[i][c1];
        conserve_qty[i][c1] -= tcc_flux;
        if (c2 < 0) mass_solutes_bc_[i] -= tcc_flux;  // if c2 is boundary, update BC
      }
      // Update or subtract (tag next) water fluxes for c1
      conserve_qty[num_components_ + 1][c1] -= dt_ * u;

    } else if (c1 >= conserve_qty.MyLength() && c2 >= 0 && c2 < conserve_qty.MyLength()) {
      // upwind cell c1 is boundary or belong to another domain owned by other processors
      // Update solute flux for c2
      for (int i = 0; i < num_advect_; i++) {
        double tcc_flux = dt_ * u * tcc_prev[i][c1];
        conserve_qty[i][c2] += tcc_flux;
      }
      // Update or add (tag next) water fluxes for c2
      conserve_qty[num_components_ + 1][c2] += dt_ * u;

    } else if (c2 < 0 && c1 >= 0 && c1 < conserve_qty.MyLength()) {
      // this case is very similar to line 1165, except c2<0 only. Why we don't update solute flux??? --PL
      conserve_qty[num_components_ + 1][c1] -= dt_ * u;

    } else if (c1 < 0 && c2 >= 0 && c2 < conserve_qty.MyLength()) {
      // why no solute update??? --PL
      conserve_qty[num_components_ + 1][c2] += dt_ * u;
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

      double u = std::abs(flux[0][f]);
      if (c2 >= 0) {
        for (int i = 0; i < ncomp; i++) {
          int k = tcc_index[i];
          if (k < num_advect_) {
            double tcc_flux = dt_ * u * values[i];
            conserve_qty[k][c2] += tcc_flux;
            mass_solutes_bc_[k] += tcc_flux;

            if (tcc_tmp_bf) (*tcc_tmp_bf)[i][bf] = values[i];
          }
        }
      }
    }
  }
  db_->WriteCellVector("cons (adv)", conserve_qty);

  // process external sources
  if (srcs_.size() != 0) {
    double time = t_physics_;
    ComputeAddSourceTerms_(time, dt_, conserve_qty, 0, num_aqueous_ - 1);
  }
  db_->WriteCellVector("cons (src)", conserve_qty);

  // recover concentration from new conservative state
  for (int c = 0; c < conserve_qty.MyLength(); c++) {
    double water_new =
      mesh_->getCellVolume(c) * phi[0][c] * (*ws_)[0][c] * (*mol_dens_)[0][c];  // next water in cell c
    double water_sink = conserve_qty[num_components_][c];  // current water in cell c (seem always 0 because conserve_qty.PutScalar(0.))
    double water_total = water_new + water_sink;           // unit: mol H20
    AMANZI_ASSERT(water_total >= water_new);
    conserve_qty[num_components_][c] = water_total;        // update current water in cell c

    for (int i = 0; i < num_advect_; i++) {
      if (water_new > water_tolerance_ && conserve_qty[i][c] > 0) {
        // there is both water and stuff present at the new time
        // this is stuff at the new time + stuff leaving through the domain coupling, divided by water of both
        tcc_next[i][c] = conserve_qty[i][c] / water_total;

        // limit tcc_next to tcc_max_, precipitating the remainder as solid residue
        if (tcc_next[i][c] > tcc_max_[i]) {
          solid_qty[i][c] += (tcc_next[i][c] - tcc_max_[i]) * water_total;
          tcc_next[i][c] = tcc_max_[i];
        }
      } else if (water_sink > water_tolerance_ && conserve_qty[i][c] > 0) {
        // there is water and stuff leaving through the domain coupling, but it all leaves (none at the new time)
        tcc_next[i][c] = 0.;
      } else {
        // there is no water leaving, and no water at the new time.  Change any stuff into solid
        solid_qty[i][c] += std::max(conserve_qty[i][c], 0.);
        conserve_qty[i][c] = 0.;
        tcc_next[i][c] = 0.;
      }
    }
  }
  db_->WriteCellVector("tcc_new", tcc_next);
  PrintSoluteExtrema(tcc_next, dt_);

  double mass_final = 0;
  for (int c = 0; c < conserve_qty.MyLength(); c++) {
    for (int i = 0; i < num_advect_; i++) { mass_final += conserve_qty[i][c]; }
  }

  tmp1 = mass_final;
  mesh_->getComm()->SumAll(&tmp1, &mass_final, 1);

  // update mass balance
  for (int i = 0; i < mass_solutes_exact_.size(); i++) {
    mass_solutes_exact_[i] += mass_solutes_source_[i] * dt_;
  }

  if (internal_tests) { CheckGEDProperty(*tcc_tmp->ViewComponent("cell"), t_physics_); }
}


/* *******************************************************************
 * We have to advance each component independently due to different
 * reconstructions. We use tcc when only owned data are needed and
 * tcc_next when owned and ghost data. This is a special routine for
 * transient flow and uses first-order time integrator.
 ****************************************************************** */
void
Transport_ATS::AdvanceSecondOrderUpwindRK1_(double dt_cycle)
{
  dt_ = dt_cycle; // overwrite the maximum stable transport step
  mass_solutes_source_.assign(num_aqueous_ + num_gaseous_, 0.0);

  // work memory
  const Epetra_Map& cmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::CELL, true);

  // distribute vector of concentrations
  tcc->ScatterMasterToGhosted("cell");
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);

  // We advect only aqueous components.
  int num_components_ = tcc_next.NumVectors();

  // populating conserved quantity
  Epetra_MultiVector& conserve_qty =
    *S_->GetW<CompositeVector>(conserve_qty_key_, tag_next_, name_).ViewComponent("cell", false);
  conserve_qty.PutScalar(0.);

  // populating solid quantity
  Epetra_MultiVector& solid_qty =
    *S_->GetW<CompositeVector>(solid_residue_mass_key_, tag_next_, name_).ViewComponent("cell", false);

  // using vectors
  S_->GetEvaluator(porosity_key_, Tags::NEXT).Update(*S_, name_);
  const Epetra_MultiVector& phi = *S_->Get<CompositeVector>(porosity_key_, Tags::NEXT).ViewComponent("cell", false);

  // prepopulate with initial water for better debugging
  for (int c = 0; c < conserve_qty.MyLength(); c++) {
    conserve_qty[num_components_ + 1][c] = mesh_->getCellVolume(c) * phi[0][c] * (*ws_prev_)[0][c] * (*mol_dens_prev_)[0][c];
  }

  for (int i = 0; i < num_advect_; i++) {
    current_component_ = i; // needed by BJ
    double T = t_physics_;
    Epetra_Vector*& component = tcc_prev(i);
    FunctionalTimeDerivative(T, *component, *(conserve_qty)(i));
  }
  db_->WriteCellVector("cons (time_deriv)", conserve_qty);

  // calculate the new conc
  for (int c = 0; c < conserve_qty.MyLength(); c++) {
    double water_old = conserve_qty[num_components_ + 1][c];
    double water_new =
      mesh_->getCellVolume(c) * phi[0][c] * (*ws_)[0][c] * (*mol_dens_)[0][c];
    double water_sink = conserve_qty[num_components_][c];
    double water_total = water_sink + water_new;
    conserve_qty[num_components_][c] = water_total;

    for (int i = 0; i != num_components_; ++i) {
      double cons_qty = (tcc_prev[i][c] + dt_ * conserve_qty[i][c]) * water_old;
      conserve_qty[i][c] = cons_qty;
      if (water_new > water_tolerance_ && cons_qty > 0) {
        // there is both water and stuff present at the new time
        // this is stuff at the new time + stuff leaving through the domain coupling, divided by water of both
        tcc_next[i][c] = cons_qty / water_total;
      } else if (water_sink > water_tolerance_ && cons_qty > 0) {
        // there is water and stuff leaving through the domain coupling, but it all leaves (none at the new time)
        tcc_next[i][c] = 0.;
      } else {
        // there is no water leaving, and no water at the new time.  Change any stuff into solid
        solid_qty[i][c] += std::max(cons_qty, 0.);
        conserve_qty[i][c] = 0.;
        tcc_next[i][c] = 0.;
      }
    }
  }
  db_->WriteCellVector("tcc_new", tcc_next);

  // update mass balance
  for (int i = 0; i < num_aqueous_ + num_gaseous_; i++) {
    mass_solutes_exact_[i] += mass_solutes_source_[i] * dt_;
  }

  if (internal_tests) { CheckGEDProperty(*tcc_tmp->ViewComponent("cell"), t_physics_); }
}


/* *******************************************************************
 * We have to advance each component independently due to different
 * reconstructions. This is a special routine for transient flow and
 * uses second-order predictor-corrector time integrator.
 ****************************************************************** */
void
Transport_ATS::AdvanceSecondOrderUpwindRK2_(double dt_cycle)
{
  // Q: Why we don't have conserve_qty calculation in this function? --PL

  dt_ = dt_cycle; // overwrite the maximum stable transport step
  mass_solutes_source_.assign(num_aqueous_ + num_gaseous_, 0.0);

  // work memory
  const Epetra_Map& cmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::CELL, true);
  Epetra_Vector f_component(cmap_wghost);

  // populating solid quantity
  Epetra_MultiVector& solid_qty =
    *S_->GetW<CompositeVector>(solid_residue_mass_key_, tag_next_, name_).ViewComponent("cell", false);

  // distribute old vector of concentrations
  tcc->ScatterMasterToGhosted("cell");
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);
  Epetra_Vector ws_ratio(Copy, *ws_prev_, 0);

  S_->GetEvaluator(porosity_key_, Tags::NEXT).Update(*S_, name_);
  const Epetra_MultiVector& phi = *S_->Get<CompositeVector>(porosity_key_, Tags::NEXT).ViewComponent("cell", false);

  for (int c = 0; c < solid_qty.MyLength(); c++) {
    if ((*ws_)[0][c] > 1e-10) {
      if ((*ws_prev_)[0][c] > 1e-10) {
        ws_ratio[c] = ((*ws_prev_)[0][c] * (*mol_dens_prev_)[0][c]) /
                      ((*ws_)[0][c] * (*mol_dens_)[0][c]);
      } else {
        ws_ratio[c] = 1;
      }
    } else
      ws_ratio[c] = 0.;
  }

  // predictor step
  for (int i = 0; i < num_advect_; i++) {
    current_component_ = i; // needed by BJ

    double T = t_physics_;
    Epetra_Vector*& component = tcc_prev(i);
    FunctionalTimeDerivative(T, *component, f_component);

    for (int c = 0; c < tcc_prev.MyLength(); c++) {
      tcc_next[i][c] = (tcc_prev[i][c] + dt_ * f_component[c]) * ws_ratio[c];
    }
  }

  tcc_tmp->ScatterMasterToGhosted("cell");

  // corrector step
  for (int i = 0; i < num_advect_; i++) {
    current_component_ = i; // needed by BJ

    double T = t_physics_;
    Epetra_Vector*& component = tcc_next(i);
    FunctionalTimeDerivative(T, *component, f_component);

    for (int c = 0; c < tcc_prev.MyLength(); c++) {
      double value = (tcc_prev[i][c] + dt_ * f_component[c]) * ws_ratio[c];
      tcc_next[i][c] = (tcc_next[i][c] + value) / 2;
      if (tcc_next[i][c] < 0) {
        double vol_phi_ws_den =
          mesh_->getCellVolume(c) * phi[0][c] * (*ws_)[0][c] * (*mol_dens_)[0][c];
        solid_qty[i][c] += std::abs(tcc_next[i][c]) * vol_phi_ws_den;
        tcc_next[i][c] = 0.;
      }
    }
  }

  // update mass balance
  for (int i = 0; i < num_aqueous_ + num_gaseous_; i++) {
    mass_solutes_exact_[i] += mass_solutes_source_[i] * dt_ / 2;
  }

  if (internal_tests) { CheckGEDProperty(*tcc_tmp->ViewComponent("cell"), t_physics_); }
}


/* ******************************************************************
* Computes source and sink terms and adds them to vector tcc.
* Returns mass rate for the tracer.
* The routine treats two cases of tcc with one and all components.
****************************************************************** */
void
Transport_ATS::ComputeAddSourceTerms_(double tp,
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

      if (c >= cons_qty.MyLength()) continue;

      if (srcs_[m]->getType() == DomainFunction_kind::COUPLING && n0 == 0) {
        cons_qty[num_vectors - 2][c] += values[num_vectors - 2];
      }

      if (convert_to_field_[srcs_[m]->getName()]) {
          std::string name = srcs_[m]->getName();
          name = Keys::cleanName(name);
          if (Keys::getDomain(name)!=domain_){
            name = Keys::getKey(domain_, name);
          }
          copyToCompositeVector(*srcs_[m],
          S_->GetW<CompositeVector>(name, Tags::NEXT, name)
          );
          changedEvaluatorPrimary(name, Tags::NEXT, *S_);
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

/* ******************************************************************
* Update source/sink terms to outflux
*   - tcc_c: Total concentration in previous step (known)
*   - total_outflux: water flux
*   - n0: lower number of components
*   - n1: upper number of components
****************************************************************** */
void
Transport_ATS::ComputeSinks2TotalOutFlux_(Epetra_MultiVector& tcc_c,
                                  std::vector<double>& total_outflux,
                                  int n0,
                                  int n1)
{
  int ncells_all =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  std::vector<double> sink_add(ncells_all, 0.0);
  //Assumption that there is only one sink per component per cell
  double t0 = S_->get_time(tag_current_);
  int num_vectors = tcc_c.NumVectors();                 // number of components (e.g. tracers)
  int nsrcs = srcs_.size();                             // number of sources / sinks

  for (int m = 0; m < nsrcs; m++) {
    srcs_[m]->Compute(t0, t0);                          // compute source term at time t0
    std::vector<int> index = srcs_[m]->tcc_index();     // get the index of the component in the global list

    for (auto it = srcs_[m]->begin(); it != srcs_[m]->end(); ++it) {
      int c = it->first;                                // cell id
      std::vector<double>& values = it->second;         // magnitude of source terms
      double val = 0;
      for (int k = 0; k < index.size(); ++k) {
        int i = index[k];
        if (i < n0 || i > n1) continue;

        int imap = i;
        if (num_vectors == 1) imap = 0;

        if ((values[k] < 0) && (tcc_c[imap][c] > 1e-16)) {
          if (srcs_[m]->getType() == DomainFunction_kind::COUPLING) {
            const Epetra_MultiVector& flux_interface_ =
              *S_->Get<CompositeVector>("surface-surface_subsurface_flux", Tags::NEXT).ViewComponent("cell", false);
            val = std::max(val, std::abs(flux_interface_[0][c]));
          }
        }
      }
      sink_add[c] = std::max(sink_add[c], val);
    }
  }

  for (int c = 0; c < ncells_all; c++) total_outflux[c] += sink_add[c];
}


/* *******************************************************************
* Populates operators' boundary data for given component.
* Returns true if at least one face was populated.
******************************************************************* */
bool
Transport_ATS::PopulateBoundaryData_(std::vector<int>& bc_model,
                                    std::vector<double>& bc_value,
                                    int component)
{
  int nfaces_all =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
  bool flag = false;

  for (int i = 0; i < bc_model.size(); i++) {
    bc_model[i] = Operators::OPERATOR_BC_NONE;
    bc_value[i] = 0.0;
  }

  // This is not efficient. We should do this only once, not every timestep.
  for (int f = 0; f < nfaces_all; f++) {
    auto cells = mesh_->getFaceCells(f);
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
Transport_ATS::IdentifyUpwindCells_()
{
  int ncells_all =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  upwind_cell_->PutValue(-1);
  downwind_cell_->PutValue(-1);

  S_->Get<CompositeVector>(flux_key_, Tags::NEXT).ScatterMasterToGhosted("face");
  const Epetra_MultiVector& flux = *S_->Get<CompositeVector>(flux_key_, Tags::NEXT).ViewComponent("face", true);

  // identify upwind and downwind cell of each face
  for (int c = 0; c < ncells_all; c++) {
    const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);

    for (int i = 0; i < faces.size(); i++) {
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


void
Transport_ATS::ChangedSolutionPK(const Tag& tag)
{
  changedEvaluatorPrimary(key_, tag, *S_);
}

/* ****************************************************************
* Find place of the given component in a multivector.
**************************************************************** */
int
Transport_ATS::FindComponentNumber_(const std::string& component_name)
{
  int ncomponents = component_names_.size();
  for (int i = 0; i < ncomponents; i++) {
    if (component_names_[i] == component_name) return i;
  }
  Errors::Message msg("TransportExplicit_PK: component \"");
  msg << component_name << "\" was requested, but this is not a known component for this PK.";
  Exceptions::amanzi_throw(msg);
  return -1;
}


} // namespace Transport
} // namespace Amanzi
