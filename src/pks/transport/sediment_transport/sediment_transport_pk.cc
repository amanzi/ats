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
#include "GMVMesh.hh"
#include "Mesh.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_Diffusion.hh"
#include "PDE_Accumulation.hh"
#include "PK_DomainFunctionFactory.hh"
#include "PK_Utils.hh"

#include "sediment_transport_pk.hh"
#include "TransportDomainFunction.hh"


namespace Amanzi {
namespace SedimentTransport {

/* ******************************************************************
* New constructor compatible with new MPC framework.
****************************************************************** */
SedimentTransport_PK::SedimentTransport_PK(Teuchos::ParameterList& pk_tree,
                                           const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                           const Teuchos::RCP<State>& S,
                                           const Teuchos::RCP<TreeVector>& soln)
  : S_(S), soln_(soln)
{
  name_ = Keys::cleanPListName(pk_tree.name());

  // Create miscaleneous lists.
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);

  tp_list_ = Teuchos::sublist(pk_list, name_, true);

  if (tp_list_->isParameter("component names")) {
    component_names_ = tp_list_->get<Teuchos::Array<std::string>>("component names").toVector();
    mol_masses_ = tp_list_->get<Teuchos::Array<double>>("component molar masses").toVector();
  } else if (glist->isSublist("Cycle Driver")) {
    if (glist->sublist("Cycle Driver").isParameter("component names")) {
      // grab the component names
      component_names_ = glist->sublist("Cycle Driver")
                           .get<Teuchos::Array<std::string>>("component names")
                           .toVector();
    } else {
      Errors::Message msg("Transport PK: parameter component names is missing.");
      Exceptions::amanzi_throw(msg);
    }
  } else {
    Errors::Message msg(
      "Transport PK: sublist Cycle Driver or parameter component names is missing.");
    Exceptions::amanzi_throw(msg);
  }

  subcycling_ = tp_list_->get<bool>("transport subcycling", false);

  // initialize io
  Teuchos::RCP<Teuchos::ParameterList> units_list = Teuchos::sublist(glist, "units");
  units_.Init(*units_list);

  vo_ = Teuchos::null;
}


/* ******************************************************************
* Define structure of this PK.
****************************************************************** */
void
SedimentTransport_PK::Setup(const Teuchos::Ptr<State>& S)
{
  passwd_ = "state"; // owner's password

  // are we subcycling internally?
  subcycling_ = plist_->get<bool>("transport subcycling", false);
  if (subcycling_) {
    tag_subcycle_current_ = Tag{ name() + "_subcycling_current" };
    tag_subcycle_next_ = Tag{ name() + "_subcycling_next" };
  } else {
    tag_subcycle_current_ = tag_current_;
    tag_subcycle_next_ = tag_next_;
  }

  domain_name_ = tp_list_->get<std::string>("domain name", "domain");

  saturation_key_ = Keys::readKey(*tp_list_, domain_name_, "saturation liquid");
  prev_saturation_key_ = Keys::readKey(*tp_list_, domain_name_, "previous saturation liquid");
  flux_key_ = Keys::readKey(*tp_list_, domain_name_, "water flux", "water_flux");
  tcc_key_ = Keys::readKey(*tp_list_, domain_name_, "concentration", "sediment");
  molar_density_key_ =
    Keys::readKey(*tp_list_, domain_name_, "molar density", "molar_density_liquid");
  solid_residue_mass_key_ =
    Keys::readKey(*tp_list_, domain_name_, "solid residue", "solid_residue_mass");
  sd_trapping_key_ = Keys::readKey(*tp_list_, domain_name_, "trapping rate", "trapping_rate");
  sd_settling_key_ = Keys::readKey(*tp_list_, domain_name_, "settling rate", "settling_rate");
  sd_erosion_key_ = Keys::readKey(*tp_list_, domain_name_, "erosion rate", "erosion_rate");
  sd_organic_key_ = Keys::readKey(*tp_list_, domain_name_, "organic rate", "organic_rate");
  horiz_mixing_key_ =
    Keys::readKey(*tp_list_, domain_name_, "horizontal mixing", "horizontal_mixing");
  elevation_increase_key_ = Keys::getKey(domain_name_, "deformation");
  porosity_key_ = Keys::getKey(domain_name_, "soil_porosity");

  water_tolerance_ = tp_list_->get<double>("water tolerance", 1e-1);
  max_tcc_ = tp_list_->get<double>("maximal concentration", 0.9);
  sediment_density_ = tp_list_->get<double>("sediment density [kg m^-3]");


  mesh_ = S->GetMesh(domain_name_);
  dim = mesh_->getSpaceDimension();

  // cross-coupling of PKs
  Teuchos::RCP<Teuchos::ParameterList> physical_models =
    Teuchos::sublist(tp_list_, "physical models and assumptions");
  bool abs_perm = physical_models->get<bool>("permeability field is required", false);


  if (!S->HasRecordSet(flux_key_)) {
    S->Require<CompositeVector, CompositeVectorSpace>(flux_key_, tag_next_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    S->RequireEvaluator(flux_key_, tag_next_);
  }

  S->Require<CompositeVector, CompositeVectorSpace>(saturation_key_, tag_next_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S->RequireEvaluator(saturation_key_, tag_next_);

  // prev_sat does not have an evaluator, this is managed by hand.  not sure why
  S->Require<CompositeVector, CompositeVectorSpace>(saturation_key_, tag_current_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S->GetRecordW(saturation_key_, tag_current_, passwd_).set_io_vis(false);

  if (!S->HasRecordSet(sd_organic_key_)) {
    S->Require<CompositeVector, CompositeVectorSpace>(sd_organic_key_, tag_next_, sd_organic_key_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireEvaluator(sd_organic_key_, tag_next_);
  }

  if (!S->HasRecordSet(sd_trapping_key_)) {
    S->Require<CompositeVector, CompositeVectorSpace>(sd_trapping_key_, tag_next_, sd_trapping_key_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireEvaluator(sd_trapping_key_, tag_next_);
  }

  if (!S->HasRecordSet(sd_settling_key_)) {
    S->Require<CompositeVector, CompositeVectorSpace>(sd_settling_key_, tag_next_, sd_settling_key_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireEvaluator(sd_settling_key_, tag_next_);
  }

  if (!S->HasRecordSet(sd_erosion_key_)) {
    S->Require<CompositeVector, CompositeVectorSpace>(sd_erosion_key_, tag_next_, sd_erosion_key_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireEvaluator(sd_erosion_key_, tag_next_);
  }

  if (!S->HasRecordSet(horiz_mixing_key_)) {
    S->Require<CompositeVector, CompositeVectorSpace>(
       horiz_mixing_key_, tag_next_, horiz_mixing_key_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireEvaluator(horiz_mixing_key_, tag_next_);
  }

  if (!S->HasRecordSet(porosity_key_)) {
    S->Require<CompositeVector, CompositeVectorSpace>(porosity_key_, tag_next_, porosity_key_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireEvaluator(porosity_key_, tag_next_);
  }


  int ncomponents = component_names_.size();
  std::vector<std::vector<std::string>> subfield_names(1);
  subfield_names[0] = component_names_;


  if (!S->HasRecordSet(solid_residue_mass_key_)) {
    S->Require<CompositeVector, CompositeVectorSpace>(solid_residue_mass_key_, tag_next_, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, ncomponents);
  }

  if (!S->HasRecordSet(molar_density_key_)) {
    S->Require<CompositeVector, CompositeVectorSpace>(
       molar_density_key_, tag_next_, molar_density_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireEvaluator(molar_density_key_, tag_next_);
  }

  // require state fields when Transport PK is on
  if (component_names_.size() == 0) {
    Errors::Message msg;
    msg << "Transport PK: list of solutes is empty.\n";
    Exceptions::amanzi_throw(msg);
  }


  S->Require<CompositeVector, CompositeVectorSpace>(tcc_key_, tag_next_, passwd_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, ncomponents);
}


/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.
****************************************************************** */
void
SedimentTransport_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  // Set initial values for transport variables.
  dt_ = dt_debug_ = t_physics_ = 0.0;
  double time = S->get_time();
  if (time >= 0.0) t_physics_ = time;

  // if (tp_list_->isSublist("initial conditions")) {
  //   S->GetW<CompositeVector>(tcc_key_, tag_current_, passwd_).Initialize(tp_list_->sublist("initial conditions"));
  // }

  diffusion_preconditioner = "identity";

  internal_tests = 0;
  //tests_tolerance = TRANSPORT_CONCENTRATION_OVERSHOOT;

  bc_scaling = 0.0;

  // Create verbosity object.
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = tp_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("TransportPK", vlist));

  Teuchos::OSTab tab = vo_->getOSTab();
  MyPID = mesh_->getComm()->MyPID();

  //   // initialize missed fields
  InitializeFields_(S);

  //create copies
  //S->RequireFieldCopy(tcc_key_, "subcycling", passwd_);
  tcc_tmp = S_->GetPtrW<CompositeVector>(tcc_key_, tag_subcycle_next_, name_);
  tcc = S_->GetPtrW<CompositeVector>(tcc_key_, tag_subcycle_current_, name_);
  *tcc_tmp = *tcc;

  ws_subcycle_start =
    S_->GetW<CompositeVector>(saturation_key_, tag_subcycle_current_, name_).ViewComponent("cell");

  ws_subcycle_end =
    S_->GetW<CompositeVector>(saturation_key_, tag_subcycle_next_, name_).ViewComponent("cell");

  // Check input parameters. Due to limited amount of checks, we can do it earlier.
  // Policy(S.ptr());

  ncells_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  ncells_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  nfaces_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  nfaces_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
  nnodes_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);

  // extract control parameters
  InitializeAll_();

  // state pre-prosessing
  Teuchos::RCP<const CompositeVector> cv;

  ws_ = S->Get<CompositeVector>(saturation_key_, tag_next_).ViewComponent("cell", false);
  ws_prev_ =
    S->GetPtr<CompositeVector>(saturation_key_, tag_current_)->ViewComponent("cell", false);

  mol_dens_ =
    S->GetPtr<CompositeVector>(molar_density_key_, tag_next_)->ViewComponent("cell", false);
  //  mol_dens_prev_ = S_->GetPtr<CompositeVector>(molar_density_key_) -> ViewComponent("cell", false);
  km_ = S->GetPtr<CompositeVector>(horiz_mixing_key_, tag_next_)->ViewComponent("cell", false);

  tcc = S->GetPtrW<CompositeVector>(tcc_key_, tag_next_, passwd_);

  flux_ = S->Get<CompositeVector>(flux_key_, tag_next_).ViewComponent("face", true);
  solid_qty_ = S->GetW<CompositeVector>(solid_residue_mass_key_, tag_next_, passwd_)
                 .ViewComponent("cell", false);

  //create vector of conserved quatities
  conserve_qty_ = Teuchos::rcp(new Epetra_MultiVector(
    *(S->Get<CompositeVector>(tcc_key_, tag_next_).ViewComponent("cell", true))));

  // memory for new components
  // tcc_tmp = Teuchos::rcp(new CompositeVector(*(S->GetPtr<CompositeVector>(tcc_key_))));
  // *tcc_tmp = *tcc;

  // upwind
  const Epetra_Map& fmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);
  upwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));
  downwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));

  IdentifyUpwindCells();

  // advection block initialization
  current_component_ = -1;

  const Epetra_Map& cmap_owned = mesh_->getMap(AmanziMesh::Entity_kind::CELL,false);

  // reconstruction initialization
  const Epetra_Map& cmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::CELL,true);
  lifting_ = Teuchos::rcp(new Operators::ReconstructionCellLinear(mesh_));

  // create boundary conditions
  if (tp_list_->isSublist("boundary conditions")) {
    // -- try tracer-type conditions
    PK_DomainFunctionFactory<TransportDomainFunction> factory(mesh_, S_);
    Teuchos::ParameterList& clist =
      tp_list_->sublist("boundary conditions").sublist("concentration");

    for (Teuchos::ParameterList::ConstIterator it = clist.begin(); it != clist.end(); ++it) {
      std::string name = it->first;
      if (clist.isSublist(name)) {
        Teuchos::ParameterList& bc_list = clist.sublist(name);
        if (name == "coupling") {
          Teuchos::ParameterList::ConstIterator it1 = bc_list.begin();
          std::string specname = it1->first;
          Teuchos::ParameterList& spec = bc_list.sublist(specname);
          Teuchos::RCP<TransportDomainFunction> bc =
            factory.Create(spec, "boundary concentration", AmanziMesh::Entity_kind::FACE, Teuchos::null);

          for (int i = 0; i < component_names_.size(); i++) {
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
          int gid = std::stoi(domain_name_.substr(last_of + 1, domain_name_.size()));
          spec.set("entity_gid_out", gid);
          Teuchos::RCP<TransportDomainFunction> bc =
            factory.Create(spec, "boundary concentration", AmanziMesh::Entity_kind::FACE, Teuchos::null);

          for (int i = 0; i < component_names_.size(); i++) {
            bc->tcc_names().push_back(component_names_[i]);
            bc->tcc_index().push_back(i);
          }


          bc->set_state(S_);
          bcs_.push_back(bc);

        } else {
          for (Teuchos::ParameterList::ConstIterator it1 = bc_list.begin(); it1 != bc_list.end();
               ++it1) {
            std::string specname = it1->first;
            Teuchos::ParameterList& spec = bc_list.sublist(specname);
            Teuchos::RCP<TransportDomainFunction> bc =
              factory.Create(spec, "boundary concentration", AmanziMesh::Entity_kind::FACE, Teuchos::null);

            std::vector<int>& tcc_index = bc->tcc_index();
            std::vector<std::string>& tcc_names = bc->tcc_names();
            bc->set_state(S_);

            tcc_names.push_back(name);
            tcc_index.push_back(0);

            bcs_.push_back(bc);
          }
        }
      }
    }
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

  // VV_CheckInfluxBC();

  // source term initialization: so far only "concentration" is available.
  if (tp_list_->isSublist("source terms")) {
    PK_DomainFunctionFactory<TransportDomainFunction> factory(mesh_, S_);
    //if (domain_name_ == "domain")  PKUtils_CalculatePermeabilityFactorInWell(S_.ptr(), Kxy);

    Teuchos::ParameterList& clist = tp_list_->sublist("source terms").sublist("concentration");
    for (Teuchos::ParameterList::ConstIterator it = clist.begin(); it != clist.end(); ++it) {
      std::string name = it->first;
      if (clist.isSublist(name)) {
        Teuchos::ParameterList& src_list = clist.sublist(name);
        if (name == "coupling") {
          Teuchos::ParameterList::ConstIterator it1 = src_list.begin();
          std::string specname = it1->first;
          Teuchos::ParameterList& spec = src_list.sublist(specname);
          Teuchos::RCP<TransportDomainFunction> src =
            factory.Create(spec, "sink", AmanziMesh::Entity_kind::CELL, Teuchos::null);

          for (int i = 0; i < component_names_.size(); i++) {
            src->tcc_names().push_back(component_names_[i]);
            src->tcc_index().push_back(i);
          }
          src->set_state(S_);
          srcs_.push_back(src);

        } else {
          for (Teuchos::ParameterList::ConstIterator it1 = src_list.begin(); it1 != src_list.end();
               ++it1) {
            std::string specname = it1->first;
            Teuchos::ParameterList& spec = src_list.sublist(specname);
            Teuchos::RCP<TransportDomainFunction> src =
              factory.Create(spec, "sink", AmanziMesh::Entity_kind::CELL, Teuchos::null);

            src->tcc_names().push_back(name);
            src->tcc_index().push_back(0);

            src->set_state(S_);
            srcs_.push_back(src);
          }
        }
      }
    }
  }

  dt_debug_ = tp_list_->get<double>("maximum time step", TRANSPORT_LARGE_TIME_STEP);


  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Number of components: " << tcc->size() << std::endl
               << "cfl=" << cfl_ << " spatial/temporal discretization: " << spatial_disc_order
               << " " << temporal_disc_order << std::endl;
    *vo_->os() << vo_->color("green") << "Initalization of PK is complete." << vo_->reset()
               << std::endl
               << std::endl;
  }
}


/* ******************************************************************
* Initalized fields left by State and other PKs.
****************************************************************** */
void
SedimentTransport_PK::InitializeFields_(const Teuchos::Ptr<State>& S)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // set popular default values when flow PK is off
  // if (S->HasRecordSet(saturation_key_)) {
  //   if (S->GetW(saturation_key_)->owner() == passwd_) {
  //     if (!S->GetW(saturation_key_, passwd_)->initialized()) {
  //       S->GetW<CompositeVector>(saturation_key_, passwd_).PutScalar(1.0);
  //       S->GetField(saturation_key_, passwd_)->set_initialized();

  //       if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
  //           *vo_->os() << "initilized saturation_liquid to value 1.0" << std::endl;
  //     }
  //     InitializeFieldFromField_(prev_saturation_key_, saturation_key_, S, true, true);
  //   }
  //   else {
  //     if (S->GetField(prev_saturation_key_)->owner() == passwd_) {
  //       if (!S->GetField(prev_saturation_key_, passwd_)->initialized()) {
  //         InitializeFieldFromField_(prev_saturation_key_, saturation_key_, S, true, true);
  //         S->GetField(prev_saturation_key_, passwd_)->set_initialized();

  //         if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
  //           *vo_->os() << "initilized prev_saturation_liquid from saturation" << std::endl;
  //       }
  //     }
  //   }
  // }

  S->GetW<CompositeVector>(solid_residue_mass_key_, tag_current_, passwd_).PutScalar(0.0);
  S->GetRecordW(solid_residue_mass_key_, passwd_).set_initialized();
}


/* ****************************************************************
* Auxiliary initialization technique.
**************************************************************** */
void
SedimentTransport_PK::InitializeFieldFromField_(const std::string& field0,
                                                const Tag& tag0,
                                                const std::string& field1,
                                                const Tag& tag1,
                                                const Teuchos::Ptr<State>& S,
                                                bool call_evaluator,
                                                bool overwrite)
{
  if (S->HasRecordSet(field0)) {
    if (S->GetRecord(field0, tag0).owner() == passwd_) {
      if ((!S->GetRecord(field0, tag0).initialized()) || (overwrite)) {
        if (call_evaluator) S->GetEvaluator(field1, tag1).Update(*S, passwd_);

        const CompositeVector& f1 = *S->GetPtr<CompositeVector>(field1, tag1);
        CompositeVector& f0 = *S->GetPtrW<CompositeVector>(field0, tag0, passwd_);

        double vmin0, vmax0, vavg0;
        double vmin1, vmax1, vavg1;

        f0 = f1;

        S->GetRecordW(field0, tag0, passwd_).set_initialized();
        if ((vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) && (!overwrite)) {
          *vo_->os() << "initiliazed " << field0 << " to " << field1 << std::endl;
        }
      }
    }
  }
}


/* ******************************************************************
* Inialization of various transport structures.
****************************************************************** */
void
SedimentTransport_PK::InitializeAll_()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // global transport parameters
  cfl_ = tp_list_->get<double>("cfl", 1.0);

  spatial_disc_order = tp_list_->get<int>("spatial discretization order", 1);
  if (spatial_disc_order < 1 || spatial_disc_order > 2) spatial_disc_order = 1;
  temporal_disc_order = tp_list_->get<int>("temporal discretization order", 1);
  if (temporal_disc_order < 1 || temporal_disc_order > 2) temporal_disc_order = 1;

  num_aqueous = tp_list_->get<int>("number of sediment components", component_names_.size());

  // mass_solutes_exact_.assign(num_aqueous + num_gaseous, 0.0);
  // mass_solutes_source_.assign(num_aqueous + num_gaseous, 0.0);
  // mass_solutes_bc_.assign(num_aqueous + num_gaseous, 0.0);
  // mass_solutes_stepstart_.assign(num_aqueous + num_gaseous, 0.0);

  if (tp_list_->isParameter("runtime diagnostics: regions")) {
    runtime_regions_ =
      tp_list_->get<Teuchos::Array<std::string>>("runtime diagnostics: regions").toVector();
  }

  internal_tests = tp_list_->get<std::string>("enable internal tests", "no") == "yes";
  // tests_tolerance = tp_list_->get<double>("internal tests tolerance", TRANSPORT_CONCENTRATION_OVERSHOOT);
  // dt_debug_ = tp_list_->get<double>("maximum time step", TRANSPORT_LARGE_TIME_STEP);
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
SedimentTransport_PK::StableTimeStep()
{
  S_->Get<CompositeVector>(flux_key_, Tags::DEFAULT).ScatterMasterToGhosted("face");

  flux_ = S_->Get<CompositeVector>(flux_key_, Tags::DEFAULT).ViewComponent("face", true);
  //*flux_copy_ = *flux_; // copy flux vector from S_next_ to S_;

  IdentifyUpwindCells();

  tcc = S_->GetPtrW<CompositeVector>(tcc_key_, tag_current_, passwd_);
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");

  // loop over faces and accumulate upwinding fluxes
  std::vector<double> total_outflux(ncells_wghost, 0.0);

  for (int f = 0; f < nfaces_wghost; f++) {
    int c = (*upwind_cell_)[f];
    if (c >= 0) total_outflux[c] += fabs((*flux_)[0][f]);
  }

  Sinks2TotalOutFlux(tcc_prev, total_outflux, 0, num_aqueous - 1);

  // loop over cells and calculate minimal time step
  double vol, outflux, dt_cell;
  vol = 0;
  dt_ = dt_cell = TRANSPORT_LARGE_TIME_STEP;
  int cmin_dt = 0;
  for (int c = 0; c < ncells_owned; c++) {
    outflux = total_outflux[c];

    if ((outflux > 0) && ((*ws_prev_)[0][c] > 1e-6) && ((*ws_)[0][c] > 1e-6)) {
      vol = mesh_->getCellVolume(c);
      dt_cell = vol * (*mol_dens_)[0][c] * std::min((*ws_prev_)[0][c], (*ws_)[0][c]) / outflux;
    }
    if (dt_cell < dt_) {
      dt_ = dt_cell;
      cmin_dt = c;
    }
  }


  if (spatial_disc_order == 2) dt_ /= 2;

  // communicate global time step
  double dt_tmp = dt_;
#ifdef HAVE_MPI
  const Epetra_Comm& comm = ws_prev_->Comm();
  comm.MinAll(&dt_tmp, &dt_, 1);
#endif

  // incorporate developers and CFL constraints
  dt_ = std::min(dt_, dt_debug_);
  dt_ *= cfl_;

  //   //print optional diagnostics using maximum cell id as the filter
  //   if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
  //     int cmin_dt_unique = (fabs(dt_tmp * cfl_ - dt_) < 1e-6 * dt_) ? cmin_dt : -1;

  // #ifdef HAVE_MPI
  //     int cmin_dt_tmp = cmin_dt_unique;
  //     comm.MaxAll(&cmin_dt_tmp, &cmin_dt_unique, 1);
  // #endif
  //     if (cmin_dt == cmin_dt_unique) {
  //       const AmanziGeometry::Point& p = mesh_->getCellCentroid(cmin_dt);

  //       Teuchos::OSTab tab = vo_->getOSTab();
  //       *vo_->os() << "cell " << cmin_dt << " has smallest dt, (" << p[0] << ", " << p[1];
  //       if (p.dim() == 3) *vo_->os() << ", " << p[2];
  //       *vo_->os() << ")" << std::endl;
  //     }
  //   }
  return dt_;
}


/* *******************************************************************
* Estimate returns last time step unless it is zero.
******************************************************************* */
double
SedimentTransport_PK::get_dt()
{
  if (subcycling_) {
    return 1e+99;
  } else {
    //  flux_ = S_next_->Get<CompositeVector>(flux_key_).ViewComponent("face", true);
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
bool
SedimentTransport_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool failed = false;
  double dt_MPC = t_new - t_old;
  Teuchos::OSTab tab = vo_->getOSTab();

  if (S_->HasEvaluator(saturation_key_, tag_current_)) {
    S_->GetEvaluator(saturation_key_, tag_next_).Update(*S_, saturation_key_);
    S_->GetEvaluator(saturation_key_, tag_current_).Update(*S_, saturation_key_);
  }
  ws_ = S_->Get<CompositeVector>(saturation_key_, tag_next_).ViewComponent("cell", false);
  ws_prev_ = S_->Get<CompositeVector>(saturation_key_, tag_current_).ViewComponent("cell", false);


  if (S_->HasEvaluator(molar_density_key_, tag_next_)) {
    S_->GetEvaluator(molar_density_key_, tag_next_).Update(*S_, molar_density_key_);
  }
  mol_dens_ = S_->Get<CompositeVector>(molar_density_key_, tag_next_).ViewComponent("cell", false);


  solid_qty_ = S_->GetW<CompositeVector>(solid_residue_mass_key_, tag_next_, passwd_)
                 .ViewComponent("cell", false);

  // We use original tcc and make a copy of it later if needed.
  tcc = S_->GetPtrW<CompositeVector>(tcc_key_, tag_current_, passwd_);
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");

  // calculate stable time step
  double dt_shift = 0.0, dt_global = dt_MPC;
  double time = t_old;
  if (time >= 0.0) {
    t_physics_ = time;
    dt_shift = time - S_->get_time(tag_current_);
    dt_global = S_->get_time(tag_next_) - S_->get_time(tag_current_);
    AMANZI_ASSERT(std::abs(dt_global - dt_MPC) < 1.e-4);
  }

  StableTimeStep();
  double dt_stable = dt_; // advance routines override dt_

  int interpolate_ws = 0; // (dt_ < dt_global) ? 1 : 0;

  interpolate_ws = (dt_ < dt_global) ? 1 : 0;

  // start subcycling
  double dt_sum = 0.0;
  double dt_cycle;
  if (interpolate_ws) {
    dt_cycle = std::min(dt_stable, dt_MPC);
    InterpolateCellVector(*ws_prev_, *ws_, dt_shift, dt_global, *ws_subcycle_start);
    InterpolateCellVector(*ws_prev_, *ws_, dt_shift + dt_cycle, dt_global, *ws_subcycle_end);
    ws_start = ws_subcycle_start;
    ws_end = ws_subcycle_end;
    mol_dens_start = mol_dens_;
    mol_dens_end = mol_dens_;
  } else {
    dt_cycle = dt_MPC;
    ws_start = ws_prev_;
    ws_end = ws_;
    mol_dens_start = mol_dens_;
    mol_dens_end = mol_dens_;
  }


  for (int c = 0; c < ncells_owned; c++) {
    double vol_ws_den;
    vol_ws_den = mesh_->getCellVolume(c) * (*ws_prev_)[0][c] * (*mol_dens_)[0][c];
    mass_sediment_stepstart_ = tcc_prev[0][c] * vol_ws_den;
  }


  int ncycles = 0, swap = 1;


  while (dt_sum < dt_MPC - 1e-5) {
    // update boundary conditions
    time = t_physics_ + dt_cycle / 2;
    for (int i = 0; i < bcs_.size(); i++) { bcs_[i]->Compute(time, time); }

    double dt_try = dt_MPC - dt_sum;
    double tol = 1e-10 * (dt_try + dt_stable);
    bool final_cycle = false;
    if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
      *vo_->os() << std::setprecision(10) << "dt_MPC " << dt_MPC << " dt_cycle " << dt_cycle
                 << " dt_sum " << dt_sum << " dt_stable " << dt_stable << " dt_try " << dt_try
                 << " " << dt_try - (dt_stable + tol) << " tol " << tol << "\n";
    }

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
        ws_start = ws_subcycle_start;
        ws_end = ws_subcycle_end;
        mol_dens_start = mol_dens_;
        mol_dens_end = mol_dens_;

        double dt_int = dt_sum + dt_shift;
        InterpolateCellVector(*ws_prev_, *ws_, dt_int, dt_global, *ws_subcycle_end);
        //InterpolateCellVector(*mol_dens_prev_, *mol_dens_, dt_int, dt_global, *mol_dens_subcycle_end);
      } else { // Initial water saturation is in 'end'.
        ws_start = ws_subcycle_end;
        ws_end = ws_subcycle_start;
        mol_dens_start = mol_dens_;
        mol_dens_end = mol_dens_;

        double dt_int = dt_sum + dt_shift;
        InterpolateCellVector(*ws_prev_, *ws_, dt_int, dt_global, *ws_subcycle_start);
        //InterpolateCellVector(*mol_dens_prev_, *mol_dens_, dt_int, dt_global, *mol_dens_subcycle_start);
      }
      swap = 1 - swap;
    }

    if (spatial_disc_order == 1) { // temporary solution (lipnikov@lanl.gov)
      AdvanceDonorUpwind(dt_cycle);
      // } else if (spatial_disc_order == 2 && temporal_disc_order == 1) {
      //   AdvanceSecondOrderUpwindRK1(dt_cycle);
      // } else if (spatial_disc_order == 2 && temporal_disc_order == 2) {
      //   AdvanceSecondOrderUpwindRK2(dt_cycle);
    }

    if (!final_cycle) { // rotate concentrations (we need new memory for tcc)
      tcc = Teuchos::RCP<CompositeVector>(new CompositeVector(*tcc_tmp));
    }

    ncycles++;
  }


  dt_ = dt_stable; // restore the original time step (just in case)

  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", false);

  Advance_Diffusion(t_old, t_new);


  // statistics output
  nsubcycles = ncycles;
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << ncycles << " sub-cycles, dt_stable=" << units_.OutputTime(dt_stable)
               << " [sec]  dt_MPC=" << units_.OutputTime(dt_MPC) << " [sec]" << std::endl;

    //VV_PrintSoluteExtrema(tcc_next, dt_MPC);
  }

  return failed;
}


void
SedimentTransport_PK ::Advance_Diffusion(double t_old, double t_new)
{
  double dt_MPC = t_new - t_old;
  // We define tracer as the species #0 as calculate some statistics.
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", false);
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");
  int num_components = tcc_prev.NumVectors();

  bool flag_diffusion(true);

  if (flag_diffusion) {
    Teuchos::ParameterList& op_list =
      tp_list_->sublist("operators").sublist("diffusion operator").sublist("matrix");

    Teuchos::RCP<Operators::BCs> bc_dummy =
      Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));

    // default boundary conditions (none inside domain and Neumann on its boundary)
    auto& bc_model = bc_dummy->bc_model();
    auto& bc_value = bc_dummy->bc_value();
    PopulateBoundaryData(bc_model, bc_value, -1);

    Operators::PDE_DiffusionFactory opfactory;
    Teuchos::RCP<Operators::PDE_Diffusion> op1 = opfactory.Create(op_list, mesh_, bc_dummy);
    op1->SetBCs(bc_dummy, bc_dummy);
    Teuchos::RCP<Operators::Operator> op = op1->global_operator();
    Teuchos::RCP<Operators::PDE_Accumulation> op2 =
      Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, op));

    const CompositeVectorSpace& cvs = op->DomainMap();
    CompositeVector sol(cvs), factor(cvs), factor0(cvs), source(cvs), zero(cvs);
    zero.PutScalar(0.0);

    // instantiate solver

    S_->GetEvaluator(horiz_mixing_key_, tag_current_).Update(*S_, horiz_mixing_key_);

    CalculateDiffusionTensor_(*km_, *ws_, *mol_dens_);

    int phase, num_itrs(0);
    bool flag_op1(true);
    double md_change, md_old(0.0), md_new, residual(0.0);

    // Disperse and diffuse aqueous components
    for (int i = 0; i < num_aqueous; i++) {
      // set initial guess
      Epetra_MultiVector& sol_cell = *sol.ViewComponent("cell");
      for (int c = 0; c < ncells_owned; c++) { sol_cell[0][c] = tcc_next[i][c]; }
      if (sol.HasComponent("face")) { sol.ViewComponent("face")->PutScalar(0.0); }

      op->Init();
      Teuchos::RCP<std::vector<WhetStone::Tensor>> Dptr = Teuchos::rcpFromRef(D_);
      op1->Setup(Dptr, Teuchos::null, Teuchos::null);
      op1->UpdateMatrices(Teuchos::null, Teuchos::null);

      // add accumulation term
      Epetra_MultiVector& fac = *factor.ViewComponent("cell");
      for (int c = 0; c < ncells_owned; c++) { fac[0][c] = (*ws_)[0][c] * (*mol_dens_)[0][c]; }
      op2->AddAccumulationDelta(sol, factor, factor, dt_MPC, "cell");

      op1->ApplyBCs(true, true, true);

      CompositeVector& rhs = *op->rhs();
      int ierr = op->ApplyInverse(rhs, sol);

      if (ierr < 0) {
        Errors::Message msg("SedimentTransport_PK solver failed with message: \"");
        msg << op->returned_code_string() << "\"";
        Exceptions::amanzi_throw(msg);
      }

      residual += op->residual();
      num_itrs += op->num_itrs();

      for (int c = 0; c < ncells_owned; c++) { tcc_next[i][c] = sol_cell[0][c]; }
    }

    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "sediment transport solver: ||r||=" << residual / num_components
                 << " itrs=" << num_itrs / num_components << std::endl;
    }
  }
}


/* *******************************************************************
* Copy the advected tcc field to the state.
******************************************************************* */
void
SedimentTransport_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  Teuchos::RCP<CompositeVector> tcc;
  tcc = S->GetPtrW<CompositeVector>(tcc_key_, tag_next_, passwd_);
  *tcc = *tcc_tmp;


  /// Not sure that it is necessary DSV
  InitializeFieldFromField_(
    prev_saturation_key_, tag_current_, saturation_key_, tag_next_, S.ptr(), false, true);

  // Copy to S_ as well
  // tcc = S_->GetPtrW<CompositeVector>(tcc_key_, passwd_);
  // *tcc = *tcc_tmp;
}


/* *******************************************************************
 * A simple first-order transport method
 ****************************************************************** */
void
SedimentTransport_PK::AdvanceDonorUpwind(double dt_cycle)
{
  dt_ = dt_cycle; // overwrite the maximum stable transport step
  mass_sediment_source_ = 0;
  mass_sediment_bc_ = 0;

  // populating next state of concentrations
  tcc->ScatterMasterToGhosted("cell");
  Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);

  // prepare conservative state in master and slave cells
  double vol_ws_den, tcc_flux;
  double mass_start = 0., tmp1, mass;

  // We advect only aqueous components.
  int num_advect = num_aqueous;

  for (int c = 0; c < ncells_owned; c++) {
    vol_ws_den = mesh_->getCellVolume(c) * (*ws_start)[0][c] * (*mol_dens_start)[0][c];
    for (int i = 0; i < num_advect; i++) {
      (*conserve_qty_)[i][c] = tcc_prev[i][c] * vol_ws_den;
      // if ((vol_ws_den > water_tolerance_) && ((*solid_qty_)[i][c] > 0 )){   // Desolve solid residual into liquid
      //   double add_mass = std::min((*solid_qty_)[i][c], max_tcc_* vol_ws_den - (*conserve_qty_)[i][c]);
      //   (*solid_qty_)[i][c] -= add_mass;
      //   (*conserve_qty_)[i][c] += add_mass;
      // }
      mass_start += (*conserve_qty_)[i][c];
    }
  }

  tmp1 = mass_start;
  mesh_->getComm()->SumAll(&tmp1, &mass_start, 1);

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    if (domain_name_ == "surface")
      *vo_->os() << std::setprecision(10) << "Surface mass start " << mass_start << "\n";
    else
      *vo_->os() << std::setprecision(10) << "Subsurface mass start " << mass_start << "\n";
  }


  // advance all components at once
  for (int f = 0; f < nfaces_wghost; f++) { // loop over master and slave faces
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];

    double u = fabs((*flux_)[0][f]);

    if (c1 >= 0 && c1 < ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      for (int i = 0; i < num_advect; i++) {
        tcc_flux = dt_ * u * tcc_prev[i][c1];
        (*conserve_qty_)[i][c1] -= tcc_flux;
        (*conserve_qty_)[i][c2] += tcc_flux;
      }

    } else if (c1 >= 0 && c1 < ncells_owned && (c2 >= ncells_owned || c2 < 0)) {
      for (int i = 0; i < num_advect; i++) {
        tcc_flux = dt_ * u * tcc_prev[i][c1];
        (*conserve_qty_)[i][c1] -= tcc_flux;
        if (c2 < 0) mass_sediment_bc_ -= tcc_flux;
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
            mass_sediment_bc_ += tcc_flux;
          }
        }
      }
    }
  }


  // process external sources
  //if (srcs_.size() != 0) {
  double time = t_physics_;
  ComputeAddSourceTerms(time, dt_, *conserve_qty_, 0, num_advect - 1);
  //}

  // recover concentration from new conservative state
  for (int c = 0; c < ncells_owned; c++) {
    vol_ws_den = mesh_->getCellVolume(c) * (*ws_end)[0][c] * (*mol_dens_end)[0][c];
    for (int i = 0; i < num_advect; i++) {
      if ((*ws_end)[0][c] > water_tolerance_ && (*conserve_qty_)[i][c] > 0) {
        tcc_next[i][c] = (*conserve_qty_)[i][c] / vol_ws_den;

      } else {
        (*solid_qty_)[i][c] += std::max((*conserve_qty_)[i][c], 0.);
        tcc_next[i][c] = 0.;
      }
    }
  }

  double mass_final = 0;
  for (int c = 0; c < ncells_owned; c++) {
    for (int i = 0; i < num_advect; i++) { mass_final += (*conserve_qty_)[i][c]; }
  }

  tmp1 = mass_final;
  mesh_->getComm()->SumAll(&tmp1, &mass_final, 1);


  // update mass balance

  mass_sediment_exact_ += mass_sediment_source_ * dt_;
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    tmp1 = mass_sediment_bc_;
    mesh_->getComm()->SumAll(&tmp1, &mass_sediment_bc_, 1);
    // *vo_->os() << "*****************\n";
    // if (domain_name_ == "surface") *vo_->os()<<"Surface mass BC "<<mass_sediment_bc_<<"\n";
    // else *vo_->os() <<"Subsurface mass BC "<<mass_sediment_bc_<<"\n";

    // tmp1 = mass_sediment_source_;
    // mesh_->getComm()->SumAll(&tmp1, &mass_sediment_source_, 1);
    // if (domain_name_ == "surface") *vo_->os()<<"Surface mass_sediment source "<<mass_sediment_source_ * dt_<<"\n";
    // else *vo_->os() << "Subsurface mass_sediment source "<<mass_sediment_source_ * dt_<<"\n";
    // *vo_->os() << "*****************\n";
  }


  // if (internal_tests) {
  //   VV_CheckGEDproperty(*tcc_tmp->ViewComponent("cell"));
  // }

  // if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH){
  //   if (domain_name_ == "surface")  *vo_->os()<<"Surface mass final "<<mass_final<<"\n";
  //   else  *vo_->os()<<"Subsurface mass final "<<mass_final<<"\n";
  // }

  // if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
  //   *vo_->os()<<"mass error "<<abs(mass_final - (mass_start + mass_sediment_bc_ + mass_sediment_source_*dt_) )<<"\n";

  //if (abs(mass_final - (mass_start + mass_sediment_bc_[0] + mass_sediment_source_[0]*dt_) )/mass_final > 1e-6) exit(-1);
}


// /* *******************************************************************
//  * We have to advance each component independently due to different
//  * reconstructions. We use tcc when only owned data are needed and
//  * tcc_next when owned and ghost data. This is a special routine for
//  * transient flow and uses first-order time integrator.
//  ****************************************************************** */
// void SedimentTransport_PK::AdvanceSecondOrderUpwindRK1(double dt_cycle)
// {
//   dt_ = dt_cycle;  // overwrite the maximum stable transport step
//   mass_sediment_source_.assign(num_aqueous + num_gaseous, 0.0);

//   // work memory
//   const Epetra_Map& cmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::CELL,true);
//   Epetra_Vector f_component(cmap_wghost);

//   // distribute vector of concentrations
//   S_->Get<CompositeVector>(tcc_key_).ScatterMasterToGhosted("cell");
//   Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
//   Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);


//   Epetra_Vector ws_ratio(Copy, *ws_start, 0);
//   for (int c = 0; c < ncells_owned; c++){
//     double vol_phi_ws_den_end = mesh_->getCellVolume(c) * (*phi_)[0][c] * (*ws_end)[0][c] * (*mol_dens_end)[0][c];
//     if (vol_phi_ws_den_end > water_tolerance_)  {
//       double vol_phi_ws_den_start = mesh_->getCellVolume(c) * (*phi_)[0][c] * (*ws_start)[0][c] * (*mol_dens_start)[0][c];
//       if (vol_phi_ws_den_start > water_tolerance_){
//         ws_ratio[c] = ( (*ws_start)[0][c] * (*mol_dens_start)[0][c] )
//                     / ( (*ws_end)[0][c]   * (*mol_dens_end)[0][c]   );
//       }else{
//         ws_ratio[c] = 1;
//       }
//     }
//     else  ws_ratio[c]=0.;
//   }


//   // We advect only aqueous components.
//   int num_advect = num_aqueous;

//   for (int i = 0; i < num_advect; i++) {
//     current_component_ = i;  // needed by BJ

//     double T = t_physics_;
//     Epetra_Vector*& component = tcc_prev(i);
//     FunctionalTimeDerivative(T, *component, f_component);

//     for (int c = 0; c < ncells_owned; c++) {
//       tcc_next[i][c] = (tcc_prev[i][c] + dt_ * f_component[c]) * ws_ratio[c];

//       if (tcc_next[i][c] < 0){
//         double vol_phi_ws_den = mesh_->getCellVolume(c) * (*phi_)[0][c] * (*ws_end)[0][c] * (*mol_dens_end)[0][c];
//         (*solid_qty_)[i][c] += abs(tcc_next[i][c])*vol_phi_ws_den;
//         tcc_next[i][c] = 0.;
//       }
//     }
//   }

//   // update mass balance
//   for (int i = 0; i < num_aqueous + num_gaseous; i++) {
//     mass_sediment_exact_[i] += mass_sediment_source_[i] * dt_;
//   }

//   if (internal_tests) {
//     VV_CheckGEDproperty(*tcc_tmp->ViewComponent("cell"));
//   }
// }


// /* *******************************************************************
//  * We have to advance each component independently due to different
//  * reconstructions. This is a special routine for transient flow and
//  * uses second-order predictor-corrector time integrator.
//  ****************************************************************** */
// void SedimentTransport_PK::AdvanceSecondOrderUpwindRK2(double dt_cycle)
// {
//   dt_ = dt_cycle;  // overwrite the maximum stable transport step
//   mass_sediment_source_.assign(num_aqueous + num_gaseous, 0.0);

//   // work memory
//   const Epetra_Map& cmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::CELL,true);
//   Epetra_Vector f_component(cmap_wghost);//,  f_component2(cmap_wghost);

//   // distribute old vector of concentrations
//   S_->Get<CompositeVector>(tcc_key_).ScatterMasterToGhosted("cell");
//   Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
//   Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);

//   Epetra_Vector ws_ratio(Copy, *ws_start, 0);
//   for (int c = 0; c < ncells_owned; c++){
//     if ((*ws_end)[0][c] > 1e-10)  {
//       if ((*ws_start)[0][c] > 1e-10){
//         ws_ratio[c] = ( (*ws_start)[0][c] * (*mol_dens_start)[0][c] )
//                     / ( (*ws_end)[0][c]   * (*mol_dens_end)[0][c]   );
//       }else{
//         ws_ratio[c] = 1;
//       }
//     }
//     else  ws_ratio[c]=0.;
//   }

//   // We advect only aqueous components.
//   int num_advect = num_aqueous;

//   // predictor step
//   for (int i = 0; i < num_advect; i++) {
//     current_component_ = i;  // needed by BJ

//     double T = t_physics_;
//     Epetra_Vector*& component = tcc_prev(i);
//     FunctionalTimeDerivative(T, *component, f_component);

//     for (int c = 0; c < ncells_owned; c++) {
//       tcc_next[i][c] = (tcc_prev[i][c] + dt_ * f_component[c]) * ws_ratio[c];
//       //if (tcc_next[i][c] < 0) tcc_next[i][c] = 0.;

//     }
//   }

//   tcc_tmp->ScatterMasterToGhosted("cell");

//   //if (domain_name_ == "surface") {
//   //*vo_->os()<<"after predictor ToTaL "<<domain_name_<<" :"<<std::setprecision(10)<<ComputeSolute( tcc_next, 0)<<"\n";
//   //}

//   // corrector step
//   for (int i = 0; i < num_advect; i++) {
//     current_component_ = i;  // needed by BJ

//     double T = t_physics_;
//     Epetra_Vector*& component = tcc_next(i);
//     FunctionalTimeDerivative(T, *component, f_component);

//     for (int c = 0; c < ncells_owned; c++) {
//       double value = (tcc_prev[i][c] + dt_ * f_component[c]) * ws_ratio[c];
//       tcc_next[i][c] = (tcc_next[i][c] + value) / 2;
//       if (tcc_next[i][c] < 0){
//         double vol_phi_ws_den = mesh_->getCellVolume(c) * (*phi_)[0][c] * (*ws_end)[0][c] * (*mol_dens_end)[0][c];
//         (*solid_qty_)[i][c] += abs(tcc_next[i][c])*vol_phi_ws_den;
//         tcc_next[i][c] = 0.;
//       }

//     }
//   }

//   // f_component2.Update(-1, f_component, 1.);
//   // double diff_norm;
//   // f_component2.NormInf(&diff_norm);
//   // *vo_->os()<<domain_name_<<" difference "<<diff_norm<<"\n";
//   //if (domain_name_ == "surface") {
//   //*vo_->os()<<"after corrector ToTaL "<<domain_name_<<" :"<<std::setprecision(10)<<ComputeSolute( tcc_next, 0)<<"\n";
//   //}

//   // update mass balance
//   for (int i = 0; i < num_aqueous + num_gaseous; i++) {
//     mass_sediment_exact_[i] += mass_sediment_source_[i] * dt_ / 2;
//   }

//   if (internal_tests) {
//     VV_CheckGEDproperty(*tcc_tmp->ViewComponent("cell"));
//   }

// }


// /* *******************************************************************
// * Advance each component independently due to different field
// * reconstructions. This routine uses generic explicit time integrator.
// ******************************************************************* */
// // void SedimentTransport_PK::AdvanceSecondOrderUpwindRKn(double dt_cycle)
// // {
// //   dt_ = dt_cycle;  // overwrite the maximum stable transport step

// //   S_->Get<CompositeVector>("total_component_concentration").ScatterMasterToGhosted("cell");
// //   Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
// //   Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);

// //   // define time integration method
// //   auto ti_method = Explicit_TI::forward_euler;
// //   if (temporal_disc_order == 2) {
// //     ti_method = Explicit_TI::heun_euler;
// //   } else if (temporal_disc_order == 3) {
// //     ti_method = Explicit_TI::kutta_3rd_order;
// //   } else if (temporal_disc_order == 3) {
// //     ti_method = Explicit_TI::runge_kutta_4th_order;
// //   }

// //   // We interpolate ws using dt which becomes local time.
// //   double T = 0.0;
// //   // We advect only aqueous components.
// //   int ncomponents = num_aqueous;

// //   for (int i = 0; i < ncomponents; i++) {
// //     current_component_ = i;  // it is needed in BJ called inside RK:fun

// //     Epetra_Vector*& component_prev = tcc_prev(i);
// //     Epetra_Vector*& component_next = tcc_next(i);

// //     Explicit_TI::RK<Epetra_Vector> TVD_RK(*this, ti_method, *component_prev);
// //     TVD_RK.TimeStep(T, dt_, *component_prev, *component_next);
// //   }
// // }


/* ******************************************************************
* Computes source and sink terms and adds them to vector tcc.
* Returns mass rate for the tracer.
* The routine treats two cases of tcc with one and all components.
****************************************************************** */
void
SedimentTransport_PK::ComputeAddSourceTerms(double tp,
                                            double dtp,
                                            Epetra_MultiVector& tcc,
                                            int n0,
                                            int n1)
{
  int num_vectors = tcc.NumVectors();
  int nsrcs = srcs_.size();

  double mass1 = 0., mass2 = 0., add_mass = 0., tmp1;
  bool chg;

  chg = S_->GetEvaluator(sd_trapping_key_, tag_next_).Update(*S_, sd_trapping_key_);
  const Epetra_MultiVector& Q_dt =
    *S_->GetPtr<CompositeVector>(sd_trapping_key_, tag_next_)->ViewComponent("cell", false);

  chg = S_->GetEvaluator(sd_settling_key_, tag_next_).Update(*S_, sd_settling_key_);
  const Epetra_MultiVector& Q_ds =
    *S_->GetPtr<CompositeVector>(sd_settling_key_, tag_next_)->ViewComponent("cell", false);

  chg = S_->GetEvaluator(sd_erosion_key_, tag_next_).Update(*S_, sd_erosion_key_);
  const Epetra_MultiVector& Q_e =
    *S_->GetPtr<CompositeVector>(sd_erosion_key_, tag_next_)->ViewComponent("cell", false);

  chg = S_->GetEvaluator(sd_organic_key_, tag_next_).Update(*S_, sd_organic_key_);
  const Epetra_MultiVector& Q_db =
    *S_->GetPtr<CompositeVector>(sd_organic_key_, tag_next_)->ViewComponent("cell", false);

  Epetra_MultiVector& dz = *S_->GetW<CompositeVector>(elevation_increase_key_, tag_next_, "state")
                              .ViewComponent("cell", false);

  const Epetra_MultiVector& poro =
    *S_->Get<CompositeVector>(porosity_key_, tag_next_).ViewComponent("cell", false);


  for (int c = 0; c < ncells_owned; c++) {
    double value = mesh_->getCellVolume(c) * (Q_e[0][c] - Q_dt[0][c] - Q_ds[0][c]);
    tcc[0][c] += value * dtp;
    mass_sediment_source_ += value;
    dz[0][c] += ((1. / sediment_density_) * ((Q_dt[0][c] + Q_ds[0][c]) - Q_e[0][c]) + +Q_db[0][c]) *
                dtp / (1 - poro[0][c]);
  }


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
          value = mesh_->getCellVolume(c) * values[k];
        }

        //add_mass += dtp * value;
        tcc[imap][c] += dtp * value;
        mass_sediment_source_ += value;
      }
    }
  }
}

void
SedimentTransport_PK::Sinks2TotalOutFlux(Epetra_MultiVector& tcc,
                                         std::vector<double>& total_outflux,
                                         int n0,
                                         int n1)
{
  std::vector<double> sink_add(ncells_wghost, 0.0);
  //Assumption that there is only one sink per component per cell
  double t0 = S_->get_time(tag_current_);
  int num_vectors = tcc.NumVectors();
  int nsrcs = srcs_.size();

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

        if ((values[k] < 0) && (tcc[imap][c] > 0)) {
          if (srcs_[m]->name() == "domain coupling") {
            // if (values[k]<0) {
            val = std::max(val, fabs(values[k]) / tcc[imap][c]);
            //}
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
SedimentTransport_PK::PopulateBoundaryData(std::vector<int>& bc_model,
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
SedimentTransport_PK::IdentifyUpwindCells()
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

// void SedimentTransport_PK::ComputeVolumeDarcyFlux(Teuchos::RCP<const Epetra_MultiVector> flux,
//                                               Teuchos::RCP<const Epetra_MultiVector> molar_density,
//                                               Teuchos::RCP<Epetra_MultiVector>& vol_darcy_flux){

//   int nfaces_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

//   for (int f = 0; f < nfaces_wghost ; f++){
//     auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
//     double n_liq=0.;
//     for (int c=0; c<cells.size();c++) n_liq += (*molar_density)[0][c];
//     n_liq /= cells.size();
//     if (n_liq > 0) (*vol_darcy_flux)[0][f] = (*flux_)[0][f]/n_liq;
//     else (*vol_darcy_flux)[0][f] = 0.;
//   }

// }
//


/* *******************************************************************
* Interpolate linearly in time between two values v0 and v1. The time
* is measuared relative to value v0; so that v1 is at time dt. The
* interpolated data are at time dt_int.
******************************************************************* */
void
SedimentTransport_PK::InterpolateCellVector(const Epetra_MultiVector& v0,
                                            const Epetra_MultiVector& v1,
                                            double dt_int,
                                            double dt,
                                            Epetra_MultiVector& v_int)
{
  double a = dt_int / dt;
  double b = 1.0 - a;
  v_int.Update(b, v0, a, v1, 0.);
}

} // namespace SedimentTransport
} // namespace Amanzi
