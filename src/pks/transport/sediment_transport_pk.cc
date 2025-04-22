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

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

#include "BCs.hh"
#include "errors.hh"
#include "Explicit_TI_RK.hh"
#include "Evaluator.hh"
#include "TensorVector.hh"
#include "Mesh.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_Diffusion.hh"
#include "PDE_Accumulation.hh"
#include "PK_DomainFunctionFactory.hh"
#include "PK_Utils.hh"

#include "pk_helpers.hh"

#include "sediment_transport_pk.hh"
#include "TransportDomainFunction.hh"


namespace Amanzi {
namespace Transport {

/* ******************************************************************
* New constructor compatible with new MPC framework.
****************************************************************** */
SedimentTransport_PK::SedimentTransport_PK(Teuchos::ParameterList& pk_tree,
                                           const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                           const Teuchos::RCP<State>& S,
                                           const Teuchos::RCP<TreeVector>& solution)
 
  : PK(pk_tree, glist, S, solution),
    Transport_ATS(pk_tree, glist, S, solution)
{

}


void
SedimentTransport_PK::parseParameterList()
{
  PK_Physical_Default::parseParameterList();

  num_components_ = 1;
  num_aqueous_ = num_components_;
  num_gaseous_ = 0;
  component_names_.resize(num_components_);
  component_names_[0] = "sediment";

  molar_masses_ = readParameterMapByComponent(plist_->sublist("component molar masses [kg / mol C]"), 1.0);
  tcc_max_ = readParameterMapByComponent(plist_->sublist("component maximum concentration [mol C / mol H2O]"), -1.0);

  sd_trapping_key_ =  Keys::readKey(*plist_, domain_, "trapping rate", "trapping_rate");
  sd_settling_key_ =  Keys::readKey(*plist_, domain_, "settling rate", "settling_rate");
  sd_erosion_key_ =  Keys::readKey(*plist_, domain_, "erosion rate", "erosion_rate");
  sd_organic_key_ =  Keys::readKey(*plist_, domain_, "organic rate", "organic_rate");
  
  biomass_key_ =  Keys::readKey(*plist_, domain_, "biomass", "biomass");
  plant_area_key_ = Keys::getKey(domain_, "plant_area");
  stem_diameter_key_ = Keys::getKey(domain_, "stem_diameter");
  stem_height_key_ = Keys::getKey(domain_, "stem_height");
  stem_density_key_ = Keys::getKey(domain_, "stem_density");
  
  //horiz_mixing_key_ = Keys::readKey(*plist_, domain_, "horizontal mixing", "horizontal_mixing");
  elevation_increase_key_ = Keys::readKey(*plist_, domain_, "deformation", "deformation");
  porosity_key_ = Keys::readKey(*plist_, domain_, "soil porosity", "soil_porosity");
  
  has_dispersion_ = false;
  has_diffusion_ = false;

  if (plist_->isSublist("sediment diffusion coefficient [m^2 s^-1]")) {
    has_diffusion_ = true;
    molec_diff_ = readParameterMapByComponent(plist_->sublist("sediment diffusion coefficient [m^2 s^-1]"), 0.);
  }
  
    // keys, dependencies, etc
  // -- flux -- only needed at new time, evaluator controlled elsewhere
  flux_key_ = Keys::readKey(*plist_, domain_, "water flux", "water_flux");

  // -- liquid water content - need at new time, copy at current time
  lwc_key_ = Keys::readKey(*plist_, domain_, "liquid water content", "water_content");
  requireAtCurrent(lwc_key_, tag_current_, *S_, name_);

  water_src_key_ = Keys::readKey(*plist_, domain_, "water source", "water_source");
  cv_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");

  // workspace, no evaluator
  conserve_qty_key_ =
    Keys::readKey(*plist_, domain_, "conserved quantity", "total_component_quantity");
  requireAtNext(conserve_qty_key_, tag_next_, *S_, name_);

  solid_residue_mass_key_ = Keys::readKey(*plist_, domain_, "solid residue", "solid_residue_quantity");
  
    // other parameters
  // -- a small amount of water, used to define when we are going to completely dry out a grid cell
  water_tolerance_ = plist_->get<double>("water tolerance [mol H2O / m^d]", 1e-6);

  // global transport parameters
  cfl_ = plist_->get<double>("cfl", 1.0);
  dt_max_ = plist_->get<double>("maximum timestep [s]", TRANSPORT_LARGE_TIME_STEP);

  // dispersion coefficient tensor
  // dispersion_tensor_key_ = Keys::readKey(*plist_, domain_, "horizontal mixing", "horizontal_mixing");
  //  has_dispersion_ = S_->HasEvaluatorList(dispersion_tensor_key_);

  // dispersion coefficient tensor
  dispersion_tensor_key_ = Keys::readKey(*plist_, domain_, "dispersion coefficient", "dispersion_coefficient");
  has_dispersion_ = S_->HasEvaluatorList(dispersion_tensor_key_);

  
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

  sediment_density_ = plist_->get<double>("sediment density [kg m^-3]");
  
  db_ = Teuchos::rcp(new Debugger(mesh_, name_, *plist_));

}


/* ******************************************************************
* Define structure of this PK.
****************************************************************** */
void
SedimentTransport_PK::SetupTransport_()
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


  // setup sources

  requireAtNext(sd_organic_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(sd_trapping_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(sd_settling_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(sd_erosion_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(biomass_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(plant_area_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(stem_diameter_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(stem_height_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(stem_density_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  S_->Require<double>("MSL", tag_next_);
   
  requireAtNext(porosity_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  
    //
  // create boundary conditions
  // --------------------------------------------------------------------------------
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
          bc_list.set("entity_gid_out", gid);

          Teuchos::RCP<TransportDomainFunction> bc = factory.Create(
            bc_list, "boundary concentration", AmanziMesh::Entity_kind::FACE, Teuchos::null, tag_current_);

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
  }

}

/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.
****************************************************************** */
void
SedimentTransport_PK::Initialize()
{

  Teuchos::OSTab tab = vo_->getOSTab();
  PK_Physical_Default::Initialize();

  // initialize missed fields
  Transport_ATS::InitializeFields_();
  
  // can also now setup the joint diffusion/dispersion workspace tensor
  int D_rank = -1;
  int D_dim = mesh_->getSpaceDimension();
  if (has_diffusion_) D_rank = 1; // scalar

  //   // dispersion rank is 1 or 2
  //D_rank = S_->Require<TensorVector, TensorVector_Factory>(dispersion_tensor_key_, tag_next_).get_rank();

  if (D_rank >= 0) {
    CompositeVectorSpace D_space;
    D_space.SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    D_ = Teuchos::rcp(new TensorVector(D_space, D_dim, D_rank, false));
  } 

  // compute the stable dt for the initial timestep
  dt_stable_ = ComputeStableTimeStep_();
  
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {   
    *vo_->os() << "cfl=" << cfl_ << " spatial/temporal discretization: " << adv_spatial_disc_order_
               << " " << temporal_disc_order_ << std::endl << std::endl;
    *vo_->os() << vo_->color("green") << "Initalization of PK is complete." << vo_->reset()
               << std::endl
               << std::endl;
    
  }
  
}

void
SedimentTransport_PK::AddSourceTerms_(double t0,
        double t1,
        Epetra_MultiVector& conserve_qty,
        int n0,
        int n1)
{


  double mass1 = 0., mass2 = 0., add_mass = 0., tmp1;
  bool chg;

  int ncells_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  
  chg = S_->GetEvaluator(sd_trapping_key_, tag_next_).Update(*S_, sd_trapping_key_);
  const Epetra_MultiVector& Q_dt =
    *S_->GetPtr<CompositeVector>(sd_trapping_key_, tag_next_)->ViewComponent("cell", false);

  chg = S_->GetEvaluator(sd_settling_key_, tag_next_).Update(*S_, sd_settling_key_);
  const Epetra_MultiVector& Q_ds =
    *S_->GetPtr<CompositeVector>(sd_settling_key_, tag_next_)->ViewComponent("cell", false);

  chg = S_->GetEvaluator(sd_erosion_key_, tag_next_).Update(*S_, sd_erosion_key_);
  const Epetra_MultiVector& Q_e =
    *S_->GetPtr<CompositeVector>(sd_erosion_key_, tag_next_)->ViewComponent("cell", false);

  // chg = S_->GetEvaluator(sd_organic_key_, tag_next_).Update(*S_, sd_organic_key_);
  // const Epetra_MultiVector& Q_db =
  //   *S_->GetPtr<CompositeVector>(sd_organic_key_, tag_next_)->ViewComponent("cell", false);

  for (int c = 0; c < ncells_owned; c++) {
    double value = mesh_->getCellVolume(c) * (Q_e[0][c] - Q_dt[0][c] - Q_ds[0][c]); /// m^3/s
    conserve_qty[0][c] += sediment_density_ * value * (t1 - t0)/ molar_masses_["sediment"] ;
  }


    // if (S_->HasRecordSet(elevation_increase_key_)) {
    //   Epetra_MultiVector& dz = *S_->GetW<CompositeVector>(elevation_increase_key_, tag_next_, elevation_increase_key_)
    //                           .ViewComponent("cell", false);

    //   const Epetra_MultiVector& poro =
    //     *S_->Get<CompositeVector>(porosity_key_, tag_next_).ViewComponent("cell", false);

    //   for (int c = 0; c < ncells_owned; c++) {
    //     dz[0][c] += ((1. / sediment_density_) * ((Q_dt[0][c] + Q_ds[0][c]) - Q_e[0][c]) + Q_db[0][c]) *
    //       dtp / (1 - poro[0][c]);
    //   }

    // }
    
}
  

/* *******************************************************************
* MPC will call this function to advance the transport state.
* Efficient subcycling requires to calculate an intermediate state of
* saturation only once, which leads to a leap-frog-type algorithm.
******************************************************************* */
// bool
// SedimentTransport_PK::AdvanceStep(double t_old, double t_new, bool reinit)
// {
//   bool failed = false;
//   double dt_MPC = t_new - t_old;
//   Teuchos::OSTab tab = vo_->getOSTab();

//   if (S_->HasEvaluator(saturation_key_, tag_current_)) {
//     S_->GetEvaluator(saturation_key_, tag_next_).Update(*S_, saturation_key_);
//     S_->GetEvaluator(saturation_key_, tag_current_).Update(*S_, saturation_key_);
//   }
//   ws_ = S_->Get<CompositeVector>(saturation_key_, tag_next_).ViewComponent("cell", false);
//   ws_prev_ = S_->Get<CompositeVector>(saturation_key_, tag_current_).ViewComponent("cell", false);


//   if (S_->HasEvaluator(molar_density_key_, tag_next_)) {
//     S_->GetEvaluator(molar_density_key_, tag_next_).Update(*S_, molar_density_key_);
//   }
//   mol_dens_ = S_->Get<CompositeVector>(molar_density_key_, tag_next_).ViewComponent("cell", false);


//   solid_qty_ = S_->GetW<CompositeVector>(solid_residue_mass_key_, tag_next_, passwd_)
//                  .ViewComponent("cell", false);

//   // We use original tcc and make a copy of it later if needed.
//   tcc = S_->GetPtrW<CompositeVector>(tcc_key_, tag_current_, passwd_);
//   Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");

//   // calculate stable timestep
//   double dt_shift = 0.0, dt_global = dt_MPC;
//   double time = t_old;
//   if (time >= 0.0) {
//     t_physics_ = time;
//     dt_shift = time - S_->get_time(tag_current_);
//     dt_global = S_->get_time(tag_next_) - S_->get_time(tag_current_);
//     AMANZI_ASSERT(std::abs(dt_global - dt_MPC) < 1.e-4);
//   }

//   StableTimeStep_();
//   double dt_stable = dt_; // advance routines override dt_

//   int interpolate_ws = 0; // (dt_ < dt_global) ? 1 : 0;

//   interpolate_ws = (dt_ < dt_global) ? 1 : 0;

//   // start subcycling
//   double dt_sum = 0.0;
//   double dt_cycle;
//   if (interpolate_ws) {
//     dt_cycle = std::min(dt_stable, dt_MPC);
//     InterpolateCellVector(*ws_prev_, *ws_, dt_shift, dt_global, *ws_subcycle_start);
//     InterpolateCellVector(*ws_prev_, *ws_, dt_shift + dt_cycle, dt_global, *ws_subcycle_end);
//     ws_start = ws_subcycle_start;
//     ws_end = ws_subcycle_end;
//     mol_dens_start = mol_dens_;
//     mol_dens_end = mol_dens_;
//   } else {
//     dt_cycle = dt_MPC;
//     ws_start = ws_prev_;
//     ws_end = ws_;
//     mol_dens_start = mol_dens_;
//     mol_dens_end = mol_dens_;
//   }


//   for (int c = 0; c < ncells_owned; c++) {
//     double vol_ws_den;
//     vol_ws_den = mesh_->getCellVolume(c) * (*ws_prev_)[0][c] * (*mol_dens_)[0][c];
//     mass_sediment_stepstart_ = tcc_prev[0][c] * vol_ws_den;
//   }


//   int ncycles = 0, swap = 1;


//   while (dt_sum < dt_MPC - 1e-5) {
//     // update boundary conditions
//     time = t_physics_ + dt_cycle / 2;
//     for (int i = 0; i < bcs_.size(); i++) { bcs_[i]->Compute(time, time); }

//     double dt_try = dt_MPC - dt_sum;
//     double tol = 1e-10 * (dt_try + dt_stable);
//     bool final_cycle = false;
//     if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
//       *vo_->os() << std::setprecision(10) << "dt_MPC " << dt_MPC << " dt_cycle " << dt_cycle
//                  << " dt_sum " << dt_sum << " dt_stable " << dt_stable << " dt_try " << dt_try
//                  << " " << dt_try - (dt_stable + tol) << " tol " << tol << "\n";
//     }

//     if (dt_try >= 2 * dt_stable) {
//       dt_cycle = dt_stable;
//     } else if (dt_try > dt_stable + tol) {
//       dt_cycle = dt_try / 2;
//     } else {
//       dt_cycle = dt_try;
//       final_cycle = true;
//     }

//     t_physics_ += dt_cycle;
//     dt_sum += dt_cycle;

//     if (interpolate_ws) {
//       if (swap) { // Initial water saturation is in 'start'.
//         ws_start = ws_subcycle_start;
//         ws_end = ws_subcycle_end;
//         mol_dens_start = mol_dens_;
//         mol_dens_end = mol_dens_;

//         double dt_int = dt_sum + dt_shift;
//         InterpolateCellVector(*ws_prev_, *ws_, dt_int, dt_global, *ws_subcycle_end);
//         //InterpolateCellVector(*mol_dens_prev_, *mol_dens_, dt_int, dt_global, *mol_dens_subcycle_end);
//       } else { // Initial water saturation is in 'end'.
//         ws_start = ws_subcycle_end;
//         ws_end = ws_subcycle_start;
//         mol_dens_start = mol_dens_;
//         mol_dens_end = mol_dens_;

//         double dt_int = dt_sum + dt_shift;
//         InterpolateCellVector(*ws_prev_, *ws_, dt_int, dt_global, *ws_subcycle_start);
//         //InterpolateCellVector(*mol_dens_prev_, *mol_dens_, dt_int, dt_global, *mol_dens_subcycle_start);
//       }
//       swap = 1 - swap;
//     }

//     if (spatial_disc_order == 1) { // temporary solution (lipnikov@lanl.gov)
//       AdvanceDonorUpwind_(dt_cycle);
//       // } else if (spatial_disc_order == 2 && temporal_disc_order == 1) {
//       //   AdvanceSecondOrderUpwindRK1_(dt_cycle);
//       // } else if (spatial_disc_order == 2 && temporal_disc_order == 2) {
//       //   AdvanceSecondOrderUpwindRK2_(dt_cycle);
//     }

//     if (!final_cycle) { // rotate concentrations (we need new memory for tcc)
//       tcc = Teuchos::RCP<CompositeVector>(new CompositeVector(*tcc_tmp));
//     }

//     ncycles++;
//   }


//   dt_ = dt_stable; // restore the original timestep (just in case)

//   Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", false);

//   Advance_Diffusion(t_old, t_new);


//   // statistics output
//   nsubcycles = ncycles;
//   if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
//     Teuchos::OSTab tab = vo_->getOSTab();
//     *vo_->os() << ncycles << " sub-cycles, dt_stable=" << units_.OutputTime(dt_stable)
//                << " [sec]  dt_MPC=" << units_.OutputTime(dt_MPC) << " [sec]" << std::endl;

//     //PrintSoluteExtrema(tcc_next, dt_MPC);
//   }

//   return failed;
// }


// void
// SedimentTransport_PK ::Advance_Diffusion(double t_old, double t_new)
// {
//   double dt_MPC = t_new - t_old;
//   // We define tracer as the species #0 as calculate some statistics.
//   Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", false);
//   Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell");
//   int num_components_ = tcc_prev.NumVectors();

//   bool flag_diffusion(true);

//   if (flag_diffusion) {
//     Teuchos::ParameterList& op_list =
//       tp_list_->sublist("operators").sublist("diffusion operator").sublist("matrix");

//     Teuchos::RCP<Operators::BCs> bc_dummy = Teuchos::rcp(
//       new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));

//     // default boundary conditions (none inside domain and Neumann on its boundary)
//     auto& bc_model = bc_dummy->bc_model();
//     auto& bc_value = bc_dummy->bc_value();
//     PopulateBoundaryData_(bc_model, bc_value, -1);

//     Operators::PDE_DiffusionFactory opfactory;
//     Teuchos::RCP<Operators::PDE_Diffusion> op1 = opfactory.Create(op_list, mesh_, bc_dummy);
//     op1->SetBCs(bc_dummy, bc_dummy);
//     Teuchos::RCP<Operators::Operator> op = op1->global_operator();
//     Teuchos::RCP<Operators::PDE_Accumulation> op2 =
//       Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, op));

//     const CompositeVectorSpace& cvs = op->DomainMap();
//     CompositeVector sol(cvs), factor(cvs), factor0(cvs), source(cvs), zero(cvs);
//     zero.PutScalar(0.0);

//     // instantiate solver

//     S_->GetEvaluator(horiz_mixing_key_, tag_current_).Update(*S_, horiz_mixing_key_);

//     CalculateDiffusionTensor_(*km_, *ws_, *mol_dens_);

//     int phase, num_itrs(0);
//     bool flag_op1(true);
//     double md_change, md_old(0.0), md_new, residual(0.0);

//     // Disperse and diffuse aqueous components
//     for (int i = 0; i < num_aqueous_; i++) {
//       // set initial guess
//       Epetra_MultiVector& sol_cell = *sol.ViewComponent("cell");
//       for (int c = 0; c < ncells_owned; c++) { sol_cell[0][c] = tcc_next[i][c]; }
//       if (sol.HasComponent("face")) { sol.ViewComponent("face")->PutScalar(0.0); }

//       op->Init();
//       Teuchos::RCP<std::vector<WhetStone::Tensor>> Dptr = Teuchos::rcpFromRef(D_);
//       op1->Setup(Dptr, Teuchos::null, Teuchos::null);
//       op1->UpdateMatrices(Teuchos::null, Teuchos::null);

//       // add accumulation term
//       Epetra_MultiVector& fac = *factor.ViewComponent("cell");
//       for (int c = 0; c < ncells_owned; c++) { fac[0][c] = (*ws_)[0][c] * (*mol_dens_)[0][c]; }
//       op2->AddAccumulationDelta(sol, factor, factor, dt_MPC, "cell");

//       op1->ApplyBCs(true, true, true);

//       CompositeVector& rhs = *op->rhs();
//       int ierr = op->ApplyInverse(rhs, sol);

//       if (ierr < 0) {
//         Errors::Message msg("SedimentTransport_PK solver failed with message: \"");
//         msg << op->returned_code_string() << "\"";
//         Exceptions::amanzi_throw(msg);
//       }

//       residual += op->residual();
//       num_itrs += op->num_itrs();

//       for (int c = 0; c < ncells_owned; c++) { tcc_next[i][c] = sol_cell[0][c]; }
//     }

//     if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
//       Teuchos::OSTab tab = vo_->getOSTab();
//       *vo_->os() << "sediment transport solver: ||r||=" << residual / num_components_
//                  << " itrs=" << num_itrs / num_components_ << std::endl;
//     }
//   }
// }


// /* *******************************************************************
// * Copy the advected tcc field to the state.
// ******************************************************************* */
// void
// SedimentTransport_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
// {
//   Teuchos::RCP<CompositeVector> tcc;
//   tcc = S->GetPtrW<CompositeVector>(tcc_key_, tag_next_, passwd_);
//   *tcc = *tcc_tmp;


//   /// Not sure that it is necessary DSV
//   InitializeFieldFromField_(
//     prev_saturation_key_, tag_current_, saturation_key_, tag_next_, S.ptr(), false, true);

//   // Copy to S_ as well
//   // tcc = S_->GetPtrW<CompositeVector>(tcc_key_, passwd_);
//   // *tcc = *tcc_tmp;
// }


// /* *******************************************************************
//  * A simple first-order transport method
//  ****************************************************************** */
// void
// SedimentTransport_PK::AdvanceDonorUpwind_(double dt_cycle)
// {
//   dt_ = dt_cycle; // overwrite the maximum stable transport step
//   mass_sediment_source_ = 0;
//   mass_sediment_bc_ = 0;

//   // populating next state of concentrations
//   tcc->ScatterMasterToGhosted("cell");
//   Epetra_MultiVector& tcc_prev = *tcc->ViewComponent("cell", true);
//   Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", true);

//   // prepare conservative state in master and slave cells
//   double vol_ws_den, tcc_flux;
//   double mass_start = 0., tmp1, mass;

//   // We advect only aqueous components.
//   int num_advect_ = num_aqueous_;

//   for (int c = 0; c < ncells_owned; c++) {
//     vol_ws_den = mesh_->getCellVolume(c) * (*ws_start)[0][c] * (*mol_dens_start)[0][c];
//     for (int i = 0; i < num_advect_; i++) {
//       (*conserve_qty_)[i][c] = tcc_prev[i][c] * vol_ws_den;
//       // if ((vol_ws_den > water_tolerance_) && ((*solid_qty_)[i][c] > 0 )){   // Desolve solid residual into liquid
//       //   double add_mass = std::min((*solid_qty_)[i][c], max_tcc_* vol_ws_den - (*conserve_qty_)[i][c]);
//       //   (*solid_qty_)[i][c] -= add_mass;
//       //   (*conserve_qty_)[i][c] += add_mass;
//       // }
//       mass_start += (*conserve_qty_)[i][c];
//     }
//   }

//   tmp1 = mass_start;
//   mesh_->getComm()->SumAll(&tmp1, &mass_start, 1);

//   if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
//     if (domain_name_ == "surface")
//       *vo_->os() << std::setprecision(10) << "Surface mass start " << mass_start << "\n";
//     else
//       *vo_->os() << std::setprecision(10) << "Subsurface mass start " << mass_start << "\n";
//   }


//   // advance all components at once
//   for (int f = 0; f < nfaces_wghost; f++) { // loop over master and slave faces
//     int c1 = (*upwind_cell_)[f];
//     int c2 = (*downwind_cell_)[f];

//     double u = fabs((*flux_)[0][f]);

//     if (c1 >= 0 && c1 < ncells_owned && c2 >= 0 && c2 < ncells_owned) {
//       for (int i = 0; i < num_advect_; i++) {
//         tcc_flux = dt_ * u * tcc_prev[i][c1];
//         (*conserve_qty_)[i][c1] -= tcc_flux;
//         (*conserve_qty_)[i][c2] += tcc_flux;
//       }

//     } else if (c1 >= 0 && c1 < ncells_owned && (c2 >= ncells_owned || c2 < 0)) {
//       for (int i = 0; i < num_advect_; i++) {
//         tcc_flux = dt_ * u * tcc_prev[i][c1];
//         (*conserve_qty_)[i][c1] -= tcc_flux;
//         if (c2 < 0) mass_sediment_bc_ -= tcc_flux;
//       }

//     } else if (c1 >= ncells_owned && c2 >= 0 && c2 < ncells_owned) {
//       for (int i = 0; i < num_advect_; i++) {
//         tcc_flux = dt_ * u * tcc_prev[i][c1];
//         (*conserve_qty_)[i][c2] += tcc_flux;
//       }
//     }
//   }


//   // loop over exterior boundary sets
//   for (int m = 0; m < bcs_.size(); m++) {
//     std::vector<int>& tcc_index = bcs_[m]->tcc_index();
//     int ncomp = tcc_index.size();

//     for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
//       int f = it->first;
//       std::vector<double>& values = it->second;
//       int c2 = (*downwind_cell_)[f];
//       if (c2 >= 0) {
//         double u = fabs((*flux_)[0][f]);
//         for (int i = 0; i < ncomp; i++) {
//           int k = tcc_index[i];
//           if (k < num_advect_) {
//             tcc_flux = dt_ * u * values[i];
//             (*conserve_qty_)[k][c2] += tcc_flux;
//             mass_sediment_bc_ += tcc_flux;
//           }
//         }
//       }
//     }
//   }


//   // process external sources
//   //if (srcs_.size() != 0) {
//   double time = t_physics_;
//   ComputeAddSourceTerms_(time, dt_, *conserve_qty_, 0, num_advect_ - 1);
//   //}

//   // recover concentration from new conservative state
//   for (int c = 0; c < ncells_owned; c++) {
//     vol_ws_den = mesh_->getCellVolume(c) * (*ws_end)[0][c] * (*mol_dens_end)[0][c];
//     for (int i = 0; i < num_advect_; i++) {
//       if ((*ws_end)[0][c] > water_tolerance_ && (*conserve_qty_)[i][c] > 0) {
//         tcc_next[i][c] = (*conserve_qty_)[i][c] / vol_ws_den;

//       } else {
//         (*solid_qty_)[i][c] += std::max((*conserve_qty_)[i][c], 0.);
//         tcc_next[i][c] = 0.;
//       }
//     }
//   }

//   double mass_final = 0;
//   for (int c = 0; c < ncells_owned; c++) {
//     for (int i = 0; i < num_advect_; i++) { mass_final += (*conserve_qty_)[i][c]; }
//   }

//   tmp1 = mass_final;
//   mesh_->getComm()->SumAll(&tmp1, &mass_final, 1);


//   // update mass balance

//   mass_sediment_exact_ += mass_sediment_source_ * dt_;
//   if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
//     tmp1 = mass_sediment_bc_;
//     mesh_->getComm()->SumAll(&tmp1, &mass_sediment_bc_, 1);
//     // *vo_->os() << "*****************\n";
//     // if (domain_name_ == "surface") *vo_->os()<<"Surface mass BC "<<mass_sediment_bc_<<"\n";
//     // else *vo_->os() <<"Subsurface mass BC "<<mass_sediment_bc_<<"\n";

//     // tmp1 = mass_sediment_source_;
//     // mesh_->getComm()->SumAll(&tmp1, &mass_sediment_source_, 1);
//     // if (domain_name_ == "surface") *vo_->os()<<"Surface mass_sediment source "<<mass_sediment_source_ * dt_<<"\n";
//     // else *vo_->os() << "Subsurface mass_sediment source "<<mass_sediment_source_ * dt_<<"\n";
//     // *vo_->os() << "*****************\n";
//   }


//   // if (internal_tests) {
//   //   CheckGEDProperty(*tcc_tmp->ViewComponent("cell"));
//   // }

//   // if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH){
//   //   if (domain_name_ == "surface")  *vo_->os()<<"Surface mass final "<<mass_final<<"\n";
//   //   else  *vo_->os()<<"Subsurface mass final "<<mass_final<<"\n";
//   // }

//   // if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
//   //   *vo_->os()<<"mass error "<<abs(mass_final - (mass_start + mass_sediment_bc_ + mass_sediment_source_*dt_) )<<"\n";

//   //if (abs(mass_final - (mass_start + mass_sediment_bc_[0] + mass_sediment_source_[0]*dt_) )/mass_final > 1e-6) exit(-1);
// }


// /* *******************************************************************
//  * We have to advance each component independently due to different
//  * reconstructions. We use tcc when only owned data are needed and
//  * tcc_next when owned and ghost data. This is a special routine for
//  * transient flow and uses first-order time integrator.
//  ****************************************************************** */
// void SedimentTransport_PK::AdvanceSecondOrderUpwindRK1_(double dt_cycle)
// {
//   dt_ = dt_cycle;  // overwrite the maximum stable transport step
//   mass_sediment_source_.assign(num_aqueous_ + num_gaseous_, 0.0);

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
//   int num_advect_ = num_aqueous_;

//   for (int i = 0; i < num_advect_; i++) {
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
//   for (int i = 0; i < num_aqueous_ + num_gaseous_; i++) {
//     mass_sediment_exact_[i] += mass_sediment_source_[i] * dt_;
//   }

//   if (internal_tests) {
//     CheckGEDProperty(*tcc_tmp->ViewComponent("cell"));
//   }
// }


// /* *******************************************************************
//  * We have to advance each component independently due to different
//  * reconstructions. This is a special routine for transient flow and
//  * uses second-order predictor-corrector time integrator.
//  ****************************************************************** */
// void SedimentTransport_PK::AdvanceSecondOrderUpwindRK2_(double dt_cycle)
// {
//   dt_ = dt_cycle;  // overwrite the maximum stable transport step
//   mass_sediment_source_.assign(num_aqueous_ + num_gaseous_, 0.0);

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
//   int num_advect_ = num_aqueous_;

//   // predictor step
//   for (int i = 0; i < num_advect_; i++) {
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
//   for (int i = 0; i < num_advect_; i++) {
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
//   for (int i = 0; i < num_aqueous_ + num_gaseous_; i++) {
//     mass_sediment_exact_[i] += mass_sediment_source_[i] * dt_ / 2;
//   }

//   if (internal_tests) {
//     CheckGEDProperty(*tcc_tmp->ViewComponent("cell"));
//   }

// }


// /* *******************************************************************
// * Advance each component independently due to different field
// * reconstructions. This routine uses generic explicit time integrator.
// ******************************************************************* */
// // void SedimentTransport_PK::AdvanceSecondOrderUpwindRKn_(double dt_cycle)
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
// //   int ncomponents = num_aqueous_;

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
// void
// SedimentTransport_PK::ComputeAddSourceTerms_(double tp,
//                                             double dtp,
//                                             Epetra_MultiVector& tcc,
//                                             int n0,
//                                             int n1)
// {
//   int num_vectors = tcc.NumVectors();
//   int nsrcs = srcs_.size();

//   double mass1 = 0., mass2 = 0., add_mass = 0., tmp1;
//   bool chg;

//   chg = S_->GetEvaluator(sd_trapping_key_, tag_next_).Update(*S_, sd_trapping_key_);
//   const Epetra_MultiVector& Q_dt =
//     *S_->GetPtr<CompositeVector>(sd_trapping_key_, tag_next_)->ViewComponent("cell", false);

//   chg = S_->GetEvaluator(sd_settling_key_, tag_next_).Update(*S_, sd_settling_key_);
//   const Epetra_MultiVector& Q_ds =
//     *S_->GetPtr<CompositeVector>(sd_settling_key_, tag_next_)->ViewComponent("cell", false);

//   chg = S_->GetEvaluator(sd_erosion_key_, tag_next_).Update(*S_, sd_erosion_key_);
//   const Epetra_MultiVector& Q_e =
//     *S_->GetPtr<CompositeVector>(sd_erosion_key_, tag_next_)->ViewComponent("cell", false);

//   chg = S_->GetEvaluator(sd_organic_key_, tag_next_).Update(*S_, sd_organic_key_);
//   const Epetra_MultiVector& Q_db =
//     *S_->GetPtr<CompositeVector>(sd_organic_key_, tag_next_)->ViewComponent("cell", false);

//   Epetra_MultiVector& dz = *S_->GetW<CompositeVector>(elevation_increase_key_, tag_next_, "state")
//                               .ViewComponent("cell", false);

//   const Epetra_MultiVector& poro =
//     *S_->Get<CompositeVector>(porosity_key_, tag_next_).ViewComponent("cell", false);


//   for (int c = 0; c < ncells_owned; c++) {
//     double value = mesh_->getCellVolume(c) * (Q_e[0][c] - Q_dt[0][c] - Q_ds[0][c]);
//     tcc[0][c] += value * dtp;
//     mass_sediment_source_ += value;
//     dz[0][c] += ((1. / sediment_density_) * ((Q_dt[0][c] + Q_ds[0][c]) - Q_e[0][c]) + +Q_db[0][c]) *
//                 dtp / (1 - poro[0][c]);
//   }


//   for (int m = 0; m < nsrcs; m++) {
//     double t0 = tp - dtp;
//     srcs_[m]->Compute(t0, tp);

//     std::vector<int> tcc_index = srcs_[m]->tcc_index();
//     for (auto it = srcs_[m]->begin(); it != srcs_[m]->end(); ++it) {
//       int c = it->first;
//       std::vector<double>& values = it->second;

//       if (c >= ncells_owned) continue;

//       for (int k = 0; k < tcc_index.size(); ++k) {
//         int i = tcc_index[k];
//         if (i < n0 || i > n1) continue;

//         int imap = i;
//         if (num_vectors == 1) imap = 0;

//         double value;
//         if (srcs_[m]->getType() == DomainFunction_kind::COUPLING) {
//           value = values[k];
//         } else {
//           value = mesh_->getCellVolume(c) * values[k];
//         }

//         //add_mass += dtp * value;
//         tcc[imap][c] += dtp * value;
//         mass_sediment_source_ += value;
//       }
//     }
//   }
// }

// void
// SedimentTransport_PK::Sinks2TotalOutFlux(Epetra_MultiVector& tcc,
//                                          std::vector<double>& total_outflux,
//                                          int n0,
//                                          int n1)
// {
//   std::vector<double> sink_add(ncells_wghost, 0.0);
//   //Assumption that there is only one sink per component per cell
//   double t0 = S_->get_time(tag_current_);
//   int num_vectors = tcc.NumVectors();
//   int nsrcs = srcs_.size();

//   for (int m = 0; m < nsrcs; m++) {
//     srcs_[m]->Compute(t0, t0);
//     std::vector<int> index = srcs_[m]->tcc_index();

//     for (auto it = srcs_[m]->begin(); it != srcs_[m]->end(); ++it) {
//       int c = it->first;
//       std::vector<double>& values = it->second;

//       double val = 0;
//       for (int k = 0; k < index.size(); ++k) {
//         int i = index[k];
//         if (i < n0 || i > n1) continue;

//         int imap = i;
//         if (num_vectors == 1) imap = 0;

//         if ((values[k] < 0) && (tcc[imap][c] > 0)) {
//           if (srcs_[m]->getType() == DomainFunction_kind::COUPLING) {
//             // if (values[k]<0) {
//             val = std::max(val, fabs(values[k]) / tcc[imap][c]);
//             //}
//           }
//         }
//       }
//       sink_add[c] = std::max(sink_add[c], val);
//     }
//   }

//   for (int c = 0; c < ncells_wghost; c++) total_outflux[c] += sink_add[c];
// }


/* *******************************************************************
* Populates operators' boundary data for given component.
* Returns true if at least one face was populated.
******************************************************************* */
// bool
// SedimentTransport_PK::PopulateBoundaryData_(std::vector<int>& bc_model,
//                                            std::vector<double>& bc_value,
//                                            int component)
// {
//   bool flag = false;

//   for (int i = 0; i < bc_model.size(); i++) {
//     bc_model[i] = Operators::OPERATOR_BC_NONE;
//     bc_value[i] = 0.0;
//   }

//   for (int f = 0; f < nfaces_wghost; f++) {
//     auto cells = mesh_->getFaceCells(f);
//     if (cells.size() == 1) bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
//   }

//   for (int m = 0; m < bcs_.size(); m++) {
//     std::vector<int>& tcc_index = bcs_[m]->tcc_index();
//     int ncomp = tcc_index.size();

//     for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
//       int f = it->first;
//       std::vector<double>& values = it->second;

//       for (int i = 0; i < ncomp; i++) {
//         int k = tcc_index[i];
//         if (k == component) {
//           bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
//           bc_value[f] = values[i];
//           flag = true;
//         }
//       }
//     }
//   }

//   return flag;
// }


// /* *******************************************************************
// * Identify flux direction based on orientation of the face normal
// * and sign of the  Darcy velocity.
// ******************************************************************* */
// void
// SedimentTransport_PK::IdentifyUpwindCells_()
// {
//   for (int f = 0; f < nfaces_wghost; f++) {
//     (*upwind_cell_)[f] = -1; // negative value indicates boundary
//     (*downwind_cell_)[f] = -1;
//   }

//   for (int c = 0; c < ncells_wghost; c++) {
//     const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);

//     for (int i = 0; i < faces.size(); i++) {
//       int f = faces[i];
//       double tmp = (*flux_)[0][f] * dirs[i];
//       if (tmp > 0.0) {
//         (*upwind_cell_)[f] = c;
//       } else if (tmp < 0.0) {
//         (*downwind_cell_)[f] = c;
//       } else if (dirs[i] > 0) {
//         (*upwind_cell_)[f] = c;
//       } else {
//         (*downwind_cell_)[f] = c;
//       }
//     }
//   }
// }

// // void SedimentTransport_PK::ComputeVolumeDarcyFlux(Teuchos::RCP<const Epetra_MultiVector> flux,
// //                                               Teuchos::RCP<const Epetra_MultiVector> molar_density,
// //                                               Teuchos::RCP<Epetra_MultiVector>& vol_darcy_flux){

// //   int nfaces_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

// //   for (int f = 0; f < nfaces_wghost ; f++){
// //     auto cells = mesh_->getFaceCells(f);
// //     double n_liq=0.;
// //     for (int c=0; c<cells.size();c++) n_liq += (*molar_density)[0][c];
// //     n_liq /= cells.size();
// //     if (n_liq > 0) (*vol_darcy_flux)[0][f] = (*flux_)[0][f]/n_liq;
// //     else (*vol_darcy_flux)[0][f] = 0.;
// //   }

// // }
// //


// /* *******************************************************************
// * Interpolate linearly in time between two values v0 and v1. The time
// * is measuared relative to value v0; so that v1 is at time dt. The
// * interpolated data are at time dt_int.
// ******************************************************************* */
// void
// SedimentTransport_PK::InterpolateCellVector(const Epetra_MultiVector& v0,
//                                             const Epetra_MultiVector& v1,
//                                             double dt_int,
//                                             double dt,
//                                             Epetra_MultiVector& v_int)
// {
//   double a = dt_int / dt;
//   double b = 1.0 - a;
//   v_int.Update(b, v0, a, v1, 0.);
// }


// /* *******************************************************************
// * Calculate dispersive tensor from given Darcy fluxes. The flux is
// * assumed to be scaled by face area.
// ******************************************************************* */
// void
// SedimentTransport_PK::CalculateDiffusionTensor_(const Epetra_MultiVector& km,
//                                                 const Epetra_MultiVector& ws,
//                                                 const Epetra_MultiVector& mol_density)
// {
//   D_.resize(ncells_owned);
//   for (int c = 0; c < ncells_owned; c++) D_[c].Init(dim, 1);

//   //AmanziGeometry::Point velocity(dim);
//   AmanziMesh::Entity_ID_List faces;
//   WhetStone::MFD3D_Diffusion mfd3d(mesh_);

//   for (int c = 0; c < ncells_owned; ++c) {
//     double mol_den = mol_density[0][c];
//     double ponded_depth = ws[0][c];
//     D_[c].PutScalar(km[0][c] * mol_den * ponded_depth);
//   }
// }

} // namespace Transport
} // namespace Amanzi
