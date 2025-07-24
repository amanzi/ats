/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*
  Special version of transport capabilities to model sediment transport in surface flows.

*/

#include <algorithm>
#include <vector>

#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "errors.hh"
#include "PK_Helpers.hh"
#include "TensorVector.hh"

#include "sediment_transport_pk.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* New constructor compatible with new MPC framework.
****************************************************************** */
SedimentTransport_PK::SedimentTransport_PK(Teuchos::ParameterList& pk_tree,
                                           const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                           const Teuchos::RCP<State>& S,
                                           const Teuchos::RCP<TreeVector>& solution)

  : PK(pk_tree, glist, S, solution), Transport_ATS(pk_tree, glist, S, solution)
{}


void
SedimentTransport_PK::parseParameterList()
{
  PK_Physical_Default::parseParameterList();

  num_components_ = 1;
  component_names_.resize(num_components_);
  component_names_[0] = "sediment";
  num_aqueous_ = 1;

  molar_masses_ =
    readParameterMapByComponent(plist_->sublist("component molar masses [kg / mol C]"), 1.0);
  tcc_max_ = readParameterMapByComponent(
    plist_->sublist("component maximum concentration [mol C / mol H2O]"), -1.0);

  sd_trapping_key_ = Keys::readKey(*plist_, domain_, "trapping rate", "trapping_rate");
  sd_settling_key_ = Keys::readKey(*plist_, domain_, "settling rate", "settling_rate");
  sd_erosion_key_ = Keys::readKey(*plist_, domain_, "erosion rate", "erosion_rate");
  sd_organic_key_ = Keys::readKey(*plist_, domain_, "organic rate", "organic_rate");

  elevation_increase_key_ = Keys::readKey(*plist_, domain_, "deformation", "deformation");
  requireEvaluatorAtNext(elevation_increase_key_, tag_next_, *S_, name_);

  porosity_key_ = Keys::readKey(*plist_, domain_, "soil porosity", "soil_porosity");

  has_dispersion_ = false;
  has_diffusion_ = false;

  if (plist_->isSublist("sediment diffusion coefficient [m^2 s^-1]")) {
    has_diffusion_ = true;
    molec_diff_ =
      readParameterMapByComponent(plist_->sublist("sediment diffusion coefficient [m^2 s^-1]"), 0.);
    tortuosity_ = readParameterMapByPhase(plist_->sublist("tortuosity [-]"), 1.);
  }

  // keys, dependencies, etc
  // -- flux -- only needed at new time, evaluator controlled elsewhere
  water_flux_key_ = Keys::readKey(*plist_, domain_, "water flux", "water_flux");

  mass_flux_advection_key_ =
    Keys::readKey(*plist_, domain_, "mass flux advection", "mass_flux_advection");
  requireEvaluatorAtNext(mass_flux_advection_key_, tag_next_, *S_, name_);

  mass_flux_diffusion_key_ =
    Keys::readKey(*plist_, domain_, "mass flux diffusion", "mass_flux_diffusion");
  requireEvaluatorAtNext(mass_flux_diffusion_key_, tag_next_, *S_, name_);

  // -- liquid water content - need at new time, copy at current time
  lwc_key_ = Keys::readKey(*plist_, domain_, "liquid water content", "water_content");
  requireEvaluatorAtCurrent(lwc_key_, tag_current_, *S_, name_);

  water_src_key_ = Keys::readKey(*plist_, domain_, "water source", "water_source");
  cv_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");

  // workspace, no evaluator
  conserve_qty_key_ =
    Keys::readKey(*plist_, domain_, "conserved quantity", "total_component_quantity");
  requireEvaluatorAtNext(conserve_qty_key_, tag_next_, *S_, name_);

  solid_residue_mass_key_ =
    Keys::readKey(*plist_, domain_, "solid residue", "solid_residue_quantity");

  // other parameters
  // -- a small amount of water, used to define when we are going to completely dry out a grid cell
  water_tolerance_ = plist_->get<double>("water tolerance [mol H2O / m^d]", 1e-6);

  // global transport parameters
  cfl_ = plist_->get<double>("cfl", 1.0);
  dt_max_ = plist_->get<double>("maximum timestep [s]", TRANSPORT_LARGE_TIME_STEP);


  // dispersion coefficient tensor
  dispersion_tensor_key_ =
    Keys::readKey(*plist_, domain_, "sediment dispersion coefficient", "horizontal_mixing");
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
    msg << "Transport_ATS: \"temporal discretization order\" must be 1 or 2, not "
        << temporal_disc_order_;
    Exceptions::amanzi_throw(msg);
  }
}


void
SedimentTransport_PK::SetupPhysicalEvaluators_()
{
  Transport_ATS::SetupPhysicalEvaluators_();

  // setup sediment transport sources
  requireEvaluatorAtNext(sd_organic_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireEvaluatorAtNext(sd_trapping_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireEvaluatorAtNext(sd_settling_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireEvaluatorAtNext(sd_erosion_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  S_->Require<double>("sediment_density", tag_next_);

  if (!elevation_increase_key_.empty()) {
    // note, these are done at the OUTER step, hard-coded to NEXT because dz is
    // accumulated through the inner steps, then applied to the mesh at the
    // outer step.
    requireEvaluatorAtNext(elevation_increase_key_, tag_next_, *S_, name_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    requireEvaluatorAtNext(porosity_key_, tag_next_, *S_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
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

  if (D_rank >= 0) {
    CompositeVectorSpace D_space;
    D_space.SetMesh(mesh_)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    D_ = Teuchos::rcp(new TensorVector(D_space, D_dim, D_rank, false));
  }

  // compute the stable dt for the initial timestep
  dt_stable_ = ComputeStableTimeStep_();

  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    *vo_->os() << "cfl=" << cfl_ << " spatial/temporal discretization: " << adv_spatial_disc_order_
               << " " << temporal_disc_order_ << std::endl
               << std::endl;
  }
}


void
SedimentTransport_PK::AddSourceTerms_(double t0, double t1, Epetra_MultiVector& conserve_qty)
{
  int ncells_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  sediment_density_ = S_->Get<double>("sediment_density", tag_next_);

  bool change1 = S_->GetEvaluator(sd_trapping_key_, tag_next_).Update(*S_, sd_trapping_key_);
  const Epetra_MultiVector& Q_dt =
    *S_->GetPtr<CompositeVector>(sd_trapping_key_, tag_next_)->ViewComponent("cell", false);

  bool change2 = S_->GetEvaluator(sd_settling_key_, tag_next_).Update(*S_, sd_settling_key_);
  const Epetra_MultiVector& Q_ds =
    *S_->GetPtr<CompositeVector>(sd_settling_key_, tag_next_)->ViewComponent("cell", false);

  bool change3 = S_->GetEvaluator(sd_erosion_key_, tag_next_).Update(*S_, sd_erosion_key_);
  const Epetra_MultiVector& Q_e =
    *S_->GetPtr<CompositeVector>(sd_erosion_key_, tag_next_)->ViewComponent("cell", false);

  for (int c = 0; c < ncells_owned; c++) {
    double value = mesh_->getCellVolume(c) * (Q_e[0][c] - Q_dt[0][c] - Q_ds[0][c]); /// m^3/s
    conserve_qty[0][c] += sediment_density_ * value * (t1 - t0) / molar_masses_["sediment"];
  }


  if (!elevation_increase_key_.empty()) {
    {
      bool change4 = S_->GetEvaluator(sd_organic_key_, tag_next_).Update(*S_, sd_organic_key_);
      const Epetra_MultiVector& Q_db =
        *S_->GetPtr<CompositeVector>(sd_organic_key_, tag_next_)->ViewComponent("cell", false);

      // Poro and DZ are hard-coded as the NEXT tag -- this should be the outer step's NEXT tag.
      const Epetra_MultiVector& poro =
        *S_->Get<CompositeVector>(porosity_key_, tag_next_).ViewComponent("cell", false);

      // NOTE: we do note zero this out here, because it gets accumulated across the outer step size
      Epetra_MultiVector& dz = *S_->GetW<CompositeVector>(elevation_increase_key_, tag_next_, name_)
                                  .ViewComponent("cell", false);
      // dz.PutScalar(0.);

      for (int c = 0; c < ncells_owned; c++) {
        dz[0][c] +=
          ((1. / sediment_density_) * ((Q_dt[0][c] + Q_ds[0][c]) - Q_e[0][c]) + Q_db[0][c]) *
          (t1 - t0) / (1 - poro[0][c]);
      }
    }

    // do NOT mark this as changed here, as that would also recompute the mesh
    // despite no deformation in the inner step
    // changedEvaluatorPrimary(elevation_increase_key_, tag_next_, *S_);
    db_->WriteVector("dz", S_->GetPtr<CompositeVector>(elevation_increase_key_, tag_next_).ptr());
  }
}

} // namespace Transport
} // namespace Amanzi
