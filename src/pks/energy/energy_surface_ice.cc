/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Process kernel for energy equation for overland flow.
------------------------------------------------------------------------- */
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Debugger.hh"
#include "thermal_conductivity_surface_evaluator.hh"
#include "enthalpy_evaluator.hh"
#include "energy_bc_factory.hh"
#include "Function.hh"
#include "FunctionFactory.hh"
#include "EvaluatorPrimary.hh"
#include "Op.hh"

#include "pk_helpers.hh"
#include "energy_surface_ice.hh"

namespace Amanzi {
namespace Energy {

#define DEBUG_FLAG 1

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------

EnergySurfaceIce::EnergySurfaceIce(Teuchos::ParameterList& FElist,
                                   const Teuchos::RCP<Teuchos::ParameterList>& plist,
                                   const Teuchos::RCP<State>& S,
                                   const Teuchos::RCP<TreeVector>& solution)
  : PK(FElist, plist, S, solution), EnergyBase(FElist, plist, S, solution)
{
  if (!plist_->isParameter("conserved quantity key suffix"))
    plist_->set("conserved quantity key suffix", "energy");

  domain_ss_ = Keys::readDomainHint(*plist_, domain_, "surface", "subsurface");
}


// -------------------------------------------------------------
// Create the physical evaluators for energy, enthalpy, thermal
// conductivity, and any sources.
// -------------------------------------------------------------
void
EnergySurfaceIce::SetupPhysicalEvaluators_()
{
  // -- thermal conductivity
  // move evaluator from PK plist to State
  if (plist_->isSublist("thermal conductivity evaluator")) {
    auto& tcm_plist = S_->GetEvaluatorList(conductivity_key_);
    tcm_plist.setParameters(plist_->sublist("thermal conductivity evaluator"));
    tcm_plist.set("evaluator type", "surface thermal conductivity");
  }
  EnergyBase::SetupPhysicalEvaluators_();

  // -- coupling to subsurface
  coupled_to_subsurface_via_temp_ =
    plist_->get<bool>("coupled to subsurface via temperature", false);
  coupled_to_subsurface_via_flux_ = plist_->get<bool>("coupled to subsurface via flux", false);
  AMANZI_ASSERT(!(coupled_to_subsurface_via_flux_ && coupled_to_subsurface_via_temp_));

  if (coupled_to_subsurface_via_temp_ || coupled_to_subsurface_via_flux_) {
    // -- ensure mass source from subsurface exists
    Key key_ss = Keys::getKey(domain_, "surface_subsurface_flux");

    S_->Require<CompositeVector, CompositeVectorSpace>(key_ss, tag_next_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    // -- ensure enthalpy exists at the new time
    requireAtNext(enthalpy_key_, tag_next_, *S_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    // -- and on the subsurface
    requireAtNext(Keys::getKey(domain_ss_, "enthalpy"), tag_next_, *S_)
      .SetMesh(S_->GetMesh(domain_ss_))
      ->SetGhosted()
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  if (coupled_to_subsurface_via_temp_) {
    // -- energy source term from subsurface
    S_->Require<CompositeVector, CompositeVectorSpace>("surface_subsurface_energy_flux", tag_next_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->RequireEvaluator("surface_subsurface_energy_flux", tag_next_);
  }
}


// -------------------------------------------------------------
// Initialize the needed models to plug in enthalpy.
// -------------------------------------------------------------
void
EnergySurfaceIce::Initialize()
{
  // -- set the cell initial condition if it is taken from the subsurface
  if (!S_->GetRecord(key_, tag_next_).initialized()) {
    if (!plist_->isSublist("initial condition")) {
      Errors::Message message;
      message << name_ << " has no initial condition parameter list.";
      Exceptions::amanzi_throw(message);
    }
  }

  // Call the base class's initialize.
  EnergyBase::Initialize();

  // Set the cell initial condition if it is taken from the subsurface
  if (!S_->GetRecord(key_, tag_next_).initialized()) {
    // TODO: make this go away!  This should be in an MPC?
    Teuchos::ParameterList& ic_plist = plist_->sublist("initial condition");
    if (ic_plist.get<bool>("initialize surface temperature from subsurface", false)) {
      Teuchos::RCP<CompositeVector> surf_temp_cv =
        S_->GetPtrW<CompositeVector>(key_, tag_next_, name_);
      Epetra_MultiVector& surf_temp = *surf_temp_cv->ViewComponent("cell", false);

      Key key_ss = Keys::readKey(*plist_, domain_ss_, "subsurface temperature", "temperature");

      Teuchos::RCP<const CompositeVector> subsurf_temp =
        S_->GetPtr<CompositeVector>(key_ss, tag_next_);
      auto ncells_surface = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

      if (subsurf_temp->HasComponent("face")) {
        const Epetra_MultiVector& temp = *subsurf_temp->ViewComponent("face", false);
        for (unsigned int c = 0; c != ncells_surface; ++c) {
          // -- get the surface cell's equivalent subsurface face and neighboring cell
          AmanziMesh::Entity_ID f = mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
          surf_temp[0][c] = temp[0][f];
        }
      } else if (subsurf_temp->HasComponent("boundary_face")) {
        const Epetra_MultiVector& temp = *subsurf_temp->ViewComponent("boundary_face", false);
        Teuchos::RCP<const AmanziMesh::Mesh> mesh_domain = S_->GetMesh("domain");

        for (unsigned int c = 0; c != ncells_surface; ++c) {
          // -- get the surface cell's equivalent subsurface face and neighboring cell
          AmanziMesh::Entity_ID f = mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
          int bf = mesh_domain->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false).LID(mesh_domain->getMap(AmanziMesh::Entity_kind::FACE,false).GID(f));
          if (bf >= 0) surf_temp[0][c] = temp[0][bf];
        }
      }

      // -- Update faces from cells if there
      DeriveFaceValuesFromCellValues(*surf_temp_cv);

      // mark as initialized
      S_->GetRecordW(key_, tag_next_, name_).set_initialized();

    } else if (ic_plist.get<bool>("initialize surface_star temperature from surface cells",
                                  false)) {
      AMANZI_ASSERT(domain_ == "surface_star");
      Teuchos::RCP<CompositeVector> surf_temp_cv =
        S_->GetPtrW<CompositeVector>(key_, tag_next_, name_);
      Epetra_MultiVector& surf_temp = *surf_temp_cv->ViewComponent("cell", false);

      auto ncells_surface = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
      for (unsigned int c = 0; c != ncells_surface; ++c) {
        int id = mesh_->getMap(AmanziMesh::Entity_kind::CELL,false).GID(c);
        std::stringstream name;
        name << "surface_column_" << id;
        const Epetra_MultiVector& temp =
          *S_->Get<CompositeVector>(Keys::getKey(name.str(), "temperature"), tag_next_)
             .ViewComponent("cell", false);
        surf_temp[0][c] = temp[0][0];
      }
      S_->GetRecordW(key_, tag_next_, name_).set_initialized();
    }
  }
}


// -------------------------------------------------------------
// Deal with the many source terms.
// -------------------------------------------------------------
void
EnergySurfaceIce::AddSources_(const Tag& tag, const Teuchos::Ptr<CompositeVector>& g)
{
  // General source term (covers advection from mass precip and air-surface
  // conduction).
  EnergyBase::AddSources_(tag, g);

  Teuchos::OSTab tab = vo_->getOSTab();
  Epetra_MultiVector& g_c = *g->ViewComponent("cell", false);

  // coupling to subsurface
  // -- two parts -- conduction and advection
  // -- advection source
  if (coupled_to_subsurface_via_temp_ || coupled_to_subsurface_via_flux_) {
    S_->GetEvaluator(Keys::getKey(domain_ss_, "enthalpy"), tag).Update(*S_, name_);
    S_->GetEvaluator(enthalpy_key_, tag).Update(*S_, name_);

    // -- advection source
    Key key_ss = Keys::getKey(domain_, "surface_subsurface_flux");

    const Epetra_MultiVector& source1 =
      *S_->Get<CompositeVector>(key_ss, tag).ViewComponent("cell", false);
    const Epetra_MultiVector& enth_surf =
      *S_->Get<CompositeVector>(enthalpy_key_, tag).ViewComponent("cell", false);
    const Epetra_MultiVector& enth_subsurf =
      *S_->Get<CompositeVector>(Keys::getKey(domain_ss_, "enthalpy"), tag)
         .ViewComponent("cell", false);

    // not needed?
    const Epetra_MultiVector& pd =
      *S_->Get<CompositeVector>(Keys::getKey(domain_, "ponded_depth"), tag)
         .ViewComponent("cell", false);

    unsigned int ncells = g_c.MyLength();
    const auto& mesh_ss = *S_->GetMesh(domain_ss_);
    for (unsigned int c = 0; c != ncells; ++c) {
      double flux = source1[0][c]; // NOTE: this flux is in mol/s

      // upwind the enthalpy
      if (flux > 0.) { // exfiltration
        // get the subsurface's enthalpy
        AmanziMesh::Entity_ID f = mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
        auto cells = mesh_ss.getFaceCells(f, AmanziMesh::Parallel_kind::ALL);

        AMANZI_ASSERT(cells.size() == 1);
        g_c[0][c] -= flux * enth_subsurf[0][cells[0]];
        //        std::cout << "source = " << flux << " * " << enth_subsurf[0][cells[0]] << " = " << -flux * enth_subsurf[0][cells[0]] << std::endl;
        //        std::cout << "OR source = " << flux << " * " << enth_surf[0][c] << " = " << -flux * enth_surf[0][c] << std::endl;
      } else { // infiltration
        g_c[0][c] -= flux * enth_surf[0][c];
      }
    }
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Adding advection to subsurface" << std::endl;
    }
    db_->WriteVector("res (s-s adv src)", g, false);
  }

  // -- conduction source
  if (coupled_to_subsurface_via_temp_) {
    const Epetra_MultiVector& e_source1 =
      *S_->Get<CompositeVector>("surface_subsurface_energy_flux", tag).ViewComponent("cell", false);

    AmanziMesh::Entity_ID_List cells;

    unsigned int ncells = g_c.MyLength();
    for (unsigned int c = 0; c != ncells; ++c) { g_c[0][c] -= e_source1[0][cells[0]]; }
    db_->WriteVector("res (s-s adv src)", g, false);
  }
}


void
EnergySurfaceIce::AddSourcesToPrecon_(double h)
{
  // Deals with nonlinear source terms that are implemented correctly as an evaluator
  EnergyBase::AddSourcesToPrecon_(h);

  // Additionally deal with nonlinear source terms that are NOT
  // implemented correctly, as they are part of a PK (surface energy
  // balance!)
  if (is_source_term_ &&
      S_->HasEvaluator(Keys::getKey(domain_, "conducted_energy_source"), tag_next_) &&
      !S_->GetEvaluator(Keys::getKey(domain_, "conducted_energy_source"), tag_next_)
         .IsDependency(*S_, key_, tag_next_) &&
      S_->HasRecordSet(Keys::getDerivKey(Keys::getKey(domain_, "conducted_energy_source"),
                                         Keys::getKey(domain_, "temperature")))) {
    // This checks if 1, there is a source, and, 2, there is a
    // conducted component to that source, and 4, someone, somewhere
    // (i.e. SEB PK) has defined a dsource_dT, but 3, the source
    // evaluator does not think it depends upon T (because it is
    // hacked in by the PK).
    CompositeVector acc(S_->GetPtrW<CompositeVector>(
                            Keys::getKey(domain_, "conducted_energy_source"), tag_next_, name_)
                          ->Map());
    Epetra_MultiVector& acc_c = *acc.ViewComponent("cell", false);

    Epetra_MultiVector& dsource_dT =
      *S_->GetPtrW<CompositeVector>(
           Keys::getDerivKey(Keys::getKey(domain_, "conducted_energy_source"),
                             Keys::getKey(domain_, "temperature")),
           tag_next_,
           name_)
         ->ViewComponent("cell", false);
    const Epetra_MultiVector& cell_vol =
      *S_->Get<CompositeVector>(cell_vol_key_, tag_next_).ViewComponent("cell", false);
    unsigned int ncells = dsource_dT.MyLength();
    for (unsigned int c = 0; c != ncells; ++c) { acc_c[0][c] = -dsource_dT[0][c] * cell_vol[0][c]; }
    preconditioner_acc_->AddAccumulationTerm(acc, "cell");

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Adding hacked source to PC:" << std::endl;
      db_->WriteVector("de_src_dT",
                       S_->GetPtrW<CompositeVector>(
                           Keys::getDerivKey(Keys::getKey(domain_, "conducted_energy_source"),
                                             Keys::getKey(domain_, "temperature")),
                           tag_next_,
                           name_)
                         .ptr(),
                       false);
    }
  }
}


} // namespace Energy
} // namespace Amanzi
