/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Process kernel for energy equation for overland flow.
------------------------------------------------------------------------- */
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Debugger.hh"
#include "thermal_conductivity_surface_evaluator.hh"
#include "enthalpy_evaluator.hh"
#include "energy_bc_factory.hh"
#include "Function.hh"
#include "FunctionFactory.hh"
#include "independent_variable_field_evaluator.hh"
#include "Op.hh"

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
                                   const Teuchos::RCP<TreeVector>& solution) :
    PK(FElist, plist, S,  solution),
    EnergyBase(FElist, plist, S,  solution),
    standalone_mode_(false),
    is_energy_source_term_(false),
    is_water_source_term_(false),
    is_air_conductivity_(false)
{
  if(!plist_->isParameter("conserved quantity key suffix"))
    plist_->set("conserved quantity key suffix", "energy");
}



// -------------------------------------------------------------
// Create the physical evaluators for energy, enthalpy, thermal
// conductivity, and any sources.
// -------------------------------------------------------------
void EnergySurfaceIce::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {

  standalone_mode_ = S->GetMesh() == S->GetMesh(domain_);

  // Get data and evaluators needed by the PK
  // -- energy, the conserved quantity
  S->RequireField(conserved_key_)->SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(conserved_key_);

  // -- thermal conductivity
  S->RequireField(conductivity_key_)->SetMesh(mesh_)
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList& tcm_plist =
    plist_->sublist("thermal conductivity evaluator");
  tcm_plist.set("evaluator name", conductivity_key_);
  Teuchos::RCP<Energy::ThermalConductivitySurfaceEvaluator> tcm =
    Teuchos::rcp(new Energy::ThermalConductivitySurfaceEvaluator(tcm_plist));
  S->SetFieldEvaluator(conductivity_key_, tcm);

  // -- coupling to subsurface
  coupled_to_subsurface_via_temp_ =
      plist_->get<bool>("coupled to subsurface via temperature", false);
  coupled_to_subsurface_via_flux_ =
      plist_->get<bool>("coupled to subsurface via flux", false);
  AMANZI_ASSERT(! (coupled_to_subsurface_via_flux_ && coupled_to_subsurface_via_temp_));

  if (coupled_to_subsurface_via_temp_ || coupled_to_subsurface_via_flux_ ) {
    // -- ensure mass source from subsurface exists
    Key key_ss = Keys::getKey(domain_,"surface_subsurface_flux");

    S->RequireField(key_ss)
        ->SetMesh(mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);
  }

  if (coupled_to_subsurface_via_temp_) {
    S->RequireFieldEvaluator("surface_subsurface_energy_flux");
    // -- energy source term from subsurface
    S->RequireField("surface_subsurface_energy_flux")
        ->SetMesh(mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);
  }
}


// -------------------------------------------------------------
// Initialize the needed models to plug in enthalpy.
// -------------------------------------------------------------
void EnergySurfaceIce::Initialize(const Teuchos::Ptr<State>& S) {
  // -- set the cell initial condition if it is taken from the subsurface
  if (!plist_->isSublist("initial condition")) {
    std::stringstream messagestream;
    messagestream << name_ << " has no initial condition parameter list.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // Call the base class's initialize.
  EnergyBase::Initialize(S);

  // Set the cell initial condition if it is taken from the subsurface
  if (!S->GetField(key_)->initialized()) {
    // TODO: make this go away!  This should be in an MPC?
    Teuchos::ParameterList& ic_plist = plist_->sublist("initial condition");
    if (ic_plist.get<bool>("initialize surface temperature from subsurface",false)) {
      Teuchos::RCP<CompositeVector> surf_temp_cv = S->GetFieldData(key_, name_);
      Epetra_MultiVector& surf_temp = *surf_temp_cv->ViewComponent("cell",false);


      Key domain_ss = Keys::readDomainHint(*plist_, domain_, "surface", "subsurface");
      Key key_ss = Keys::readKey(*plist_, domain_ss, "subsurface temperature", "temperature");

      Teuchos::RCP<const CompositeVector> subsurf_temp = S->GetFieldData(key_ss);
      unsigned int ncells_surface = mesh_->num_entities(AmanziMesh::CELL,AmanziMesh::Parallel_type::OWNED);

      if (subsurf_temp->HasComponent("face")) {
        const Epetra_MultiVector& temp = *subsurf_temp->ViewComponent("face",false);
        for (unsigned int c=0; c!=ncells_surface; ++c) {
          // -- get the surface cell's equivalent subsurface face and neighboring cell
          AmanziMesh::Entity_ID f =
            mesh_->entity_get_parent(AmanziMesh::CELL, c);
          surf_temp[0][c] = temp[0][f];
        }
      } else if (subsurf_temp->HasComponent("boundary_face")) {
        const Epetra_MultiVector& temp = *subsurf_temp->ViewComponent("boundary_face",false);
        Teuchos::RCP<const AmanziMesh::Mesh> mesh_domain = S->GetMesh("domain");
        unsigned int ncells_sub = mesh_domain->num_entities(AmanziMesh::CELL,AmanziMesh::Parallel_type::OWNED);

        for (unsigned int c=0; c!=ncells_surface; ++c) {
          // -- get the surface cell's equivalent subsurface face and neighboring cell
          AmanziMesh::Entity_ID f = mesh_->entity_get_parent(AmanziMesh::CELL, c);
          int bf = mesh_domain->exterior_face_map(false).LID(mesh_domain->face_map(false).GID(f));
          if (bf >=0) surf_temp[0][c] = temp[0][bf];
        }
      }

      // -- Update faces from cells if there
      DeriveFaceValuesFromCellValues(*surf_temp_cv);

      // mark as initialized
      S->GetField(key_,name_)->set_initialized();

    } else if (ic_plist.get<bool>("initialize surface_star temperature from surface cells",false)) {

      AMANZI_ASSERT(domain_ == "surface_star");
      Epetra_MultiVector& surf_temp = *S->GetFieldData(key_, name_)->ViewComponent("cell",false);

      unsigned int ncells_surface = mesh_->num_entities(AmanziMesh::CELL,AmanziMesh::Parallel_type::OWNED);

      for (unsigned int c=0; c!=ncells_surface; ++c) {
        int id = mesh_->cell_map(false).GID(c);
        std::stringstream name;
        name << "surface_column_" << id;
        const Epetra_MultiVector& temp = *S->GetFieldData(Keys::getKey(name.str(),"temperature"))->ViewComponent("cell",false);
        surf_temp[0][c] = temp[0][0];
      }
      S->GetField(key_,name_)->set_initialized();
    }
  }

}


// -------------------------------------------------------------
// Deal with the many source terms.
// -------------------------------------------------------------
void EnergySurfaceIce::AddSources_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {
  // General source term (covers advection from mass precip and air-surface
  // conduction).
  EnergyBase::AddSources_(S,g);

  Teuchos::OSTab tab = vo_->getOSTab();
  Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);

  // coupling to subsurface
  // -- two parts -- conduction and advection
  // -- advection source
  if (coupled_to_subsurface_via_temp_ || coupled_to_subsurface_via_flux_) {

    // FIXME -- remove this in favor of adding Keys::readDomainHint calls!
    Key domain_ss;
    if (Keys::starts_with(domain_, "surface")  && domain_.find("column") != std::string::npos) {
      domain_ss = plist_->get<std::string>("subsurface domain name",
              domain_.substr(8,domain_.size()));
    } else {
      domain_ss = plist_->get<std::string>("subsurface domain name", "domain");
    }

    S->GetFieldEvaluator(Keys::getKey(domain_ss,"enthalpy"))->HasFieldChanged(S.ptr(), name_);
    S->GetFieldEvaluator(enthalpy_key_)->HasFieldChanged(S.ptr(), name_);

    // -- advection source
    Key key_ss = Keys::getKey(domain_,"surface_subsurface_flux");

    const Epetra_MultiVector& source1 =
      *S->GetFieldData(key_ss)->ViewComponent("cell",false);
    const Epetra_MultiVector& enth_surf =
      *S->GetFieldData(enthalpy_key_)->ViewComponent("cell",false);
    const Epetra_MultiVector& enth_subsurf =
      *S->GetFieldData(Keys::getKey(domain_ss,"enthalpy"))->ViewComponent("cell",false);

    const Epetra_MultiVector& pd =
      *S->GetFieldData(Keys::getKey(domain_,"ponded_depth"))->ViewComponent("cell",false);

    AmanziMesh::Entity_ID_List cells;

    unsigned int ncells = g_c.MyLength();

    for (unsigned int c=0; c!=ncells; ++c) {
      double flux = source1[0][c]; // NOTE: this flux is in mol/s

      // upwind the enthalpy
      if (flux > 0.) { // exfiltration
        // get the subsurface's enthalpy
        AmanziMesh::Entity_ID f = mesh_->entity_get_parent(AmanziMesh::CELL, c);
        S->GetMesh(domain_ss)->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

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
    db_->WriteVector("res (src post surf-subsurf adv)", g, false);
  }

  // -- conduction source
  if (coupled_to_subsurface_via_temp_) {
    const Epetra_MultiVector& e_source1 =
        *S->GetFieldData("surface_subsurface_energy_flux")
        ->ViewComponent("cell",false);

    AmanziMesh::Entity_ID_List cells;

    unsigned int ncells = g_c.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      g_c[0][c] -= e_source1[0][cells[0]];
    }
    db_->WriteVector("res (src post surf-subsurf diff)", g, false);
  }

}


void EnergySurfaceIce::AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h) {
  // Deals with nonlinear source terms that are implemented correctly as an evaluator
  EnergyBase::AddSourcesToPrecon_(S,h);

  // Additionally deal with nonlinear source terms that are NOT
  // implemented correctly, as they are part of a PK (surface energy
  // balance!)
  if (is_source_term_ &&
      S->HasFieldEvaluator(Keys::getKey(domain_,"conducted_energy_source")) &&
      !S->GetFieldEvaluator(Keys::getKey(domain_,"conducted_energy_source"))->IsDependency(S, key_) &&
      S->HasField(Keys::getDerivKey(Keys::getKey(domain_,"conducted_energy_source"), Keys::getKey(domain_,"temperature")))) {

    // This checks if 1, there is a source, and, 2, there is a
    // conducted component to that source, and 4, someone, somewhere
    // (i.e. SEB PK) has defined a dsource_dT, but 3, the source
    // evaluator does not think it depends upon T (because it is
    // hacked in by the PK).
    CompositeVector acc(S->GetFieldData(Keys::getKey(domain_,"conducted_energy_source"))->Map());
    Epetra_MultiVector& acc_c = *acc.ViewComponent("cell", false);

    const Epetra_MultiVector& dsource_dT =
      *S->GetFieldData(Keys::getDerivKey(Keys::getKey(domain_,"conducted_energy_source"), Keys::getKey(domain_,"temperature")))->ViewComponent("cell",false);
    const Epetra_MultiVector& cell_vol = *S->GetFieldData(cell_vol_key_)->ViewComponent("cell",false);
    unsigned int ncells = dsource_dT.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      acc_c[0][c] = -dsource_dT[0][c] * cell_vol[0][c];
    }
    preconditioner_acc_->AddAccumulationTerm(acc, "cell");

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Adding hacked source to PC:" << std::endl;
      db_->WriteVector("de_src_dT", S->GetFieldData(Keys::getDerivKey(Keys::getKey(domain_,"conducted_energy_source"), Keys::getKey(domain_,"temperature"))).ptr(), false);
    }

  }
}


} // namespace
} // namespace
