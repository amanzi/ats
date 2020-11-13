/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Author: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "EpetraExt_MultiVectorOut.h"
#include "Epetra_MultiVector.h"

#include "flow_bc_factory.hh"
#include "Mesh.hh"
#include "Point.hh"
#include "Op.hh"

#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"
#include "independent_variable_field_evaluator.hh"
#include "primary_variable_field_evaluator.hh"


#include "upwind_potential_difference.hh"
#include "upwind_cell_centered.hh"
#include "upwind_total_flux.hh"
#include "pres_elev_evaluator.hh"
#include "elevation_evaluator.hh"
#include "meshed_elevation_evaluator.hh"
#include "elevation_evaluator_column.hh"
#include "standalone_elevation_evaluator.hh"
#include "overland_conductivity_evaluator.hh"
#include "overland_conductivity_model.hh"
#include "overland_pressure_water_content_evaluator.hh"
#include "height_model.hh"
#include "height_evaluator.hh"
//#include "overland_source_from_subsurface_flux_evaluator.hh"

#include "UpwindFluxFactory.hh"

#include "PDE_DiffusionFactory.hh"

#include "overland_pressure.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 0
#define DEBUG_RES_FLAG 0

OverlandPressureFlow::OverlandPressureFlow(Teuchos::ParameterList& pk_tree,
                                           const Teuchos::RCP<Teuchos::ParameterList>& plist,
                                           const Teuchos::RCP<State>& S,                        
                                           const Teuchos::RCP<TreeVector>& solution) :   
    PK_PhysicalBDF_Default(pk_tree, plist, S, solution),
    PK(pk_tree, plist, S, solution),
    standalone_mode_(false),
    is_source_term_(false),
    coupled_to_subsurface_via_head_(false),
    coupled_to_subsurface_via_flux_(false),
    perm_update_required_(true),
    update_flux_(UPDATE_FLUX_ITERATION),
    niter_(0),
    source_only_if_unfrozen_(false),
    precon_used_(true),
    jacobian_(false),
    jacobian_lag_(0),
    iter_(0),
    iter_counter_time_(0.)
{
  if(!plist_->isParameter("conserved quantity key suffix"))
    plist_->set("conserved quantity key suffix", "water_content");

  // set a default absolute tolerance
  if (!plist_->isParameter("absolute error tolerance"))
    plist_->set("absolute error tolerance", 0.01 * 55000.0); // h * nl
}


// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
void OverlandPressureFlow::Setup(const Teuchos::Ptr<State>& S) {

  // set up the meshes
  standalone_mode_ = S->GetMesh() == S->GetMesh(domain_);

  PK_PhysicalBDF_Default::Setup(S);
  
  // -- water content
  S->RequireField(Keys::getKey(domain_,"water_content"))->SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(Keys::getKey(domain_,"water_content"));


  // add _bar evaluators
  Teuchos::ParameterList pd_bar_list = S->FEList().sublist(Keys::getKey(domain_, "ponded_depth"));
  pd_bar_list.set("allow negative ponded depth", true);
  pd_bar_list.setName(Keys::getKey(domain_, "ponded_depth_bar"));
  S->FEList().set(Keys::getKey(domain_, "ponded_depth_bar"), pd_bar_list);

  Teuchos::ParameterList wc_bar_list = S->FEList().sublist(Keys::getKey(domain_, "water_content"));
  wc_bar_list.set("allow negative water content", true);
  wc_bar_list.setName(Keys::getKey(domain_, "water_content_bar"));
  S->FEList().set(Keys::getKey(domain_, "water_content_bar"), wc_bar_list);
  
  // this pk uses density
  S->RequireField(Keys::getKey(domain_,"molar_density_liquid"))->SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(Keys::getKey(domain_,"molar_density_liquid"));

  SetupOverlandFlow_(S);
  SetupPhysicalEvaluators_(S);

}


void OverlandPressureFlow::SetupOverlandFlow_(const Teuchos::Ptr<State>& S) {
  // -- cell volume and evaluator
  S->RequireFieldEvaluator(Keys::getKey(domain_,"cell_volume"));

  S->RequireGravity();
  S->RequireScalar("atmospheric_pressure");

  // Set up Operators
  // -- boundary conditions
  Teuchos::ParameterList bc_plist = plist_->sublist("boundary conditions", true);
  FlowBCFactory bc_factory(mesh_, bc_plist);
  bc_head_ = bc_factory.CreateHead();
  bc_pressure_ = bc_factory.CreatePressure();
  bc_zero_gradient_ = bc_factory.CreateZeroGradient();
  bc_flux_ = bc_factory.CreateMassFlux();
  bc_level_ = bc_factory.CreateFixedLevel();
  bc_seepage_head_ = bc_factory.CreateSeepageFaceHead();
  bc_seepage_pressure_ = bc_factory.CreateSeepageFacePressure();
  bc_critical_depth_ = bc_factory.CreateCriticalDepth();

  bc_dynamic_ = bc_factory.CreateDynamic();
  bc_tidal_ = bc_factory.CreateTidalHead();

  bc_level_flux_lvl_ = bc_factory.CreateFixedLevelFlux_Level();
  bc_level_flux_vel_ = bc_factory.CreateFixedLevelFlux_Velocity();


  
  if (bc_plist.isParameter("seepage face")) {
    // old style! DEPRECATED
    Errors::Message message;
    message << name_ << ": DEPRECATION: \"seepage face\" boundary condition parameter list names have changed to \"seepage face head\" and \"seepage face pressure\".";
    Exceptions::amanzi_throw(message);
  }
  
  // -- nonlinear coefficients/upwinding
  Teuchos::ParameterList cond_plist;// = plist_->sublist("overland conductivity evaluator");
  Operators::UpwindFluxFactory upwfactory;

  upwinding_ = upwfactory.Create(cond_plist, name_,
       Keys::getKey(domain_,"overland_conductivity"), Keys::getKey(domain_,"upwind_overland_conductivity"),
                                 Keys::getKey(domain_,"mass_flux_direction"));

  // -- require the data on appropriate locations
  std::string coef_location = upwinding_->CoefficientLocation();
  if (coef_location == "upwind: face") {  

    S->RequireField(Keys::getKey(domain_,"upwind_overland_conductivity"), name_)->SetMesh(mesh_)
      ->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
  } else if (coef_location == "standard: cell") {
    S->RequireField(Keys::getKey(domain_,"upwind_overland_conductivity"), name_)->SetMesh(mesh_)
      ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
  } else {
    Errors::Message message;
    message << name_ << ": Unknown upwind coefficient location in overland flow.";
    Exceptions::amanzi_throw(message);
  }

  S->GetField(Keys::getKey(domain_,"upwind_overland_conductivity"),name_)->set_io_vis(false);

  // -- create the forward operator for the diffusion term
  // DEPRECATED OPTIONS
  if (plist_->isParameter("Diffusion") ||
      plist_->isParameter("Diffusion PC")) {
    Errors::Message message("Richards PK: DEPRECATION: Discretization lists \"Diffusion\" and \"Diffusion PC\" have been renamed \"diffusion\" and \"diffusion preconditioner\", respectively.");
    Exceptions::amanzi_throw(message);
  }

  Teuchos::ParameterList& mfd_plist = plist_->sublist("diffusion");
  mfd_plist.set("nonlinear coefficient", coef_location);
  if (!mfd_plist.isParameter("scaled constraint equation"))
    mfd_plist.set("scaled constraint equation", true);
  if (mfd_plist.isParameter("Newton correction")) {
    Errors::Message message;
    message << name_ << ": The forward operator for Diffusion should not set a "
            << "\"Newton correction\" term, perhaps you meant to put this in a "
            << "\"Diffusion PC\" sublist.";
    Exceptions::amanzi_throw(message);
  }    
  
  Operators::PDE_DiffusionFactory opfactory;
  matrix_diff_ = opfactory.Create(mfd_plist, mesh_, bc_);
  matrix_diff_->SetTensorCoefficient(Teuchos::null);
  matrix_ = matrix_diff_->global_operator();
  
  // -- create the operator, data for flux directions
  Teuchos::ParameterList face_diff_list(mfd_plist);
  face_diff_list.set("nonlinear coefficient", "none");
  face_matrix_diff_ = opfactory.Create(face_diff_list, mesh_, bc_);
  face_matrix_diff_->SetTensorCoefficient(Teuchos::null);
  face_matrix_diff_->SetScalarCoefficient(Teuchos::null, Teuchos::null);
  face_matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);


  S->RequireField(Keys::getKey(domain_,"mass_flux_direction"), name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("face", AmanziMesh::FACE, 1);
  
  // -- create the operators for the preconditioner
  //    diffusion
  Teuchos::ParameterList& mfd_pc_plist = plist_->sublist("diffusion preconditioner");
  mfd_pc_plist.set("nonlinear coefficient", coef_location);
  mfd_pc_plist.set("scaled constraint equation",
                   mfd_plist.get<bool>("scaled constraint equation"));
  mfd_pc_plist.set("constraint equation scaling cutoff",
                   mfd_plist.get<double>("constraint equation scaling cutoff", 1.0));
  if (!mfd_pc_plist.isParameter("discretization primary"))
    mfd_pc_plist.set("discretization primary",
                     mfd_plist.get<std::string>("discretization primary"));
  if (!mfd_pc_plist.isParameter("discretization secondary") &&
      mfd_plist.isParameter("discretization secondary"))
    mfd_pc_plist.set("discretization secondary",
                     mfd_plist.get<std::string>("discretization secondary"));
  if (!mfd_pc_plist.isParameter("schema") && mfd_plist.isParameter("schema"))
    mfd_pc_plist.set("schema",
                     mfd_plist.get<Teuchos::Array<std::string> >("schema"));
  if (mfd_pc_plist.get<bool>("include Newton correction", false)) {
    if (mfd_pc_plist.get<std::string>("discretization primary") == "fv: default") {
      mfd_pc_plist.set("Newton correction", "true Jacobian");
    } else {
      mfd_pc_plist.set("Newton correction", "approximate Jacobian");
    }
  }

  //    this needs to get some better logic now that symbolic doesn't use this. --etc
  precon_used_ = plist_->isSublist("preconditioner") ||
    plist_->isSublist("inverse") ||
    plist_->isSublist("linear solver");

  if (precon_used_) {
    auto& inv_list = mfd_pc_plist.sublist("inverse");
    inv_list.setParameters(plist_->sublist("inverse"));
    // old style... deprecate me!
    inv_list.setParameters(plist_->sublist("preconditioner"));
    inv_list.setParameters(plist_->sublist("linear solver"));
  }

  preconditioner_diff_ = opfactory.Create(mfd_pc_plist, mesh_, bc_);
  preconditioner_diff_->SetTensorCoefficient(Teuchos::null);
  preconditioner_ = preconditioner_diff_->global_operator();

  // If using approximate Jacobian for the preconditioner, we also need derivative information.
  jacobian_ = (mfd_pc_plist.get<std::string>("Newton correction", "none") != "none");
  if (jacobian_) {
    jacobian_lag_ = mfd_pc_plist.get<int>("Newton correction lag", 0);
    
    if (preconditioner_->RangeMap().HasComponent("face")) {
      // MFD -- upwind required
      S->RequireField(Keys::getDerivKey(Keys::getKey(domain_,"upwind_overland_conductivity"),Keys::getKey(domain_,"ponded_depth")), name_)
        ->SetMesh(mesh_)->SetGhosted()
        ->SetComponent("face", AmanziMesh::FACE, 1);

      upwinding_dkdp_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
                                    Keys::getDerivKey(Keys::getKey(domain_,"overland_conductivity"),Keys::getKey(domain_,"ponded_depth")),
                                    Keys::getDerivKey(Keys::getKey(domain_,"upwind_overland_conductivity"),Keys::getKey(domain_,"ponded_depth")),
                                    Keys::getKey(domain_,"mass_flux_direction"),1.e-12));
    }
  }
  
  // -- coupling to subsurface
  coupled_to_subsurface_via_flux_ =
      plist_->get<bool>("coupled to subsurface via flux", false);
  coupled_to_subsurface_via_head_ =
      plist_->get<bool>("coupled to subsurface via head", false);
  AMANZI_ASSERT(!(coupled_to_subsurface_via_flux_ && coupled_to_subsurface_via_head_));

  if (coupled_to_subsurface_via_head_) {
    // -- source term from subsurface, filled in by evaluator,
    //    which picks the fluxes from "mass_flux" field.
    S->RequireFieldEvaluator("surface_subsurface_flux");
    S->RequireField("surface_subsurface_flux")
        ->SetMesh(mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  
  //    accumulation
  Teuchos::ParameterList& acc_pc_plist = plist_->sublist("accumulation preconditioner");
  acc_pc_plist.set("entity kind", "cell");
  preconditioner_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(acc_pc_plist, preconditioner_));
  
  // primary variable
  S->RequireField(key_, name_)->Update(matrix_->RangeMap())->SetGhosted()
      ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S->RequireField(Keys::getKey(domain_,"pressure"))->Update(matrix_->RangeMap())->SetGhosted();

  // fluxes
  S->RequireField(Keys::getKey(domain_,"mass_flux"), name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("face", AmanziMesh::FACE, 1);

  S->RequireField(Keys::getKey(domain_,"velocity"), name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 3);

  // limiters
  p_limit_ = plist_->get<double>("limit correction to pressure change [Pa]", -1.);
  patm_limit_ = plist_->get<double>("limit correction when crossing atmospheric pressure [Pa]", -1.);
  patm_hard_limit_ = plist_->get<bool>("allow no negative ponded depths", false);
  min_vel_ponded_depth_ = plist_->get<double>("min ponded depth for velocity calculation", 1e-2);
  min_tidal_bc_ponded_depth_ = plist_->get<double>("min ponded depth for tidal bc", 0.02);  
  
  
};


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void OverlandPressureFlow::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2, 1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";
  
  // -- evaluator for surface geometry.
  Key elev_key = Keys::readKey(*plist_, domain_, "elevation", "elevation");
  S->RequireField(elev_key)->SetMesh(S->GetMesh(domain_))->SetGhosted()
      ->AddComponents(names2, locations2, num_dofs2);

  Key slope_key = Keys::readKey(*plist_, domain_, "slope magnitude", "slope_magnitude");
  S->RequireField(slope_key)->SetMesh(S->GetMesh(domain_))
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // TODO: move these to state
  Teuchos::RCP<Flow::ElevationEvaluator> elev_evaluator;
  if (S->FEList().isSublist(elev_key)) {
    S->RequireFieldEvaluator(elev_key);
    S->RequireFieldEvaluator(slope_key);
  } else {
    if (standalone_mode_) {
      AMANZI_ASSERT(plist_->isSublist("elevation evaluator"));
      Teuchos::ParameterList elev_plist = plist_->sublist("elevation evaluator");
      elev_plist.set("evaluator name", Keys::getKey(domain_, "elevation"));
      elev_evaluator = Teuchos::rcp(new Flow::StandaloneElevationEvaluator(elev_plist));
    } else {
      Teuchos::ParameterList elev_plist = plist_->sublist("elevation evaluator");
      elev_plist.set("evaluator name", Keys::getKey(domain_, "elevation"));
      elev_evaluator = Teuchos::rcp(new Flow::MeshedElevationEvaluator(elev_plist));
    }

    S->SetFieldEvaluator(elev_key, elev_evaluator);
    S->SetFieldEvaluator(slope_key, elev_evaluator);
  }

  // TODO: move this to state
  // -- evaluator for potential field, h + z
  S->RequireField(Keys::getKey(domain_,"pres_elev"))->Update(matrix_->RangeMap())->SetGhosted();
  Teuchos::ParameterList pres_elev_plist = plist_->sublist("potential evaluator");
  pres_elev_plist.set("evaluator name", Keys::getKey(domain_,"pres_elev"));

  Teuchos::RCP<Flow::PresElevEvaluator> pres_elev_eval =
      Teuchos::rcp(new Flow::PresElevEvaluator(pres_elev_plist));
  S->SetFieldEvaluator(Keys::getKey(domain_,"pres_elev"), pres_elev_eval);


  // -- evaluator for source term
  is_source_term_ = plist_->get<bool>("source term");
  if (is_source_term_) {
    if (source_key_.empty()) {
      source_key_ = Keys::readKey(*plist_, domain_, "source", "mass_source");
    }
    source_in_meters_ = plist_->get<bool>("mass source in meters", true);

    S->RequireField(source_key_)->SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(source_key_);

    if (source_in_meters_){
      // density of incoming water [mol/m^3]
      S->RequireField(Keys::getKey(domain_,"source_molar_density"))->SetMesh(mesh_)
          ->AddComponent("cell", AmanziMesh::CELL, 1);
      S->RequireFieldEvaluator(Keys::getKey(domain_,"source_molar_density"));
    }
  }

  // -- water content bar (can be negative)
  S->RequireField(Keys::getKey(domain_,"water_content_bar"))->SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(Keys::getKey(domain_,"water_content_bar"));

  // -- ponded depth
  S->RequireField(Keys::getKey(domain_,"ponded_depth"))->Update(matrix_->RangeMap())->SetGhosted();
  S->RequireFieldEvaluator(Keys::getKey(domain_,"ponded_depth"));

  // -- ponded depth bar (can be negative)
  S->RequireField(Keys::getKey(domain_,"ponded_depth_bar"))->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(Keys::getKey(domain_,"ponded_depth_bar"));
  
  // -- conductivity evaluator
  S->RequireField(Keys::getKey(domain_,"overland_conductivity"))->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);

  if (plist_->isSublist("overland conductivity evaluator")) {
    Teuchos::ParameterList cond_plist = plist_->sublist("overland conductivity evaluator");
    cond_plist.set("evaluator name", Keys::getKey(domain_, "overland_conductivity"));

    Teuchos::RCP<Flow::OverlandConductivityEvaluator> cond_evaluator =
        Teuchos::rcp(new Flow::OverlandConductivityEvaluator(cond_plist));

    S->SetFieldEvaluator(Keys::getKey(domain_,"overland_conductivity"), cond_evaluator);
  } else {
    S->RequireFieldEvaluator(Keys::getKey(domain_,"overland_conductivity"));
  }
}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void OverlandPressureFlow::Initialize(const Teuchos::Ptr<State>& S) {



#if DEBUG_RES_FLAG
  for (int i=1; i!=23; ++i) {
    std::stringstream namestream;
    namestream << "flow_residual_" << i;
    S->GetFieldData(namestream.str(),name_)->PutScalar(0.);
    S->GetField(namestream.str(),name_)->set_initialized();

    std::stringstream solnstream;
    solnstream << "flow_solution_" << i;
    S->GetFieldData(solnstream.str(),name_)->PutScalar(0.);
    S->GetField(solnstream.str(),name_)->set_initialized();
  }
#endif


  Teuchos::RCP<CompositeVector> pres_cv = S->GetFieldData(key_, name_);

  // initial condition is tricky
  if (!S->GetField(key_)->initialized()) {
    if (!plist_->isSublist("initial condition")) {
      std::stringstream messagestream;
      messagestream << name_ << " has no initial condition parameter list.";
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
    }
    pres_cv->PutScalar(0.);
  }

  // Initialize BDF stuff and physical domain stuff.
  PK_PhysicalBDF_Default::Initialize(S);
 
  if (!S->GetField(key_)->initialized()) {
    // TODO: can this be deprecated?  Shouldn't this get handled by the MPC?    
    // -- set the cell initial condition if it is taken from the subsurface
    Teuchos::ParameterList ic_plist = plist_->sublist("initial condition");
    if (ic_plist.get<bool>("initialize surface head from subsurface",false)) {
      Epetra_MultiVector& pres = *pres_cv->ViewComponent("cell",false);
      Key key_ss;

      if (boost::starts_with(domain_, "surface") && domain_.find("column") != std::string::npos) {
        Key domain_ss;
        if (domain_ == "surface") domain_ss = "domain";
        else domain_ss = domain_.substr(8,domain_.size());
        key_ss = ic_plist.get<std::string>("subsurface pressure key",
                Keys::getKey(domain_ss, "pressure"));
      } else {
        key_ss = ic_plist.get<std::string>("subsurface pressure key", "pressure");
      }
      
    Teuchos::RCP<const CompositeVector> subsurf_pres = S->GetFieldData(key_ss);
    unsigned int ncells_surface = mesh_->num_entities(AmanziMesh::CELL,AmanziMesh::Parallel_type::OWNED);
      if (subsurf_pres->HasComponent("face")) {

        const Epetra_MultiVector& subsurf_pres = *S->GetFieldData(key_ss)
                                                 ->ViewComponent("face",false);
        unsigned int ncells_surface = mesh_->num_entities(AmanziMesh::CELL,AmanziMesh::Parallel_type::OWNED);
        for (unsigned int c=0; c!=ncells_surface; ++c) {
          // -- get the surface cell's equivalent subsurface face and neighboring cell
          AmanziMesh::Entity_ID f =
              mesh_->entity_get_parent(AmanziMesh::CELL, c);
          pres[0][c] = subsurf_pres[0][f];
        }
          
      } else if (subsurf_pres->HasComponent("boundary_face")) {

        const Epetra_MultiVector& subsurf_pres_vec = *subsurf_pres->ViewComponent("boundary_face",false);
        Teuchos::RCP<const AmanziMesh::Mesh> mesh_domain = S->GetMesh("domain");
        unsigned int ncells_sub = mesh_domain->num_entities(AmanziMesh::CELL,AmanziMesh::Parallel_type::OWNED);
          
        for (unsigned int c=0; c!=ncells_surface; ++c) {
          // -- get the surface cell's equivalent subsurface face and neighboring cell
          AmanziMesh::Entity_ID f = mesh_->entity_get_parent(AmanziMesh::CELL, c);
          int bf = mesh_domain->exterior_face_map(false).LID(mesh_domain->face_map(false).GID(f));
          if (bf >=0)   pres[0][c] = subsurf_pres_vec[0][bf];
        }
      }

    }

    // -- Update faces from cells if there
    DeriveFaceValuesFromCellValues_(pres_cv.ptr());

    // mark as initialized
    if (ic_plist.get<bool>("initialize surface head from subsurface",false)) {
      S->GetField(key_,name_)->set_initialized();
    }
    else if (ic_plist.get<bool>("initialize surface_star head from surface cells",false)) {
      // TODO: can't this move into an MPC?
      AMANZI_ASSERT(domain_ == "surface_star");
      Epetra_MultiVector& pres_star = *pres_cv->ViewComponent("cell",false);
    
      unsigned int ncells_surface = mesh_->num_entities(AmanziMesh::CELL,AmanziMesh::Parallel_type::OWNED);
      for (unsigned int c=0; c!=ncells_surface; ++c) {
        int id = mesh_->cell_map(false).GID(c);
        
        std::stringstream name;
        name << "surface_column_"<< id;
        
        const Epetra_MultiVector& pres = *S->GetFieldData(Keys::getKey(name.str(),"pressure"))->ViewComponent("cell",false);
        
        // -- get the surface cell's equivalent subsurface face and neighboring cell
        if (pres[0][0] > 101325.)
          pres_star[0][c] = pres[0][0];
        else
          pres_star[0][c] = 101325.0;
      }
      
      // mark as initialized
      S->GetField(key_,name_)->set_initialized();
    }
    
  }
  
  // Initialize BC values
  bc_head_->Compute(S->time());
  bc_pressure_->Compute(S->time());
  bc_zero_gradient_->Compute(S->time());
  bc_flux_->Compute(S->time());
  bc_level_->Compute(S->time());
  bc_level_flux_lvl_->Compute(S->time());
  bc_level_flux_vel_->Compute(S->time());  
  
  bc_seepage_head_->Compute(S->time());
  bc_seepage_pressure_->Compute(S->time());
  bc_critical_depth_->Compute(S->time());
  bc_dynamic_->Compute(S->time());
  bc_tidal_->Compute(S->time());  
  
  // Set extra fields as initialized -- these don't currently have evaluators.
  S->GetFieldData(Keys::getKey(domain_,"upwind_overland_conductivity"),name_)->PutScalar(0.0);
  S->GetField(Keys::getKey(domain_,"upwind_overland_conductivity"),name_)->set_initialized();

  if (jacobian_ && preconditioner_->RangeMap().HasComponent("face")) {
    S->GetFieldData(Keys::getDerivKey(Keys::getKey(domain_,"upwind_overland_conductivity"),Keys::getKey(domain_,"ponded_depth")),name_)->PutScalar(1.0);
    S->GetField(Keys::getDerivKey(Keys::getKey(domain_,"upwind_overland_conductivity"),Keys::getKey(domain_,"ponded_depth")),name_)->set_initialized();
  }

  S->GetField(Keys::getKey(domain_,"mass_flux"), name_)->set_initialized();
  S->GetFieldData(Keys::getKey(domain_,"mass_flux_direction"), name_)->PutScalar(0.);
  S->GetField(Keys::getKey(domain_,"mass_flux_direction"), name_)->set_initialized();
  S->GetFieldData(Keys::getKey(domain_,"velocity"), name_)->PutScalar(0.);
  S->GetField(Keys::getKey(domain_,"velocity"), name_)->set_initialized();
};



// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void OverlandPressureFlow::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {
  niter_ = 0;

  double dt = t_new - t_old;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Commiting state." << std::endl;

  PK_PhysicalBDF_Default::CommitStep(t_old, t_new, S);

  // update boundary conditions
  bc_head_->Compute(S->time());
  bc_pressure_->Compute(S->time());
  bc_flux_->Compute(S->time());
  bc_level_->Compute(S->time());
  bc_level_flux_lvl_->Compute(S->time());
  bc_level_flux_vel_->Compute(S->time());  
  
  bc_seepage_head_->Compute(S->time());
  bc_seepage_pressure_->Compute(S->time());
  bc_critical_depth_->Compute(S->time());
  bc_dynamic_->Compute(S->time());
  bc_tidal_->Compute(S->time());
  UpdateBoundaryConditions_(S.ptr());

  // Update flux if rel perm or h + Z has changed.
  bool update = UpdatePermeabilityData_(S.ptr());
  update |= S->GetFieldEvaluator(Keys::getKey(domain_,"pres_elev"))->HasFieldChanged(S.ptr(), name_);

  // update the stiffness matrix with the new rel perm
  auto cond = S->GetFieldData(Keys::getKey(domain_,"upwind_overland_conductivity"));

  // update the stiffness matrix
  matrix_->Init();
  matrix_diff_->SetScalarCoefficient(cond, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  FixBCsForOperator_(S.ptr(), matrix_diff_.ptr()); // deals with zero gradient case
  matrix_diff_->ApplyBCs(true, true, true);
  
  // derive the fluxes
    
  Teuchos::RCP<const CompositeVector> potential = S->GetFieldData(Keys::getKey(domain_,"pres_elev"));
  Teuchos::RCP<CompositeVector> flux = S->GetFieldData(Keys::getKey(domain_,"mass_flux"), name_);
  matrix_diff_->UpdateFlux(potential.ptr(), flux.ptr());
  
};


// -----------------------------------------------------------------------------
// Update diagnostics -- used prior to vis.
// -----------------------------------------------------------------------------
void OverlandPressureFlow::CalculateDiagnostics(const Teuchos::RCP<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Calculating diagnostic variables." << std::endl;

  // update the cell velocities
  UpdateBoundaryConditions_(S.ptr());

  Teuchos::RCP<const CompositeVector> conductivity = S->GetFieldData(Keys::getKey(domain_,"upwind_overland_conductivity"));

  // update the stiffness matrix
  matrix_diff_->SetScalarCoefficient(conductivity, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  FixBCsForOperator_(S.ptr(), matrix_diff_.ptr()); // deals with zero gradient case
  matrix_diff_->ApplyBCs(true, true, true);

  // derive fluxes
  Teuchos::RCP<const CompositeVector> potential = S->GetFieldData(Keys::getKey(domain_,"pres_elev"));
  Teuchos::RCP<CompositeVector> flux = S->GetFieldData(Keys::getKey(domain_,"mass_flux"), name_);
  matrix_diff_->UpdateFlux(potential.ptr(), flux.ptr());

  // update velocity
  Epetra_MultiVector& velocity = *S->GetFieldData(Keys::getKey(domain_,"velocity"), name_)
      ->ViewComponent("cell", true);
  flux->ScatterMasterToGhosted("face");
  const Epetra_MultiVector& flux_f = *flux->ViewComponent("face",true);
  const Epetra_MultiVector& nliq_c = *S->GetFieldData(Keys::getKey(domain_,"molar_density_liquid"))
    ->ViewComponent("cell");
  const Epetra_MultiVector& pd_c = *S->GetFieldData(Keys::getKey(domain_,"ponded_depth"))
    ->ViewComponent("cell");
  
  int d(mesh_->space_dimension());
  AmanziGeometry::Point local_velocity(d);

  Teuchos::LAPACK<int, double> lapack;
  Teuchos::SerialDenseMatrix<int, double> matrix(d, d);
  double rhs[d];

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  AmanziMesh::Entity_ID_List faces;
  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int i=0; i!=d; ++i) rhs[i] = 0.0;
    matrix.putScalar(0.0);

    for (int n=0; n!=nfaces; ++n) {  // populate least-square matrix
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double area = mesh_->face_area(f);

      for (int i=0; i!=d; ++i) {
        rhs[i] += normal[i] * flux_f[0][f];
        matrix(i, i) += normal[i] * normal[i];
        for (int j = i+1; j < d; ++j) {
          matrix(j, i) = matrix(i, j) += normal[i] * normal[j];
        }
      }
    }

    int info;
    lapack.POSV('U', d, 1, matrix.values(), d, rhs, d, &info);

    // NOTE this is probably wrong in the frozen case?  pd --> uf*pd?
    for (int i=0; i!=d; ++i) velocity[i][c] = pd_c[0][c] > min_vel_ponded_depth_ ? rhs[i] / (nliq_c[0][c] * pd_c[0][c]) : 0.;
  }


};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool OverlandPressureFlow::UpdatePermeabilityData_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating permeability?";

  bool update_perm = S->GetFieldEvaluator(Keys::getKey(domain_,"ponded_depth"))->HasFieldChanged(S, name_);
  
  // this is an ugly hack to get boundary conditions into conductivities
  Teuchos::RCP<CompositeVector> pd = S->GetFieldData(Keys::getKey(domain_,"ponded_depth"), Keys::getKey(domain_,"ponded_depth"));
  Teuchos::RCP<const CompositeVector> elev = S->GetFieldData(Keys::getKey(domain_,"elevation"));
  ApplyBoundaryConditions_(pd.ptr(), elev.ptr());

  update_perm |= S->GetFieldEvaluator(Keys::getKey(domain_,"pres_elev"))->HasFieldChanged(S, name_);
  update_perm |= S->GetFieldEvaluator(Keys::getKey(domain_,"overland_conductivity"))->HasFieldChanged(S, name_);
  update_perm |= perm_update_required_;

  if (update_perm) {
    // Update the perm only if needed.
    perm_update_required_ = false;
    
    // update the direction of the flux -- note this is NOT the flux, but grad(h+z).
    Teuchos::RCP<CompositeVector> flux_dir = S->GetFieldData(Keys::getKey(domain_,"mass_flux_direction"), name_);
    Teuchos::RCP<const CompositeVector> pres_elev = S->GetFieldData(Keys::getKey(domain_,"pres_elev"));

    // -- note the matrices are constant coefficient, and do not change, but
    //    the boundary local matrices have been overwritten with 0 by the last
    //    call to ApplyBCs().  Recover the shadow matrices.
    face_matrix_diff_->local_op()->CopyShadowToMaster();

    // -- need to apply BCs to get boundary flux directions correct
    FixBCsForOperator_(S.ptr(), face_matrix_diff_.ptr()); // deals with zero gradient condition
    face_matrix_diff_->ApplyBCs(true, true, true);

    // -- now we can calculate the flux direction
    face_matrix_diff_->UpdateFlux(pres_elev.ptr(), flux_dir.ptr());

    // upwind
    // -- get upwind conductivity data
    Teuchos::RCP<CompositeVector> uw_cond = S->GetFieldData(Keys::getKey(domain_,"upwind_overland_conductivity"), name_);

    // -- get conductivity data
    Teuchos::RCP<const CompositeVector> cond = S->GetFieldData(Keys::getKey(domain_,"overland_conductivity"));
    
    // -- Move rel perm on boundary_faces into uw_rel_perm on faces
    {
      const Epetra_Import& vandelay = mesh_->exterior_face_importer();
      const Epetra_MultiVector& cond_bf = *cond->ViewComponent("boundary_face",false);
      Epetra_MultiVector& uw_cond_f = *uw_cond->ViewComponent("face",false);
      uw_cond_f.Export(cond_bf, vandelay, Insert);
    }

    // -- upwind
    upwinding_->Update(S);
    uw_cond->ScatterMasterToGhosted("face");
  }

  if (update_perm && vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << " TRUE." << std::endl;
  return update_perm;
}


// -----------------------------------------------------------------------------
// Derivatives of the overland conductivity, upwinded.
// -----------------------------------------------------------------------------
bool OverlandPressureFlow::UpdatePermeabilityDerivativeData_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating permeability derivatives?";


  bool update_perm = S->GetFieldEvaluator(Keys::getKey(domain_,"overland_conductivity"))
    ->HasFieldDerivativeChanged(S, name_, Keys::getKey(domain_,"ponded_depth"));
  Teuchos::RCP<const CompositeVector> dcond =
    S->GetFieldData(Keys::getDerivKey(Keys::getKey(domain_,"overland_conductivity"), Keys::getKey(domain_,"ponded_depth")));

  if (update_perm) {
    if (preconditioner_->RangeMap().HasComponent("face")) {
      // get upwind conductivity data
      Teuchos::RCP<CompositeVector> duw_cond =
        S->GetFieldData(Keys::getDerivKey(Keys::getKey(domain_,"upwind_overland_conductivity"),Keys::getKey(domain_,"ponded_depth")), name_);

      duw_cond->PutScalar(0.);
    
      // Then upwind.  This overwrites the boundary if upwinding says so.
      upwinding_dkdp_->Update(S);
      duw_cond->ScatterMasterToGhosted("face");
    } else {
      dcond->ScatterMasterToGhosted("cell");
    }
  }

  if (update_perm && vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << " TRUE." << std::endl;
  return update_perm;
}


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void OverlandPressureFlow::UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating BCs." << std::endl;

  auto& markers = bc_markers();
  auto& values = bc_values();

  AmanziMesh::Entity_ID_List cells;

  const Epetra_MultiVector& elevation = *S->GetFieldData(Keys::getKey(domain_,"elevation"))
      ->ViewComponent("face",false);

  // initialize all as null
  for (unsigned int n=0; n!=markers.size(); ++n) {
    markers[n] = Operators::OPERATOR_BC_NONE;
    values[n] = 0.0;
  }

  // Head BCs are standard Dirichlet, plus elevation
  for (Functions::BoundaryFunction::Iterator bc=bc_head_->begin();
       bc!=bc_head_->end(); ++bc) {
    int f = bc->first;
    markers[f] = Operators::OPERATOR_BC_DIRICHLET;
    values[f] = bc->second + elevation[0][f];
  }

  // Pressure BCs require a change in coordinates from pressure to head
  if (bc_pressure_->size() > 0) {

    S->GetFieldEvaluator(Keys::getKey(domain_,"ponded_depth"))->HasFieldChanged(S.ptr(), name_);

    const Epetra_MultiVector& h_cells = *S->GetFieldData(Keys::getKey(domain_,"ponded_depth"))->ViewComponent("cell");
    const Epetra_MultiVector& elevation_cells = *S->GetFieldData(Keys::getKey(domain_,"elevation"))->ViewComponent("cell");
    const Epetra_MultiVector& rho_l = *S->GetFieldData(Keys::getKey(domain_,"mass_density_liquid"))->ViewComponent("cell");
    double gz = -(*S->GetConstantVectorData("gravity"))[2];
    const double& p_atm = *S->GetScalarData("atmospheric_pressure");

    if (S->HasFieldEvaluator(Keys::getKey(domain_,"mass_density_ice"))) {
      // thermal model of height
      const Epetra_MultiVector& eta = *S->GetFieldData(Keys::getKey(domain_,"unfrozen_fraction"))->ViewComponent("cell");
      const Epetra_MultiVector& rho_i = *S->GetFieldData(Keys::getKey(domain_,"mass_density_ice"))->ViewComponent("cell");

      for (Functions::BoundaryFunction::Iterator bc = bc_pressure_->begin(); 
           bc != bc_pressure_->end(); ++bc) {
        int f = bc->first;
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        int c = cells[0];

        double p0 = bc->second > p_atm ? bc->second : p_atm;
        double h0 = (p0 - p_atm) / ((eta[0][c]*rho_l[0][c] + (1.-eta[0][c])*rho_i[0][c]) * gz);

        markers[f] = Operators::OPERATOR_BC_DIRICHLET;
        values[f] = h0 + elevation[0][f];
      }

    } else {
      // non-thermal model
      for (Functions::BoundaryFunction::Iterator bc = bc_pressure_->begin(); 
           bc != bc_pressure_->end(); ++bc) {
        int f = bc->first;
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        int c = cells[0];

        double p0 = bc->second > p_atm ? bc->second : p_atm;
        double h0 = (p0 - p_atm) / (rho_l[0][c] * gz);

        markers[f] = Operators::OPERATOR_BC_DIRICHLET;
        values[f] = h0 + elevation[0][f];
      }
    }
  }
  
  // Head BCs for fixed water level
  for (Functions::BoundaryFunction::Iterator bc=bc_level_->begin();
       bc!=bc_level_->end(); ++bc) {
    int f = bc->first;
    markers[f] = Operators::OPERATOR_BC_DIRICHLET;
    double val = bc->second;

    if (elevation[0][f] > val) values[f] = 0;
    else values[f] = val;
  }

  //for (auto bc : bc_dynamic_) {
  if (bc_dynamic_ -> size() > 0){
    double time = S->time();
    int id = bc_dynamic_->Func_ID(time);
    for (Functions::DynamicBoundaryFunction::Iterator bc=bc_dynamic_->GetFunction(id) -> begin();
         bc!=bc_dynamic_-> GetFunction(id) -> end(); ++bc){
      int f = bc->first;
      if (id == 0){
        markers[f] = Operators::OPERATOR_BC_DIRICHLET;
        values[f] = bc->second + elevation[0][f];
      }else if (id==1){
        markers[f] = Operators::OPERATOR_BC_NEUMANN;
        values[f] = bc->second;
      }
    }
  }

  // Standard Neumann data for flux
  for (Functions::BoundaryFunction::Iterator bc=bc_flux_->begin();
       bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    markers[f] = Operators::OPERATOR_BC_NEUMANN;
    values[f] = bc->second;
  }

  AMANZI_ASSERT(bc_level_flux_lvl_->size()==bc_level_flux_vel_->size());

  for (auto bc_lvl=bc_level_flux_lvl_->begin(), bc_vel=bc_level_flux_vel_->begin();
       bc_lvl != bc_level_flux_lvl_->end(); ++bc_lvl, ++bc_vel){

    int f = bc_lvl->first;
    markers[f] = Operators::OPERATOR_BC_NEUMANN;
    double val = bc_lvl->second;
    if (elevation[0][f] > val) values[f] = 0;
    else {
      values[f] = val * bc_vel->second;
    }
  }
  
  // zero gradient: grad h = 0 implies that q = -k grad z
  // -- cannot be done yet as rel perm update is done after this and is needed.
  // -- Instead zero gradient BCs are done in FixBCs methods.

  // Seepage face head boundary condition
  if (bc_seepage_head_->size() > 0) {

    S->GetFieldEvaluator(Keys::getKey(domain_,"ponded_depth"))->HasFieldChanged(S.ptr(), name_);

    const Epetra_MultiVector& h_c = *S->GetFieldData(Keys::getKey(domain_,"ponded_depth"))->ViewComponent("cell");
    const Epetra_MultiVector& elevation_c = *S->GetFieldData(Keys::getKey(domain_,"elevation"))->ViewComponent("cell");


    for (Functions::BoundaryFunction::Iterator bc = bc_seepage_head_->begin(); 
         bc != bc_seepage_head_->end(); ++bc) {
      int f = bc->first;
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int c = cells[0];

      double hz_f = bc->second + elevation[0][f];
      double hz_c = h_c[0][c] + elevation_c[0][c];

      if (hz_f >= hz_c) {
        markers[f] = Operators::OPERATOR_BC_NEUMANN;
        values[f] = 0.0;
      } else {
        markers[f] = Operators::OPERATOR_BC_DIRICHLET;
        values[f] = hz_f;
      }
    }
  }

  // Seepage face pressure boundary condition
  if (bc_seepage_pressure_->size() > 0) {
    S->GetFieldEvaluator(Keys::getKey(domain_,"ponded_depth"))->HasFieldChanged(S.ptr(), name_);

    const Epetra_MultiVector& h_cells = *S->GetFieldData(Keys::getKey(domain_,"ponded_depth"))->ViewComponent("cell");
    const Epetra_MultiVector& elevation_cells = *S->GetFieldData(Keys::getKey(domain_,"elevation"))->ViewComponent("cell");
    const Epetra_MultiVector& rho_l = *S->GetFieldData(Keys::getKey(domain_,"mass_density_liquid"))->ViewComponent("cell");
    double gz = -(*S->GetConstantVectorData("gravity"))[2];
    const double& p_atm = *S->GetScalarData("atmospheric_pressure");

    if (S->HasFieldEvaluator(Keys::getKey(domain_,"mass_density_ice"))) {
      // thermal model of height
      const Epetra_MultiVector& eta = *S->GetFieldData(Keys::getKey(domain_,"unfrozen_fraction"))->ViewComponent("cell");
      const Epetra_MultiVector& rho_i = *S->GetFieldData(Keys::getKey(domain_,"mass_density_ice"))->ViewComponent("cell");

      for (Functions::BoundaryFunction::Iterator bc = bc_seepage_pressure_->begin(); 
           bc != bc_seepage_pressure_->end(); ++bc) {
        int f = bc->first;
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        int c = cells[0];

        double p0 = bc->second > p_atm ? bc->second : p_atm;
        double h0 = (p0 - p_atm) / ((eta[0][c]*rho_l[0][c] + (1.-eta[0][c])*rho_i[0][c]) * gz);
        double dz = elevation_cells[0][c] - elevation[0][f];

        if (h_cells[0][c] + dz < h0) {
          markers[f] = Operators::OPERATOR_BC_NEUMANN;
          values[f] = 0.0;
        } else {
          markers[f] = Operators::OPERATOR_BC_DIRICHLET;
          values[f] = h0 + elevation[0][f];
        }
      }

    } else {
      // non-thermal model
      for (Functions::BoundaryFunction::Iterator bc = bc_seepage_pressure_->begin(); 
           bc != bc_seepage_pressure_->end(); ++bc) {
        int f = bc->first;
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        int c = cells[0];

        double p0 = bc->second > p_atm ? bc->second : p_atm;
        double h0 = (p0 - p_atm) / (rho_l[0][c] * gz);
        double dz = elevation_cells[0][c] - elevation[0][f];

        if (h_cells[0][c] + dz < h0) {
          markers[f] = Operators::OPERATOR_BC_NEUMANN;
          values[f] = 0.0;
        } else {
          markers[f] = Operators::OPERATOR_BC_DIRICHLET;
          values[f] = h0 + elevation[0][f];
        }
      }
    }
  }

  // Critical depth boundary condition
  if (bc_critical_depth_->size() > 0) {
    S->GetFieldEvaluator(Keys::getKey(domain_,"ponded_depth"))->HasFieldChanged(S.ptr(), name_);
    
    const Epetra_MultiVector& h_c = *S->GetFieldData(Keys::getKey(domain_,"ponded_depth"))->ViewComponent("cell");
    const Epetra_MultiVector& nliq_c = *S->GetFieldData(Keys::getKey(domain_,"molar_density_liquid"))
                                       ->ViewComponent("cell");
    double gz = -(*S->GetConstantVectorData("gravity"))[2];
    
    for (Functions::BoundaryFunction::Iterator bc = bc_critical_depth_->begin();
         bc != bc_critical_depth_->end(); ++bc) {
      int f = bc->first;
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int c = cells[0];
      
      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      values[f] = sqrt(gz)*std::pow(h_c[0][c], 1.5)*nliq_c[0][c];
    }
  }


  //Tidal BC

  Teuchos::RCP<const CompositeVector> elev = S->GetFieldData(Keys::getKey(domain_,"elevation"));
  Teuchos::RCP<const CompositeVector> ponded_depth = S->GetFieldData(Keys::getKey(domain_,"ponded_depth"));

  elev->ScatterMasterToGhosted();
  ponded_depth->ScatterMasterToGhosted();
  
  const Epetra_MultiVector& elevation_f = *elev -> ViewComponent("face",false);
  const Epetra_MultiVector& elevation_c = *elev -> ViewComponent("cell",false);
  const Epetra_MultiVector& ponded_c = *ponded_depth -> ViewComponent("cell",false);
  Teuchos::RCP<const CompositeVector> flux =
       S->GetFieldData(Keys::getKey(domain_,"mass_flux"));
  const Epetra_MultiVector& flux_f = *flux -> ViewComponent("face",false);  
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  for (Functions::BoundaryFunction::Iterator bc=bc_tidal_->begin();
       bc!=bc_tidal_->end(); ++bc) {

    int f = bc->first;

    AmanziMesh::Entity_ID_List cells, faces;
    std::vector<int> fdirs;
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    AMANZI_ASSERT(cells.size() == 1);
    AmanziMesh::Entity_ID c = cells[0];
      
    if (f < nfaces_owned) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
      int j=0;
      for (j=0; j<faces.size(); j++){
        if (faces[j] == f) break;
      }
      if (j >= faces.size()) {
        Errors::Message message("Overland_pressure PK: boundary face is not found in the boundary cell.");
        Exceptions::amanzi_throw(message);
      }
      double h0 = bc->second;

      if (vo_->os_OK(Teuchos::VERB_HIGH)){
        *vo_->os() << "Tidal BC: f="<<f<< "h0 "<<h0<<" h0 - elevation_f[0][f] "<<h0 - elevation_f[0][f]<<"\n";
      }
      
      if ((h0 - elevation_f[0][f]  < min_tidal_bc_ponded_depth_) ) {      
          
        markers[f] = Operators::OPERATOR_BC_NEUMANN;
        values[f] = 0.;
        
      }else{        
        markers[f] = Operators::OPERATOR_BC_DIRICHLET;
        values[f] = h0;
        
      }

      if (vo_->os_OK(Teuchos::VERB_HIGH)){
           *vo_->os() << "Tidal BC: f="<<f<<" type "<<markers[f]<<" val "<<values[f]<<"\n";
      }
      
    }

  }

  // check that there are no internal faces and mark all remaining boundary conditions as zero flux conditions

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    if ((markers[f] != Operators::OPERATOR_BC_NONE) && (ncells == 2)) {
      std::stringstream mesg_stream;
      mesg_stream << "Tried to set a boundary condition on an internal face ";
      Errors::Message mesg(mesg_stream.str());
      amanzi_throw(mesg);
    }
    if ((markers[f] == Operators::OPERATOR_BC_NONE) && (ncells == 1)) {
        markers[f] = Operators::OPERATOR_BC_NEUMANN;
        values[f] = 0.0;
    }
  }
  

}


// -----------------------------------------------------------------------------
// Add a boundary marker to owned faces.
// -----------------------------------------------------------------------------
void
OverlandPressureFlow::ApplyBoundaryConditions_(const Teuchos::Ptr<CompositeVector>& u,
        const Teuchos::Ptr<const CompositeVector>& elev) {

  auto& markers = bc_markers();
  auto& values = bc_values();

  if (u->HasComponent("face")) {
    const Epetra_MultiVector& elevation = *elev->ViewComponent("face");

    Epetra_MultiVector& u_f = *u->ViewComponent("face",false);
    unsigned int nfaces = u_f.MyLength();
    for (unsigned int f=0; f!=nfaces; ++f) {
      if (markers[f] == Operators::OPERATOR_BC_DIRICHLET) {
        u_f[0][f] = (values[f] - elevation[0][f]);
      }
    }
  } else if (u->HasComponent("boundary_face")) {
    const Epetra_MultiVector& elevation = *elev->ViewComponent("face");
    const Epetra_Map& vandalay_map = mesh_->exterior_face_map(false);
    const Epetra_Map& face_map = mesh_->face_map(false);

    const Epetra_MultiVector& u_c = *u->ViewComponent("cell",false);
    Epetra_MultiVector& u_bf = *u->ViewComponent("boundary_face",false);
    unsigned int nfaces = u_bf.MyLength();
    for (unsigned int bf=0; bf!=nfaces; ++bf) {
      AmanziMesh::Entity_ID f = face_map.LID(vandalay_map.GID(bf));
      if (markers[f] == Operators::OPERATOR_BC_DIRICHLET) {
        u_bf[0][bf] = (values[f] - elevation[0][f]);
      } else if (markers[f] == Operators::OPERATOR_BC_NEUMANN) {
        // NOTE: this is not really the correct value, but it is close.  Doing
        //  this simply makes the upwinded conductivities make more sense, and
        //  changes no answers as boundary faces and their resulting
        //  conductivity are not used in Neumann conditions.
        AmanziMesh::Entity_ID_List cells;
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        u_bf[0][bf] = u_c[0][cells[0]];
      }
    }
  }    
};


void OverlandPressureFlow::FixBCsForOperator_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<Operators::PDE_Diffusion>& diff_op) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "    Tweaking BCs for the Operator." << std::endl;

  auto& markers = bc_markers();
  auto& values = bc_values();

  // Now we can safely calculate q = -k grad z for zero-gradient problems
  Teuchos::RCP<const CompositeVector> elev = S->GetFieldData(Keys::getKey(domain_,"elevation"));
  Teuchos::RCP<const CompositeVector> ponded_depth = S->GetFieldData(Keys::getKey(domain_,"ponded_depth"));  
  Teuchos::RCP<const CompositeVector> upw_conductivity =
    S->GetFieldData(Keys::getKey(domain_,"upwind_overland_conductivity"));
  Teuchos::RCP<const CompositeVector> conductivity =
       S->GetFieldData(Keys::getKey(domain_,"overland_conductivity"));
  Teuchos::RCP<const CompositeVector> flux =
       S->GetFieldData(Keys::getKey(domain_,"mass_flux"));

  const Epetra_MultiVector& elevation_f = *elev -> ViewComponent("face",false);
  const Epetra_MultiVector& elevation_c = *elev -> ViewComponent("cell",false);
  const Epetra_MultiVector& ponded_c = *ponded_depth -> ViewComponent("cell",false);
  const Epetra_MultiVector& upw_cond_f = *upw_conductivity -> ViewComponent("face",false);
  const Epetra_MultiVector& cond_c = *conductivity -> ViewComponent("cell",false);
  const Epetra_MultiVector& flux_f = *flux -> ViewComponent("face",false);  

  std::vector<WhetStone::DenseMatrix>& Aff = diff_op->local_op()->matrices;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
       bc!=bc_zero_gradient_->end(); ++bc) {
    int f = bc->first;

    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    AMANZI_ASSERT(cells.size() == 1);
    AmanziMesh::Entity_ID c = cells[0];

    if (f < nfaces_owned) {
      double dp = elevation_f[0][f] - elevation_c[0][c];
      double bc_val = -dp * Aff[f](0,0);

      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      values[f] = bc_val / mesh_->face_area(f);        

    }
  }

  const double& p_atm = *S->GetScalarData("atmospheric_pressure");
  const Epetra_MultiVector& rho_l = *S->GetFieldData(Keys::getKey(domain_,"mass_density_liquid"))->ViewComponent("cell");
  double gz = -(*S->GetConstantVectorData("gravity"))[2];  

  for (Functions::BoundaryFunction::Iterator bc=bc_tidal_->begin();
       bc!=bc_tidal_->end(); ++bc) {

    int f = bc->first;
    AmanziMesh::Entity_ID_List cells, faces;
    std::vector<int> fdirs;
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    AMANZI_ASSERT(cells.size() == 1);
    AmanziMesh::Entity_ID c = cells[0];
                                               
    if (f < nfaces_owned) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
      int j=0;
      for (j=0; j<faces.size(); j++){
        if (faces[j] == f) break;
      }
      if (j >= faces.size()) {
        Errors::Message message("Overland_pressure PK: boundary face is not found in the boundary cell.");
        Exceptions::amanzi_throw(message);
      }
      double h0 = bc->second;

      
      if ((h0 - elevation_f[0][f] < min_tidal_bc_ponded_depth_)) {
        double dp = elevation_f[0][f] - elevation_c[0][c];
        if (dp < 0){
          double bc_val = -dp * Aff[f](0,0);         
          markers[f] = Operators::OPERATOR_BC_NEUMANN;
          values[f] = bc_val / mesh_->face_area(f);
        }else{
          markers[f] = Operators::OPERATOR_BC_NEUMANN;
          values[f] = 0.;
        }
      }else{        
        markers[f] = Operators::OPERATOR_BC_DIRICHLET;
        values[f] = h0;
      }

      if (vo_->os_OK(Teuchos::VERB_HIGH)){
        *vo_->os() << "Tidal BC2: f = "<<f<<" type "<<markers[f]<<" val "<<values[f]<<"\n";
      }

    }
  }

  
};


void OverlandPressureFlow::FixBCsForPrecon_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "    Tweaking BCs for the PC." << std::endl;

  // // Attempt of a hack to deal with zero rel perm
  // double eps = 1.e-30;
  // Teuchos::RCP<CompositeVector> relperm =
  //     S->GetFieldData("upwind_overland_conductivity", name_);
  // for (unsigned int f=0; f!=relperm->size("face"); ++f) {
  //   if ((*relperm)("face",f) < eps) {
  //     if (bc_markers_[f] == Operators::OPERATOR_BC_NEUMANN) {
  //       bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
  //     } else if (bc_markers_[f] == Operators::OPERATOR_BC_NONE) {
  //       bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
  //     }
  //   }
  // }
};

// void OverlandPressureFlow::FixBCsForConsistentFaces_(const Teuchos::Ptr<State>& S) {
//   Teuchos::OSTab tab = vo_->getOSTab();
//   if (vo_->os_OK(Teuchos::VERB_EXTREME))
//     *vo_->os() << "    Tweaking BCs for calculation of consistent faces." << std::endl;

//   // // If the rel perm is 0, the face value drops out and is unconstrained.
//   // // Therefore we set it to Dirichlet to eliminate it from the system.
//   // double eps = 1.e-30;
//   // const Epetra_MultiVector& elevation = *S->GetFieldData("elevation")
//   //     ->ViewComponent("face",false);
//   // Teuchos::RCP<CompositeVector> relperm =
//   //     S->GetFieldData("upwind_overland_conductivity", name_);

//   // for (unsigned int f=0; f!=relperm->size("face"); ++f) {
//   //   if ((*relperm)("face",f) < eps) {
//   //     if (bc_markers_[f] == Operators::OPERATOR_BC_NEUMANN) {
//   //       bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
//   //       values[f] =  elevation[0][f];
//   //     } else if (bc_markers_[f] == Operators::OPERATOR_BC_NONE) {
//   //       bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
//   //       values[f] =  elevation[0][f];
//   //     }
//   //   }
//   // }

//   // Now we can safely calculate q = -k grad z for zero-gradient problems
//   Teuchos::RCP<const CompositeVector> elev = S->GetFieldData("elevation");
//   elev->ScatterMasterToGhosted();
//   const Epetra_MultiVector& elevation_f = *elev->ViewComponent("face",false);
//   const Epetra_MultiVector& elevation_c = *elev->ViewComponent("cell",false);

//   std::vector<WhetStone::DenseMatrix>& Aff_cells =
//       matrix_diff_->local_op()->matrices;
//   Epetra_MultiVector& rhs_f = *matrix_->rhs()->ViewComponent("face",false);

//   int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
//   for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
//        bc!=bc_zero_gradient_->end(); ++bc) {
//     int f = bc->first;

//     AmanziMesh::Entity_ID_List cells;
//     mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
//     AMANZI_ASSERT(cells.size() == 1);
//     AmanziMesh::Entity_ID c = cells[0];

//     if (c < ncells_owned) {
//       AmanziMesh::Entity_ID_List faces;
//       mesh_->cell_get_faces(c, &faces);

//       std::vector<double> dp(faces.size());
//       for (unsigned int n=0; n!=faces.size(); ++n) {
//         dp[n] = elevation_f[0][faces[n]] - elevation_c[0][c];
//       }
//       unsigned int my_n = std::find(faces.begin(), faces.end(), f) - faces.begin();
//       AMANZI_ASSERT(my_n !=faces.size());

//       double bc_val = 0.;
//       for (unsigned int m=0; m!=faces.size(); ++m) {
//         bc_val -= Aff_cells[c](my_n,m) * dp[m];
//       }

//       // Apply the BC to the matrix
//       rhs_f[0][f] -= bc_val;
//     }
//   }
// };


bool OverlandPressureFlow::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
                                       Teuchos::RCP<TreeVector> u) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Modifying predictor:" << std::endl;

  return false;
};


void OverlandPressureFlow::CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u) {
  // // VerboseObject stuff.
  // Teuchos::OSTab tab = vo_->getOSTab();

  // // update the rel perm according to the scheme of choice
  // ChangedSolution();
  // UpdatePermeabilityData_(S_next_.ptr());

  // // update boundary conditions
  // bc_head_->Compute(S_next_->time());
  // bc_flux_->Compute(S_next_->time());
  // bc_seepage_head_->Compute(S_next_->time());
  // bc_seepage_pressure_->Compute(S_next_->time());
  // UpdateBoundaryConditions_(S_next_.ptr());

  // // update the stiffness matrix
  // Teuchos::RCP<const CompositeVector> cond =
  //   S_next_->GetFieldData("upwind_overland_conductivity", name_);
  // matrix_->CreateMFDstiffnessMatrices(cond.ptr());
  // matrix_->CreateMFDrhsVectors();

  // // Patch up BCs in the case of zero conductivity
  // FixBCsForConsistentFaces_(S_next_.ptr());

  // // Grab needed data.
  // S_next_->GetFieldEvaluator("pres_elev")->HasFieldChanged(S_next_.ptr(), name_);
  // Teuchos::RCP<CompositeVector> pres_elev = S_next_->GetFieldData("pres_elev","pres_elev");

  // // Update the preconditioner with darcy and gravity fluxes
  // // skip accumulation terms, they're not needed
  // // Assemble
  // matrix_->ApplyBoundaryConditions(bc_markers_, values);

  // // derive the consistent faces, involves a solve
  // matrix_->UpdateConsistentFaceConstraints(pres_elev.ptr());

  // // back out heights from pres_elev
  // const Epetra_MultiVector& elevation = *S_next_->GetFieldData("elevation")
  //     ->ViewComponent("face",false);
  // u->ViewComponent("face",false)->Update(1., *pres_elev->ViewComponent("face",false),
  //         -1., elevation, 0.);
}


AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
OverlandPressureFlow::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                 Teuchos::RCP<const TreeVector> u,
                 Teuchos::RCP<TreeVector> du) {
  Teuchos::OSTab tab = vo_->getOSTab();

  // if the primary variable has boundary face, this is for upwinding rel
  // perms and is never actually used.  Make sure it does not go to undefined
  // pressures.
  if (du->Data()->HasComponent("boundary_face")) {
    du->Data()->ViewComponent("boundary_face")->PutScalar(0.);
  }

  // debugging -- remove me! --etc
  for (CompositeVector::name_iterator comp=du->Data()->begin();
       comp!=du->Data()->end(); ++comp) {
    Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp,false);
    double max, l2;
    du_c.NormInf(&max);
    du_c.Norm2(&l2);
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Linf, L2 pressure correction (" << *comp << ") = " << max << ", " << l2 << std::endl;
    }
  }

  // limit by capping corrections when they cross atmospheric pressure
  // (where pressure derivatives are discontinuous)
  int my_limited = 0;
  int n_limited_spurt = 0;
  if (patm_limit_ > 0.) {
    double patm = *S_next_->GetScalarData("atmospheric_pressure");

    Epetra_MultiVector& du_c = *du->Data()->ViewComponent("cell",false);
    const Epetra_MultiVector& u_c = *u->Data()->ViewComponent("cell",false);

    for (int c=0; c!=du_c.MyLength(); ++c) {
      if ((u_c[0][c] < patm) &&
          (u_c[0][c] - du_c[0][c] > patm + patm_limit_)) {
        du_c[0][c] = u_c[0][c] - (patm + patm_limit_);          
        my_limited++;
      } else if ((u_c[0][c] > patm) &&
                 (u_c[0][c] - du_c[0][c] < patm - patm_limit_)) {
        du_c[0][c] = u_c[0][c] - (patm - patm_limit_);          
        my_limited++;
      }
    }
    mesh_->get_comm()->MaxAll(&my_limited, &n_limited_spurt, 1);
  }

  if (patm_hard_limit_) {
    double patm = *S_next_->GetScalarData("atmospheric_pressure");

    Epetra_MultiVector& du_c = *du->Data()->ViewComponent("cell",false);
    const Epetra_MultiVector& u_c = *u->Data()->ViewComponent("cell",false);

    for (int c=0; c!=du_c.MyLength(); ++c) {
      if (u_c[0][c] - du_c[0][c] < patm) {
        du_c[0][c] = u_c[0][c] - patm;          
        my_limited++;
      }
    }
    mesh_->get_comm()->MaxAll(&my_limited, &n_limited_spurt, 1);
  }

  if (n_limited_spurt > 0) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "  limiting the spurt." << std::endl;
    }
  }

  // debugging -- remove me! --etc
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    for (CompositeVector::name_iterator comp=du->Data()->begin();
         comp!=du->Data()->end(); ++comp) {
      Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp,false);
      double max, l2;
      du_c.NormInf(&max);
      du_c.Norm2(&l2);
      *vo_->os() << "Linf, L2 pressure correction (" << *comp << ") = " << max << ", " << l2 << std::endl;
    }
  }
  
  // Limit based on a max pressure change
  my_limited = 0;
  int n_limited_change = 0;

  if (p_limit_ > 0.) {
    for (CompositeVector::name_iterator comp=du->Data()->begin();
         comp!=du->Data()->end(); ++comp) {
      Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp,false);

      double max;
      du_c.NormInf(&max);
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "Max overland pressure correction (" << *comp << ") = " << max << std::endl;
      }
      
      for (int c=0; c!=du_c.MyLength(); ++c) {
        if (std::abs(du_c[0][c]) > p_limit_) {
          du_c[0][c] = ((du_c[0][c] > 0) - (du_c[0][c] < 0)) * p_limit_;
          my_limited++;
        }
      }
    }

    mesh_->get_comm()->SumAll(&my_limited, &n_limited_change, 1);
  }

  if (n_limited_change > 0) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "  limited by pressure." << std::endl;
    }
  }
  
  // debugging -- remove me! --etc
  for (CompositeVector::name_iterator comp=du->Data()->begin();
       comp!=du->Data()->end(); ++comp) {
    Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp,false);
    double max, l2;
    du_c.NormInf(&max);
    du_c.Norm2(&l2);
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Linf, L2 pressure correction (" << *comp << ") = " << max << ", " << l2 << std::endl;
    }
  }

  if (n_limited_spurt > 0) {
    return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED_LAG_BACKTRACKING;
  } else if (n_limited_change > 0) {

    return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
  }
  return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}
  
} // namespace
} // namespace

