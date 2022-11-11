/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code.
License: BSD
Authors: Daniil Svyatsky(dasvyat@gmail.com)
------------------------------------------------------------------------- */
#include "boost/math/special_functions/fpclassify.hpp"

#include "boost/algorithm/string/predicate.hpp"

#include "Epetra_Import.h"

#include "flow_bc_factory.hh"

#include "Point.hh"

#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"

#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"

#include "predictor_delegate_bc_flux.hh"
#include "wrm_evaluator.hh"
#include "rel_perm_evaluator.hh"
#include "richards_water_content_evaluator.hh"
#include "OperatorDefs.hh"
#include "BoundaryFlux.hh"
#include "pk_helpers.hh"

#include "preferential.hh"

#define DEBUG_RES_FLAG 0


namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------

Preferential::Preferential(Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& glist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution) :
    PK(pk_tree, glist,  S, solution),
    Richards(pk_tree, glist,  S, solution)
{
    coef_grav_key_ = Keys::readKey(*plist_, domain_, "gravity conductivity", "gravity_relative_permeability");
    // set up an additional primary variable evaluator for flux
    Teuchos::ParameterList& pv_sublist = S->GetEvaluatorList(flux_key_);
    pv_sublist.set("field evaluator type", "primary variable");
}

// -------------------------------------------------------------
// Setup data
// -------------------------------------------------------------
void Preferential::Setup()
{
  PK_PhysicalBDF_Default::Setup();
  SetupPreferentialFlow_();
  SetupDiscretization_();
  SetupPhysicalEvaluators_();
};


// -------------------------------------------------------------
// Pieces of the construction process that are common to all
// Preferential-like PKs.
// -------------------------------------------------------------
void Preferential::SetupPreferentialFlow_()
{
  // Get data for special-case entities.
  S_->Require<CompositeVector,CompositeVectorSpace>(cell_vol_key_, tag_next_).SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(cell_vol_key_, tag_next_);


  // Set up Operators
  // -- boundary conditions
  Teuchos::ParameterList bc_plist = plist_->sublist("boundary conditions", true);
  FlowBCFactory bc_factory(mesh_, bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_head_ = bc_factory.CreateHead();
  bc_level_ = bc_factory.CreateFixedLevel();
  bc_flux_ = bc_factory.CreateMassFlux();
  bc_seepage_ = bc_factory.CreateSeepageFacePressure();
  bc_seepage_->Compute(0.); // compute at t=0 to set up
  std::tie(bc_seepage_infilt_explicit_, bc_seepage_infilt_) =
    bc_factory.CreateSeepageFacePressureWithInfiltration();
  bc_seepage_infilt_->Compute(0.); // compute at t=0 to set up
  bc_rho_water_ = bc_plist.get<double>("hydrostatic water density [kg m^-3]",1000.);

  // scaling for permeability
  perm_scale_ = plist_->get<double>("permeability rescaling", 1.e7);

  // permeability type - scalar or tensor?
  Teuchos::ParameterList& perm_list = S_->GetEvaluatorList(perm_key_);
  std::string perm_type = perm_list.get<std::string>("permeability type", "scalar");
  if (perm_type == "scalar") {
    perm_tensor_rank_ = 1;
    num_perm_vals_ = 1;
  } else if (perm_type == "horizontal and vertical") {
    perm_tensor_rank_ = 2;
    num_perm_vals_ = 2;
  } else if (perm_type == "diagonal tensor") {
    perm_tensor_rank_ = 2;
    num_perm_vals_ = mesh_->space_dimension();
  } else if (perm_type == "full tensor") {
    perm_tensor_rank_ = 2;
    num_perm_vals_ = (mesh_->space_dimension() == 3) ? 6 : 3;
  } else {
    Errors::Message message("`permeability type` must be one of the following: \"scalar\", \"diagonal tensor\", \"full tensor\", or \"horizontal and vertical\".");
    Exceptions::amanzi_throw(message);
  }

  // -- linear tensor coefficients
  unsigned int c_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  K_ = Teuchos::rcp(new std::vector<WhetStone::Tensor>(c_owned));
  for (unsigned int c=0; c!=c_owned; ++c) {
    (*K_)[c].Init(mesh_->space_dimension(),perm_tensor_rank_);
  }

  // -- nonlinear coefficients/upwinding
  if (plist_->isParameter("surface rel perm strategy")) {
    clobber_policy_ = plist_->get<std::string>("surface rel perm strategy");
  } else if (plist_->get<bool>("clobber surface rel perm", false)) {
    clobber_policy_ = "clobber";
  } else if (plist_->get<bool>("max surface rel perm", false)) {
    clobber_policy_ = "max";
  } else if (plist_->get<bool>("unsaturated clobber surface rel perm", false)) {
    clobber_policy_ = "unsaturated";
  } else {
    clobber_policy_ = "none";
  }
  clobber_boundary_flux_dir_ = plist_->get<bool>("clobber boundary flux direction for upwinding", false);

  std::string method_name = plist_->get<std::string>("relative permeability method", "upwind with Darcy flux");
  if (method_name == "upwind with gravity") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindGravityFlux(name_,
                                                               tag_next_, K_));
    Krel_method_ = Operators::UPWIND_METHOD_GRAVITY;
  } else if (method_name == "cell centered") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindCellCentered(name_, tag_next_));
    Krel_method_ = Operators::UPWIND_METHOD_CENTERED;
  } else if (method_name == "upwind with Darcy flux") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
            tag_next_, flux_dir_key_, 1.e-5));
    Krel_method_ = Operators::UPWIND_METHOD_TOTAL_FLUX;
  } else if (method_name == "arithmetic mean") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindArithmeticMean(name_,
            tag_next_));
    Krel_method_ = Operators::UPWIND_METHOD_ARITHMETIC_MEAN;
  } else {
    std::stringstream messagestream;
    messagestream << "Preferential PK has no upwinding method named: " << method_name;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // -- require the data on appropriate locations
  std::string coef_location = upwinding_->CoefficientLocation();
  if (coef_location == "upwind: face") {
    std::vector<AmanziMesh::Entity_kind> locations2(2);
    std::vector<std::string> names2(2);
    std::vector<int> num_dofs2(2,1);
    locations2[0] = AmanziMesh::FACE;
    locations2[1] = AmanziMesh::FACE;
    names2[0] = "face";
    names2[1] = "grav";    
    S_->Require<CompositeVector,CompositeVectorSpace>(uw_coef_key_, tag_next_, name_).SetMesh(mesh_)
      ->SetGhosted()->SetComponents(names2, locations2, num_dofs2);
  } else if (coef_location == "standard: cell") {
    std::vector<AmanziMesh::Entity_kind> locations2(2);
    std::vector<std::string> names2(2);
    std::vector<int> num_dofs2(2,1);
    locations2[0] = AmanziMesh::CELL;
    locations2[1] = AmanziMesh::FACE;
    names2[0] = "cell";
    names2[1] = "grav";        
    S_->Require<CompositeVector,CompositeVectorSpace>(uw_coef_key_, tag_next_, name_).SetMesh(mesh_)
      ->SetGhosted()->SetComponents(names2, locations2, num_dofs2);
  } else {
    Errors::Message message("Unknown upwind coefficient location in Preferential flow.");
    Exceptions::amanzi_throw(message);
  }
  S_->GetRecordW(uw_coef_key_, tag_next_, name_).set_io_vis(false);


  Teuchos::ParameterList& mfd_plist = plist_->sublist("diffusion");
  mfd_plist.set("nonlinear coefficient", coef_location);
  mfd_plist.set("gravity", true);  
  Teuchos::ParameterList& mfd_pc_plist = plist_->sublist("diffusion preconditioner");
  mfd_pc_plist.set("nonlinear coefficient", coef_location);
  mfd_pc_plist.set("gravity", true);
  if (!mfd_pc_plist.isParameter("discretization primary"))
    mfd_pc_plist.set("discretization primary", mfd_plist.get<std::string>("discretization primary"));
  if (!mfd_pc_plist.isParameter("discretization secondary") && mfd_plist.isParameter("discretization secondary"))
    mfd_pc_plist.set("discretization secondary", mfd_plist.get<std::string>("discretization secondary"));
  if (!mfd_pc_plist.isParameter("schema") && mfd_plist.isParameter("schema"))
    mfd_pc_plist.set("schema", mfd_plist.get<Teuchos::Array<std::string> >("schema"));
  if (mfd_pc_plist.get<bool>("include Newton correction", false)) {
    if (mfd_pc_plist.get<std::string>("discretization primary") == "fv: default") {
      mfd_pc_plist.set("Newton correction", "true Jacobian");
    } else {
      mfd_pc_plist.set("Newton correction", "approximate Jacobian");
    }
  }

  precon_used_ = plist_->isSublist("preconditioner") ||
    plist_->isSublist("inverse") ||
    plist_->isSublist("linear solver");
  if (precon_used_) {
    mfd_pc_plist.set("inverse", plist_->sublist("inverse"));
    // old style... deprecate me!
    mfd_pc_plist.sublist("inverse").setParameters(plist_->sublist("preconditioner"));
    mfd_pc_plist.sublist("inverse").setParameters(plist_->sublist("linear solver"));
  }

  Operators::PDE_DiffusionFactory opfactory;
  preconditioner_diff_ = opfactory.CreateWithGravity(mfd_pc_plist, mesh_, bc_);
  preconditioner_ = preconditioner_diff_->global_operator();
  
  //    If using approximate Jacobian for the preconditioner, we also need
  //    derivative information.  For now this means upwinding the derivative.
  jacobian_ = mfd_pc_plist.get<std::string>("Newton correction", "none") != "none";
  if (jacobian_) {
    jacobian_lag_ = mfd_pc_plist.get<int>("Newton correction lag", 0);

    if (mfd_pc_plist.get<std::string>("discretization primary") != "fv: default"){
      // MFD -- upwind required
      //      dcoef_key_ = Keys::getDerivKey(coef_key_, key_);
      //dcoef_grav_key_ = Keys::getDerivKey(coef_grav_key_, key_);
      duw_coef_key_ = Keys::getDerivKey(uw_coef_key_, key_);
      std::vector<AmanziMesh::Entity_kind> locations2(2);
      std::vector<std::string> names2(2);
      std::vector<int> num_dofs2(2,1);
      locations2[0] = AmanziMesh::FACE;
      locations2[1] = AmanziMesh::FACE;
      names2[0] = "face";
      names2[1] = "grav";        
      S_->Require<CompositeVector,CompositeVectorSpace>(duw_coef_key_, tag_next_, name_).SetMesh(mesh_)
        ->SetGhosted()->SetComponents(names2, locations2, num_dofs2);

      upwinding_deriv_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
                                            tag_next_, flux_dir_key_, 1.e-8));
    } else {
      // FV -- no upwinding
      //dcoef_key_ = Keys::getDerivKey(coef_key_, key_);
      //dcoef_grav_key_ = Keys::getDerivKey(coef_grav_key_, key_);
      duw_coef_key_ = std::string();
    }
  }
  
  // Globalization and other timestep control flags
  // -- predictors
  modify_predictor_with_consistent_faces_ =
    plist_->get<bool>("modify predictor with consistent faces", false);
  modify_predictor_bc_flux_ =
    plist_->get<bool>("modify predictor for flux BCs", false);
  modify_predictor_first_bc_flux_ =
    plist_->get<bool>("modify predictor for initial flux BCs", false);
  modify_predictor_wc_ =
    plist_->get<bool>("modify predictor via water content", false);

  // -- correctors
  p_limit_ = plist_->get<double>("limit correction to pressure change [Pa]", -1.);
  patm_limit_ = plist_->get<double>("limit correction to pressure change when crossing atmospheric [Pa]", -1.);

  // -- valid step controls
  sat_change_limit_ = plist_->get<double>("max valid change in saturation in a time step [-]", -1.);
  sat_ice_change_limit_ = plist_->get<double>("max valid change in ice saturation in a time step [-]", -1.);
}


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void Preferential::SetupPhysicalEvaluators_() {
  // -- Absolute permeability.
  //       For now, we assume scalar permeability.  This will change.
  S_->Require<CompositeVector,CompositeVectorSpace>(perm_key_, tag_next_).SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, num_perm_vals_);
  S_->RequireEvaluator(perm_key_, tag_next_);

  // -- water content, and evaluator
  S_->Require<CompositeVector,CompositeVectorSpace>(conserved_key_, tag_next_).SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireDerivative<CompositeVector,CompositeVectorSpace>(conserved_key_,
          tag_next_, key_, tag_next_);

  //    and at the current time, where it is a copy evaluator
  requireAtCurrent(conserved_key_, tag_current_, *S_, name_);

  // -- Water retention evaluators
  // This deals with deprecated location for the WRM list (in the PK).  Move it
  // to state.
  // if (plist_->isSublist("water retention evaluator")) {
  //   auto& wrm_plist = S_->GetEvaluatorList(sat_key_);
  //   wrm_plist.setParameters(plist_->sublist("water retention evaluator"));
  //   wrm_plist.set("field evaluator type", "WRM");
  // }
    if (plist_->isSublist("water retention evaluator")) {
    auto& wrm_plist = S_->GetEvaluatorList(sat_key_);
    wrm_plist.setParameters(plist_->sublist("water retention evaluator"));
    wrm_plist.set("evaluator type", "WRM");
  }
  if (!S_->HasEvaluator(coef_key_, tag_next_) &&
      (S_->GetEvaluatorList(coef_key_).numParams() == 0)) {
    Teuchos::ParameterList& kr_plist = S_->GetEvaluatorList(coef_key_);
    kr_plist.setParameters(S_->GetEvaluatorList(sat_key_));
    kr_plist.set<std::string>("evaluator type", "WRM rel perm");
  }


    // -- Water retention evaluators for gravity term
  // This deals with deprecated location for the WRM list (in the PK).  Move it
  // if (plist_->isSublist("water retention evaluator for gravity term")) {
  //   auto& wrm_plist = S_->GetEvaluatorList(sat_key_);
  //   wrm_plist.setParameters(plist_->sublist("water retention evaluator for gravity term"));
  //   wrm_plist.set("field evaluator type", "WRM");
  // }

  if (plist_->isSublist("water retention evaluator for gravity term")) {
    auto& wrm_plist = S_->GetEvaluatorList(sat_key_);
    wrm_plist.setParameters(plist_->sublist("water retention evaluator for gravity term"));
    wrm_plist.set("evaluator type", "WRM");
  }
  if (!S_->HasEvaluator(coef_grav_key_, tag_next_) &&
      (S_->GetEvaluatorList(coef_grav_key_).numParams() == 0)) {
    Teuchos::ParameterList& kr_plist = S_->GetEvaluatorList(coef_grav_key_);
    kr_plist.setParameters(S_->GetEvaluatorList(sat_key_));
    kr_plist.set<std::string>("evaluator type", "WRM rel perm");
  }
  
  // -- saturation
  requireAtNext(sat_key_, tag_next_, *S_)
    .SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(sat_gas_key_, tag_next_, *S_)
    .SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  auto& wrm = S_->RequireEvaluator(sat_key_, tag_next_);

  //    and at the current time, where it is a copy evaluator
  requireAtCurrent(sat_key_, tag_current_, *S_, name_);

  // -- rel perm
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::BOUNDARY_FACE;
  names2[0] = "cell";
  names2[1] = "boundary_face";
  S_->Require<CompositeVector,CompositeVectorSpace>(coef_key_, tag_next_)
    .SetMesh(mesh_)->SetGhosted()->AddComponents(names2, locations2, num_dofs2);
  S_->GetEvaluatorList(coef_key_).set<double>("permeability rescaling", perm_scale_);
  S_->RequireEvaluator(coef_key_, tag_next_);

  // -- rel perm grav
  S_->Require<CompositeVector,CompositeVectorSpace>(coef_grav_key_, tag_next_).SetMesh(mesh_)->SetGhosted()
      ->AddComponents(names2, locations2, num_dofs2);
  S_->GetEvaluatorList(coef_grav_key_).set<double>("permeability rescaling", perm_scale_);
  S_->RequireEvaluator(coef_grav_key_, tag_next_);

  // -- get the WRM models
  auto wrm_eval = dynamic_cast<Flow::WRMEvaluator*>(&wrm);
  AMANZI_ASSERT(wrm_eval != nullptr);
  wrms_ = wrm_eval->get_WRMs();

  // -- Liquid density and viscosity for the transmissivity.
  S_->Require<CompositeVector,CompositeVectorSpace>(molar_dens_key_, tag_next_).SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(molar_dens_key_, tag_next_);

  // -- liquid mass density for the gravity fluxes
  S_->Require<CompositeVector,CompositeVectorSpace>(mass_dens_key_, tag_next_).SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(mass_dens_key_, tag_next_); // simply picks up the molar density one.

}




// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool Preferential::UpdatePermeabilityData_(const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating permeability?";

  Teuchos::RCP<const CompositeVector> rel_perm = S_->GetPtr<CompositeVector>(coef_key_, tag);
  Teuchos::RCP<const CompositeVector> rel_perm_grav = S_->GetPtr<CompositeVector>(coef_grav_key_, tag);
  bool update_perm = S_->GetEvaluator(coef_key_, tag).Update(*S_, name_);
  update_perm |= S_->GetEvaluator(coef_grav_key_, tag).Update(*S_, name_);

  // requirements due to the upwinding method
  if (Krel_method_ == Operators::UPWIND_METHOD_TOTAL_FLUX) {
    bool update_dir = S_->GetEvaluator(mass_dens_key_, tag).Update(*S_, name_);
    update_dir |= S_->GetEvaluator(key_, tag).Update(*S_, name_);

    if (update_dir) {
      // update the direction of the flux -- note this is NOT the flux
      Teuchos::RCP<const CompositeVector> rho = S_->GetPtr<CompositeVector>(mass_dens_key_, tag);
      Teuchos::RCP<CompositeVector> flux_dir = S_->GetPtrW<CompositeVector>(flux_dir_key_, tag, name_);
      Teuchos::RCP<const CompositeVector> pres = S_->GetPtr<CompositeVector>(key_, tag);

      if (!deform_key_.empty() && S_->GetEvaluator(deform_key_, tag_next_).Update(*S_, name_+" flux dir"))
        face_matrix_diff_->SetTensorCoefficient(K_);
      
      face_matrix_diff_->SetDensity(rho);
      face_matrix_diff_->UpdateMatrices(Teuchos::null, pres.ptr());
      //if (!pres->HasComponent("face"))
      face_matrix_diff_->ApplyBCs(true, true, true);
      face_matrix_diff_->UpdateFlux(pres.ptr(), flux_dir.ptr());

      if (clobber_boundary_flux_dir_) {
        Epetra_MultiVector& flux_dir_f = *flux_dir->ViewComponent("face",false);

        auto& markers = bc_markers();
        auto& values = bc_values();

        for (int f=0; f!=markers.size(); ++f) {
          if (markers[f] == Operators::OPERATOR_BC_NEUMANN) {
            AmanziMesh::Entity_ID_List cells;
            mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
            AMANZI_ASSERT(cells.size() == 1);
            int c = cells[0];
            AmanziMesh::Entity_ID_List faces;
            std::vector<int> dirs;
            mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
            int i = std::find(faces.begin(), faces.end(), f) - faces.begin();

            flux_dir_f[0][f] = values[f]*dirs[i];
          }
        }
      }
    }
      
    update_perm |= update_dir;
  }

  if (update_perm) {
    Teuchos::RCP<CompositeVector> uw_rel_perm = S_->GetPtrW<CompositeVector>(uw_coef_key_, tag, name_);

    // // Move rel perm on boundary_faces into uw_rel_perm on faces
    const Epetra_Import& vandelay = mesh_->exterior_face_importer();
    const Epetra_MultiVector& rel_perm_bf =
         *rel_perm->ViewComponent("boundary_face",false);
    const Epetra_MultiVector& rel_perm_grav_bf =
         *rel_perm_grav->ViewComponent("boundary_face",false);

    Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face",false);
    uw_rel_perm_f.Export(rel_perm_bf, vandelay, Insert);
    Epetra_MultiVector& uw_rel_perm_grav = *uw_rel_perm->ViewComponent("grav",false);
    uw_rel_perm_grav.Export(rel_perm_grav_bf, vandelay, Insert);
    
    // Upwind, only overwriting boundary faces if the wind says to do so.
    upwinding_->Update(*rel_perm, "cell", *uw_rel_perm, "face", *S_);
    upwinding_->Update(*rel_perm, "cell", *uw_rel_perm, "grav", *S_);

    // Epetra_MultiVector& test_face = *uw_rel_perm->ViewComponent("face",false);
    // Epetra_MultiVector& test_grav = *uw_rel_perm->ViewComponent("grav",false);
    // const Epetra_MultiVector& test_coef = *S_->GetFieldData(coef_key_)->ViewComponent("cell",false);
    // const Epetra_MultiVector& test_coef_grav = *S_->GetFieldData(coef_grav_key_)->ViewComponent("cell",false);


    // for (int f=0; f<test_face.MyLength(); f++){
    //   if (abs(test_face[0][f] - test_grav[0][f]) > 1e-10){
    //     //std::cout<<"coef "<<test_coef[0][0]<<" coef_grav "<<test_coef_grav[0][0]<<"\n";
    //     std::cout<<"face "<<f<<": "<<test_face[0][f]<<" "<<test_grav[0][f]<<"\n";
    //     exit(0);
    //   }
    // }
    
    if (clobber_policy_ == "clobber") {
      Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face",false);
      uw_rel_perm_f.Export(rel_perm_bf, vandelay, Insert);
    } else if (clobber_policy_ == "max") {
      Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face",false);
      const auto& fmap = mesh_->face_map(true);
      const auto& bfmap = mesh_->exterior_face_map(true);
      for (int bf=0; bf!=rel_perm_bf.MyLength(); ++bf) {
        auto f = fmap.LID(bfmap.GID(bf));
        if (rel_perm_bf[0][bf] > uw_rel_perm_f[0][f]) {
          uw_rel_perm_f[0][f] = rel_perm_bf[0][bf];
        }
      }
    } else if (clobber_policy_ == "unsaturated") {
      // clobber only when the interior cell is unsaturated
      Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face",false);
      const Epetra_MultiVector& pres = *S_->Get<CompositeVector>(key_, tag).ViewComponent("cell",false);
      const auto& fmap = mesh_->face_map(true);
      const auto& bfmap = mesh_->exterior_face_map(true);
      for (int bf=0; bf!=rel_perm_bf.MyLength(); ++bf) {
        auto f = fmap.LID(bfmap.GID(bf));
        AmanziMesh::Entity_ID_List fcells;
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &fcells);
        AMANZI_ASSERT(fcells.size() == 1);
        if (pres[0][fcells[0]] < 101225.) {
          uw_rel_perm_f[0][f] = rel_perm_bf[0][bf];
        } else if (pres[0][fcells[0]] < 101325.) {
          double frac = (101325. - pres[0][fcells[0]])/100.;
          uw_rel_perm_f[0][f] = rel_perm_bf[0][bf] * frac + uw_rel_perm_f[0][f] * (1-frac);
        }
      }
    }

    if (uw_rel_perm->HasComponent("face"))
      uw_rel_perm->ScatterMasterToGhosted("face");
    if (uw_rel_perm->HasComponent("grav"))
      uw_rel_perm->ScatterMasterToGhosted("grav");
  }

  // Teuchos::RCP<CompositeVector> uw_rel_perm = S_->GetFieldData(uw_coef_key_, name_);
  // Epetra_MultiVector& test_face = *uw_rel_perm->ViewComponent("face",false);
  // Epetra_MultiVector& test_grav = *uw_rel_perm->ViewComponent("grav",false);
  // const Epetra_MultiVector& test_coef = *S_->GetFieldData(coef_key_)->ViewComponent("cell",false);
  // const Epetra_MultiVector& test_coef_grav = *S_->GetFieldData(coef_grav_key_)->ViewComponent("cell",false);
  // for (int f=0; f<test_face.MyLength(); f++){
  //   if (abs(test_face[0][f] - test_grav[0][f]) > 1e-10){
  //     //std::cout<<"coef "<<test_coef[0][0]<<" coef_grav "<<test_coef_grav[0][0]<<"\n";
  //     std::cout<<"face "<<f<<": "<<test_face[0][f]<<" "<<test_grav[0][f]<<"\n";
  //     exit(0);
  //   }
  // }

  
  // debugging
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << " " << update_perm << std::endl;
  }

  return update_perm;
};


bool Preferential::UpdatePermeabilityDerivativeData_(const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating permeability derivatives?";

  bool update_perm = S_->GetEvaluator(coef_key_, tag).UpdateDerivative(*S_, name_, key_, tag);  
  update_perm |= S_->GetEvaluator(coef_grav_key_, tag).UpdateDerivative(*S_, name_, key_, tag);;  
  if (update_perm) {

    const CompositeVector& drel_perm = S_->GetDerivative<CompositeVector>(coef_key_, tag, key_, tag);
    const CompositeVector& drel_grav_perm = S_->GetDerivative<CompositeVector>(coef_grav_key_, tag, key_, tag);

    if (!duw_coef_key_.empty()) {
      CompositeVector& duw_rel_perm = S_->GetW<CompositeVector>(duw_coef_key_, tag, name_);
      duw_rel_perm.PutScalar(0.);

      // Upwind, only overwriting boundary faces if the wind says to do so.
      upwinding_deriv_->Update(drel_perm, "cell", duw_rel_perm, "face", *S_);
      upwinding_deriv_->Update(drel_grav_perm, "cell", duw_rel_perm, "grav", *S_);

      duw_rel_perm.ScatterMasterToGhosted("face");
      duw_rel_perm.ScatterMasterToGhosted("grav");
    } else {
      drel_perm.ScatterMasterToGhosted("cell");
      drel_grav_perm.ScatterMasterToGhosted("cell");
    }
  }

  // Teuchos::RCP<CompositeVector> duw_rel_perm = S_->GetFieldData(duw_coef_key_, name_);
  // Epetra_MultiVector& test_face = *duw_rel_perm->ViewComponent("face",false);
  // Epetra_MultiVector& test_grav = *duw_rel_perm->ViewComponent("grav",false);
  // for (int f=0; f<test_face.MyLength(); f++){
  //   if (abs(test_face[0][f] - test_grav[0][f]) > 1e-10){
  //     //std::cout<<"coef "<<test_coef[0][0]<<" coef_grav "<<test_coef_grav[0][0]<<"\n";
  //     std::cout<<"duw face "<<f<<": "<<test_face[0][f]<<" "<<test_grav[0][f]<<"\n";
  //     exit(0);
  //   }
  // }
  
  // debugging
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << " " << update_perm << std::endl;
  }
  return update_perm;
};

} // namespace
} // namespace
