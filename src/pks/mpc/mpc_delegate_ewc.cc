/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for EWC, a helper class that does projections and preconditioners in
energy/water-content space instead of temperature/pressure space.
------------------------------------------------------------------------- */
#include "FieldEvaluator.hh"
#include "ewc_model.hh"
#include "mpc_delegate_ewc.hh"

namespace Amanzi {

#define DEBUG_FLAG 1

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
MPCDelegateEWC::MPCDelegateEWC(Teuchos::ParameterList& plist) :
    plist_(Teuchos::rcpFromRef(plist)) {
  // set up the VerboseObject
  std::string name = plist_->get<std::string>("PK name")+std::string(" EWC");
  vo_ = Teuchos::rcp(new VerboseObject(name, *plist_));
}

// -----------------------------------------------------------------------------
// Allocate any data or models required.
// -----------------------------------------------------------------------------
void MPCDelegateEWC::setup(const Teuchos::Ptr<State>& S) {
  // Verbosity
  std::string name = plist_->get<std::string>("PK name")+std::string(" EWC");

  // Get the mesh
  Key domain = plist_->get<std::string>("domain name", "");

  if (domain.size() != 0) {
    mesh_ = S->GetMesh(domain);
  } else {
    mesh_ = S->GetMesh("domain");
  }

  // set up a debugger
  db_ = Teuchos::rcp(new Debugger(mesh_, name, *plist_));

  // Process the parameter list for data Keys
  pres_key_ = Keys::readKey(*plist_, domain, "pressure", "pressure");
  temp_key_ = Keys::readKey(*plist_, domain, "temperature", "temperature");

  e_key_ = Keys::readKey(*plist_, domain, "energy", "energy");
  wc_key_ = Keys::readKey(*plist_, domain, "water content", "water_content");
  cv_key_ = Keys::readKey(*plist_, domain, "cell volume", "cell_volume");

  // Process the parameter list for methods
  std::string precon_string = plist_->get<std::string>("preconditioner type", "none");
  if (precon_string == "none") {
    precon_type_ = PRECON_NONE;
  } else if (precon_string == "ewc") {
    precon_type_ = PRECON_EWC;
    AMANZI_ASSERT(0);
  } else if (precon_string == "smart ewc") {
    precon_type_ = PRECON_SMART_EWC;
  } else {
    Errors::Message message;
    message << "EWC Delegate: invalid preconditioner string: \"" << precon_string
            << "\", valid are \"none\", \"ewc\", \"smart ewc\"";
    Exceptions::amanzi_throw(message);
  }

  // select the method used for nonlinear prediction
  std::string predictor_string = plist_->get<std::string>("predictor type", "none");
  if (predictor_string == "none") {
    predictor_type_ = PREDICTOR_NONE;
  } else if (predictor_string == "ewc") {
    predictor_type_ = PREDICTOR_EWC;
    AMANZI_ASSERT(0);
  } else if (predictor_string == "smart ewc") {
    predictor_type_ = PREDICTOR_SMART_EWC;
  } else {
    Errors::Message message;
    message << "EWC Delegate: invalid predictor string: \"" << predictor_string
            << "\", valid are \"none\", \"ewc\", \"smart ewc\"";
    Exceptions::amanzi_throw(message);
  }

  // Smart EWC uses a heuristic to guess when we need the EWC instead of using
  // it blindly.
  if (predictor_type_ == PREDICTOR_SMART_EWC || precon_type_ == PRECON_SMART_EWC) {
    if (plist_->isParameter("freeze-thaw cusp width [K]")) {
      cusp_size_T_freezing_ = plist_->get<double>("freeze-thaw cusp width [K]");
      cusp_size_T_thawing_ = cusp_size_T_freezing_;
    } else {
      cusp_size_T_freezing_ = plist_->get<double>("freeze-thaw cusp width (freezing) [K]", 0.);
      cusp_size_T_thawing_ = plist_->get<double>("freeze-thaw cusp width (thawing) [K]", 0.);
    }
  }
}


// -----------------------------------------------------------------------------
// Initialize any data required in setup.
// -----------------------------------------------------------------------------
void MPCDelegateEWC::initialize(const Teuchos::Ptr<State>& S) {
  // Create and initialize old stored data for previous steps.

  if (predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_SMART_EWC) {
    const Epetra_MultiVector& wc = *S->GetFieldData(wc_key_)
        ->ViewComponent("cell",false);
    const Epetra_MultiVector& e = *S->GetFieldData(e_key_)
        ->ViewComponent("cell",false);

    wc_prev2_ = Teuchos::rcp(new Epetra_MultiVector(wc));
    e_prev2_ = Teuchos::rcp(new Epetra_MultiVector(e));
    wc_prev2_->PutScalar(0.);
    e_prev2_->PutScalar(0.);

    time_prev2_ = S->time();
  }

  // initialize the Jacobian
  if (precon_type_ == PRECON_EWC || precon_type_ == PRECON_SMART_EWC) {
    int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    jac_.resize(ncells, WhetStone::Tensor(2,2));
  }

  // initialize the model, which grabs all needed models from state
  model_->InitializeModel(S, *plist_);
}


// -----------------------------------------------------------------------------
// Set state pointers.
// -----------------------------------------------------------------------------
void MPCDelegateEWC::set_states(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next) {
  S_= S;
  S_inter_ = S_inter;
  S_next_ = S_next;
}


// -----------------------------------------------------------------------------
// Save info from previous iterations if needed.
// -----------------------------------------------------------------------------
void MPCDelegateEWC::commit_state(double dt, const Teuchos::RCP<State>& S) {
  if ((predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_SMART_EWC)
      && S_inter_ != Teuchos::null) {
    // stash water content and energy in S_work.
    *wc_prev2_ = *S_inter_->GetFieldData(wc_key_)->ViewComponent("cell",false);
    *e_prev2_ = *S_inter_->GetFieldData(e_key_)->ViewComponent("cell",false);
    time_prev2_ = S_inter_->time();
  }
}


// -----------------------------------------------------------------------------
// Modify the prediction from linearization of the time integration.
// -----------------------------------------------------------------------------
bool MPCDelegateEWC::ModifyPredictor(double h, Teuchos::RCP<TreeVector> up) {
  bool modified = false;
  double dt_prev = S_inter_->time() - time_prev2_;

  if (predictor_type_ == PREDICTOR_EWC) {
    if (dt_prev > 0.) {
      AMANZI_ASSERT(0);
      //      modified = modify_predictor_ewc_(h,up);
    }
  } else if (predictor_type_ == PREDICTOR_SMART_EWC) {
    if (dt_prev > 0.) {
      modified = modify_predictor_smart_ewc_(h,up);
    }
  }
  return modified;
}


// -----------------------------------------------------------------------------
// Update of the preconditioner
// -----------------------------------------------------------------------------
void MPCDelegateEWC::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
  if (precon_type_ == PRECON_EWC || precon_type_ == PRECON_SMART_EWC) {
    update_precon_ewc_(t,up,h);
  }
}


// -----------------------------------------------------------------------------
// Application of the preconditioner.
// -----------------------------------------------------------------------------
int MPCDelegateEWC::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  int ierr = 0;
  if ((precon_type_ == PRECON_EWC) || (precon_type_ == PRECON_SMART_EWC)) {
    precon_ewc_(u,Pu);
  }
  else ierr = 1;
  
  return ierr;
}



void MPCDelegateEWC::update_precon_ewc_(double t, Teuchos::RCP<const TreeVector> up, double h) {
  Key dedT_key = std::string("d")+e_key_+std::string("_d")+temp_key_;
  S_next_->GetFieldEvaluator(e_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), "ewc", temp_key_);
  const Epetra_MultiVector& dedT = *S_next_->GetFieldData(dedT_key)
      ->ViewComponent("cell",false);

  Key dedp_key = std::string("d")+e_key_+std::string("_d")+pres_key_;
  S_next_->GetFieldEvaluator(e_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), "ewc", pres_key_);
  const Epetra_MultiVector& dedp = *S_next_->GetFieldData(dedp_key)
      ->ViewComponent("cell",false);

  Key dwcdT_key = std::string("d")+wc_key_+std::string("_d")+temp_key_;
  S_next_->GetFieldEvaluator(wc_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), "ewc", temp_key_);
  const Epetra_MultiVector& dwcdT = *S_next_->GetFieldData(dwcdT_key)
      ->ViewComponent("cell",false);

  Key dwcdp_key = std::string("d")+wc_key_+std::string("_d")+pres_key_;
  S_next_->GetFieldEvaluator(wc_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), "ewc", pres_key_);
  const Epetra_MultiVector& dwcdp = *S_next_->GetFieldData(dwcdp_key)
      ->ViewComponent("cell",false);

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c=0; c!=ncells; ++c) {
    jac_[c](0,0) = dwcdp[0][c];
    jac_[c](0,1) = dwcdT[0][c];
    jac_[c](1,0) = dedp[0][c];
    jac_[c](1,1) = dedT[0][c];
  }
}

} // namespace
