/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Default base with default implementations of methods for a physical PK.
   ------------------------------------------------------------------------- */
#include "StateDefs.hh"
#include "pk_physical_default.hh"

namespace Amanzi {

PK_Physical_Default::PK_Physical_Default(Teuchos::ParameterList& pk_tree,
                                         const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                         const Teuchos::RCP<State>& S,
                                         const Teuchos::RCP<TreeVector>& solution) :
    PK(pk_tree, glist, S, solution),
    PK_Physical(pk_tree, glist, S, solution)
{
  domain_ = plist_->get<std::string>("domain name", "domain");
  key_ = Keys::readKey(*plist_, domain_, "primary variable");
  
  // set up the primary variable solution, and its evaluator
  Teuchos::ParameterList& FElist = S->FEList();
  Teuchos::ParameterList& pv_sublist = FElist.sublist(key_);
  pv_sublist.set("evaluator name", key_);
  pv_sublist.set("field evaluator type", "primary variable");

  // primary variable max change
  max_valid_change_ = plist_->get<double>("max valid change", -1.0);

  // verbose object
  if (plist_->isSublist(name_ + " verbose object")) 
    plist_->set("verbose object", plist_->sublist(name_ + " verbose object"));
  vo_ = Teuchos::rcp(new VerboseObject(*S->GetMesh(domain_)->get_comm(), name_, *plist_));
}

// -----------------------------------------------------------------------------
// Construction of data.
// -----------------------------------------------------------------------------

void PK_Physical_Default::Setup(const Teuchos::Ptr<State>& S) {
  //PKDefaultBase::setup(S);

  // get the mesh
  mesh_ = S->GetMesh(domain_);

  // set up the debugger
  db_ = Teuchos::rcp(new Debugger(mesh_, name_, *plist_));

  // require primary variable evaluator
  S->RequireFieldEvaluator(key_);
  Teuchos::RCP<FieldEvaluator> fm = S->GetFieldEvaluator(key_);
  solution_evaluator_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  AMANZI_ASSERT(solution_evaluator_ != Teuchos::null);
};


// -----------------------------------------------------------------------------
// Transfer operators -- ONLY COPIES POINTERS
// -----------------------------------------------------------------------------
void PK_Physical_Default::State_to_Solution(const Teuchos::RCP<State>& S,
        TreeVector& solution) {
  solution.SetData(S->GetFieldData(key_, name_));
};


// -----------------------------------------------------------------------------
// Transfer operators -- ONLY COPIES POINTERS
// -----------------------------------------------------------------------------
void PK_Physical_Default::Solution_to_State(TreeVector& solution,
        const Teuchos::RCP<State>& S) {
  AMANZI_ASSERT(solution.Data() == S->GetFieldData(key_));
  //  S->SetData(key_, name_, solution->Data());
  //  solution_evaluator_->SetFieldAsChanged();
};


void PK_Physical_Default::Solution_to_State(const TreeVector& solution,
        const Teuchos::RCP<State>& S) {
  AMANZI_ASSERT(solution.Data() == S->GetFieldData(key_));
  //  TreeVector* soln_nc_ptr = const_cast<TreeVector*>(&solution);
  //  Solution_to_State(*soln_nc_ptr, S);
};


// -----------------------------------------------------------------------------
// Experimental approach -- we must pull out S_next_'s solution_evaluator_ to
// stay current for ChangedSolution()
// -----------------------------------------------------------------------------
void PK_Physical_Default::set_states(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next) {
  //  PKDefaultBase::set_states(S, S_inter, S_next);

  S_ = S;
  S_inter_ = S_inter;
  S_next_ = S_next;

  // Get the FE and mark it as changed.
  // Note that this is necessary because we need this to point at the
  // FE in S_next_, not the one which we created in S_.
  Teuchos::RCP<FieldEvaluator> fm = S_next->GetFieldEvaluator(key_);
#if ENABLE_DBC
  solution_evaluator_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  AMANZI_ASSERT(solution_evaluator_ != Teuchos::null);
#else
  solution_evaluator_ = Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(fm);
#endif

  solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
};


// -----------------------------------------------------------------------------
// Ensures the step size is smaller than max_valid_change
// -----------------------------------------------------------------------------
bool PK_Physical_Default::ValidStep() {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Validating time step." << std::endl;

  if (max_valid_change_ > 0.0) {
    const CompositeVector& var_new = *S_next_->GetFieldData(key_);
    const CompositeVector& var_old = *S_inter_->GetFieldData(key_);
    CompositeVector dvar(var_new);
    dvar.Update(-1., var_old, 1.);
    double change = 0.;
    dvar.NormInf(&change);
    if (change > max_valid_change_) {
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << "Invalid time step, max primary variable change="
                   << change << " > limit=" << max_valid_change_ << std::endl;
      return false;
    }
  }
  return true;
}


// -----------------------------------------------------------------------------
//  Marks as changed
// -----------------------------------------------------------------------------
void PK_Physical_Default::ChangedSolutionPK(const Teuchos::Ptr<State>& S) {
  assert(S != Teuchos::null);

  Teuchos::RCP<FieldEvaluator> fm = S->GetFieldEvaluator(key_);

  Teuchos::RCP<PrimaryVariableFieldEvaluator> solution_evaluator =
    Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  AMANZI_ASSERT(solution_evaluator != Teuchos::null);
  solution_evaluator->SetFieldAsChanged(S);
}


// -----------------------------------------------------------------------------
// Initialization of the PK data.
// -----------------------------------------------------------------------------
void PK_Physical_Default::Initialize(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<Field> field = S->GetField(key_, name_);

  if (!field->initialized()) {
    // initial conditions
    // -- Get the IC function plist.
    if (!plist_->isSublist("initial condition")) {
      std::stringstream messagestream;
      messagestream << name_ << " has no initial condition parameter list.";
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
    }

    // -- Calculate the IC.
    Teuchos::ParameterList ic_plist = plist_->sublist("initial condition");
    field->Initialize(ic_plist);

    // communicate just to make sure values are initialized for valgrind's sake
    field->GetFieldData()->ScatterMasterToGhosted();
    solution_evaluator_->SetFieldAsChanged(S);
  }

  // -- Push the data into the solution.
  solution_->SetData(field->GetFieldData());
};

} // namespace
