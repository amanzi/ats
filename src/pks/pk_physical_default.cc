/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Default base with default implementations of methods for a physical PK.
   ------------------------------------------------------------------------- */

#include "EvaluatorPrimary.hh"
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
  key_ = Keys::readKey(*plist_, domain_, "primary variable");

  // primary variable max change
  max_valid_change_ = plist_->get<double>("max valid change", -1.0);

  // verbose object
  if (plist_->isSublist(name() + " verbose object"))
    plist_->set("verbose object", plist_->sublist(name() + " verbose object"));
  vo_ = Teuchos::rcp(new VerboseObject(*S->GetMesh(domain_)->get_comm(), name(), *plist_));
}

// -----------------------------------------------------------------------------
// Construction of data.
// -----------------------------------------------------------------------------

void PK_Physical_Default::Setup()
{
  // get the mesh
  mesh_ = S_->GetMesh(domain_);

  // set up the debugger
  db_ = Teuchos::rcp(new Debugger(mesh_, name_, *plist_));

  // require primary variable evaluators
  S_->Require<CompositeVector, CompositeVectorSpace>(key_, tag_next_, name_);
  RequireEvaluatorPrimary_(key_, tag_next_);
  S_->Require<CompositeVector, CompositeVectorSpace>(key_, tag_current_, name_);
  RequireEvaluatorPrimary_(key_, tag_current_);
};


void PK_Physical_Default::CommitStep(double t_old, double t_new, const Tag& tag)
{
  S_->GetW<CompositeVector>(key_, tag_current_, name_) =
    S_->Get<CompositeVector>(key_, tag_next_);
  ChangedEvaluatorPrimary_(key_, tag_current_);
}

void PK_Physical_Default::FailStep(double t_old, double t_new, const Tag& tag)
{
  S_->GetW<CompositeVector>(key_, tag_next_, name_) =
    S_->Get<CompositeVector>(key_, tag_current_);
  ChangedEvaluatorPrimary_(key_, tag_next_);
}


// -----------------------------------------------------------------------------
// Ensures the step size is smaller than max_valid_change
// -----------------------------------------------------------------------------
bool PK_Physical_Default::ValidStep()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Validating time step." << std::endl;

  if (max_valid_change_ > 0.0) {
    const CompositeVector& var_new = S_->Get<CompositeVector>(key_, tag_next_);
    const CompositeVector& var_old = S_->Get<CompositeVector>(key_, tag_current_);
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
void PK_Physical_Default::ChangedSolutionPK(const Tag& tag)
{
  ChangedEvaluatorPrimary_(key_, tag);
}


// -----------------------------------------------------------------------------
// Initialization of the PK data.
// -----------------------------------------------------------------------------
void PK_Physical_Default::Initialize()
{
  // Get the record
  Record& record = S_->GetRecordW(key_, tag_current_, name());

  // Initialize the data
  if (!record.initialized()) {
    // initial conditions
    // -- Get the IC function plist.
    if (!plist_->isSublist("initial condition")) {
      Errors::Message message;
      message << name() << " has no initial condition parameter list.";
      Exceptions::amanzi_throw(message);
    }

    // -- Calculate the IC.
    Teuchos::ParameterList ic_plist = plist_->sublist("initial condition");
    record.Initialize(ic_plist);

    // communicate just to make sure values are initialized for valgrind's sake
    if (record.Get<CompositeVector>().Ghosted())
      record.Get<CompositeVector>().ScatterMasterToGhosted();
    ChangedSolutionPK(tag_next_);
  }

  // Push the data into the solution.
  solution_->SetData(record.GetPtrW<CompositeVector>(name()));
};

Teuchos::RCP<EvaluatorPrimaryCV>
PK_Physical_Default::RequireEvaluatorPrimary_(const Key& key,
        const Tag& tag)
{
  // first check, is there one already
  if (S_->HasEvaluator(key, tag)) {
    // if so, make sure it is primary
    Teuchos::RCP<Evaluator> eval = S_->GetEvaluatorPtr(key, tag);
    Teuchos::RCP<EvaluatorPrimaryCV> eval_pv =
      Teuchos::rcp_dynamic_cast<EvaluatorPrimaryCV>(eval);
    if (eval_pv == Teuchos::null) {
      Errors::Message msg;
      msg << "PK " << name() << " expected primary variable evaluator for "
          << key << " @ " << tag.get();
      Exceptions::amanzi_throw(msg);
    }
    return eval_pv;
  }

  // if not, create one, only at this tag, not to be shared across tags.  By
  // this, we mean we don't stick the "type" = "primary" back into the
  // evaluator list -- this allows "copy evaluators" e.g. "water content at the
  // old tag" to differ from the standard evalulator, e.g. "water content at
  // the new tag" which is likely a secondary variable evaluator.
  Teuchos::ParameterList plist(key);
  plist.set("evaluator type", "primary variable");
  plist.set("tag", tag.get());
  auto eval_pv = Teuchos::rcp(new EvaluatorPrimaryCV(plist));
  S_->SetEvaluator(key, tag, eval_pv);
  return eval_pv;
}

void
PK_Physical_Default::ChangedEvaluatorPrimary_(const Key& key,
        const Tag& tag)
{
  Teuchos::RCP<Evaluator> eval = S_->GetEvaluatorPtr(key, tag);
  Teuchos::RCP<EvaluatorPrimaryCV> eval_pv =
    Teuchos::rcp_dynamic_cast<EvaluatorPrimaryCV>(eval);
  if (eval_pv == Teuchos::null) {
    Errors::Message msg;
    msg << "PK " << name() << " expected primary variable evaluator for "
        << key << " @ " << tag.get();
    Exceptions::amanzi_throw(msg);
  }
  eval_pv->SetChanged();
}


} // namespace
