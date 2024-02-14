/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   ATS

   Default base with default implementations of methods for a physical PK.
   ------------------------------------------------------------------------- */

#include "EvaluatorPrimary.hh"
#include "StateDefs.hh"
#include "pk_helpers.hh"
#include "pk_physical_default.hh"

namespace Amanzi {

PK_Physical_Default::PK_Physical_Default(const Comm_ptr_type& comm,
                                         Teuchos::ParameterList& pk_tree,
                                         const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                         const Teuchos::RCP<State>& S)
  : PK(comm, pk_tree, glist, S), PK_Physical(comm, pk_tree, glist, S)
{}

void
PK_Physical_Default::ParseParameterList_()
{
  domain_ = plist_->get<std::string>("domain name", "domain");
  mesh_ = S_->GetMesh(domain_);

  key_ = Keys::readKey(*plist_, domain_, "primary variable");
  PK_Physical::ParseParameterList_();

  // primary variable max change
  max_valid_change_ = plist_->get<double>("max valid change", -1.0);
}


// -----------------------------------------------------------------------------
// Construction of data.
// -----------------------------------------------------------------------------
void
PK_Physical_Default::Setup()
{
  // get the mesh
  mesh_ = S_->GetMesh(domain_);

  // set up the debugger
  Teuchos::RCP<Teuchos::ParameterList> vo_plist = plist_;
  if (plist_->isSublist(name_ + " verbose object")) {
    vo_plist = Teuchos::rcp(new Teuchos::ParameterList(*plist_));
    vo_plist->set("verbose object", plist_->sublist(name_ + " verbose object"));
  }

  db_ = Teuchos::rcp(new Debugger(mesh_, name_, *vo_plist));

  // require primary variable evaluators
  PKHelpers::requireAtNext(key_, tag_next_, *S_, name_);
  PKHelpers::requireAtCurrent(key_, tag_current_, *S_, name_);
};


void
PK_Physical_Default::CommitStep(double t_old, double t_new, const Tag& tag_next)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Commiting state @ " << tag_next << std::endl;

  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;
  PKHelpers::assign(key_, tag_current, tag_next, *S_);
}


void
PK_Physical_Default::FailStep(double t_old, double t_new, const Tag& tag_next)
{
  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;
  PKHelpers::assign(key_, tag_next, tag_current, *S_);
}


// -----------------------------------------------------------------------------
// Ensures the step size is smaller than max_valid_change
// -----------------------------------------------------------------------------
bool
PK_Physical_Default::ValidStep()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Validating time step." << std::endl;

  if (max_valid_change_ > 0.0) {
    const CompositeVector& var_new = S_->Get<CompositeVector>(key_, tag_next_);
    const CompositeVector& var_old = S_->Get<CompositeVector>(key_, tag_current_);
    CompositeVector dvar(var_new);
    dvar.update(-1., var_old, 1.);
    double change = dvar.normInf();
    if (change > max_valid_change_) {
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << "Invalid time step, max primary variable change=" << change
                   << " > limit=" << max_valid_change_ << std::endl;
      return false;
    }
  }
  return true;
}


// -----------------------------------------------------------------------------
//  Marks as changed
// -----------------------------------------------------------------------------
void
PK_Physical_Default::ChangedSolutionPK(const Tag& tag)
{
  PKHelpers::changedEvaluatorPrimary(key_, tag, *S_);
}


// -----------------------------------------------------------------------------
// Initialization of the PK data.
// -----------------------------------------------------------------------------
void
PK_Physical_Default::Initialize()
{
  // Get the record
  Record& record = S_->GetRecordW(key_, tag_next_, getName());

  // Initialize the data
  if (!record.initialized()) {
    // initial conditions
    // -- Get the IC function plist.
    if (!plist_->isSublist("initial condition")) {
      Errors::Message message;
      message << getName() << " has no initial condition parameter list.";
      Exceptions::amanzi_throw(message);
    }

    // -- Calculate the IC.
    Teuchos::ParameterList ic_plist = plist_->sublist("initial condition");
    ic_plist.setName(key_);
    record.Initialize(ic_plist);

    // debug
    db_->WriteVector(key_ + " IC", record.GetPtr<CompositeVector>().ptr(), true);

    // communicate just to make sure values are initialized for valgrind's sake
    record.Get<CompositeVector>().scatterMasterToGhosted();
    ChangedSolutionPK(tag_next_);
  }
};


} // namespace Amanzi
