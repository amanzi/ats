/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with default implementations of methods for a PK integrated using
BDF.
------------------------------------------------------------------------- */

#include "Teuchos_TimeMonitor.hpp"
#include "BDF1_TI.hh"
#include "pk_bdf_default.hh"
#include "State.hh"

namespace Amanzi {


// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
void PK_BDF_Default::Setup(const Teuchos::Ptr<State>& S)
{
  // initial timestep
  dt_ = plist_->get<double>("initial time step", 1.);

  // preconditioner assembly
  assemble_preconditioner_ = plist_->get<bool>("assemble preconditioner", true);

  if (!plist_->get<bool>("strongly coupled PK", false)) {
    Teuchos::ParameterList& bdf_plist = plist_->sublist("time integrator");
    // -- check if continuation method
    // -- ETC Note this needs fixed if more than one continuation method used
    if (bdf_plist.isSublist("continuation parameters")) {
      S->RequireScalar("continuation_parameter", name_);
    }
  }
};


// -----------------------------------------------------------------------------
// Initialization of timestepper.
// -----------------------------------------------------------------------------
void PK_BDF_Default::Initialize(const Teuchos::Ptr<State>& S)
{
  // set up the timestepping algorithm
  if (!plist_->get<bool>("strongly coupled PK", false)) {
    // -- instantiate time stepper
    Teuchos::ParameterList& bdf_plist = plist_->sublist("time integrator");
    bdf_plist.set("initial time", S->time());
    if (!bdf_plist.isSublist("verbose object"))
      bdf_plist.set("verbose object", plist_->sublist("verbose object"));
    time_stepper_ = Teuchos::rcp(new BDF1_TI<TreeVector,TreeVectorSpace>(*this, bdf_plist, solution_));

    // initialize continuation parameter if needed.
    if (bdf_plist.isSublist("continuation parameters")) {
      *S->GetScalarData("continuation_parameter", name_) = 1.;
      S->GetField("continuation_parameter", name_)->set_initialized();
    }

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp(new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->SetInitialState(S->time(), solution_, solution_dot);
  }
};

void PK_BDF_Default::ResetTimeStepper(double time)
{
  // -- initialize time derivative
  Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp(new TreeVector(*solution_));
  solution_dot->PutScalar(0.0);

  // -- set initial state
  time_stepper_->SetInitialState(time, solution_, solution_dot);
  return;
}

// -----------------------------------------------------------------------------
// Initialization of timestepper.
// -----------------------------------------------------------------------------
double PK_BDF_Default::get_dt() { return dt_; }

void PK_BDF_Default::set_dt(double dt) { dt_ = dt; }

// -- Commit any secondary (dependent) variables.
void PK_BDF_Default::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) \
{
  double dt = t_new - t_old;
  if (time_stepper_ != Teuchos::null) {
    if (dt <= 0) {
      ResetTimeStepper(t_old);
    } else {
      time_stepper_->CommitSolution(dt, solution_, true);
    }
  }
}

void PK_BDF_Default::set_states(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next)
{
  S_ = S;
  S_inter_ = S_inter;
  S_next_ = S_next;
}


// -----------------------------------------------------------------------------
// Advance from state S to state S_next at time S.time + dt.
// -----------------------------------------------------------------------------
bool PK_BDF_Default::AdvanceStep(double t_old, double t_new, bool reinit)
{
  double dt = t_new -t_old;
  Teuchos::OSTab out = vo_->getOSTab();

  if (vo_->os_OK(Teuchos::VERB_LOW))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << S_inter_->time()
               << " t1 = " << S_next_->time() << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;

  State_to_Solution(S_next_, *solution_);

  // take a bdf timestep
  double dt_solver;
  bool fail;
  if (true) { // this is here simply to create a context for timer,
              // which stops the clock when it is destroyed at the
              // closing brace.
    fail = time_stepper_->TimeStep(dt, dt_solver, solution_);
  }

  if (!fail) {
    // check step validity
    bool valid = ValidStep();
    if (valid) {
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << "successful timestep" << std::endl;
      // update the timestep size
      if (dt_solver < dt_ && dt_solver >= dt) {
        // We took a smaller step than we recommended, and it worked fine (not
        // suprisingly).  Likely this was due to constraints from other PKs or
        // vis.  Do not reduce our recommendation.
      } else {
        dt_ = dt_solver;
      }
    } else {
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << "successful advance, but not valid" << std::endl;
      time_stepper_->CommitSolution(dt_, solution_, valid);
      dt_ = 0.5*dt_;
    }
  } else {
    if (vo_->os_OK(Teuchos::VERB_LOW))
      *vo_->os() << "unsuccessful timestep" << std::endl;
    // take the decreased timestep size
    dt_ = dt_solver;
  }

  return fail;
};


// update the continuation parameter
void PK_BDF_Default::UpdateContinuationParameter(double lambda)
{
  *S_next_->GetScalarData("continuation_parameter", name_) = lambda;
  ChangedSolution();
}

} // namespace
