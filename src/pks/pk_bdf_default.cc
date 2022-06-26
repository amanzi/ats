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
void PK_BDF_Default::Setup()
{
  // initial timestep
  dt_ = plist_->get<double>("initial time step", 1.);
  dt_ = plist_->get<double>("initial time step [s]", dt_);

  // preconditioner assembly
  assemble_preconditioner_ = plist_->get<bool>("assemble preconditioner", true);

  if (!plist_->get<bool>("strongly coupled PK", false)) {
    Teuchos::ParameterList& bdf_plist = plist_->sublist("time integrator");
    // -- check if continuation method
    // -- ETC Note this needs fixed if more than one continuation method used
    if (bdf_plist.isSublist("continuation parameters")) {
      S_->Require<double>("continuation_parameter", Tags::DEFAULT, name_);
    }
  }
};


// -----------------------------------------------------------------------------
// Initialization of timestepper.
// -----------------------------------------------------------------------------
void PK_BDF_Default::Initialize()
{
  // set up the timestepping algorithm
  if (!plist_->get<bool>("strongly coupled PK", false)) {
    // -- instantiate time stepper
    Teuchos::ParameterList& bdf_plist = plist_->sublist("time integrator");
    bdf_plist.set("initial time", S_->get_time());
    if (!bdf_plist.isSublist("verbose object"))
      bdf_plist.set("verbose object", plist_->sublist("verbose object"));
    time_stepper_ = Teuchos::rcp(new BDF1_TI<TreeVector,TreeVectorSpace>(*this,
            bdf_plist, solution_));

    // initialize continuation parameter if needed.
    if (bdf_plist.isSublist("continuation parameters")) {
      S_->GetW<double>("continuation_parameter", Tags::DEFAULT, name_) = 1.;
      S_->GetRecordW("continuation_parameter", Tags::DEFAULT, name_).set_initialized();
    }

    // -- initialize time derivative
    auto solution_dot = Teuchos::rcp(new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->SetInitialState(S_->get_time(), solution_, solution_dot);
  }
};

void PK_BDF_Default::ResetTimeStepper(double time)
{
  // -- initialize time derivative
  auto solution_dot = Teuchos::rcp(new TreeVector(*solution_));
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
void PK_BDF_Default::CommitStep(double t_old, double t_new, const Tag& tag) \
{
  if (tag == tag_next_) {
    double dt = t_new - t_old;
    if (time_stepper_ != Teuchos::null) {
      if (dt <= 0) {
        ResetTimeStepper(t_old);
      } else {
        time_stepper_->CommitSolution(dt, solution_, true);
      }
    }
  }
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
               << "Advancing: t0 = " << S_->get_time(tag_current_)
               << " t1 = " << S_->get_time(tag_next_) << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;

  State_to_Solution(Tags::NEXT, *solution_);

  // take a bdf timestep
  double dt_solver;
  bool fail = time_stepper_->TimeStep(dt, dt_solver, solution_);

  if (!fail) {
    // check step validity
    bool valid = ValidStep();
    if (valid) {
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << "successful advance" << std::endl;
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
      // when including Valid here, make fail = true refs #110
    }
  } else {
    if (vo_->os_OK(Teuchos::VERB_LOW))
      *vo_->os() << "unsuccessful advance" << std::endl;
    // take the decreased timestep size
    dt_ = dt_solver;
  }

  return fail;
};


// update the continuation parameter
void PK_BDF_Default::UpdateContinuationParameter(double lambda)
{
  S_->GetW<double>("continuation_parameter", Tags::DEFAULT, name()) = lambda;
  ChangedSolution();
}

} // namespace
