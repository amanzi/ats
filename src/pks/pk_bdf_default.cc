/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Default base with default implementations of methods for a PK integrated using
BDF.
------------------------------------------------------------------------- */

#include "Teuchos_TimeMonitor.hpp"

#include "Event.hh"
#include "State.hh"
#include "BDF1_TI.hh"
#include "pk_bdf_default.hh"

namespace Amanzi {


// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
void
PK_BDF_Default::Setup()
{
  // preconditioner assembly
  assemble_preconditioner_ = plist_->get<bool>("assemble preconditioner", true);
  strongly_coupled_ = plist_->get<bool>("strongly coupled PK", false);


  if (!strongly_coupled_) {
    Teuchos::ParameterList& bdf_plist = plist_->sublist("time integrator");

    // check if continuation method and require continuation parameter
    // -- ETC Note this needs fixed if more than one continuation method used
    if (bdf_plist.isSublist("continuation parameters")) {
      S_->Require<double>(
        Keys::cleanName(name_, true) + "_continuation_parameter", Tags::DEFAULT, name_);
    }
  }
};


// -----------------------------------------------------------------------------
// Initialization of timestepper.
// -----------------------------------------------------------------------------
void
PK_BDF_Default::Initialize()
{
  if (!strongly_coupled_) {
    // -- initialize continuation parameter if needed.
    if (S_->HasRecord(Keys::cleanName(name_, true) + "_continuation_parameter", Tags::DEFAULT)) {
      S_->Assign(
        Keys::cleanName(name_, true) + "_continuation_parameter", Tags::DEFAULT, name_, (double)1.);
      S_->GetRecordW(Keys::cleanName(name_, true) + "_continuation_parameter", Tags::DEFAULT, name_)
        .set_initialized();
    }

    // -- initialize time derivative
    auto solution_dot = Teuchos::rcp(new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // set up the timestepping algorithm -- note this is done now because the
    // solution space is not known until after Setup() is complete.
    // -- construct the time integrator
    Teuchos::ParameterList& bdf_plist = plist_->sublist("time integrator");
    time_stepper_ = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(
      name() + "_TI", bdf_plist, *this, solution_->get_map(), S_));

    // -- set initial state
    time_stepper_->SetInitialState(S_->get_time(), solution_, solution_dot);
  }
};


// -----------------------------------------------------------------------------
// Initialization of timestepper.
// -----------------------------------------------------------------------------
double
PK_BDF_Default::get_dt()
{
  if (dt_next_ < 0.) dt_next_ = time_stepper_->initial_timestep();
  return dt_next_;
}

void
PK_BDF_Default::set_dt(double dt)
{
  if (!strongly_coupled_) dt_next_ = dt;
}

// -- Commit any secondary (dependent) variables.
void
PK_BDF_Default::CommitStep(double t_old, double t_new, const Tag& tag)
{
  if (tag == tag_next_) {
    double dt = t_new - t_old;
    if (time_stepper_ != Teuchos::null && dt > 0) {
      time_stepper_->CommitSolution(dt, solution_);
    }
  }
}


// -----------------------------------------------------------------------------
// Advance from state S to state S_next at time S.time + dt.
// -----------------------------------------------------------------------------
bool
PK_BDF_Default::AdvanceStep(double t_old, double t_new, bool reinit)
{
  double dt = t_new - t_old;
  Teuchos::OSTab out = vo_->getOSTab();

  if (vo_->os_OK(Teuchos::VERB_LOW))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << t_old << " t1 = " << t_new << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;

  AMANZI_ASSERT(std::abs(S_->get_time(tag_current_) - t_old) < 1.e-4);
  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t_new) < 1.e-4);
  State_to_Solution(Tags::NEXT, *solution_);

  // take a bdf timestep
  bool fail = false;
  try {
    fail = time_stepper_->AdvanceStep(dt, dt_next_, solution_);
  } catch (Errors::TimestepCrash& e) {
    // inject more information into the crash message
    std::stringstream msg_str;
    msg_str << "TimestepCrash in PK: \"" << name() << "\"" << std::endl
            << "  at t = " << t_old << " with dt = " << dt << std::endl
            << "  error message: " << std::endl
            << std::endl
            << e.what() << std::endl
            << std::endl;
    Errors::TimestepCrash msg(msg_str.str());
    Exceptions::amanzi_throw(msg);
  }
  return fail;
};


// update the continuation parameter
void
PK_BDF_Default::UpdateContinuationParameter(double lambda)
{
  S_->Assign(
    Keys::cleanName(name_, true) + "_continuation_parameter", Tags::DEFAULT, name_, lambda);
  ChangedSolution();
}

} // namespace Amanzi
