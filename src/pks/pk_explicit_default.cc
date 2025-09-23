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
Explicit.
------------------------------------------------------------------------- */

#include "Teuchos_TimeMonitor.hpp"
#include "PK.hh"
#include "State.hh"
#include "pk_explicit_default.hh"

namespace Amanzi {
namespace ATS_Physics {    

// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
void
PK_Explicit_Default::Setup()
{
  // initial timestep
  dt_ = plist_->get<double>("initial timestep", 1.);
};


// -----------------------------------------------------------------------------
// Initialization of timestepper.
// -----------------------------------------------------------------------------
void
PK_Explicit_Default::Initialize()
{
  // set up the timestepping algorithm
  if (!plist_->get<bool>("strongly coupled PK", false)) {
    // -- instantiate timestepper
    Teuchos::ParameterList& ti_plist = plist_->sublist("time integrator");
    ti_plist.set("initial time", S_->get_time());

    time_stepper_ = Teuchos::rcp(new Explicit_TI::RK<Amanzi::TreeVector>(*this, ti_plist, *solution_));
    solution_old_ = Teuchos::rcp(new Amanzi::TreeVector(*solution_));
    
  }
};


// -----------------------------------------------------------------------------
// Advance from state S to state S_next at time S.time + dt.
// -----------------------------------------------------------------------------
bool
PK_Explicit_Default::AdvanceStep(double t_old, double t_new, bool reinit)
{
  double dt = t_new - t_old;

  Teuchos::OSTab out = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << S_->get_time(tag_current_)
               << " t1 = " << S_->get_time(tag_next_) << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;

  State_to_Solution(tag_current_, *solution_old_);
  State_to_Solution(tag_next_, *solution_);

  // take a timestep
  time_stepper_->TimeStep(S_->get_time(tag_current_), dt, *solution_old_, *solution_);
  return false;
};

} // namespace ATS_Physics      
} // namespace Amanzi
