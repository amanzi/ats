/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the derived WeakMPC class.  Provides only the advance()
method missing from MPC.hh.  In weak coupling, we simply loop over the
sub-PKs, calling their advance() methods and returning failure if any fail.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#include "weak_mpc.hh"

namespace Amanzi {

WeakMPC::WeakMPC(Teuchos::ParameterList& FElist,
                 const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& solution)
    : PK(FElist, plist, S, solution),
      MPC<PK>(FElist, plist, S, solution) {
  MPC<PK>::init_(S, solution_->Comm());
};


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double WeakMPC::get_dt() {
  double dt = 1.0e99;
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    dt = std::min<double>(dt, (*pk)->get_dt());
  }
  return dt;
};

// -----------------------------------------------------------------------------
// Set timestep for sub PKs 
// -----------------------------------------------------------------------------
void WeakMPC::set_dt( double dt) {
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_dt(dt);
  }

};

// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool WeakMPC::AdvanceStep(double t_old, double t_new, bool reinit) {
  bool fail = false;
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    fail = (*pk)->AdvanceStep(t_old, t_new, reinit);
    if (fail) {
      return fail;
    }
  }
  return fail;
};


} // namespace
