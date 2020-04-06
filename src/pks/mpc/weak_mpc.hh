/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Multi process coupler for sequential coupling.

/*!

Noniterative sequential coupling simply calls each PK's AdvanceStep() method in
order.

.. _weak-mpc-spec:
.. admonition:: weak-mpc-spec

    INCLUDES:

    - ``[mpc-spec]`` *Is a* MPC_.

*/


#ifndef PKS_MPC_WEAKMPC_HH_
#define PKS_MPC_WEAKMPC_HH_

#include "PK.hh"
#include "mpc.hh"

namespace Amanzi {

class WeakMPC : public MPC<PK> {

public:

  WeakMPC(Teuchos::ParameterList& FElist,
          const Teuchos::RCP<Teuchos::ParameterList>& plist,
          const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~WeakMPC() = default;

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);

  virtual void set_dt(double dt);

private:
  // factory registration
  static RegisteredPKFactory<WeakMPC> reg_;


};
} // close namespace Amanzi

#endif
