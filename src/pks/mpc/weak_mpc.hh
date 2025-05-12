/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

Noniterative sequential coupling simply calls each PK's AdvanceStep() method in
order.

`"PK type`" = `"weak MPC`"

.. _pk-weak-mpc-spec:
.. admonition:: pk-weak-mpc-spec

   INCLUDES:

   - ``[mpc-spec]`` *Is a* :ref:`MPC`.

*/


#pragma once

#include "PK.hh"
#include "mpc.hh"

namespace Amanzi {

class WeakMPC : public MPC<PK> {
 public:
  WeakMPC(Teuchos::ParameterList& pk_tree,
          const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
          const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& solution);

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt() override;

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override;

  virtual void set_dt(double dt) override;

 private:
  // factory registration
  static RegisteredPKFactory<WeakMPC> reg_;
};
} // namespace Amanzi
