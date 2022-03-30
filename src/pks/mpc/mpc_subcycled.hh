/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

/*
  MPC for subcycling one PK relative to another.

  * `"subcycling`" ``[bool]`` **True** Whether to subcycle or not.
  * `"subcycled PK index`" ``[int]`` **1** Index, in the PK list, of the
    subcyled PK.
  * `"minimum subcycled relative dt`" ``[double]`` **1.e-5** Sets the minimum
    time step size of the subcycled PK, as a multiple of the standard PK's
    timestep size.
    

*/

#ifndef ATS_AMANZI_SUBCYCLED_MPC_HH_
#define ATS_AMANZI_SUBCYCLED_MPC_HH_

#include "PK.hh"
#include "mpc.hh"

namespace Amanzi {

class MPCSubcycled : public MPC<PK> {

public:
  MPCSubcycled(Teuchos::ParameterList& pk_tree_or_fe_list,
                      const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  virtual void set_dt(double dt);

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);

 protected:
  int standard_;
  int subcycled_;
  double standard_dt_;
  double subcycled_dt_;
  double min_dt_;
  bool subcycling_;

  // states
  Teuchos::RCP<State> S_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCSubcycled> reg_;
};

}  // namespace Amanzi

#endif
