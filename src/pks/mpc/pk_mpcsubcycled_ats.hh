/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Class for subcycling a slave step within a master step.
  Assumes that intermediate_time() can be used (i.e. this is not nestable?)

  See additional documentation in the base class src/pks/mpc_pk/PK_MPC.hh
*/

#ifndef ATS_AMANZI_SUBCYCLED_MPC_HH_
#define ATS_AMANZI_SUBCYCLED_MPC_HH_

#include "PK.hh"
#include "mpc.hh"

namespace Amanzi {

class PK_MPCSubcycled_ATS : public MPC<PK> {

public:
  PK_MPCSubcycled_ATS(Teuchos::ParameterList& pk_tree_or_fe_list,
                      const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& soln);

 // Virtual destructor
  virtual ~PK_MPCSubcycled_ATS() = default;

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  virtual void set_dt(double dt) {};

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);

  virtual std::string name() {return "pk_mpcsubcycled_ats";};

 protected:
  int master_;
  int slave_;
  double master_dt_;
  double slave_dt_;
  double min_dt_;
  bool subcycling;

  // states
  Teuchos::RCP<State> S_;

 private:
  // factory registration
  static RegisteredPKFactory<PK_MPCSubcycled_ATS> reg_;
};

}  // namespace Amanzi

#endif
