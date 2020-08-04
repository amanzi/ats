/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  PK for coupling of Flow PK with Transport_PK and Chemestry_PK
*/

#ifndef ATS_AMANZI_FLOWREACTIVETRANSPORT_PK_HH_
#define ATS_AMANZI_FLOWREACTIVETRANSPORT_PK_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "PK.hh"
#include "mpc.hh"
#include "pk_mpcsubcycled_ats.hh"


namespace Amanzi {

class FlowReactiveTransport_PK_ATS : public PK_MPCSubcycled_ATS {
 public:
  FlowReactiveTransport_PK_ATS(Teuchos::ParameterList& pk_tree_or_fe_list,
                               const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                               const Teuchos::RCP<State>& S,
                               const Teuchos::RCP<TreeVector>& soln);

  ~FlowReactiveTransport_PK_ATS() {
  }

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  virtual void set_dt(double dt);

  // -- advance each sub pk from t_old to t_new.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);

  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);

  std::string name() { return "flow reactive transport ATS";} 

 private:

  Teuchos::RCP<Teuchos::Time> flow_timer_;
  Teuchos::RCP<Teuchos::Time> rt_timer_;

  
  // factory registration
  static RegisteredPKFactory<FlowReactiveTransport_PK_ATS> reg_;
};

}  // namespace Amanzi
#endif
