/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon, Adam Atchley, Satish Karra

   DOCUMENT ME
   Surface Energy Balance for Snow Surface and Ground Surface
   Calculates Energy flux, rate or water, and water temperature
   entering through the surface skin.  Snow surface energy balance
   is calculated at equilibrium with ground/surface water and Air.

   ------------------------------------------------------------------------- */

#ifndef PK_SURFACE_BALANCE_CLM_HH_
#define PK_SURFACE_BALANCE_CLM_HH_

#include "PK_Factory.hh"
#include "pk_physical_bdf_default.hh"

namespace Amanzi {
namespace SurfaceBalance {

class SurfaceBalanceCLM : public PK_Physical_Default {

public:

  SurfaceBalanceCLM(Teuchos::ParameterList& pk_tree,
                         const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                         const Teuchos::RCP<State>& S,
                         const Teuchos::RCP<TreeVector>& solution);

  // main methods
  // -- Setup data.
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S);


  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new,  const Teuchos::RCP<State>& S) {}

  // -- Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) {}

  virtual void set_dt(double dt) {}
  virtual double get_dt() {return dt_;}

  // Advance PK from time t_old to time t_new. True value of the last 
  // parameter indicates drastic change of boundary and/or source terms
  // that may need PK's attention. 
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);
  

 protected:
  Key domain_ss_;
  Teuchos::RCP<const AmanziMesh::Mesh> subsurf_mesh_;

  Teuchos::RCP<EvaluatorPrimary> pvfe_wsource_;
  Teuchos::RCP<EvaluatorPrimary> pvfe_w_sub_source_;

  double dt_;
  double my_next_time_;
  
 private:
  // factory registration
  static RegisteredPKFactory<SurfaceBalanceCLM> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
