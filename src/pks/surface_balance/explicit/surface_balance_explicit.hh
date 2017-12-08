/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

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

#ifndef PK_SURFACE_BALANCE_EXPLICIT_HH_
#define PK_SURFACE_BALANCE_EXPLICIT_HH_

#include "pk_factory.hh"
#include "pk_physical_base.hh"

namespace Amanzi {
namespace SurfaceBalance {

class SurfaceBalanceExplicit : public PKPhysicalBase {

public:
  SurfaceBalanceExplicit(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                         Teuchos::ParameterList& FElist,
                         const Teuchos::RCP<TreeVector>& solution);

  // main methods
  // -- Setup data.
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- provide a timestep size
  virtual double get_dt() {
    return dt_;
  }

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {}

  // -- Calculate any diagnostics prior to doing vis
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}

  // -- advance via one of a few methods
  virtual bool advance(double dt);

 protected:
  // multiple primary variables
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_esource_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_wsource_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_w_v_source_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_wtemp_;

  double dt_;
  double min_wind_speed_;
  double snow_ground_trans_;
  double min_snow_trans_;
  double snow_lambda_fresh_;
  double snow_lambda_exponent_;

  Teuchos::RCP<const AmanziMesh::Mesh> subsurf_mesh_;

 private:
  // factory registration
  static RegisteredPKFactory<SurfaceBalanceExplicit> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
