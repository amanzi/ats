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

#ifndef PK_SURFACE_BALANCE_IMPLICIT_HH_
#define PK_SURFACE_BALANCE_IMPLICIT_HH_

#include "pk_factory.hh"
#include "pk_physical_bdf_base.hh"

namespace Amanzi {
namespace SurfaceBalance {

class SurfaceBalanceImplicit : public PKPhysicalBDFBase {

public:
  SurfaceBalanceImplicit(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                         Teuchos::ParameterList& FElist,
                         const Teuchos::RCP<TreeVector>& solution);

  // main methods
  // -- Setup data.
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  virtual void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  // error monitor
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du);

  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {}

  // -- Calculate any diagnostics prior to doing vis
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}


 protected:
  // multiple primary variables
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_esource_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_wsource_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_w_v_source_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_wtemp_;

  bool eval_derivatives_;
  bool implicit_snow_;
  bool modify_predictor_advance_;
  bool longwave_input_;

  double min_wind_speed_;
  double wind_speed_ref_ht_;
  double snow_ground_trans_;
  double min_snow_trans_;
  double snow_lambda_fresh_;
  double snow_lambda_exponent_;

  Teuchos::RCP<const AmanziMesh::Mesh> subsurf_mesh_;

 private:
  // factory registration
  static RegisteredPKFactory<SurfaceBalanceImplicit> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
