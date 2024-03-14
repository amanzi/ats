/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon, Adam Atchley, Satish Karra
*/

/* -------------------------------------------------------------------------
   ATS

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
  virtual void Setup() override;

  // -- Initialize owned (dependent) variables.
  virtual void Initialize() override;


  // // -- Commit any secondary (dependent) variables.
  // virtual void CommitStep(double t_old, double t_new,  );

  // -- Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Tag& tag) override {}

  virtual void set_dt(double dt) override { AMANZI_ASSERT(std::abs(dt - dt_) < 1.e-4); }
  virtual double get_dt() override { return dt_; }

  // Advance PK from time t_old to time t_new. True value of the last
  // parameter indicates drastic change of boundary and/or source terms
  // that may need PK's attention.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override;

 protected:
  void SetupDependencies_(const Tag& tag);
  void InitializeCLM_(const Tag& tag);
  void InitializePrimaryVariables_(const Tag& tag);

 protected:
  Key domain_ss_;
  Key domain_snow_;
  Key domain_can_;

  double dt_;

  Key surf_water_src_key_;
  Key ss_water_src_key_;

  Key met_sw_key_;
  Key met_lw_key_;
  Key met_air_temp_key_;
  Key met_vp_air_key_;
  Key met_wind_speed_key_;
  Key met_prain_key_;
  Key met_psnow_key_;

  Key qE_lh_key_;
  Key qE_sh_key_;
  Key qE_lw_out_key_;
  Key qE_cond_key_;

  Key snow_swe_key_;
  Key can_wc_key_;
  Key surf_temp_key_;
  Key soil_temp_key_;
  Key can_temp_key_;

  Key pres_key_;
  Key poro_key_;
  Key sl_key_;

  Key sand_frac_key_;
  Key silt_frac_key_;
  Key clay_frac_key_;
  Key color_index_key_;
  Key pft_index_key_;

 private:
  // factory registration
  static RegisteredPKFactory<SurfaceBalanceCLM> reg_;
};

} // namespace SurfaceBalance
} // namespace Amanzi

#endif
