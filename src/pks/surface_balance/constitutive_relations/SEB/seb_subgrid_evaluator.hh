/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon (coonet @ ornl.gov)
   
 ------------------------------------------------------------------------- */

//! SubgridEvaluator: evaluates the Surface Energy Balance model on subgrid units.

/*!

Sets up a collection of patches, for portions of the column covered in snow,
ponded water, and vegetated/bare ground.  The surface energy balance on these
area weighted patches are individually calculated then averaged to form the
total quantities.  All down- and up-scaling of relevant quantities are done
through the area weighting, which is calculated by a minimum threshold in snow
and a depression depth/geometry-based approach for water.  All snow is assumed
to first cover water (likely ice), then cover land, as both water and snow
prefer low-lying depressions due to gravity- and wind-driven redistributions,
respectively.

* 

*/


#ifndef SEB_SUBGRID_EVALUATOR_HH_
#define SEB_SUBGRID_EVALUATOR_HH_

#include "Factory.hh"
#include "Debugger.hh"
#include "secondary_variables_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

class SubgridEvaluator : public SecondaryVariablesFieldEvaluator {
 public:
  explicit
  SubgridEvaluator(Teuchos::ParameterList& plist);
  SubgridEvaluator(const SubgridEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new SubgridEvaluator(*this));
  }

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);
  
 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results);

  // this is non-standard practice.  Implementing UpdateFieldDerivative_ to
  // override the default chain rule behavior, instead doing a numerical
  // finite difference
  virtual void UpdateFieldDerivative_(const Teuchos::Ptr<State>& S, Key wrt_key);
  
 protected:
  Key mass_source_key_, energy_source_key_;
  Key ss_mass_source_key_, ss_energy_source_key_;
  Key snow_source_key_, new_snow_key_;
  Key met_sw_key_, met_lw_key_, met_air_temp_key_, met_rel_hum_key_;
  Key met_wind_speed_key_, met_prain_key_, met_psnow_key_;
  Key snow_depth_key_, snow_dens_key_, snow_death_rate_key_;
  Key ponded_depth_key_, unfrozen_fraction_key_;
  Key sg_albedo_key_, sg_emissivity_key_, area_frac_key_;
  Key surf_temp_key_, surf_pres_key_;
  Key sat_gas_key_, poro_key_,ss_pres_key_;

  Key melt_key_, evap_key_;
  Key snow_temp_key_;
  Key qE_sh_key_, qE_lh_key_, qE_sm_key_, qE_lw_out_key_, qE_cond_key_;
  Key albedo_key_;

  Key domain_;
  Key domain_ss_;
  Key domain_snow_;

  double min_wind_speed_;       // wind speed of 0, under this model, would have 0 latent or sensible heat?
  double wind_speed_ref_ht_;    // reference height of the met data
  double roughness_bare_ground_;
  double roughness_snow_covered_ground_; // fetch lengths? Or elevation differences?  Or some other smoothness measure? [m]

  double snow_ground_trans_;    // snow depth at which soil starts to appear
  double min_snow_trans_;       // snow depth at which snow reaches no area coverage

  double dessicated_zone_thickness_; // max thickness of the zone over which
                                     // evaporation dessicates the soil, and
                                     // therefore vapor diffusion must act to
                                     // bring evaporated water to the surface.
                                     // A limiter on evaporation as the water
                                     // table drops below the surface.
  bool ss_topcell_based_evap_;
  bool diagnostics_;
  Teuchos::RCP<Debugger> db_;
  Teuchos::RCP<Debugger> db_ss_;
  Teuchos::ParameterList plist_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator,SubgridEvaluator> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
