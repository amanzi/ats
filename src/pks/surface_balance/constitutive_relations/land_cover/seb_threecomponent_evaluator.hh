/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Calculates source terms for surface fluxes to and from the atmosphere and a ground surface.
/*!

Calculates source terms for surface fluxes to and from the atmosphere and a
ground surface characterized by three components -- snow, water-covered ground,
and vegetated/bare soil.

The surface energy balance on these area weighted patches are individually
calculated then averaged to form the total quantities.  All down- and
up-scaling of relevant quantities are done through the area weighting, which is
calculated by a minimum threshold in snow and a depression depth/geometry-based
approach for water.  All snow is assumed to first cover water (likely ice),
then cover land, as both water and snow prefer low-lying depressions due to
gravity- and wind-driven redistributions, respectively.

`"evaluator type`" = `"surface energy balance, two components`"

.. _seb-threecomponent-evaluator-spec:
.. admonition:: seb-threecomponent-evaluator-spec

   * `"wind speed reference height [m]`" ``[double]`` **2.0** Reference height at which
     wind speed is measured.
   * `"minimum wind speed [m s^-1]`" ``[double]`` **1.0** Sets a floor on wind speed for
     potential wierd data.  Models have trouble with no wind.

   * `"save diagnostic data`" ``[bool]`` **false** Saves a suite of diagnostic variables to vis.

   * `"surface domain name`" ``[string]`` **DEFAULT** Default set by parameterlist name.
   * `"subsurface domain name`" ``[string]`` **DEFAULT** Default set relative to surface domain name.
   * `"snow domain name`" ``[string]`` **DEFAULT** Default set relative to surface domain name.

   KEYS:

   - `"surface water source`" **DOMAIN-water_source**  [m s^-1]
   - `"surface energy source`" **DOMAIN-total_energy_source** [MW m^-2]
   - `"subsurface water source`" **DOMAIN-water_source**  [mol s^-1]
   - `"subsurface energy source`" **DOMAIN-total_energy_source** [MW m^-3]
   - `"snow mass source - sink`" **DOMAIN-source_sink** [m_SWE s^-1]
   - `"new snow source`" **DOMAIN-source** [m_SWE s^-1]

   - `"albedo`" **DOMAIN-albedo** [-] A single variate diagnostic of the final albedo.
   - `"snowmelt`" **DOMAIN_SNOW-melt** [m_SWE s^-1]
   - `"evaporation`" **DOMAIN-evaporative_flux** [m s^-1]
   - `"snow temperature`" **DOMAIN_SNOW-temperature** [K]
   - `"sensible heat flux`" **DOMAIN-qE_sensible_heat** [W m^-2]
   - `"latent heat of evaporation`" **DOMAIN-qE_latent_heat** [W m^-2]
   - `"latent heat of snowmelt`" **DOMAIN-qE_snowmelt** [W m^-2]
   - `"outgoing longwave radiation`" **DOMAIN-qE_lw_out** [W m^-2]
   - `"conducted energy flux`" **DOMAIN-qE_conducted** [W m^-2]

   DEPENDENCIES:

   - `"incoming shortwave radiation`" **DOMAIN-incoming_shortwave_radiation** [W m^-2]
   - `"incoming longwave radiation`" **DOMAIN-incoming_longwave_radiation** [W m^-2]
   - `"air temperature`" **DOMAIN-air_temperature** [K]
   - `"vapor pressure air`" **DOMAIN-vapor_pressure_air** [Pa]
   - `"wind speed`" **DOMAIN-wind_speed** [m s^-1]
   - `"precipitation rain`" **DOMAIN-precipitation_rain** [m s^-1]
   - `"precipitation snow`" **DOMAIN_SNOW-precipitation** [m_SWE s^-1]

   - `"snow depth`" **DOMAIN_SNOW-depth** [m]
   - `"snow density`" **DOMAIN_SNOW-density** [kg m^-3]
   - `"snow death rate`" **DOMAIN_SNOW-death_rate** [m s^-1]  Snow "death" refers to the last bit of snowmelt that we want to remove discretely.
   - `"ponded depth`" **DOMAIN-ponded_depth** [m]
   - `"unfrozen fraction`" **DOMAIN-unfrozen_fraction** [-]  1 --> all surface water, 0 --> all surface ice
   - `"subgrid albedos`" **DOMAIN-albedos** [-] Dimension 2 field of (no-snow, snow) albedos.
   - `"subgrid emissivity`" **DOMAIN-emissivities** [-] Dimension 2 field of (no-snow, snow) emissivities.
   - `"area fractions`" **DOMAIN-fractional_areas** Dimension 2 field of (no-snow, snow) area fractions (sum to 1).

   - `"temperature`" **DOMAIN-temperature**  [K] surface skin temperature.
   - `"pressure`" **DOMAIN-pressure** [Pa] surface skin pressure.
   - `"gas saturation`" **DOMAIN_SS-saturation_gas** [-] subsurface gas saturation
   - `"liquid saturation`" **DOMAIN_SS-saturation_liquid** [-] subsurface liquid saturation
   - `"porosity`" [-] subsurface porosity
   - `"subsurface pressure`" **DOMAIN_SS-pressure** [Pa]
   - `"molar density liquid`" **DOMAIN-molar_density_liquid** [mol m^-3]
   - `"mass density liquid`" **DOMAIN-mass_density_liquid** [kg m^-3]


.. note:

   This also depends upon multiple parameters from the LandCover_ types:

   - `"roughness length of bare ground [m]`" ``[double]`` **0.04** Defines a fetch controlling
     latent and sensible heat fluxes.
   - `"roughness length of snow-covered ground [m]`" ``[double]`` **0.004** Defines a
     fetch controlling latent and sensible heat fluxes.
   - `"dessicated zone thickness [m]`" ``[double]`` Thickness of the immediate surface
     layer over which vapor pressure diffusion must move water to evaporate
     from dry soil.  More implies less evaporation.

*/

#pragma once

#include "Factory.hh"
#include "Debugger.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class SEBThreeComponentEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit SEBThreeComponentEvaluator(Teuchos::ParameterList& plist);
  SEBThreeComponentEvaluator(const SEBThreeComponentEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new SEBThreeComponentEvaluator(*this));
  }

 protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // make sure the structure are set up on all diagnostic variables
  virtual void EnsureCompatibility_Structure_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  // this is non-standard practice.  Implementing UpdateFieldDerivative_ to
  // override the default chain rule behavior, instead doing a numerical
  // finite difference
  virtual void UpdateFieldDerivative_(const State& S, const Key& wrt_key, const Tag& wrt_tag);

 protected:
  Key water_source_key_, energy_source_key_;
  Key ss_water_source_key_, ss_energy_source_key_;
  Key snow_source_key_, new_snow_key_;
  Key met_sw_key_, met_lw_key_, met_air_temp_key_, met_vp_air_key_;
  Key met_wind_speed_key_, met_prain_key_, met_psnow_key_;
  Key snow_depth_key_, snow_dens_key_, snow_death_rate_key_;
  Key ponded_depth_key_, unfrozen_fraction_key_;
  Key sg_albedo_key_, sg_emissivity_key_, area_frac_key_;
  Key surf_temp_key_, surf_pres_key_;
  Key sat_gas_key_, sat_liq_key_, poro_key_, ss_pres_key_;
  Key mol_dens_key_, mass_dens_key_;

  Key melt_key_, evap_key_;
  Key snow_temp_key_;
  Key qE_sh_key_, qE_lh_key_, qE_sm_key_, qE_lw_out_key_, qE_cond_key_;
  Key albedo_key_;

  Key domain_;
  Key domain_ss_;
  Key domain_snow_;

  double
    min_wind_speed_; // wind speed of 0, under this model, would have 0 latent or sensible heat?
  double wind_speed_ref_ht_; // reference height of the met data

  LandCoverMap land_cover_;

  bool compatible_;
  bool diagnostics_;
  Teuchos::RCP<Debugger> db_;
  Teuchos::RCP<Debugger> db_ss_;
  Teuchos::ParameterList plist_;

 private:
  static Utils::RegisteredFactory<Evaluator, SEBThreeComponentEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
