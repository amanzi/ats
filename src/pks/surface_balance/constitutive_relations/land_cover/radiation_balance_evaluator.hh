/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Evaluates a net radiation balance for ground and canopy.
/*!

Here the net radiation is positive for energy inputs to the layer.  Note that
ground is based on the two-channel (land + snow) while canopy is assumed to be
a simple, single layer.

Requires the use of LandCover types, for canopy albedo and emissivity.

.. _radiation-balance-evaluator-spec:
.. admonition:: radiation-balance-evaluator-spec

   * `"albedo ice [-]`" ``[double]`` **0.44**
   * `"albedo water [-]`" ``[double]`` **0.1168**

   * `"emissivity ice [-]`" ``[double]`` **0.98**
   * `"emissivity water [-]`" ``[double]`` **0.995**
   * `"emissivity snow [-]`" ``[double]`` **0.98**

   KEYS:
   - `"surface albedos`" **SURFACE_DOMAIN-albedos**
   - `"surface emissivities`" **SURFACE_DOMAIN-emissivities**
   - `"incoming shortwave radiation`" **SURFACE_DOMAIN-incoming_shortwave_radiation**
   - `"incoming longwave radiation`" **SURFACE_DOMAIN-incoming_longwave_radiation**
   - `"surface temperature`" **SURFACE_DOMAIN-temperature**
   - `"snow temperature`" **SNOW_DOMAIN-temperature**
   - `"canopy temperature`" **CANOPY_DOMAIN-temperature**
   - `"leaf area index`" **CANOPY_DOMAIN-leaf_area_index**
   - `"area fractions`" **SURFACE_DOMAIN-area_fractions**

*/

#pragma once

#include "Factory.hh"
#include "secondary_variables_field_evaluator.hh"
#include "LandCover.hh"
#include "seb_physics_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class RadiationBalanceEvaluator : public SecondaryVariablesFieldEvaluator {
 public:
  explicit
  RadiationBalanceEvaluator(Teuchos::ParameterList& plist);
  RadiationBalanceEvaluator(const RadiationBalanceEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const override {
    return Teuchos::rcp(new RadiationBalanceEvaluator(*this));
  }

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;

 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results) override;

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results) override;

 protected:
  Key domain_surf_;
  Key domain_snow_;
  Key domain_canopy_;

  Key rad_bal_surf_key_;
  Key rad_bal_snow_key_;
  Key rad_bal_can_key_;

  Key albedo_surf_key_, emissivity_surf_key_;
  Key sw_in_key_, lw_in_key_;
  Key temp_surf_key_, temp_canopy_key_, temp_snow_key_;
  Key area_frac_key_;
  Key lai_key_;

  bool compatible_;

  // this is horrid, because this cannot yet live in state
  // bring on new state!
  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,RadiationBalanceEvaluator> reg_;
};

}  // namespace Relations
}  // namespace SurfaceBalance
}  // namespace Amanzi
