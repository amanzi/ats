/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
//! Evaluates albedos and emissivities in a two-component subgrid model.
/*!

Evaluates the albedo and emissivity as an interpolation on the surface
properties and cover.  This allows for two components -- snow and not snow
(water/ice/land).  Note this internally calculates albedo of snow based upon
snow density.

Components are indexed by: 0 = land/ice/water, 1 = snow.

Requires the use of LandCover types, for ground albedo and emissivity.

.. _albedo-evaluator-spec:
.. admonition:: albedo-evaluator-spec

   * `"albedo ice [-]`" ``[double]`` **0.44**
   * `"albedo water [-]`" ``[double]`` **0.1168**

   * `"emissivity ice [-]`" ``[double]`` **0.98**
   * `"emissivity water [-]`" ``[double]`` **0.995**
   * `"emissivity snow [-]`" ``[double]`` **0.98**

   KEYS:
   - `"subgrid albedos`" **DOMAIN-subgrid_albedos**
   - `"subgrid emissivities`" **DOMAIN-subgrid_emissivities**

   DEPENDENCIES:
   - `"snow density`" **SNOW_DOMAIN-density**
   - `"ponded depth`" **DOMAIN-ponded_depth**
   - `"unfrozen fraction`" **DOMAIN-unfrozen_fraction**

*/

#pragma once

#include "Factory.hh"
#include "secondary_variables_field_evaluator.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class AlbedoTwoComponentEvaluator : public SecondaryVariablesFieldEvaluator {
 public:
  explicit
  AlbedoTwoComponentEvaluator(Teuchos::ParameterList& plist);
  AlbedoTwoComponentEvaluator(const AlbedoTwoComponentEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const override {
    return Teuchos::rcp(new AlbedoTwoComponentEvaluator(*this));
  }

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;

 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results) override;

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results) override;

 protected:
  Key domain_;
  Key domain_snow_;

  Key albedo_key_, emissivity_key_;
  Key snow_dens_key_, ponded_depth_key_, unfrozen_fraction_key_;

  double a_ice_, a_water_;
  double e_snow_, e_ice_, e_water_;

  // this is horrid, because this cannot yet live in state
  // bring on new state!
  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,AlbedoTwoComponentEvaluator> reg_;
};

}  // namespace Relations
}  // namespace SurfaceBalance
}  // namespace Amanzi
