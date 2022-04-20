/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Evaluates albedos and emissivities in a three-component subgrid model.
/*!

Evaluates the albedo and emissivity as an interpolation on the surface
properties and cover.  This allows for three components -- water/ice, land, and
snow.  Note this internally calculates albedo of snow based upon snow density.

Components are: 0 = land, 1 = ice/water, 2 = snow.

Requires the use of LandCover types, for ground albedo and emissivity.

.. _albedo-evaluator-subgrid-spec:
.. admonition:: albedo-evaluator-subgrid-spec

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
   - `"unfrozen fraction`" **DOMAIN-unfrozen_fraction**

*/

#pragma once

#include "Factory.hh"
#include "secondary_variables_field_evaluator.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class AlbedoThreeComponentEvaluator : public SecondaryVariablesFieldEvaluator {
 public:
  explicit
  AlbedoThreeComponentEvaluator(Teuchos::ParameterList& plist);
  AlbedoThreeComponentEvaluator(const AlbedoThreeComponentEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new AlbedoThreeComponentEvaluator(*this));
  }

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results);

 protected:
  Key domain_;
  Key domain_snow_;

  Key albedo_key_, emissivity_key_;
  Key snow_dens_key_;
  Key unfrozen_fraction_key_;

  double a_water_, a_ice_;
  double e_water_, e_ice_, e_snow_;

  // this is horrid, because this cannot yet live in state
  // bring on new state!
  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,AlbedoThreeComponentEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi


