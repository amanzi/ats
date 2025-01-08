/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Evaluates albedos and emissivities in a three-component subgrid model.
/*!

Evaluates the albedo and emissivity as an interpolation on the surface
properties and cover.  This allows for three components -- water/ice, land, and
snow.  Note this internally calculates albedo of snow based upon snow density.

Components are: 0 = land, 1 = ice/water, 2 = snow.

Requires the use of LandCover types, for ground albedo and emissivity.

.. _albedo_evaluator_subgrid-spec:
.. admonition:: albedo_evaluator_subgrid-spec

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
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class AlbedoThreeComponentEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit AlbedoThreeComponentEvaluator(Teuchos::ParameterList& plist);
  AlbedoThreeComponentEvaluator(const AlbedoThreeComponentEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new AlbedoThreeComponentEvaluator(*this));
  }

 protected:
  // custom EC used to set subfield names
  virtual void EnsureCompatibility_Structure_(State& S) override;

  // custom EC used because deps have 1 component not 3
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  Key domain_;
  Key domain_snow_;

  Key albedo_key_, emissivity_key_;
  Key snow_dens_key_;
  Key unfrozen_fraction_key_;

  bool is_constant_snow_albedo_;
  double a_water_, a_ice_, a_snow_;
  double e_water_, e_ice_, e_snow_;

  // this is horrid, because this cannot yet live in state
  // bring on new state!
  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<Evaluator, AlbedoThreeComponentEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
