/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
/*!

Evaluates the albedo and emissivity as an interpolation on the surface
properties and cover.  This allows for two components -- snow and not snow
(water/ice/land).  Note this internally calculates albedo of snow based upon
snow density.

Components are indexed by: 0 = land/ice/water, 1 = snow.

`"evaluator type`" = `"subgrid albedos, two components`"

.. _evaluator-subgrid-albedos-two-components-spec:
.. admonition:: evaluator-subgrid-albedos-two-components-spec

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


.. note::

   This evaluator also uses the :ref:`Land Cover` types.  From that struct, it
   requires the value of the following parameters:

   - `"emissivity of bare ground [-]`"
   - `"albedo of bare ground [-]`"

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class AlbedoTwoComponentEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit AlbedoTwoComponentEvaluator(Teuchos::ParameterList& plist);
  AlbedoTwoComponentEvaluator(const AlbedoTwoComponentEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new AlbedoTwoComponentEvaluator(*this));
  }

 protected:
  // custom EC used to set subfield names
  virtual void EnsureCompatibility_Structure_(State& S) override;

  // custom EC used because deps have 1 component not 2
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
  Key snow_dens_key_, ponded_depth_key_, unfrozen_fraction_key_;

  bool is_constant_snow_albedo_;
  double a_ice_, a_water_, a_snow_;
  double e_snow_, e_ice_, e_water_;

  // this is horrid, because this cannot yet live in state
  // bring on new state!
  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<Evaluator, AlbedoTwoComponentEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
