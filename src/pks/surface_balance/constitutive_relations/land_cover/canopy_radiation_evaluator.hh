/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates the canopy radiation balance, providing canopy net and radiation to the snow/surface.
/*!

Computes and sums the downward radiation terms, determining the total radiation
sent down to the surface from the canopy and above.

Requires the use of LandCover types, for albedo and emissivity of the canopy
itself, along with Beer's law coefficients.

Computes:

1. canopy-downward_shortwave_radiation -- transmitted shortwave.  Note that
   incoming shortwave is attenuated by Beer's law, and partially transmitted
   without attenuation when there are gaps (e.g. LAI < 1) in the canopy.

2. canopy-downward_longwave_radiation -- transmitted longwave (see above,
   noting that Beer's law coefficients should be used that absorb most if not
   all the longwave radiation), along with longwave emitted by the canopy
   computed using a canopy leaf temperature and a Bolzmann equation.

3. canopy-downward_net_radiation -- this is a partial computation of the net
   radiation experienced by the canopy.  It includes the portion of shortwave
   and longwave from the atmosphere that are absorbed via Beer's law, minus the
   outgoing longwave emitted from the canopy (see downward above).  It does NOT
   include upward longwave radiation emitted by the snow or surface.  It also
   does not include any secondary bounces (e.g. reflected
   atmosphere->canopy->cloud->back to canopy, or transmitted by the canopy,
   reflected by snow).

Here the net radiation is positive for energy added to the canopy, while the
other two are positive for energy sent to the layer below.

In the canopy-downward_net_radiation, we cannot include the upward terms YET,
because these are a function of snow and surface temperature, which in turn
depend upon the downward radiation computed here.  So we choose to break the
loop here, by computing downard terms first, then iterating to compute snow
temperature, then compute upward terms.  The alternative would be to have a
formal snow energy PK that computed snow temperature, at which point we would
solve all of these balances to convergence simultaneously.

`"evaluator type`" = `"canopy radiation balance from above`"

.. _canopy-radiation-evaluator-spec:
.. admonition:: canopy-radiation-evaluator-spec

   KEYS:
   - `"incoming shortwave radiation`" **SURFACE_DOMAIN-incoming_shortwave_radiation**
   - `"incoming longwave radiation`" **SURFACE_DOMAIN-incoming_longwave_radiation**
   - `"canopy temperature`" **CANOPY_DOMAIN-temperature**
   - `"leaf area index`" **CANOPY_DOMAIN-leaf_area_index**

Note that this is a subset of the physics in the "radiation balance evaluator,"
and is therefore mutually exclusive with that model.

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"
#include "seb_physics_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class CanopyRadiationEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit CanopyRadiationEvaluator(Teuchos::ParameterList& plist);
  CanopyRadiationEvaluator(const CanopyRadiationEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new CanopyRadiationEvaluator(*this));
  }

 protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override;
  virtual void EnsureCompatibility_Structure_(State& S) override {
    EnsureCompatibility_StructureSame_(S);
  }

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  Key domain_canopy_;

  Key can_down_sw_key_;
  Key can_down_lw_key_;
  Key rad_bal_can_key_;

  Key sw_in_key_, lw_in_key_;
  Key temp_canopy_key_;
  Key lai_key_;

  bool compatible_;

  // this is horrid, because this cannot yet live in state
  // bring on new state!
  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<Evaluator, CanopyRadiationEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
