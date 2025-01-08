/*
  Copyright 2010-202x held jointly by participating institutions.
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

This evaluator requires that the surface temperature, snow temperature, and
canopy temperature are known, or at least being solved for.

Requires the use of LandCover types, for albedo and Beer's law coefficients.

This is combination of CLM v4.5 Tech Note and Beer's law for attenuation of
radiation absorption.  In particular, long-wave is exactly as Figure 4.1c in CLM
4.5 Tech Note.  The main difference comes in how absorptivity (which is equal
to emissivity, epsilon in that document) is defined.  Here we use Beer's law
which is an exponential decay with LAI.

Unlike CLM 4.5, here we do not split shortwave into direct and diffuse light.

Computes:

1. "surface radiation balance" -- Net radiation seen by the bare soil/ponded
   water, this includes radiation transmitted to the surface through the
   canopy, longwave emitted by the canopy, and less the longwave emitted by the
   surface itself.  [W m^-2] of actual area -- this does NOT include the
   surface area fraction factor which would be required to compute a total
   energy flux in W.

2. "snow radiation balance" -- Net radiation seen by the snow.  See surface
   above -- all are the same except using snow properties. [W m^-2]

3. "canopy radiation balance" -- this is a compute computation of the net
   radiation experienced by the canopy.  It includes the portion of shortwave
   and longwave from the atmosphere that are absorbed via Beer's law, minus the
   outgoing longwave emitted from the canopy, plus upward longwave radiation
   emitted by the snow and surface.  It also does not include any secondary
   bounces (e.g. reflected atmosphere->canopy->cloud->back to canopy, or
   transmitted by the canopy, reflected by snow/surface.

Requires the use of LandCover types, for canopy albedo and Beer's law
coefficients.

`"evaluator type`" = `"radiation balance, surface and canopy`"

.. _radiation_balance_evaluator-spec:
.. admonition:: radiation_balance_evaluator-spec

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

Note that this is a superset of the physics in the "canopy radiation
evaluator," and is therefore mutually exclusive with that model.

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"
#include "seb_physics_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class RadiationBalanceEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit RadiationBalanceEvaluator(Teuchos::ParameterList& plist);
  RadiationBalanceEvaluator(const RadiationBalanceEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new RadiationBalanceEvaluator(*this));
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
  static Utils::RegisteredFactory<Evaluator, RadiationBalanceEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
