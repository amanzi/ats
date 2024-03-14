/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Distributes and downregulates potential transpiration to the rooting zone.
/*!

The transpiration distribution evaluator looks to take a potential
evapotranspiration and distribute it across the vertical column based on water
availability and rooting depths.  It also potentially limits the transpiration
to avoid taking water where it is not available (thereby crashing the code).

This model is based off of versions from both CLM 4.5 and PRMS.  It requires:

1. A root distribution profile.
2. A plant wilting factor (e.g. how water stressed is the plant?)
3. A potential transpiration, typically calculated from data or a potential
   difference based on a latent heat calculation.

Note this also requires columnar meshes -- meaning that the subsurface mesh
must have `"build columns from set`" provided.

A normalized fraction of water is calculated through multiplying the water
factor by the root distribution factor, integrating across the column, and
dividing by the integral.  This gives a factor which sums to 1 and can be used
to distribute the potential ET throughout the soil column.

Then, this potential ET is down-regulated and multiplied by the plant wilting
factor.  If there is no water locally, it cannot be taken.  Note that almost
always, if there is no water, this did not contribute (much) to the integral
and so is already small.  But if the entire root zone is dried out, this might
have been a bunch of small numbers which then got normalized to 1, meaning they
are all now significant.

Finally, transpiration may be turned off for winter -- relative to time zero,
parameters `"leaf on doy`" and `"leaf off doy`" are used to control when ET
is zero.  By default these are set to 0 and 366 days, ensuring that
transpiration is never turned off and the potential ET must include this
factor.  This is the right choice for, e.g. ELM output, or eddy covariance flux
tower data (where leaf on and off are already included in the potential
calculation).  It is the wrong choice for, e.g. Priestly-Taylor or
Penmann-Montief models, which are happy to predict transpiration all winter
long.  Good choices for those models depend upon the local climate, but may be
something like Julian day 101 for leaf on and Julian day 254 for leaf off (PRMS
defaults for US temperate forests).

Note that `"leaf on doy`" and `"leaf off doy`" are relative to the simulation's
zero time, not the start time.  Typically these are Julian day of the year, but
this assumes that the 0 time of the simulation (not the "start time", but time
0!) is Jan 1.  This leaf on/off cycle is modulo the `"year duration`"
(typically 1 noleap).  Note if `"leaf off doy`" < `"leaf on time`" is ok too --
this is the case if simulation time zero is mid-summer.  These parameters come
from the LandCover type.

.. _transpiration-distribution-evaluator-spec:
.. admonition:: transpiration-distribution-evaluator-spec
   * `"year duration`" ``[double]`` **1**
   * `"year duration units`" ``[string]`` **noleap**

   * `"water limiter function`" ``[function-spec]`` **optional** If provided,
     limit the total water sink as a function of the integral of the water
     potential * rooting fraction.

   KEYS:

   - `"plant wilting factor`" **DOMAIN-plant_wilting_factor**
   - `"rooting depth fraction`" **DOMAIN-rooting_depth_fraction**
   - `"potential transpiration`" **DOMAIN_SURF-potential_transpiration**
   - `"cell volume`" **DOMAIN-cell_volume**
   - `"surface cell volume`" **DOMAIN_SURF-cell_volume**

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {

class Function;

namespace SurfaceBalance {
namespace Relations {

class TranspirationDistributionEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit TranspirationDistributionEvaluator(Teuchos::ParameterList& plist);
  TranspirationDistributionEvaluator(const TranspirationDistributionEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    // calculate of derivatives of this is a tricky thing to do, with
    // non-cell-local terms due to rescaling.  Just turn off derivatives
    // instead.
    return false;
  }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  // need a custom EnsureCompatibility as some vectors cross meshes.
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  void InitializeFromPlist_();
  bool TranspirationPeriod_(double time, double leaf_on_doy, double leaf_off_doy);

 protected:
  Key domain_surf_;
  Key domain_sub_;

  Key f_wp_key_;
  Key f_root_key_;
  Key potential_trans_key_;
  Key cv_key_;
  Key surf_cv_key_;

  double year_duration_;

  LandCoverMap land_cover_;

  bool limiter_local_;
  Teuchos::RCP<Function> limiter_;

 private:
  static Utils::RegisteredFactory<Evaluator, TranspirationDistributionEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
