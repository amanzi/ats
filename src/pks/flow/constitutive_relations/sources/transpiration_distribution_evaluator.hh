/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! TranspirationDistributionEvaluator distributes potential ET across the rooting zone.
/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

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
parameters `"leaf on time`" and `"leaf off time`" are used to control when ET
is zero.  By default these are set to 0 and 366 days, ensuring that
transpiration is never turned off and the potential ET must include this
factor.  This is the right choice for, e.g. ELM output, or eddy covariance flux
tower data (where leaf on and off are already included in the potential
calculation).  It is the wrong choice for, e.g. Priestly-Taylor or
Penmann-Montief models, which are happy to predict transpiration all winter
long.  Good choices for those models depend upon the local climate, but may be
something like Julian day 101 for leaf on and Julian day 254 for leaf off (PRMS
defaults for US temperate forests).

Note that `"leaf on time`" and `"leaf off time`" are relative to the simulation
start time.  If the run starts on Jan 1, these are Julain day.  This leaf
on/off cycle is modulo the `"year duration`" (typically either 1 noleap or 1
year).  If the simulation is run across a water year (e.g. starts October 1),
then `"leaf on time`" and `"leaf off time`" must be adjusted relative to this
start date.  Note if `"leaf off time`" < `"leaf on time`" is ok too -- this is
the case if the simulation starts mid-summer.

* `"year duration`" ``[double]`` **1**
* `"year duration units`" ``[string]`` **noleap**
* `"leaf on time`" ``[double]`` **0**
* `"leaf on time units`" ``[string]`` **d**
* `"leaf off time`" ``[double]`` **366**
* `"leaf off time units`" ``[string]`` **d**

WARNING: there is some odd code here in which some quantities are PFT based and
some are not.  It isn't clear that those are correct for all usages.  And there
is currently no accumulation on the result of this evaluator (which is
PFT-based) to collapse it into a single sink for use in the water balance.
DON'T USE THIS WITH MULTIPLE PFTS without really understanding what you are
doing! (NOTE, error added for NPFTS > 1)  --etc


* `"number of PFTs`" ``[int]`` **1** NOTE: must be 1 currently.

* `"water limiter function`" ``[function-spec]`` **optional** If provided,
  limit the total water sink as a function of the integral of the water
  potential * rooting fraction.


Keys:

* `"plant wilting factor`"
* `"rooting depth fraction`"
* `"potential transpiration`"
* `"cell volume`"
* `"surface cell volume`"

*/

#ifndef AMANZI_FLOW_TRANSPIRATION_DISTRIBUTION_EVALUATOR_HH_
#define AMANZI_FLOW_TRANSPIRATION_DISTRIBUTION_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {

class Function;

namespace Flow {
namespace Relations {

class TranspirationDistributionEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  TranspirationDistributionEvaluator(Teuchos::ParameterList& plist);
  TranspirationDistributionEvaluator(const TranspirationDistributionEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  // need a custom EnsureCompatibility as some vectors cross meshes.
  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);
  
 protected:
  void InitializeFromPlist_();
  bool TranspirationPeriod_(double time);

  Key f_wp_key_;
  Key f_root_key_;
  Key potential_trans_key_;
  Key cv_key_;
  Key surf_cv_key_;
  int npfts_;

  double leaf_on_time_, leaf_off_time_, year_duration_;
  bool limiter_local_;
  Teuchos::RCP<Function> limiter_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator,TranspirationDistributionEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
