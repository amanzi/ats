/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The elevation evaluator gets the surface elevation, slope, and updates pres + elev.

*/

#include "snow_skin_potential_evaluator.hh"

namespace Amanzi {
namespace Flow {

SnowSkinPotentialEvaluator::SnowSkinPotentialEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  Key surf_domain = Keys::readDomainHint(plist_, domain, "snow", "surface");

  pd_key_ = Keys::readKey(plist_, surf_domain, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ pd_key_, tag });

  sd_key_ = Keys::readKey(plist_, domain, "snow depth", "depth");
  dependencies_.insert(KeyTag{ sd_key_, tag });

  precip_key_ = Keys::readKey(plist_, domain, "snow precipitation", "precipitation");
  dependencies_.insert(KeyTag{ precip_key_, tag });

  elev_key_ = Keys::readKey(plist_, surf_domain, "elevation", "elevation");
  dependencies_.insert(KeyTag{ elev_key_, tag });

  factor_ = plist_.get<double>("dt factor [s]", -1.0);
}


Teuchos::RCP<Evaluator>
SnowSkinPotentialEvaluator::Clone() const
{
  return Teuchos::rcp(new SnowSkinPotentialEvaluator(*this));
}


void
SnowSkinPotentialEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  // update pressure + elevation
  Teuchos::RCP<const CompositeVector> pd = S.GetPtr<CompositeVector>(pd_key_, tag);
  Teuchos::RCP<const CompositeVector> sd = S.GetPtr<CompositeVector>(sd_key_, tag);
  Teuchos::RCP<const CompositeVector> precip = S.GetPtr<CompositeVector>(precip_key_, tag);
  Teuchos::RCP<const CompositeVector> elev = S.GetPtr<CompositeVector>(elev_key_, tag);

  // note factor of 10 accounts for change from precip in m SWE to actual m.
  result[0]->Update(1.0, *elev, 1.0, *pd, 0.0);
  if (factor_ > 0.) {
    result[0]->Update(1.0, *sd, 10 * factor_, *precip, 1.0);
  } else {
    double dt = S.get_time() - S.last_time();
    result[0]->Update(1.0, *sd, 10 * dt, *precip, 1.0);
  }
}


// This is hopefully never called?
void
SnowSkinPotentialEvaluator::EvaluatePartialDerivative_(const State& S,
                                                       const Key& wrt_key,
                                                       const Tag& wrt_tag,
                                                       const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(0);
  result[0]->PutScalar(1.0);
}

} // namespace Flow
} // namespace Amanzi
