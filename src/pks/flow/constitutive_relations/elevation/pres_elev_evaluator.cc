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

#include "pres_elev_evaluator.hh"

namespace Amanzi {
namespace Flow {

PresElevEvaluator::PresElevEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  pres_key_ = Keys::readKey(plist_, domain, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ pres_key_, tag });
  elev_key_ = Keys::readKey(plist_, domain, "elevation", "elevation");
  dependencies_.insert(KeyTag{ elev_key_, tag });
}


Teuchos::RCP<Evaluator>
PresElevEvaluator::Clone() const
{
  return Teuchos::rcp(new PresElevEvaluator(*this));
}


void
PresElevEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  // update pressure + elevation
  Teuchos::RCP<const CompositeVector> pres = S.GetPtr<CompositeVector>(pres_key_, tag);
  Teuchos::RCP<const CompositeVector> elev = S.GetPtr<CompositeVector>(elev_key_, tag);
  result[0]->Update(1.0, *elev, 1.0, *pres, 0.0);
}


// This is hopefully never called?
void
PresElevEvaluator::EvaluatePartialDerivative_(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag,
                                              const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(1.0);
}

} // namespace Flow
} // namespace Amanzi
