/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Determining the molar fraction of a gas component within a gas mixture.

*/

#include "organic_matter_evaluator.hh"
#include "boost/math/constants/constants.hpp"

namespace Amanzi {

OrganicMatterRateEvaluator::OrganicMatterRateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;
  Key domain_name = Keys::getDomain(my_keys_.front().first);

  biomass_key_ = Keys::readKey(plist_, domain_name, "biomass", "biomass");
  dependencies_.insert(KeyTag{ biomass_key_, tag });

  // please put units on all of these!  --etc
  Bmax_ = 1.0 / plist_.get<double>("maximum biomass");
  Q_db0_ = plist_.get<double>("empirical coefficient");
  Q_on_Bmax_ = Q_db0_ / Bmax_;
}


Teuchos::RCP<Evaluator>
OrganicMatterRateEvaluator::Clone() const
{
  return Teuchos::rcp(new OrganicMatterRateEvaluator(*this));
}


void
OrganicMatterRateEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const Epetra_MultiVector& bio = *S.Get<CompositeVector>(biomass_key_, tag).ViewComponent("cell");
  Epetra_MultiVector& result_c = *result[0]->ViewComponent("cell");

  result_c.PutScalar(0.);
  for (int c = 0; c < result_c.MyLength(); c++) {
    for (int j = 0; j < bio.NumVectors(); j++) { result_c[0][c] += Q_on_Bmax_ * bio[j][c]; }
  }
}


void
OrganicMatterRateEvaluator::EvaluatePartialDerivative_(const State& S,
                                                       const Key& wrt_key,
                                                       const Tag& wrt_tag,
                                                       const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(0);
}

} // namespace Amanzi
