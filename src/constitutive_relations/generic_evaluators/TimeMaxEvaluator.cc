/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "TimeMaxEvaluator.hh"

namespace Amanzi {
namespace Relations {

TimeMaxEvaluator::TimeMaxEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist), evaluated_once_(false)
{
  operator_ = plist.get<std::string>("operator", "max");
  if (operator_ != "max" && operator_ != "min") {
    Errors::Message msg;
    msg << "TimeMaxEvaluator: " << my_keys_.front().first << ": invalid operator \"" << operator_
        << "\"";
    Exceptions::amanzi_throw(msg);
  }

  if (dependencies_.size() == 0) {
    Tag tag = my_keys_.front().second;
    Key key = my_keys_.front().first;
    AMANZI_ASSERT(Keys::starts_with(key, operator_ + "_"));
    dependencies_.insert(KeyTag{ key.substr(4, key.size()), tag });
  }
}


Teuchos::RCP<Evaluator>
TimeMaxEvaluator::Clone() const
{
  return Teuchos::rcp(new TimeMaxEvaluator(*this));
}


void
TimeMaxEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  if (!evaluated_once_) {
    if (operator_ == "max") {
      result[0]->PutScalar(-1.e16);
    } else {
      result[0]->PutScalar(1.e16);
    }
    evaluated_once_ = true;
  }

  Epetra_MultiVector& res = *result[0]->ViewComponent("cell", false);
  const Epetra_MultiVector& dep =
    *S.Get<CompositeVector>(dependencies_.front().first, dependencies_.front().second)
       .ViewComponent("cell", false);

  if (operator_ == "max") {
    for (int c = 0; c != res.MyLength(); ++c) {
      res[0][c] = std::max(res[0][c], dep[0][c]);
    }
  }
}


} // namespace Relations
} // namespace Amanzi
