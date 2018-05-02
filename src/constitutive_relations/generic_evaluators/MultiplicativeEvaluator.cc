/*
  MultiplicativeEvaluator is the generic evaluator for multipying two vectors.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "MultiplicativeEvaluator.hh"

namespace Amanzi {
namespace Relations {

MultiplicativeEvaluator::MultiplicativeEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondary(plist) {
  coef_ = plist_.get<double>("coefficient", 1.0);
}

MultiplicativeEvaluator::MultiplicativeEvaluator(const MultiplicativeEvaluator& other) :
    EvaluatorSecondary(other),
    coef_(other.coef_)
{}

Teuchos::RCP<Evaluator>
MultiplicativeEvaluator::Clone() const
{
  return Teuchos::rcp(new MultiplicativeEvaluator(*this));
}


// Required methods from EvaluatorSecondary
void
MultiplicativeEvaluator::Evaluate_(const State& S,
                                   CompositeVector& result)
{
  ASSERT(dependencies_.size() > 1);
  auto  key = dependencies_.begin();
  result = S.Get<CompositeVector>((*key).first, (*key).second);
  result.Scale(coef_);
  key++;

  for (; key!=dependencies_.end(); ++key) {
    const CompositeVector& dep = S.Get<CompositeVector>((*key).first, (*key).second);
    for (CompositeVector::name_iterator lcv=result.begin(); lcv!=result.end(); ++lcv) {
      Epetra_MultiVector& res_c = *result.ViewComponent(*lcv, false);
      const Epetra_MultiVector& dep_c = *dep.ViewComponent(*lcv, false);

      for (int c=0; c!=res_c.MyLength(); ++c) res_c[0][c] *= dep_c[0][c];
    }
  }
}

void
MultiplicativeEvaluator::EvaluatePartialDerivative_(const State& S,
                                                    const Key& wrt_key,
                                                    const Key& wrt_tag,
                                                    CompositeVector& result)
{
  ASSERT(dependencies_.size() > 1);

  auto  key = dependencies_.begin();
  while ((*key).first == wrt_key) key++;
  result = S.Get<CompositeVector>((*key).first, (*key).second);;
  result.Scale(coef_);
  key++;

  for (; key!=dependencies_.end(); ++key) {
    if ((*key).first != wrt_key) {
      const CompositeVector& dep = S.Get<CompositeVector>((*key).first, (*key).second);
      for (CompositeVector::name_iterator lcv=result.begin(); lcv!=result.end(); ++lcv) {
        Epetra_MultiVector& res_c = *result.ViewComponent(*lcv, false);
        const Epetra_MultiVector& dep_c = *dep.ViewComponent(*lcv, false);

        for (int c=0; c!=res_c.MyLength(); ++c) res_c[0][c] *= dep_c[0][c];
      }
    }
  }
}


} // namespace
} // namespace

