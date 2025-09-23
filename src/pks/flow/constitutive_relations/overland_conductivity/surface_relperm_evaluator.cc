/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates the conductivity of surface flow.

*/

#include "surface_relperm_evaluator.hh"
#include "surface_relperm_model_factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

SurfaceRelPermEvaluator::SurfaceRelPermEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;
  Key domain = Keys::getDomain(my_keys_.front().first);

  // create the model
  SurfaceRelPermModelFactory fac;
  model_ = fac.createModel(plist_.sublist("surface rel perm model"));

  // set up the height dependency
  h_key_ = Keys::readKey(plist_, domain, "pressure key", "pressure");
  dependencies_.insert(KeyTag{ h_key_, tag });

  // set up the temperature dependency
  is_temp_ = model_->TemperatureDependent();
  if (is_temp_) {
    uf_key_ = Keys::readKey(plist_, domain, "unfrozen fraction", "unfrozen_fraction");
    dependencies_.insert(KeyTag{ uf_key_, tag });
  }
}


Teuchos::RCP<Evaluator>
SurfaceRelPermEvaluator::Clone() const
{
  return Teuchos::rcp(new SurfaceRelPermEvaluator(*this));
}


// Required methods from EvaluatorSecondaryMonotypeCV
void
SurfaceRelPermEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  if (is_temp_) {
    Teuchos::RCP<const CompositeVector> uf = S.GetPtr<CompositeVector>(uf_key_, tag);
    Teuchos::RCP<const CompositeVector> h = S.GetPtr<CompositeVector>(h_key_, tag);

    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& uf_v = *uf->ViewComponent(*comp, false);
      const Epetra_MultiVector& h_v = *h->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->SurfaceRelPerm(uf_v[0][i], h_v[0][i]);
      }
    }

  } else {
    Teuchos::RCP<const CompositeVector> h = S.GetPtr<CompositeVector>(h_key_, tag);

    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& h_v = *h->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->SurfaceRelPerm(0., h_v[0][i]);
      }
    }
  }
}


void
SurfaceRelPermEvaluator::EvaluatePartialDerivative_(const State& S,
                                                    const Key& wrt_key,
                                                    const Tag& wrt_tag,
                                                    const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(0);
}


} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
