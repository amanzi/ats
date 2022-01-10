/*
  Interception/throughfall rates.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "interception_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

InterceptionEvaluator::InterceptionEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  std::string domain = plist_.get<std::string>("layer name");
  
  // Set up my dependencies.
  ai_key_ = plist_.get<std::string>("area index key", Keys::getKey(domain, "area_index"));
  dependencies_.insert(ai_key_);

  source_key_ = plist_.get<std::string>("source key", Keys::getKey(domain, "source"));
  dependencies_.insert(source_key_);
  
  source_in_meters_ = plist_.get<bool>("source in meters", false);
  n_liq_ = plist_.get<double>("density of liquid water [mol/m^3]", 1000. / 0.0180153);
  throughfall_coef_ = plist_.get<double>("throughfall coefficient [m^2 / m^2 biomass area]", 0.0);
};


InterceptionEvaluator::InterceptionEvaluator(const InterceptionEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    ai_key_(other.ai_key_),
    source_key_(other.source_key_),
    source_in_meters_(other.source_in_meters_),
    n_liq_(other.n_liq_),
    throughfall_coef_(other.throughfall_coef_)
{}


Teuchos::RCP<FieldEvaluator> InterceptionEvaluator::Clone() const {
  return Teuchos::rcp(new InterceptionEvaluator(*this));
}


void InterceptionEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const Teuchos::Ptr<CompositeVector>& result) {
  // Pull dependencies out of state.
  const Epetra_MultiVector& ai =
      *S->Get<CompositeVector>(ai_key_).ViewComponent("cell",false);
  const Epetra_MultiVector& source =
      *S->Get<CompositeVector>(source_key_).ViewComponent("cell",false);  
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);

  double coef = source_in_meters_ ? n_liq_ : 1.;
  
  // evaluate the model
  for (int c=0; c!=res_c.MyLength(); ++c) {
    res_c[0][c] = std::min(ai[0][c] * throughfall_coef_, 1.0) * source[0][c] * coef;
  }
}


void InterceptionEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result) {

  // Pull dependencies out of state.
  const Epetra_MultiVector& ai =
      *S->Get<CompositeVector>(ai_key_).ViewComponent("cell",false);
  const Epetra_MultiVector& source =
      *S->Get<CompositeVector>(source_key_).ViewComponent("cell",false);  
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);

  double coef = source_in_meters_ ? n_liq_ : 1.;

  if (wrt_key == ai_key_) {
    for (int c=0; c!=res_c.MyLength(); ++c) {
      res_c[0][c] = ai[0][c] * throughfall_coef_ > 1. ? source[0][c] * coef : throughfall_coef_ * source[0][c] * coef;
    }

  } else if (wrt_key == source_key_) {
    for (int c=0; c!=res_c.MyLength(); ++c) {
      res_c[0][c] = std::min(ai[0][c] * throughfall_coef_, 1.0) * coef;
    }
    
  } else {
    AMANZI_ASSERT(0);
  }  
}

} // namespace
} // namespace
} // namespace
