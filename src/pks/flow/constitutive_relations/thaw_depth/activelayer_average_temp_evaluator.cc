/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The column average temperature evaluator gets the subsurface temperature and number of cells (related to depth), and returns the average column temperature.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "activelayer_average_temp_evaluator.hh"

namespace Amanzi {
namespace Flow {



ActiveLayerAverageTempEvaluator::ActiveLayerAverageTempEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotypeCV(plist)
{
  domain_ = Keys::getDomain(my_keys_.front().first); //surface_column domain
  Tag tag = my_keys_.front().second;

  Key domain_ss = Keys::readDomainHint(plist_, domain_, "surface", "subsurface");
  temp_key_ = Keys::readKey(plist, domain_ss,"temperature", "temperature");
  dependencies_.insert(KeyTag{temp_key_, tag});

  trans_width_ =  plist_.get<double>("transition width [K]", 0.2);
}


Teuchos::RCP<Evaluator>
ActiveLayerAverageTempEvaluator::Clone() const
{
  return Teuchos::rcp(new ActiveLayerAverageTempEvaluator(*this));
}


void
ActiveLayerAverageTempEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  double trans_temp = 273.15 + 0.5*trans_width_;
  Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell",false);
  const auto& temp_c = *S.Get<CompositeVector>(temp_key_, tag).ViewComponent("cell", false);

  AMANZI_ASSERT(res_c.MyLength() == 1);
  int col_cells = temp_c.MyLength();
  double temp_sum = 0;
  int count = 0 ;

  for (int i=0; i!=col_cells; ++i) {
    if (temp_c[0][i] >= trans_temp) {
      temp_sum += temp_c[0][i];
      count += 1;
    }
  }
  res_c[0][0] = count > 0 ? temp_sum/count : 0.0;
}

void
ActiveLayerAverageTempEvaluator::EvaluatePartialDerivative_(const State& S,
               const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{}

// Custom EnsureCompatibility forces this to be updated once.
bool
ActiveLayerAverageTempEvaluator::Update(State& S, const Key& request)
{
  bool changed = EvaluatorSecondaryMonotypeCV::Update(S,request);

  if (!updated_once_) {
    Update_(S);
    updated_once_ = true;
    return true;
  }
  return changed;
}

void
ActiveLayerAverageTempEvaluator::EnsureCompatibility(State& S)
{
  // note, no derivs are valid here
  EnsureCompatibility_ClaimOwnership_(S);
  EnsureCompatibility_Flags_(S);
  EnsureCompatibility_DepEvals_(S);

  CompositeVectorSpace fac;
  fac.AddComponent("cell", AmanziMesh::CELL, 1);
  EnsureCompatibility_DepsFromFac_(S, fac, true);

  EnsureCompatibility_DepEnsureCompatibility_(S);
}


} //namespace
} //namespace
