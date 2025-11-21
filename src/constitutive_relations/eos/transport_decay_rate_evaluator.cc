/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Bo Gao (gaob@ornl.gov)
*/

/*
  The elevation evaluator calculates temperature-dependent decay rate.

*/

#include "transport_decay_rate_evaluator.hh"

namespace Amanzi {
namespace Relations {

TransportDecayRateEvaluator::TransportDecayRateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;
  domain_ = Keys::getDomain(my_keys_.front().first); // column, domain

  temp_key_ = Keys::readKey(plist, domain_, "temperature", "temperature");
  dependencies_.insert(KeyTag{ temp_key_, tag });

  q10_ = plist_.get<double>("Q10 [-]", 2.0);
  temp_ref_ = plist_.get<double>("reference temperature [K]", 278.15);
  decay_ref_ = plist_.get<double>("reference decay rate [s^-1]", 1.925e-4);
}


Teuchos::RCP<Evaluator>
TransportDecayRateEvaluator::Clone() const
{
  return Teuchos::rcp(new TransportDecayRateEvaluator(*this));
}


void
TransportDecayRateEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell", false);

  const auto& temp_c = *S.Get<CompositeVector>(temp_key_, tag).ViewComponent("cell", false);
  const auto& mesh = *S.GetMesh(domain_);

  for (int col = 0; col != mesh.columns.num_columns_owned; ++col) {
    const auto& col_cells = mesh.columns.getCells(col);
    for (int c = 0; c != col_cells.size(); ++c) {
      double f_temp = Func_Temp(temp_c[0][col_cells[c]], temp_ref_, q10_);
      res_c[0][col_cells[c]] = f_temp * decay_ref_; // 1/(s)
    }
  }
}


double
TransportDecayRateEvaluator::Func_Temp(double temp, double temp_ref, double q10) const
{
  return std::pow(q10, (temp - temp_ref) / 10.0);
}

} // namespace Relations
} // namespace Amanzi
