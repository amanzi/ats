/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for water drainage through a pipe network.
  This depends on:
  1) flow depth above surface elevation
  2) hydraulic head in the pipe network

  Authors: Giacomo Capodaglio (gcapodaglio@lanl.gov)
*/

#include "pipe_drain_evaluator.hh"

namespace Amanzi {
namespace Flow {


PipeDrainEvaluator::PipeDrainEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist)
{

  manhole_radius_ = plist_.get<double>("manhole radius", 0.24);
  energy_losses_coeff_ = plist.get<double>("energy losses coeff", 0.1);
  //H_max_ = plist.get<double>("pipe max height", 0.1);
  //H_ = plist.get<double>("bed flume height", 0.478);

  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // my dependencies
  head_key_ = Keys::readKey(plist_, domain_name, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{head_key_, tag});

  mark_key_ = Keys::readKey(plist_, domain_name, "manhole locations", "manhole_locations");
  dependencies_.insert(KeyTag{mark_key_, tag});
  
}


Teuchos::RCP<Evaluator>
PipeDrainEvaluator::Clone() const {
  return Teuchos::rcp(new PipeDrainEvaluator(*this));
}


void PipeDrainEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Epetra_MultiVector& res = *result[0]->ViewComponent("cell",false);

  const Epetra_MultiVector& head = *S.GetPtr<CompositeVector>(head_key_, tag)
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& mark = *S.GetPtr<CompositeVector>(mark_key_, tag)
      ->ViewComponent("cell",false);

  const auto& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double gz = -gravity[2];  // check this
  double pi = 3.14159265359;
  double manhole_area_ = pi * manhole_radius_ * manhole_radius_;

  int ncells = res.MyLength();
  for (int c=0; c!=ncells; ++c) {
    //if (hp < H_) {
     res[0][c] = - mark[0][c] *  4.0 / 3.0 * energy_losses_coeff_ * pi * manhole_radius_ * sqrt(2.0 * gz) * pow(head[0][c],3.0/2.0) / manhole_area_; 
          // } else if (H_ < hp[c] && hp < (H_ + head[0][c]) ){
          //res[0][c] = - energy_losses_coeff_ * manhole_area_ * sqrt(2.0 * gz) * sqrt(head[0][c] + H_ - hp[c]);   
          //} else if (hp > (H_ + head[0][c])) {
          //res[0][c] = energy_losses_coeff_ * manhole_area_ * sqrt(2.0 * gz) * sqrt(hp[c] - H_ - head[0][c]);
       //}
  }

}

void PipeDrainEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  AMANZI_ASSERT(wrt_key == head_key_);

  Epetra_MultiVector& res = *result[0]->ViewComponent("cell",false);
  
  const Epetra_MultiVector& head = *S.GetPtr<CompositeVector>(head_key_, tag)
       ->ViewComponent("cell",false);

    const Epetra_MultiVector& mark = *S.GetPtr<CompositeVector>(mark_key_, tag)
      ->ViewComponent("cell",false);

  const auto& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double gz = -gravity[2];  // check this
  double pi = 3.14159265359;
  double manhole_area_ = pi * manhole_radius_ * manhole_radius_;

  if (wrt_key == head_key_) {
     int ncells = res.MyLength();
     for (int c=0; c!=ncells; ++c) {
       res[0][c] = - mark[0][c] *  2.0 * energy_losses_coeff_ * pi * manhole_radius_ * sqrt(2.0 * gz) * sqrt(head[0][c]) / manhole_area_;
     }
     
  }
}


} //namespace
} //namespace
