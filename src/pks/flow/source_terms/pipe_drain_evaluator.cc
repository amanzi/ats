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
  H_max_ = plist.get<double>("pipe max height", 0.1);
  H_ = plist.get<double>("bed flume height", 0.478);
// need to get an input array with info on where the manhole is

  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // TODO
  // my dependencies
  //pres_key_ = Keys::readKey(plist_, domain_name, "pressure", "pressure");
  //dependencies_.insert(KeyTag{pres_key_, tag});
  //cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
  //dependencies_.insert(KeyTag{cv_key_, tag});
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

//change these two to what we need, let's call then hp and h for now
//////
  const Epetra_MultiVector& pres = *S.GetPtr<CompositeVector>(pres_key_, tag)
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& cv = *S.GetPtr<CompositeVector>(cv_key_, tag)
      ->ViewComponent("cell",false);
/////

  const auto& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double gz = -gravity[2];  // check this
  //do we also have pi somewhere? let's assume we do

  //let's assume res for us is also defined on the cells
  int ncells = res.MyLength();
// we need a characteristic function to apply the source term ONLY where the manhole is
  for (int c=0; c!=ncells; ++c) {
       if (hp < H_) {
          res[0][c] = - 4.0 / 3.0 * energy_losses_coeff_ * pi * manhole_radius * sqrt(2.0 * gz) * pow(h[c],3.0/2.0); 
       } else if (H_ < hp[c] && hp < (H_ + h[c]) ){
          res[0][c] = - energy_losses_coeff_ * pi * manhole_radius * manhole_radius * sqrt(2.0 * gz) * sqrt(h[c] + H_ - hp[c]);   
       } else if (hp > (H_ + h[c])) {
          res[0][c] = energy_losses_coeff_ * pi * manhole_radius * manhole_radius * sqrt(2.0 * gz) * sqrt(hp[c] - H_ - h[c]);
       }
  }

}

// do we need a derivative?
void OverlandPressureWaterContentEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  AMANZI_ASSERT(wrt_key == pres_key_);

  Epetra_MultiVector& res = *result[0]->ViewComponent("cell",false);
  const Epetra_MultiVector& pres = *S.GetPtr<CompositeVector>(pres_key_, tag)
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& cv = *S.GetPtr<CompositeVector>(cv_key_, tag)
      ->ViewComponent("cell",false);

  const double& p_atm = S.Get<double>("atmospheric_pressure", Tags::DEFAULT);
  const auto& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double gz = -gravity[2];  // check this

  if (wrt_key == pres_key_) {
    int ncells = res.MyLength();
    if (bar_) {
      for (int c=0; c!=ncells; ++c) {
        res[0][c] = cv[0][c] / (gz * M_);
      }
    } else if (rollover_ > 0.) {
      for (int c=0; c!=ncells; ++c) {
        double dp = pres[0][c] - p_atm;
        double ddp_eff = dp < 0. ? 0. :
          dp < rollover_ ? dp/rollover_ : 1.;
        res[0][c] = cv[0][c] * ddp_eff / (gz * M_);
      }
    } else {
      for (int c=0; c!=ncells; ++c) {
        res[0][c] = pres[0][c] < p_atm ? 0. :
          cv[0][c] / (gz * M_);
      }
    }
  } else {
    res.PutScalar(0.);
  }
}


} //namespace
} //namespace
