/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for water drainage through a pipe network.
  This depends on:
  1) flow depth above surface elevation (surface_depth_key_)
  2) hydraulic head in the pipe network (pipe_depth_key_)

  Authors: Giacomo Capodaglio (gcapodaglio@lanl.gov)
*/

#include "pipe_drain_evaluator.hh"

namespace Amanzi {
namespace Flow {


PipeDrainEvaluator::PipeDrainEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist)
{

  manhole_radius_ = plist_.get<double>("manhole radius", 0.24);
  energ_loss_coeff_ = plist.get<double>("energy losses coeff", 0.1);
  drain_length_ = plist.get<double>("drain length", 0.478);

  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // my dependencies
  surface_depth_key_ = Keys::readKey(plist_, domain_name, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{surface_depth_key_, tag});

  pipe_depth_key_ = Keys::readKey(plist_, domain_name, "water depth", "water_depth");
  dependencies_.insert(KeyTag{pipe_depth_key_, tag});

  mask_key_ = Keys::readKey(plist_, domain_name, "manhole locations", "manhole_locations");
  dependencies_.insert(KeyTag{mask_key_, tag});
  
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

  const Epetra_MultiVector& srfcDepth = *S.GetPtr<CompositeVector>(surface_depth_key_, tag)
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& pipeDepth = *S.GetPtr<CompositeVector>(pipe_depth_key_, tag)
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& mnhMask = *S.GetPtr<CompositeVector>(mask_key_, tag)
      ->ViewComponent("cell",false);

  const auto& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double g = -gravity[1]; 
  double pi = 3.14159265359;
  double mnhArea = pi * manhole_radius_ * manhole_radius_;
  double mnhPerimeter = 2.0 * pi * manhole_radius_;   
  double sqrtTwoG = sqrt(2.0 * g); 

  int ncells = res.MyLength();
  for (int c=0; c!=ncells; ++c) {
    if (pipeDepth[0][c] < drain_length_) {
       res[0][c] = - mnhMask[0][c] *  2.0 / 3.0 * energ_loss_coeff_ * mnhPerimeter * sqrtTwoG * pow(srfcDepth[0][c],3.0/2.0);
    } 
    else if (drain_length_ < pipeDepth[0][c] && pipeDepth[0][c] < (drain_length_ + srfcDepth[0][c]) ){
       res[0][c] = - mnhMask[0][c] * energ_loss_coeff_ * mnhArea * sqrtTwoG * sqrt(srfcDepth[0][c] + drain_length_ - pipeDepth[0][c]);   
    } 
    else if (pipeDepth[0][c] > (drain_length_ + srfcDepth[0][c])) {
       res[0][c] = mnhMask[0][c] * energ_loss_coeff_ * mnhArea * sqrtTwoG * sqrt(pipeDepth[0][c] - drain_length_ - srfcDepth[0][c]);
    }
  }

}

void PipeDrainEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  Epetra_MultiVector& res = *result[0]->ViewComponent("cell",false);

  const Epetra_MultiVector& srfcDepth = *S.GetPtr<CompositeVector>(surface_depth_key_, tag)
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& pipeDepth = *S.GetPtr<CompositeVector>(pipe_depth_key_, tag)
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& mnhMask = *S.GetPtr<CompositeVector>(mask_key_, tag)
      ->ViewComponent("cell",false);

  const auto& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double g = -gravity[1];
  double pi = 3.14159265359;
  double mnhArea = pi * manhole_radius_ * manhole_radius_;
  double mnhPerimeter = 2.0 * pi * manhole_radius_;
  double sqrtTwoG = sqrt(2.0 * g);

  int ncells = res.MyLength();
  if (wrt_key == surface_depth_key_) {
     for (int c=0; c!=ncells; ++c) {
        if (pipeDepth[0][c] < drain_length_) {
           res[0][c] = - mnhMask[0][c] * energ_loss_coeff_ * mnhPerimeter* sqrtTwoG * sqrt(srfcDepth[0][c]);
        }
        else if (drain_length_ < pipeDepth[0][c] && pipeDepth[0][c] < (drain_length_ + srfcDepth[0][c]) ){
           res[0][c] = - 0.5 * mnhMask[0][c] * energ_loss_coeff_ * mnhArea * sqrtTwoG / sqrt(srfcDepth[0][c] + drain_length_ - pipeDepth[0][c]);
        }
        else if (pipeDepth[0][c] > (drain_length_ + srfcDepth[0][c])) {
           res[0][c] = - 0.5 * mnhMask[0][c] * energ_loss_coeff_ * mnhArea * sqrtTwoG / sqrt(pipeDepth[0][c] - drain_length_ - srfcDepth[0][c]);
        }
     }   
  }
  else if (wrt_key == pipe_depth_key_) {
     for (int c=0; c!=ncells; ++c) {
        if (pipeDepth[0][c] < drain_length_) {
           res[0][c] = 0.0; 
        }
        else if (drain_length_ < pipeDepth[0][c] && pipeDepth[0][c] < (drain_length_ + srfcDepth[0][c]) ){
           res[0][c] = 0.5 * mnhMask[0][c] * energ_loss_coeff_ * mnhArea * sqrtTwoG / sqrt(srfcDepth[0][c] + drain_length_ - pipeDepth[0][c]);
        }
        else if (pipeDepth[0][c] > (drain_length_ + srfcDepth[0][c])) {
           res[0][c] = 0.5 * mnhMask[0][c] * energ_loss_coeff_ * mnhArea * sqrtTwoG / sqrt(pipeDepth[0][c] - drain_length_ - srfcDepth[0][c]);
        }
     }
  }
  else {
    AMANZI_ASSERT(0);
  }
}


} //namespace
} //namespace
