/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for water drainage through a pipe network.
  This depends on:
  1) flow depth above surface elevation (surface_depth_key_)
  2) hydraulic head in the pipe network (pressure_head_key_)

  Authors: Giacomo Capodaglio (gcapodaglio@lanl.gov)
*/

#include "pipe_drain_evaluator.hh"

namespace Amanzi {
namespace Flow {


PipeDrainEvaluator::PipeDrainEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist)
{

  manhole_radius_ = plist_.get<double>("manhole radius", 0.24);
  energ_loss_coeff_weir_ = plist.get<double>("energy losses coeff weir", 0.54);
  energ_loss_coeff_subweir_ = plist.get<double>("energy losses coeff submerged weir", 0.056);
  energ_loss_coeff_orifice_ = plist.get<double>("energy losses coeff orifice", 0.167);
  drain_length_ = plist.get<double>("drain length", 0.478);
  sw_domain_name_ = plist.get<std::string>("sw domain name", "surface"); 
  pipe_domain_name_ = plist.get<std::string>("pipe domain name", ""); 
  Tag tag = my_keys_.front().second;

  // my dependencies
  surface_depth_key_ = Keys::readKey(plist_, sw_domain_name_, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{surface_depth_key_, tag});

  if(!pipe_domain_name_.empty()){
     pressure_head_key_ = Keys::readKey(plist_, pipe_domain_name_, "pressure head", "pressure_head");
     dependencies_.insert(KeyTag{pressure_head_key_, tag});
  }

  mask_key_ = Keys::readKey(plist_, sw_domain_name_, "manhole locations", "manhole_locations"); 
  dependencies_.insert(KeyTag{mask_key_, tag});


  manhole_map_key_ = Keys::readKey(plist_, sw_domain_name_, "manhole map", "manhole_map"); 
  dependencies_.insert(KeyTag{manhole_map_key_, tag});
  
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

  const Epetra_MultiVector& mnhMask = *S.GetPtr<CompositeVector>(mask_key_, tag)
      ->ViewComponent("cell",false);

  double g = norm(S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT));
  double pi = 3.14159265359;
  double mnhArea = pi * manhole_radius_ * manhole_radius_;
  double mnhPerimeter = 2.0 * pi * manhole_radius_;
  double sqrtTwoG = sqrt(2.0 * g);

  int ncells = res.MyLength();


  // manhole cell map
  const Epetra_MultiVector& mhmap_c = *S.GetPtr<CompositeVector>(manhole_map_key_, tag)->ViewComponent("cell",false);
  /*
  if(manhole_map_key_.empty()){
    for (int c = 0; c < ncells; c++) {
      mhmap_c[0][c] = 1.0;
    }
  }
  */

  if(!pipe_domain_name_.empty()){

     const Epetra_MultiVector& pressHead = *S.GetPtr<CompositeVector>(pressure_head_key_, tag)
         ->ViewComponent("cell",false);

     for (int c=0; c!=ncells; ++c) {

        int c1 = mhmap_c[0][c];

        if (pressHead[0][c] < drain_length_) {
           res[0][c] = - mnhMask[0][c] *  2.0 / 3.0 * energ_loss_coeff_weir_ * mnhPerimeter * sqrtTwoG * pow(srfcDepth[0][c1],3.0/2.0);
        } 
        else if (drain_length_ < pressHead[0][c] && pressHead[0][c] < (drain_length_ + srfcDepth[0][c1]) ){
           res[0][c] = - mnhMask[0][c] * energ_loss_coeff_subweir_ * mnhArea * sqrtTwoG 
                       * sqrt(srfcDepth[0][c1] + drain_length_ - pressHead[0][c]);   
        } 
        else if (pressHead[0][c] > (drain_length_ + srfcDepth[0][c1])) {
           res[0][c] = mnhMask[0][c] * energ_loss_coeff_orifice_ * mnhArea * sqrtTwoG * sqrt(pressHead[0][c] - drain_length_ - srfcDepth[0][c1]);
        }
     }
  }
  else {
     for (int c=0; c!=ncells; ++c) {
           res[0][c] = - mnhMask[0][c] *  2.0 / 3.0 * energ_loss_coeff_weir_ * mnhPerimeter * sqrtTwoG * pow(srfcDepth[0][c],3.0/2.0);
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

  const Epetra_MultiVector& mnhMask = *S.GetPtr<CompositeVector>(mask_key_, tag)
      ->ViewComponent("cell",false);

  double g = norm(S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT));
  double pi = 3.14159265359;
  double mnhArea = pi * manhole_radius_ * manhole_radius_;
  double mnhPerimeter = 2.0 * pi * manhole_radius_;
  double sqrtTwoG = sqrt(2.0 * g);

  int ncells = res.MyLength();


  if(!pipe_domain_name_.empty()){

     const Epetra_MultiVector& pressHead = *S.GetPtr<CompositeVector>(pressure_head_key_, tag)
        ->ViewComponent("cell",false);

     if (wrt_key == surface_depth_key_) {
        for (int c=0; c!=ncells; ++c) {
           if (pressHead[0][c] < drain_length_) {
              res[0][c] = - mnhMask[0][c] * energ_loss_coeff_weir_ * mnhPerimeter* sqrtTwoG * sqrt(srfcDepth[0][c]);
           }
           else if (drain_length_ < pressHead[0][c] && pressHead[0][c] < (drain_length_ + srfcDepth[0][c]) ){
              res[0][c] = - 0.5 * mnhMask[0][c] * energ_loss_coeff_subweir_ * mnhArea * sqrtTwoG 
                          / sqrt(srfcDepth[0][c] + drain_length_ - pressHead[0][c]);
           }
           else if (pressHead[0][c] > (drain_length_ + srfcDepth[0][c])) {
              res[0][c] = - 0.5 * mnhMask[0][c] * energ_loss_coeff_orifice_ * mnhArea * sqrtTwoG 
                          / sqrt(pressHead[0][c] - drain_length_ - srfcDepth[0][c]);
           }
        }   
     }
     else if (wrt_key == pressure_head_key_) {
        for (int c=0; c!=ncells; ++c) {
           if (pressHead[0][c] < drain_length_) {
              res[0][c] = 0.0; 
           }
           else if (drain_length_ < pressHead[0][c] && pressHead[0][c] < (drain_length_ + srfcDepth[0][c]) ){
              res[0][c] = 0.5 * mnhMask[0][c] * energ_loss_coeff_subweir_ * mnhArea * sqrtTwoG 
                          / sqrt(srfcDepth[0][c] + drain_length_ - pressHead[0][c]);
           }
           else if (pressHead[0][c] > (drain_length_ + srfcDepth[0][c])) {
              res[0][c] = 0.5 * mnhMask[0][c] * energ_loss_coeff_orifice_ * mnhArea * sqrtTwoG 
                          / sqrt(pressHead[0][c] - drain_length_ - srfcDepth[0][c]);
           }
        }
     }
     else {
       AMANZI_ASSERT(0);
     }
  }
  else {
     if (wrt_key == surface_depth_key_) {
        for (int c=0; c!=ncells; ++c) {
              res[0][c] = - mnhMask[0][c] * energ_loss_coeff_weir_ * mnhPerimeter* sqrtTwoG * sqrt(srfcDepth[0][c]);
        }
     }
     else {
       AMANZI_ASSERT(0);
     }
  }
}

void PipeDrainEvaluator::isManhole(AmanziGeometry::Point xc)
{
  double x = xc[0], y = xc[1];
}

} //namespace
} //namespace
