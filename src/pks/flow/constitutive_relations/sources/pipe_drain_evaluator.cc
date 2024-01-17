/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for water drainage through a pipe network.
  This depends on:
  1) flow depth above surface elevation (surface_depth_key_)
  2) hydraulic head in the pipe network (pressure_head_key_)

  Authors: Giacomo Capodaglio (gcapodaglio@lanl.gov)
           Naren Vohra (vohra@lanl.gov)
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
  sink_source_coeff_ = plist.get<double>("sink-source coefficient", 1.0);
  Tag tag = my_keys_.front().second;

  // my dependencies
  surface_depth_key_ = Keys::readKey(plist_, sw_domain_name_, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{surface_depth_key_, tag});

  if(!pipe_domain_name_.empty()){
     pressure_head_key_ = Keys::readKey(plist_, pipe_domain_name_, "pressure head", "pressure_head");
     dependencies_.insert(KeyTag{pressure_head_key_, tag});
  }

  if (my_keys_.front().first == "network-source_drain") {
    pipe_flag = true;
    sw_flag = false;

    mask_key_ = Keys::readKey(plist_, pipe_domain_name_, "manhole locations", "manhole_locations"); 
    dependencies_.insert(KeyTag{mask_key_, tag});

    manhole_map_key_ = Keys::readKey(plist_, pipe_domain_name_, "manhole map", "manhole_map"); 
    dependencies_.insert(KeyTag{manhole_map_key_, tag});

  }
  
  else if (my_keys_.front().first == "surface-source_drain") {
    pipe_flag = false;
    sw_flag = true;

    mask_key_ = Keys::readKey(plist_, sw_domain_name_, "manhole locations", "manhole_locations"); 
    dependencies_.insert(KeyTag{mask_key_, tag});

    manhole_map_key_ = Keys::readKey(plist_, sw_domain_name_, "manhole map", "manhole_map"); 
    dependencies_.insert(KeyTag{manhole_map_key_, tag});
  }

  cell_map_flag = 0.0;

}


Teuchos::RCP<Evaluator>
PipeDrainEvaluator::Clone() const {
  
  return Teuchos::rcp(new PipeDrainEvaluator(*this));
}


void PipeDrainEvaluator::CreateCellMap(const State& S) 
{
  // grab the meshes 
  Teuchos::RCP<const AmanziMesh::Mesh> pipe_mesh = S.GetMesh(pipe_domain_name_);
  Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh = S.GetMesh(sw_domain_name_);
  
 
  // get the manhole locations/ cell maps
  Tag tag = my_keys_.front().second;

  Key pipe_mask_key_ = Keys::readKey(plist_, pipe_domain_name_, "manhole locations", "manhole_locations"); 
  Key sw_mask_key_ = Keys::readKey(plist_, sw_domain_name_, "manhole locations", "manhole_locations"); 

  Key pipe_manhole_map_key_ = Keys::readKey(plist_, pipe_domain_name_, "manhole map", "manhole_map"); 
  Key sw_manhole_map_key_ = Keys::readKey(plist_, sw_domain_name_, "manhole map", "manhole_map"); 

  const Epetra_MultiVector& mnhMask_pipe = *S.GetPtr<CompositeVector>(pipe_mask_key_, tag)->ViewComponent("cell",false);
  const Epetra_MultiVector& mnhMask_sw = *S.GetPtr<CompositeVector>(sw_mask_key_, tag)->ViewComponent("cell",false);

  const Epetra_MultiVector& mnhmap_pipe = *S.GetPtr<CompositeVector>(pipe_manhole_map_key_, tag)->ViewComponent("cell",false);
  const Epetra_MultiVector& mnhmap_sw = *S.GetPtr<CompositeVector>(sw_manhole_map_key_, tag)->ViewComponent("cell",false);
  
  

  // loop over mesh cells and create map
  
  int ncells_pipe = pipe_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_sw = surface_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (int c_sw = 0; c_sw < ncells_sw; ++c_sw) {
    const Amanzi::AmanziGeometry::Point &xc_sw = surface_mesh->cell_centroid(c_sw);
    for (int c_pipe = 0; c_pipe < ncells_pipe; ++c_pipe) {
      const Amanzi::AmanziGeometry::Point& xc_pipe = pipe_mesh->cell_centroid(c_pipe);
      if (std::abs(mnhMask_sw[0][c_sw] - 1.0) < 1.e-12 && (std::abs(xc_sw[0] - xc_pipe[0]) < 1.e-12 ) && (std::abs(xc_sw[1] - xc_pipe[1]) < 1.e-12 ) ) {
        mnhMask_sw[0][c_sw] = c_pipe; 
      }
    }
  }

  

  for (int c_pipe = 0; c_pipe < ncells_pipe; ++c_pipe) {
    const Amanzi::AmanziGeometry::Point &xc_pipe = pipe_mesh->cell_centroid(c_pipe);
    for (int c_sw = 0; c_sw < ncells_sw; ++c_sw) {
      const Amanzi::AmanziGeometry::Point& xc_sw = surface_mesh->cell_centroid(c_sw);
      if (std::abs(mnhMask_sw[0][c_pipe] - 1.0) < 1.e-12 && (std::abs(xc_sw[0] - xc_pipe[0]) < 1.e-12 ) && (std::abs(xc_sw[1] - xc_pipe[1]) < 1.e-12 ) ) {
        mnhMask_pipe[0][c_pipe] = c_sw; 
      }
    }
  }

  cell_map_flag += 1.0;

}


void PipeDrainEvaluator::EnsureCompatibility_ToDeps_(State& S, const CompositeVectorSpace& fac)
{                                         

  if (my_keys_.front().first == "network-source_drain") {            
    for (const auto& dep : dependencies_) {
      auto domain = Keys::getDomain(dep.first);
      if (pipe_domain_name_ == domain) {
        auto& dep_fac = S.Require<CompositeVector,CompositeVectorSpace>(dep.first, dep.second);    
        dep_fac.Update(fac);
      }
    }
  }
  
  else if (my_keys_.front().first == "surface-source_drain") {            
    for (const auto& dep : dependencies_) {
      auto domain = Keys::getDomain(dep.first);
      if (sw_domain_name_ == domain) {
        auto& dep_fac = S.Require<CompositeVector,CompositeVectorSpace>(dep.first, dep.second);    
        dep_fac.Update(fac);
      }
    }
  }

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

  if (std::abs(cell_map_flag) < 1.0) {
    CreateCellMap(S);
  }

  // manhole cell map
    
  const Epetra_MultiVector& mhmap_c = *S.GetPtr<CompositeVector>(manhole_map_key_, tag)->ViewComponent("cell",false);

  if(!pipe_domain_name_.empty()){

     const Epetra_MultiVector& pressHead = *S.GetPtr<CompositeVector>(pressure_head_key_, tag)
         ->ViewComponent("cell",false);
     int c_pipe, c_sw;
     for (int c=0; c!=ncells; ++c) {

        int c1 = mhmap_c[0][c];
        //int c1 = c;
        if (pipe_flag == true) {
          c_pipe = c;
          c_sw = mhmap_c[0][c];
        } else if (sw_flag == true) {
          c_pipe = mhmap_c[0][c];
          c_sw = c;
        }

        if (pressHead[0][c_pipe] < drain_length_) {
           res[0][c] = - mnhMask[0][c] *  2.0 / 3.0 * energ_loss_coeff_weir_ * mnhPerimeter * sqrtTwoG * pow(srfcDepth[0][c_sw],3.0/2.0);
        } 
        else if (drain_length_ < pressHead[0][c_pipe] && pressHead[0][c_pipe] < (drain_length_ + srfcDepth[0][c_sw]) ){
           res[0][c] = - mnhMask[0][c] * energ_loss_coeff_subweir_ * mnhArea * sqrtTwoG 
                       * sqrt(srfcDepth[0][c_sw] + drain_length_ - pressHead[0][c_pipe]);   
        } 
        else if (pressHead[0][c_pipe] > (drain_length_ + srfcDepth[0][c_sw])) {
           res[0][c] = mnhMask[0][c] * energ_loss_coeff_orifice_ * mnhArea * sqrtTwoG * sqrt(pressHead[0][c_pipe] - drain_length_ - srfcDepth[0][c_sw]);
        }
     res[0][c] *= sink_source_coeff_;
     }
  }
  else {
     for (int c=0; c!=ncells; ++c) {
           res[0][c] = - mnhMask[0][c] *  2.0 / 3.0 * energ_loss_coeff_weir_ * mnhPerimeter * sqrtTwoG * pow(srfcDepth[0][c],3.0/2.0);
           res[0][c] *= sink_source_coeff_;
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


} //namespace
} //namespace
