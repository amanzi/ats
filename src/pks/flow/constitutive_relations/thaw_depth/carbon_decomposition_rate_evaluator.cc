/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the subsurface temperature and computes moisture content 
  over time.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "carbon_decomposition_rate_evaluator.hh"

namespace Amanzi {
namespace Flow {


CarbonDecomposeRateEvaluator::CarbonDecomposeRateEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotypeCV(plist)
{
  Key dset_name = plist.get<std::string>("domain set name", "column");
  
  domain_ = Keys::getDomain(my_key_); //surface_column domain
  int col_id = Keys::getDomainSetIndex<int>(domain_);
  
  Key domain_ss = Keys::getDomainInSet(dset_name, col_id);
  
  temp_key_ =  Keys::readKey(plist, domain_ss, "temperature","temperature");
  dependencies_.insert(temp_key_);
  
  pres_key_ =  Keys::readKey(plist, domain_ss, "pressure", "pressure");
  dependencies_.insert(pres_key_);
  
  por_key_ =  Keys::readKey(plist, domain_ss, "porosity", "porosity");
  dependencies_.insert(por_key_);

  cv_key_ =  Keys::readKey(plist, domain_ss, "cell volume","cell_volume");
  dependencies_.insert(cv_key_);
  
  //trans_width_ =  plist_.get<double>("transition width [K]", 0.2);
  q10_ =  plist_.get<double>("Q10 [-]", 2.0);

}
  

CarbonDecomposeRateEvaluator::CarbonDecomposeRateEvaluator(const CarbonDecomposeRateEvaluator& other)
  : EvaluatorSecondaryMonotypeCV(other),
    temp_key_(other.temp_key_),
    pres_key_(other.pres_key_),
    por_key_(other.por_key_),
    cv_key_(other.cv_key_),
    q10_(other.q10_)
{}
  
Teuchos::RCP<Evaluator>
CarbonDecomposeRateEvaluator::Clone() const
{
  return Teuchos::rcp(new CarbonDecomposeRateEvaluator(*this));
}


void
CarbonDecomposeRateEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{ 
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);

  const auto& temp_c = *S->Get<CompositeVector>(temp_key_).ViewComponent("cell", false);
  const auto& pres_c = *S->Get<CompositeVector>(pres_key_).ViewComponent("cell", false);
  const auto& por_c = *S->Get<CompositeVector>(por_key_).ViewComponent("cell", false);
  const auto& vol_c = *S->Get<CompositeVector>(cv_key_).ViewComponent("cell", false);

  std::string domain_ss = Keys::getDomain(temp_key_);
  const auto& top_z_centroid = S->GetMesh(domain_ss)->face_centroid(0);
  AmanziGeometry::Point z_up_centroid(top_z_centroid);
  AmanziGeometry::Point z_down_centroid(top_z_centroid);
  AmanziGeometry::Point z_depth(top_z_centroid);
    
  int col_cells = temp_c.MyLength();
  double col_sum = 0;

  for (int i=0; i!=col_cells; ++i) {
    if (temp_c[0][i] >= 273.15) {
      
      z_up_centroid = S->GetMesh(domain_ss)->face_centroid(i);
      z_down_centroid = S->GetMesh(domain_ss)->face_centroid(i+1);
      
      double dz = z_up_centroid[2] - z_down_centroid[2];
      double depth = top_z_centroid[2] - z_down_centroid[2];

      double f_temp = std::pow(q10_,(temp_c[0][i]-293.15)/10.0);
      
      double f_depth = depth < 1.0 ? 1 : std::exp(-0.5*(depth-1));

      double f_pres_temp = Func_TempPres(temp_c[0][i],pres_c[0][i]);

      double soil = (1-por_c[0][i]) * vol_c[0][i];
      col_sum += (f_temp * f_depth * f_pres_temp * dz) * soil;
    }

  }
  
  res_c[0][0] = col_sum;
}
  
void
CarbonDecomposeRateEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{}

 
// Custom EnsureCompatibility forces this to be updated once.
bool
CarbonDecomposeRateEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
        Key request)
{
  bool changed = EvaluatorSecondaryMonotypeCV::HasFieldChanged(S,request);

  if (!updated_once_) {
    UpdateField_(S);
    updated_once_ = true;
    return true;
  }
  return changed;
}

void
CarbonDecomposeRateEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{

  AMANZI_ASSERT(my_key_ != std::string(""));
   
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->Require<CompositeVector,CompositeVectorSpace>(my_key_, Tags::NEXT,  my_key_);
  
  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);
  
  if (my_fac->Mesh() != Teuchos::null) {
    // Recurse into the tree to propagate info to leaves.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      S->RequireEvaluator(*key)->EnsureCompatibility(S);
    }
  }
}

double
CarbonDecomposeRateEvaluator::Func_TempPres(double temp, double pres) {
  double p_min = -1.0e7;
  double p_max = -1.0e4;
  double p_atm = 101325.;
  
  double pn_star = pres - p_atm;
  double pn = std::min(0.0, std::max(pn_star,p_min));
  
  double f = -1.e18;
  
  if (pn >= p_max) {
    f = pn / p_max;
  }
  else {
    AMANZI_ASSERT(pn != 0.);
    f = std::log(p_min/pn) / std::log(p_min/p_max);
  }
  
  return f; 
  
}

} //namespace
} //namespace
