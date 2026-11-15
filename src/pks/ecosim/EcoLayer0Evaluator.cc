#include "EcoLayer0Evaluator.hh"
#include "State.hh"
#include "Key.hh"
#include "errors.hh"
#include "exceptions.hh"
#include <cmath>
extern "C" void ecosim_layer0_retention(const double*,double,double,double*,double*,int*);
namespace Amanzi { namespace EcoSIM {
EcoLayer0Evaluator::EcoLayer0Evaluator(Teuchos::ParameterList& plist)
 : EvaluatorSecondaryMonotypeCV(plist),energy_(plist.get<std::string>("quantity")=="energy") {
 const auto domain=Keys::getDomain(my_keys_.front().first);
 const auto tag=my_keys_.front().second;
 base_=plist.get<std::string>("mobile storage key");
 material_=Keys::getKey(domain,"ecosim_layer0_material");
 pressure_=Keys::getKey(domain,"pressure");area_=Keys::getKey(domain,"cell_volume");
 temperature_=Keys::getKey(domain,"temperature");fraction_=Keys::getKey(domain,"unfrozen_fraction");
 liquid_=Keys::getKey(domain,"internal_energy_liquid");ice_=Keys::getKey(domain,"internal_energy_ice");
 for(const auto& k:{base_,material_,pressure_,area_})dependencies_.insert({k,tag});
 if(energy_)for(const auto& k:{temperature_,fraction_,liquid_,ice_})dependencies_.insert({k,tag});
}
void EcoLayer0Evaluator::EnsureCompatibility_ToDeps_(State& S) {
 auto tag=my_keys_.front().second;
 const auto& fac=S.Require<CompositeVector,CompositeVectorSpace>(my_keys_.front().first,tag);
 if(fac.Mesh()!=Teuchos::null)for(const auto& dep:dependencies_) {
  if(dep.first==material_)S.Require<CompositeVector,CompositeVectorSpace>(dep.first,dep.second)
    .SetMesh(fac.Mesh())->SetGhosted()->AddComponent("cell",AmanziMesh::CELL,14);
  else S.Require<CompositeVector,CompositeVectorSpace>(dep.first,dep.second).Update(fac);
 }
 S.Require<double>("atmospheric_pressure",Tags::DEFAULT);
}
void EcoLayer0Evaluator::Evaluate_(const State& S,const std::vector<CompositeVector*>& r) {
 Calculate_(S,"",r);
}
void EcoLayer0Evaluator::EvaluatePartialDerivative_(const State& S,const Key& wrt,const Tag&,
 const std::vector<CompositeVector*>& r) { Calculate_(S,wrt,r); }
void EcoLayer0Evaluator::Calculate_(const State& S,const Key& wrt,
 const std::vector<CompositeVector*>& result) const {
 const auto tag=my_keys_.front().second;
 auto field=[&](const Key& k)->const Epetra_MultiVector& {
  return *S.Get<CompositeVector>(k,tag).ViewComponent("cell",false);
 };
 const auto& base=field(base_);const auto& p=field(pressure_);const auto& area=field(area_);
 const auto& mat=field(material_);auto& out=*result[0]->ViewComponent("cell",false);
 const double patm=S.Get<double>("atmospheric_pressure",Tags::DEFAULT);
 constexpr double water_moles_per_m3=1000./0.0180153;
 for(int c=0;c<out.MyLength();++c){
  double props[14];for(int j=0;j<14;++j)props[j]=mat[j][c];
  double h,dh;int status;ecosim_layer0_retention(props,p[0][c],patm,&h,&dh,&status);
  if(status || !std::isfinite(h) || !std::isfinite(dh)) {
   Errors::Message msg("EcoSIM layer zero: invalid native litter retention parameters");
   Exceptions::amanzi_throw(msg);
  }
  double retained=water_moles_per_m3*h;
  if(wrt.empty()){
   if(energy_){
    double f=field(fraction_)[0][c],u=f*field(liquid_)[0][c]+(1-f)*field(ice_)[0][c];
    out[0][c]=base[0][c]+area[0][c]*(retained*u+props[1]*(field(temperature_)[0][c]-273.15));
   }else out[0][c]=base[0][c]+area[0][c]*retained;
  }else if(wrt==base_)out[0][c]=1.;
  else if(!energy_){
   out[0][c]=wrt==pressure_?area[0][c]*water_moles_per_m3*dh:wrt==area_?retained:0.;
  }else {
   double f=field(fraction_)[0][c],ul=field(liquid_)[0][c],ui=field(ice_)[0][c];
   if(wrt==pressure_)out[0][c]=area[0][c]*water_moles_per_m3*dh*(f*ul+(1-f)*ui);
   else if(wrt==temperature_)out[0][c]=area[0][c]*props[1];
   else if(wrt==fraction_)out[0][c]=area[0][c]*retained*(ul-ui);
   else if(wrt==liquid_)out[0][c]=area[0][c]*retained*f;
   else if(wrt==ice_)out[0][c]=area[0][c]*retained*(1-f);
   else if(wrt==area_)out[0][c]=retained*(f*ul+(1-f)*ui)+props[1]*(field(temperature_)[0][c]-273.15);
   else out[0][c]=0.; // native material is primary and frozen during an ATS solve
  }
 }
}
} }
