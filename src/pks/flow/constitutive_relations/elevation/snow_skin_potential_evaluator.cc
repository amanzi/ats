/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the surface elevation, slope, and updates pres + elev.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "boost/algorithm/string/predicate.hpp"
#include "snow_skin_potential_evaluator.hh"

namespace Amanzi {
namespace Flow {

SnowSkinPotentialEvaluator::SnowSkinPotentialEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  Key domain = Keys::getDomain(my_key_);
  Key surf_domain;
  if (domain == "snow") {
    surf_domain = "surface";
    surf_domain = plist_.get<std::string>("surface domain", surf_domain);
  } else if (boost::starts_with(domain, "snow")) {
    surf_domain = Key("surface")+domain.substr(4,domain.size());
  } else {
    surf_domain = "surface";
  }
  
  pd_key_ = Keys::readKey(plist_, surf_domain, "ponded depth", "ponded_depth");
  dependencies_.insert(pd_key_);

  sd_key_ = Keys::readKey(plist_, domain, "snow depth", "depth");
  dependencies_.insert(sd_key_);

  precip_key_ = Keys::readKey(plist_, domain, "snow precipitation", "precipitation");
  dependencies_.insert(precip_key_);

  elev_key_ = Keys::readKey(plist_, surf_domain, "elevation", "elevation");
  dependencies_.insert(elev_key_);

  factor_ = plist_.get<double>("dt factor [s]", -1.0);
}


SnowSkinPotentialEvaluator::SnowSkinPotentialEvaluator(const SnowSkinPotentialEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    elev_key_(other.elev_key_),
    pd_key_(other.pd_key_),
    sd_key_(other.sd_key_),
    precip_key_(other.precip_key_),
    factor_(other.factor_)
{};


Teuchos::RCP<FieldEvaluator>
SnowSkinPotentialEvaluator::Clone() const {
  return Teuchos::rcp(new SnowSkinPotentialEvaluator(*this));
}


void SnowSkinPotentialEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  // update pressure + elevation
  Teuchos::RCP<const CompositeVector> pd = S->GetFieldData(pd_key_);
  Teuchos::RCP<const CompositeVector> sd = S->GetFieldData(sd_key_);
  Teuchos::RCP<const CompositeVector> precip = S->GetFieldData(precip_key_);
  Teuchos::RCP<const CompositeVector> elev = S->GetFieldData(elev_key_);


  // note factor of 10 accounts for change from precip in m SWE to actual m.
  result->Update(1.0, *elev, 1.0, *pd, 0.0);
  if (factor_ > 0.) {
    result->Update(1.0, *sd, 10*factor_, *precip, 1.0);
  } else {
    double dt = S->time() - S->last_time();
    result->Update(1.0, *sd, 10*dt, *precip, 1.0);
  }
}


// This is hopefully never called?
void SnowSkinPotentialEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  AMANZI_ASSERT(0);
  result->PutScalar(1.0);
}

} //namespace
} //namespace
