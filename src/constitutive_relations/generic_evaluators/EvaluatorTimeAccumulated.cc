/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "EvaluatorTimeAccumulated.hh"

namespace Amanzi {
namespace Relations {

EvaluatorTimeAccumulated::EvaluatorTimeAccumulated(Teuchos::ParameterList& plist)
  : EvaluatorSecondary(plist)
{
  // my_keys_[0] was populated by EvaluatorSecondary from the plist name/tag.
  const auto& my_key = my_keys_.front().first;
  const auto& my_tag = my_keys_.front().second;

  accumulation_type_ = plist_.get<std::string>("accumulation type", "integral");
  if (accumulation_type_ != "integral" && accumulation_type_ != "min" &&
      accumulation_type_ != "max") {
    Errors::Message msg;
    msg << "EvaluatorTimeAccumulated for \"" << my_key << "\": invalid \"accumulation type\" \""
        << accumulation_type_ << "\", must be one of: integral, min, max";
    Exceptions::amanzi_throw(msg);
  }

  auto domain = Keys::getDomain(my_key);
  accumulated_key_ = Keys::readKey(plist_, domain, "accumulated", Keys::getVarName(my_key));

  if (!plist_.isParameter("accumulated tag")) {
    Errors::Message msg;
    msg << "EvaluatorTimeAccumulated for \"" << my_key
        << "\": missing required parameter \"accumulated tag\"";
    Exceptions::amanzi_throw(msg);
  }
  accumulated_tag_ = Tag(plist_.get<std::string>("accumulated tag"));

  // my_keys_[1]: the accumulated-dt scalar, same tag as the CV result
  Key acc_dt_key = my_key + "_accumulated_dt";
  my_keys_.emplace_back(KeyTag{ acc_dt_key, my_tag });

  // accumulated_key@accumulated_tag changes each inner step — triggers Update_()
  dependencies_.insert(KeyTag(accumulated_key_, accumulated_tag_));

  // no evaluator for dt, so can't add it to dependencies
  // dt@accumulated_tag provides the inner timestep size for weighting
  //dependencies_.insert(KeyTag("dt", accumulated_tag_));
}


Teuchos::RCP<Evaluator>
EvaluatorTimeAccumulated::Clone() const
{
  return Teuchos::rcp(new EvaluatorTimeAccumulated(*this));
}


void
EvaluatorTimeAccumulated::EnsureCompatibility(State& S)
{
  // my_keys_[0]: CV result — claim ownership; structure flows from client requirements
  const auto& [cv_key, cv_tag] = my_keys_[0];
  auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(cv_key, cv_tag, cv_key);

  // my_keys_[1]: accumulated-dt scalar — claim ownership
  const auto& [dt_key, dt_tag] = my_keys_[1];
  S.Require<double>(dt_key, dt_tag, dt_key);

  // wire dependency evaluators into the graph
  if (my_fac.Mesh() != Teuchos::null) {
    for (const auto& dep : dependencies_)
      // dependency fac includes my fac
      S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second)
        .Update(my_fac);
  }

  EvaluatorSecondary::EnsureCompatibility_DepEnsureCompatibility_(S);
  EvaluatorSecondary::EnsureCompatibility_Flags_(S);
}


void
EvaluatorTimeAccumulated::Update_(State& S)
{
  const auto& [cv_key, cv_tag] = my_keys_[0];
  const auto& [dt_key, dt_tag] = my_keys_[1];

  double dt_inner = S.Get<double>("dt", accumulated_tag_);
  const auto& source = S.Get<CompositeVector>(accumulated_key_, accumulated_tag_);
  double acc_dt = S.Get<double>(dt_key, dt_tag);
  if (dt_inner > 0) {
    // cv is read and overwritten in place — valid because we own it
    auto& cv = S.GetW<CompositeVector>(cv_key, cv_tag, cv_key);

    if (accumulation_type_ == "integral") {
      double total_dt = acc_dt + dt_inner;
      for (const auto& comp : cv) {
        auto& res = *cv.ViewComponent(comp, false);
        const auto& src = *source.ViewComponent(comp, false);
        for (int j = 0; j != res.NumVectors(); ++j)
          for (int i = 0; i != res.MyLength(); ++i)
            res[j][i] = (res[j][i] * acc_dt + src[j][i] * dt_inner) / total_dt;
      }
    } else if (accumulation_type_ == "max") {
      for (const auto& comp : cv) {
        auto& res = *cv.ViewComponent(comp, false);
        const auto& src = *source.ViewComponent(comp, false);
        for (int j = 0; j != res.NumVectors(); ++j)
          for (int i = 0; i != res.MyLength(); ++i)
            res[j][i] = std::max(res[j][i], src[j][i]);
      }
    } else { // min
      for (const auto& comp : cv) {
        auto& res = *cv.ViewComponent(comp, false);
        const auto& src = *source.ViewComponent(comp, false);
        for (int j = 0; j != res.NumVectors(); ++j)
          for (int i = 0; i != res.MyLength(); ++i)
            res[j][i] = std::min(res[j][i], src[j][i]);
      }
    }

    S.GetW<double>(dt_key, dt_tag, dt_key) = acc_dt + dt_inner;

    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      Teuchos::OSTab tab = vo_.getOSTab();
      *vo_.os() << "Updated time accumulation time from " << acc_dt << " to " << acc_dt + dt_inner << std::endl;
      cv.Print(std::cout);
    }
  }
}


void
EvaluatorTimeAccumulated::Reset(State& S)
{
  const auto& [cv_key, cv_tag] = my_keys_[0];
  const auto& [dt_key, dt_tag] = my_keys_[1];

  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo_.getOSTab();
    *vo_.os() << "Resetting time accumulation" << std::endl;
  }

  auto& cv = S.GetW<CompositeVector>(cv_key, cv_tag, cv_key);
  if (accumulation_type_ == "max") {
    cv.PutScalar(-1.e16);
  } else if (accumulation_type_ == "min") {
    cv.PutScalar(1.e16);
  } else {
    cv.PutScalar(0.);
  }

  S.GetW<double>(dt_key, dt_tag, dt_key) = 0.;

  // invalidate lazy cache so the next Update() call always runs Update_()
  requests_.clear();
}


} // namespace Relations
} // namespace Amanzi
