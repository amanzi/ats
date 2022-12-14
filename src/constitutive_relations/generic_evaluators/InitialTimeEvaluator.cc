/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "InitialTimeEvaluator.hh"
#include "Units.hh"

namespace Amanzi {
namespace Relations {


InitialTimeEvaluator::InitialTimeEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist), evaluated_once_(false)
{
  auto name = my_keys_.front().first;
  auto tag = my_keys_.front().second;

  // get the default
  auto domain = Keys::getDomain(name);
  auto varname = Keys::getVarName(name);
  Key dep_key;
  if (Keys::starts_with(varname, "init_")) {
    varname = varname.substr(std::string("init_").size(), std::string::npos);
    dep_key = Keys::readKey(plist, domain, "dependency", varname);
  } else if (Keys::starts_with(varname, "initial_")) {
    varname = varname.substr(std::string("initial_").size(), std::string::npos);
    dep_key = Keys::readKey(plist, domain, "dependency", varname);
  } else {
    dep_key = Keys::readKey(plist, domain, "dependency");
  }
  dependencies_.insert(KeyTag(dep_key, tag));

  time_ = plist_.get<double>("initial time", 0.);
  auto units_str = plist_.get<std::string>("initial time units", "s");

  Utils::Units units;
  bool success(false);
  time_ = units.ConvertTime(time_, units_str, "s", success);
  if (!success) {
    Errors::Message msg;
    msg << "InitialTimeEvaluator for \"" << name << "@" << tag.get()
        << "\" provided invalid time unit \"" << units_str << ", valid are "
        << units.ValidTimeStrings();
    Exceptions::amanzi_throw(msg);
  }
}


Teuchos::RCP<Evaluator>
InitialTimeEvaluator::Clone() const
{
  return Teuchos::rcp(new InitialTimeEvaluator(*this));
}


// evaluate once, at the right time
void
InitialTimeEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  auto dep = dependencies_.front();
  if (S.get_time(tag) == time_ && !evaluated_once_) {
    *result[0] = S.Get<CompositeVector>(dep.first, dep.second);
    evaluated_once_ = true;
  }
}


// must be checkpointed
void
InitialTimeEvaluator::EnsureCompatibility_Flags_(State& S)
{
  EvaluatorSecondaryMonotypeCV::EnsureCompatibility_Flags_(S);
  auto keytag = my_keys_.front();
  S.GetRecordW(keytag.first, keytag.second, keytag.first).set_io_checkpoint(true);
}


} // namespace Relations
} // namespace Amanzi
