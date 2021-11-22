/*
  AdditiveEvaluator is the generic evaluator for adding N other fields.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "AdditiveEvaluator.hh"

namespace Amanzi {
namespace Relations {

AdditiveEvaluator::AdditiveEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::Array<std::string> names;
  if (!plist.isParameter("evaluator dependencies")) {
    if (plist.isParameter("evaluator dependency suffixes")) {
      names = plist_.get<Teuchos::Array<std::string> >("evaluator dependency suffixes");
      Key domain = Keys::getDomain(my_key_);
      for (auto name : names) {
        Key varname = Keys::getKey(domain, name);
        dependencies_.insert(varname);
        Key pname = name + std::string(" coefficient");
        coefs_[varname] = plist.get<double>(pname, 1.0);
      }
    } else {
      Errors::Message msg;
      msg << "AdditiveEvaluator for: \"" << my_key_ << "\" has no dependencies.";
      Exceptions::amanzi_throw(msg);
    }
  } else {
    names = plist.get<Teuchos::Array<std::string> >("evaluator dependencies");
    for (auto name : names) {
      Key pname = name + std::string(" coefficient");
      coefs_[name] = plist.get<double>(pname, 1.0);
    }  
  }

  shift_ = plist.get<double>("constant shift", 0.);
}


Teuchos::RCP<FieldEvaluator>
AdditiveEvaluator::Clone() const
{
  return Teuchos::rcp(new AdditiveEvaluator(*this));
}

// Required methods from SecondaryVariableFieldEvaluator
void
AdditiveEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  result->PutScalar(shift_);

  for (std::map<Key, double>::const_iterator it=coefs_.begin();
       it!=coefs_.end(); ++it) {
    Teuchos::RCP<const CompositeVector> dep = S->GetFieldData(it->first);
    result->Update(it->second, *dep, 1.0);
  }
}

void
AdditiveEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  result->PutScalar(coefs_[wrt_key]);
}


} // namespace
} // namespace

