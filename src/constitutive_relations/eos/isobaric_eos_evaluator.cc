/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EOSFieldEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_factory.hh"
#include "isobaric_eos_evaluator.hh"

namespace Amanzi {
namespace Relations {

IsobaricEOSEvaluator::IsobaricEOSEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariablesFieldEvaluator(plist) {

  // Process the list for my provided field.
  std::string mode = plist_.get<std::string>("EOS basis", "molar");
  if (mode == "molar") {
    mode_ = EOS_MODE_MOLAR;
  } else if (mode == "mass") {
    mode_ = EOS_MODE_MASS;
  } else if (mode == "both") {
    mode_ = EOS_MODE_BOTH;
  } else {
    AMANZI_ASSERT(0);
  }

  // my keys
  if (mode_ == EOS_MODE_MOLAR || mode_ == EOS_MODE_BOTH) {
    a_key_ = plist_.get<std::string>("molar density key");
    my_keys_.push_back(a_key_);
  }

  if (mode_ == EOS_MODE_MASS || mode_ == EOS_MODE_BOTH) {
    a_key_ = plist_.get<std::string>("mass density key");
    my_keys_.push_back(a_key_);
  }

  // Set up my dependencies.
  std::size_t end = a_key_.find_first_of("_");
  std::string domain_name = a_key_.substr(0,end);
  if (domain_name == std::string("density") ||
      domain_name == std::string("molar") ||
      domain_name == std::string("mass")) {
    domain_name = std::string("");
  } else {
    domain_name = domain_name+std::string("_");
  }

  // -- dependence
  dep_key_ = plist_.get<std::string>("dependence key",
          domain_name+std::string("temperature"));
  dependencies_.insert(dep_key_);

  // -- pressure
  pres_key_ = plist_.get<std::string>("pressure key", "atmospheric_pressure");

  // Construct my EOS model
  AMANZI_ASSERT(plist_.isSublist("EOS parameters"));
  EOSFactory eos_fac;
  eos_ = eos_fac.createEOS(plist_.sublist("EOS parameters"));
};


IsobaricEOSEvaluator::IsobaricEOSEvaluator(const IsobaricEOSEvaluator& other) :
    SecondaryVariablesFieldEvaluator(other),
    eos_(other.eos_),
    mode_(other.mode_),
    dep_key_(other.dep_key_),
    pres_key_(other.pres_key_) {}


Teuchos::RCP<FieldEvaluator> IsobaricEOSEvaluator::Clone() const {
  return Teuchos::rcp(new IsobaricEOSEvaluator(*this));
}


void IsobaricEOSEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> dep_cv = S->GetPtr<CompositeVector>(dep_key_);
  Teuchos::RCP<const double> pres = S->GetScalarData(pres_key_);
  std::vector<double> eos_params(2);

  int index = 0; // index to the results list
  if (mode_ == EOS_MODE_MOLAR || mode_ == EOS_MODE_BOTH) {
    // evaluate MolarDensity()
    Teuchos::Ptr<CompositeVector> result = results[index];
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& dep_v = *(dep_cv->ViewComponent(*comp,false));
      Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

      int count = result->size(*comp);
      for (int id=0; id!=count; ++id) {
        eos_params[0] = dep_v[0][id];
        eos_params[1] = *pres;               
        result_v[0][id] = eos_->MolarDensity(eos_params);
      }
    }
    index++;
  }

  if (mode_ == EOS_MODE_BOTH && eos_->IsConstantMolarMass()) {
    // calculate MassDensity from MolarDensity and molar mass.
    double M = eos_->MolarMass();
    results[1]->Update(M, *(results[0]), 0.0);
  } else if (mode_ == EOS_MODE_MASS || mode_ == EOS_MODE_BOTH) {
    // evaluate MassDensity()
    Teuchos::Ptr<CompositeVector> result = results[index];
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& dep_v = *(dep_cv->ViewComponent(*comp,false));
      Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

      int count = result->size(*comp);
      for (int id=0; id!=count; ++id) {
        eos_params[0] = dep_v[0][id];
        eos_params[1] = *pres;        
        result_v[0][id] = eos_->MassDensity(eos_params);
      }
    }
  }
}


void IsobaricEOSEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results) {

  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> dep_cv = S->GetPtr<CompositeVector>(dep_key_);
  Teuchos::RCP<const double> pres = S->GetScalarData(pres_key_);
  std::vector<double> eos_params(2);  

  if (wrt_key == dep_key_) {

    int index = 0; // index to the results list
    if (mode_ == EOS_MODE_MOLAR || mode_ == EOS_MODE_BOTH) {
      // evaluate DMolarDensityDT()
      Teuchos::Ptr<CompositeVector> result = results[index];
      for (CompositeVector::name_iterator comp=result->begin();
           comp!=result->end(); ++comp) {
        const Epetra_MultiVector& dep_v = *(dep_cv->ViewComponent(*comp,false));
        Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

        int count = result->size(*comp);
        for (int id=0; id!=count; ++id) {
          eos_params[0] = dep_v[0][id];
          eos_params[1] = *pres;       
          result_v[0][id] = eos_->DMolarDensityDT(eos_params);
        }
      }
      index++;
    }

    if (mode_ == EOS_MODE_BOTH && eos_->IsConstantMolarMass()) {
      // calculate DMassDensityDT from DMolarDensityDT and molar mass.
      double M = eos_->MolarMass();
      results[1]->Update(M, *results[0], 0.0);
    } else if (mode_ == EOS_MODE_MASS || mode_ == EOS_MODE_BOTH) {
      // evaluate DMassDensityDT()
      Teuchos::Ptr<CompositeVector> result = results[index];
      for (CompositeVector::name_iterator comp=result->begin();
           comp!=result->end(); ++comp) {
        const Epetra_MultiVector& dep_v = *(dep_cv->ViewComponent(*comp,false));
        Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

        int count = result->size(*comp);
        for (int id=0; id!=count; ++id) {
          eos_params[0] = dep_v[0][id];
          eos_params[1] = *pres;                 
          result_v[0][id] = eos_->DMassDensityDT(eos_params);
        }
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}

} // namespace
} // namespace
