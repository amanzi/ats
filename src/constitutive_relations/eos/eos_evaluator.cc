/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EOSFieldEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_factory.hh"
#include "eos_evaluator.hh"

namespace Amanzi {
namespace Relations {

EOSEvaluator::EOSEvaluator(Teuchos::ParameterList& plist) :
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
  Key name = plist_.get<std::string>("evaluator name");
  if (mode_ == EOS_MODE_MOLAR || mode_ == EOS_MODE_BOTH) {
    std::size_t molar_pos = name.find("molar");
    if (molar_pos != std::string::npos) {
      Key molar_key = plist_.get<std::string>("molar density key", name);
      my_keys_.push_back(molar_key);
    } else {
      std::size_t mass_pos = name.find("mass");
      if (mass_pos != std::string::npos) {
        Key molar_key = name.substr(0,mass_pos)+"molar"+name.substr(mass_pos+4, name.size());
        molar_key = plist_.get<std::string>("molar density key", molar_key);
        my_keys_.push_back(molar_key);
      } else {
        Key molar_key = plist_.get<std::string>("molar density key");
        my_keys_.push_back(molar_key);
      }
    }
  }

  if (mode_ == EOS_MODE_MASS || mode_ == EOS_MODE_BOTH) {
    std::size_t mass_pos = name.find("mass");
    if (mass_pos != std::string::npos) {
      Key mass_key = plist_.get<std::string>("mass density key", name);
      my_keys_.push_back(mass_key);
    } else {
      std::size_t molar_pos = name.find("molar");
      if (molar_pos != std::string::npos) {
        Key mass_key = name.substr(0,molar_pos)+"mass"+name.substr(molar_pos+5, name.size());
        mass_key = plist_.get<std::string>("mass density key", mass_key);
        my_keys_.push_back(mass_key);
      } else {
        Key mass_key = plist_.get<std::string>("mass density key");
        my_keys_.push_back(mass_key);
      }
    }
  }

  // Set up my dependencies.
  Key domain_name = Keys::getDomain(name);

  // -- logging
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    for (KeySet::const_iterator dep=dependencies_.begin();
         dep!=dependencies_.end(); ++dep) {
      *vo_->os() << " dep: " << *dep << std::endl;
    }
  }

  // Construct my EOS model
  AMANZI_ASSERT(plist_.isSublist("EOS parameters"));
  EOSFactory eos_fac;
  eos_ = eos_fac.createEOS(plist_.sublist("EOS parameters"));
};


EOSEvaluator::EOSEvaluator(const EOSEvaluator& other) :
    SecondaryVariablesFieldEvaluator(other),
    eos_(other.eos_),
    mode_(other.mode_)
 {}


Teuchos::RCP<FieldEvaluator> EOSEvaluator::Clone() const {
  return Teuchos::rcp(new EOSEvaluator(*this));
}


void EOSEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  
  int num_dep = dependencies_.size();  
  std::vector<double> eos_params(num_dep);
  std::vector< Teuchos::RCP<const CompositeVector> > dep_cv(num_dep);
  std::vector< Teuchos::RCP<const Epetra_MultiVector> > dep_vec(num_dep);

  // Pull dependencies out of state.  
  KeySet::const_iterator dep;
  int k;
  for (k=0, dep=dependencies_.begin();
       dep!=dependencies_.end(); ++k, ++dep)  {
    dep_cv[k] = S->GetFieldData(*dep);
  }

  Teuchos::Ptr<CompositeVector> molar_dens, mass_dens;
  if (mode_ == EOS_MODE_MOLAR) {
    molar_dens = results[0];
  } else if (mode_ == EOS_MODE_MASS) {
    mass_dens = results[0];
  } else {
    molar_dens = results[0];
    mass_dens = results[1];
  }

  if (molar_dens != Teuchos::null) {
    // evaluate MolarDensity()
    for (CompositeVector::name_iterator comp=molar_dens->begin();
         comp!=molar_dens->end(); ++comp) {
      for (k=0; k<num_dep; k++){
        dep_vec[k] = dep_cv[k]->ViewComponent(*comp,false);
      }
      
      Epetra_MultiVector& dens_v = *(molar_dens->ViewComponent(*comp,false));
      int count = dens_v.MyLength();
      for (int id=0; id!=count; ++id) {
        for (k=0; k<num_dep; k++) {
          eos_params[k] = (*dep_vec[k])[0][id];
        }
        dens_v[0][id] = eos_->MolarDensity(eos_params);
      }
    }
  }

  if (mass_dens != Teuchos::null) {
    for (CompositeVector::name_iterator comp=mass_dens->begin();
         comp!=mass_dens->end(); ++comp) {
      if (mode_ == EOS_MODE_BOTH && eos_->IsConstantMolarMass() &&
          molar_dens->HasComponent(*comp)) {
        // calculate MassDensity from MolarDensity and molar mass.
        double M = eos_->MolarMass();
        mass_dens->ViewComponent(*comp,false)->Update(M,
                *molar_dens->ViewComponent(*comp,false), 0.);
      } else {
        // evaluate MassDensity() directly
        for (k=0; k<num_dep; k++){
          dep_vec[k] = dep_cv[k]->ViewComponent(*comp,false);
        }
   
        Epetra_MultiVector& dens_v = *(mass_dens->ViewComponent(*comp,false));
        int count = dens_v.MyLength();
        for (int id=0; id!=count; ++id) {
          for (k=0; k<num_dep; k++) eos_params[k] = (*dep_vec[k])[0][id];  
          dens_v[0][id] = eos_->MassDensity(eos_params);
        }
      }
    }
  }


#ifdef ENABLE_DBC
  for (const auto& vec : results) {
    double min_val = 1.;
    vec->MinValue(&min_val);
    if (min_val <= 0) {
      Errors::Message msg("EOSEvaluator: input data resulted in density calculation out of range.  If this is at the start of a run, perhaps you forgot to initialize a temperature or other value to non-zero (e.g. on \"boundary_face\" or \"face\" components?)");
      Exceptions::amanzi_throw(msg);
    }
  }
#endif
  
}


void EOSEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results)
{
  Errors::Message msg("Derivative computation is missing for these EOSEvaluator");
  Exceptions::amanzi_throw(msg);
}


} // namespace
} // namespace
