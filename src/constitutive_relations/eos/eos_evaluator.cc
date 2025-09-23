/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOSEvaluator is the interface between state/data and the model, an EOS.

*/

#include "eos_factory.hh"
#include "eos_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

void
EOSEvaluator::ParsePlistKeys_()
{
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
  AMANZI_ASSERT(my_keys_.size() > 0);
  KeyTag key_tag = my_keys_.front();
  my_keys_.clear();
  Key key = key_tag.first;
  Tag tag = key_tag.second;
  Key domain = Keys::getDomain(key);
  Key varname = Keys::getVarName(key);

  if (mode_ == EOS_MODE_MOLAR || mode_ == EOS_MODE_BOTH) {
    std::size_t molar_pos = varname.find("molar");
    if (molar_pos != std::string::npos) {
      Key molar_key = Keys::readKey(plist_, domain, "molar density", varname);
      my_keys_.emplace_back(KeyTag{ molar_key, tag });
    } else {
      std::size_t mass_pos = varname.find("mass");
      if (mass_pos != std::string::npos) {
        Key molar_key =
          varname.substr(0, mass_pos) + "molar" + varname.substr(mass_pos + 4, varname.size());
        molar_key = Keys::readKey(plist_, domain, "molar density", molar_key);
        my_keys_.emplace_back(KeyTag{ molar_key, tag });
      } else {
        Key molar_key = Keys::readKey(plist_, domain, "molar density");
        my_keys_.emplace_back(KeyTag{ molar_key, tag });
      }
    }
  }

  if (mode_ == EOS_MODE_MASS || mode_ == EOS_MODE_BOTH) {
    std::size_t mass_pos = varname.find("mass");
    if (mass_pos != std::string::npos) {
      Key mass_key = Keys::readKey(plist_, domain, "mass density", varname);
      my_keys_.emplace_back(KeyTag{ mass_key, tag });
    } else {
      std::size_t molar_pos = varname.find("molar");
      if (molar_pos != std::string::npos) {
        Key mass_key =
          varname.substr(0, molar_pos) + "mass" + varname.substr(molar_pos + 5, varname.size());
        mass_key = Keys::readKey(plist_, domain, "mass density", mass_key);
        my_keys_.emplace_back(KeyTag{ mass_key, tag });
      } else {
        Key mass_key = Keys::readKey(plist_, domain, "mass density");
        my_keys_.emplace_back(KeyTag{ mass_key, tag });
      }
    }
  }
}

void
EOSEvaluator::ParsePlistTemp_()
{
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // -- temperature
  temp_key_ = Keys::readKeyTag(plist_, domain_name, "temperature", "temperature", tag);
  dependencies_.insert(temp_key_);
}


void
EOSEvaluator::ParsePlistPres_()
{
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // -- pressure
  pres_key_ = Keys::readKeyTag(plist_, domain_name, "pressure", "pressure", tag);
  dependencies_.insert(pres_key_);
}


void
EOSEvaluator::ParsePlistConc_()
{
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // -- concentration
  conc_key_ = Keys::readKeyTag(plist_, domain_name, "mole fraction", "mole_fraction", tag);
  dependencies_.insert(conc_key_);
}


EOSEvaluator::EOSEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  ParsePlistKeys_();

  // Construct my EOS model
  EOSFactory eos_fac;
  eos_ = eos_fac.createEOS(plist_.sublist("EOS parameters"));

  if (eos_->IsMoleFraction() ) ParsePlistConc_();
  if (eos_->IsTemperature() ) ParsePlistTemp_();
  if (eos_->IsPressure() ) ParsePlistPres_();

  // -- logging
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo_.getOSTab();
    for (const auto& dep : dependencies_) {
      *vo_.os() << " dep: " << dep.first << "@" << dep.second << std::endl;
    }
  }
};


Teuchos::RCP<Evaluator>
EOSEvaluator::Clone() const
{
  return Teuchos::rcp(new EOSEvaluator(*this));
}


void
EOSEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  int num_dep = dependencies_.size();
  std::vector<double> eos_params(num_dep);
  std::vector<const CompositeVector*> dep_cv;
  std::vector<const Epetra_MultiVector*> dep_vec(num_dep, nullptr);

  // Pull dependencies out of state.
  auto tag = my_keys_.front().second;
  if (eos_->IsMoleFraction())
    dep_cv.emplace_back(S.GetPtr<CompositeVector>(conc_key_.first, conc_key_.second).get());
  if (eos_->IsTemperature())
    dep_cv.emplace_back(S.GetPtr<CompositeVector>(temp_key_.first, temp_key_.second).get());
  if (eos_->IsPressure())
    dep_cv.emplace_back(S.GetPtr<CompositeVector>(pres_key_.first, pres_key_.second).get());

  CompositeVector* molar_dens(nullptr);
  CompositeVector* mass_dens(nullptr);
  if (mode_ == EOS_MODE_MOLAR) {
    molar_dens = results[0];
  } else if (mode_ == EOS_MODE_MASS) {
    mass_dens = results[0];
  } else {
    molar_dens = results[0];
    mass_dens = results[1];
  }

  if (molar_dens != nullptr) {
    // evaluate MolarDensity()
    for (CompositeVector::name_iterator comp = molar_dens->begin() ; comp != molar_dens->end();
         ++comp) {
      for (int k = 0; k < num_dep; k++) {
        dep_vec[k] = dep_cv[k]->ViewComponent(*comp, false).get();
      }

      auto& dens_v = *(molar_dens->ViewComponent(*comp, false));
      int count = dens_v.MyLength();
      for (int id = 0; id != count; ++id) {
        for (int k = 0; k < num_dep; k++) {
          eos_params[k] = (*dep_vec[k])[0][id];
        }
        dens_v[0][id] = eos_->MolarDensity(eos_params);
        AMANZI_ASSERT(dens_v[0][id] > 0);
      }
    }
  }

  if (mass_dens != nullptr) {
    for (CompositeVector::name_iterator comp = mass_dens->begin() ; comp != mass_dens->end();
         ++comp) {
      if (mode_ == EOS_MODE_BOTH && eos_->IsConstantMolarMass() &&
          molar_dens->HasComponent(*comp)) {
        // calculate MassDensity from MolarDensity and molar mass.
        double M = eos_->MolarMass();
        mass_dens->ViewComponent(*comp, false)
          ->Update(M, *molar_dens->ViewComponent(*comp, false), 0.);
      } else {
        // evaluate MassDensity() directly
        for (int k = 0; k < num_dep; k++) {
          dep_vec[k] = dep_cv[k]->ViewComponent(*comp, false).get();
        }

        auto& dens_v = *(mass_dens->ViewComponent(*comp, false));
        int count = dens_v.MyLength();
        for (int id = 0; id != count; ++id) {
          for (int k = 0; k < num_dep; k++) eos_params[k] = (*dep_vec[k])[0][id];
          dens_v[0][id] = eos_->MassDensity(eos_params);
          AMANZI_ASSERT(dens_v[0][id] > 0);
        }
      }
    }
  }

#ifdef ENABLE_DBC
  for (const auto& vec : results) {
    double min_val = 1.;
    vec->MinValue(&min_val);
    if (min_val <= 0) {
      Errors::Message msg(
        "EOSEvaluator: input data resulted in density calculation out of range.  If this is at the "
        "start of a run, perhaps you forgot to initialize a temperature or other value to non-zero "
        "(e.g. on \"boundary_face\" or \"face\" components?)");
      Exceptions::amanzi_throw(msg);
    }
  }
#endif
}


void
EOSEvaluator::EvaluatePartialDerivative_(const State& S,
                                         const Key& wrt_key,
                                         const Tag& wrt_tag,
                                         const std::vector<CompositeVector*>& results)
{
  int num_dep = dependencies_.size();
  std::vector<double> eos_params(num_dep);
  std::vector<const CompositeVector*> dep_cv;
  std::vector<const Epetra_MultiVector*> dep_vec(num_dep, nullptr);
  KeyTag wrt{ wrt_key, wrt_tag };

  // Pull dependencies out of state.
  auto tag = my_keys_.front().second;
  if (eos_->IsMoleFraction())
    dep_cv.emplace_back(S.GetPtr<CompositeVector>(conc_key_.first, conc_key_.second).get());
  if (eos_->IsTemperature())
    dep_cv.emplace_back(S.GetPtr<CompositeVector>(temp_key_.first, temp_key_.second).get());
  if (eos_->IsPressure())
    dep_cv.emplace_back(S.GetPtr<CompositeVector>(pres_key_.first, pres_key_.second).get());

  CompositeVector* molar_dens(nullptr);
  CompositeVector* mass_dens(nullptr);
  if (mode_ == EOS_MODE_MOLAR) {
    molar_dens = results[0];
  } else if (mode_ == EOS_MODE_MASS) {
    mass_dens = results[0];
  } else {
    molar_dens = results[0];
    mass_dens = results[1];
  }

  if (molar_dens != nullptr) {
    // evaluate MolarDensity()
    for (CompositeVector::name_iterator comp = molar_dens->begin() ; comp != molar_dens->end();
         ++comp) {
      for (int k = 0; k < num_dep; k++) {
        dep_vec[k] = dep_cv[k]->ViewComponent(*comp, false).get();
      }

      auto& dens_v = *(molar_dens->ViewComponent(*comp, false));
      int count = dens_v.MyLength();

      if (wrt == conc_key_) {
        for (int id = 0; id != count; ++id) {
          for (int k = 0; k < num_dep; k++) {
            eos_params[k] = (*dep_vec[k])[0][id];
          }
          dens_v[0][id] = eos_->DMolarDensityDMoleFraction(eos_params);
        }
      } else if (wrt == pres_key_) {
        for (int id = 0; id != count; ++id) {
          for (int k = 0; k < num_dep; k++) {
            eos_params[k] = (*dep_vec[k])[0][id];
          }
          dens_v[0][id] = eos_->DMolarDensityDp(eos_params);
        }
      } else if (wrt == temp_key_) {
        for (int id = 0; id != count; ++id) {
          for (int k = 0; k < num_dep; k++) {
            eos_params[k] = (*dep_vec[k])[0][id];
          }
          dens_v[0][id] = eos_->DMolarDensityDT(eos_params);
        }
      } else {
        AMANZI_ASSERT(false);
      }
    }
  }

  if (mass_dens != nullptr) {
    for (CompositeVector::name_iterator comp = mass_dens->begin() ; comp != mass_dens->end();
         ++comp) {
      if (mode_ == EOS_MODE_BOTH && eos_->IsConstantMolarMass() &&
          molar_dens->HasComponent(*comp)) {
        // calculate MassDensity from MolarDensity and molar mass.
        double M = eos_->MolarMass();
        mass_dens->ViewComponent(*comp, false)
          ->Update(M, *molar_dens->ViewComponent(*comp, false), 0.);
      } else {
        // evaluate DMassDensity() directly
        for (int k = 0; k < num_dep; k++) {
          dep_vec[k] = dep_cv[k]->ViewComponent(*comp, false).get();
        }

        auto& dens_v = *(mass_dens->ViewComponent(*comp, false));
        int count = dens_v.MyLength();

        if (wrt == conc_key_) {
          for (int id = 0; id != count; ++id) {
            for (int k = 0; k < num_dep; k++) eos_params[k] = (*dep_vec[k])[0][id];
            dens_v[0][id] = eos_->DMassDensityDMoleFraction(eos_params);
          }
        } else if (wrt == pres_key_) {
          for (int id = 0; id != count; ++id) {
            for (int k = 0; k < num_dep; k++) eos_params[k] = (*dep_vec[k])[0][id];
            dens_v[0][id] = eos_->DMassDensityDp(eos_params);
          }
        } else if (wrt == temp_key_) {
          for (int id = 0; id != count; ++id) {
            for (int k = 0; k < num_dep; k++) eos_params[k] = (*dep_vec[k])[0][id];
            dens_v[0][id] = eos_->DMassDensityDT(eos_params);
          }
        } else {
          AMANZI_ASSERT(false);
        }
      }
    }
  }
}


} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi
