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
    EvaluatorSecondaries(plist) {

  tag_ = plist.get<std::string>("tag", "");

  // Process the list for my provided field.
  std::string mode = plist_.get<std::string>("EOS basis", "molar");
  if (mode == "molar") {
    mode_ = EOS_MODE_MOLAR;
  } else if (mode == "mass") {
    mode_ = EOS_MODE_MASS;
  } else if (mode == "both") {
    mode_ = EOS_MODE_BOTH;
  } else {
    ASSERT(0);
  }

  // my keys
  Key name = plist_.get<std::string>("evaluator name");
  if (mode_ == EOS_MODE_MOLAR || mode_ == EOS_MODE_BOTH) {
    std::size_t molar_pos = name.find("molar");
    if (molar_pos != std::string::npos) {
      Key molar_key = plist_.get<std::string>("molar density key", name);
      my_keys_.emplace_back(std::make_pair(molar_key, tag_));
    } else {
      std::size_t mass_pos = name.find("mass");
      if (mass_pos != std::string::npos) {
        Key molar_key = name.substr(0,mass_pos)+"molar"+name.substr(mass_pos+4, name.size());
        molar_key = plist_.get<std::string>("molar density key", molar_key);
        my_keys_.emplace_back(std::make_pair(molar_key, tag_));
      } else {
        Key molar_key = plist_.get<std::string>("molar density key");
        my_keys_.emplace_back(std::make_pair(molar_key, tag_));
      }
    }
  }

  if (mode_ == EOS_MODE_MASS || mode_ == EOS_MODE_BOTH) {
    std::size_t mass_pos = name.find("mass");
    if (mass_pos != std::string::npos) {
      Key mass_key = plist_.get<std::string>("mass density key", name);
      my_keys_.emplace_back(std::make_pair(mass_key, tag_));
    } else {
      std::size_t molar_pos = name.find("molar");
      if (molar_pos != std::string::npos) {
        Key mass_key = name.substr(0,molar_pos)+"mass"+name.substr(molar_pos+5, name.size());
        mass_key = plist_.get<std::string>("mass density key", mass_key);
        my_keys_.emplace_back(std::make_pair(mass_key, tag_));
      } else {
        Key mass_key = plist_.get<std::string>("mass density key");
        my_keys_.emplace_back(std::make_pair(mass_key, tag_));
      }
    }
  }

  // Set up my dependencies.
  Key domain_name = Keys::getDomain(name);

  // -- temperature
  temp_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.emplace_back(std::make_pair(temp_key_, tag_));


  // -- pressure
  pres_key_ = Keys::readKey(plist_, domain_name, "pressure", "effective_pressure");
  dependencies_.emplace_back(std::make_pair(pres_key_, tag_));

  // -- logging
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo_.getOSTab();
    for (auto dep : dependencies_) {
      *vo_.os() << " dep: " << dep.first << std::endl;
    }
  }

  // Construct my EOS model
  ASSERT(plist_.isSublist("EOS parameters"));
  EOSFactory eos_fac;
  eos_ = eos_fac.createEOS(plist_.sublist("EOS parameters"));
};


EOSEvaluator::EOSEvaluator(const EOSEvaluator& other) :
    EvaluatorSecondaries(other),
    eos_(other.eos_),
    mode_(other.mode_),
    temp_key_(other.temp_key_),
    pres_key_(other.pres_key_) {}


Teuchos::RCP<Evaluator> EOSEvaluator::Clone() const {
  //EOSEvaluator test(*this);
  return Teuchos::rcp(new EOSEvaluator(*this));
}


void EOSEvaluator::Update_(State& S) {
  // Pull dependencies out of state.
  // Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  // Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);

  // Teuchos::Ptr<CompositeVector> molar_dens, mass_dens;
  // if (mode_ == EOS_MODE_MOLAR) {
  //   molar_dens = results[0];
  // } else if (mode_ == EOS_MODE_MASS) {
  //   mass_dens = results[0];
  // } else {
  //   molar_dens = results[0];
  //   mass_dens = results[1];
  // }


  // if (molar_dens != Teuchos::null) {
  //   // evaluate MolarDensity()
  //   for (CompositeVector::name_iterator comp=molar_dens->begin();
  //        comp!=molar_dens->end(); ++comp) {
  //     const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
  //     const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
  //     Epetra_MultiVector& dens_v = *(molar_dens->ViewComponent(*comp,false));

  //     int count = dens_v.MyLength();
  //     for (int i=0; i!=count; ++i) {
  //       dens_v[0][i] = eos_->MolarDensity(temp_v[0][i], pres_v[0][i]);
  //       ASSERT(dens_v[0][i] > 0.);
  //     }
  //   }
  // }

  // if (mass_dens != Teuchos::null) {
  //   for (CompositeVector::name_iterator comp=mass_dens->begin();
  //        comp!=mass_dens->end(); ++comp) {
  //     if (mode_ == EOS_MODE_BOTH && eos_->IsConstantMolarMass() &&
  //         molar_dens->HasComponent(*comp)) {
  //       // calculate MassDensity from MolarDensity and molar mass.
  //       double M = eos_->MolarMass();

  //       mass_dens->ViewComponent(*comp,false)->Update(M,
  //               *molar_dens->ViewComponent(*comp,false), 0.);
  //     } else {
  //       // evaluate MassDensity() directly
  //       const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
  //       const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
  //       Epetra_MultiVector& dens_v = *(mass_dens->ViewComponent(*comp,false));

  //       int count = dens_v.MyLength();
  //       for (int i=0; i!=count; ++i) {
  //         dens_v[0][i] = eos_->MassDensity(temp_v[0][i], pres_v[0][i]);
  //         ASSERT(dens_v[0][i] > 0.);
  //       }
  //     }
  //   }
  // }
}


void EOSEvaluator::UpdateDerivative_(State& S, const Key& wrt_key, const Key& wrt_tag){

  // Pull dependencies out of state.
  // Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  // Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);

  // Teuchos::Ptr<CompositeVector> molar_dens, mass_dens;
  // if (mode_ == EOS_MODE_MOLAR) {
  //   molar_dens = results[0];
  // } else if (mode_ == EOS_MODE_MASS) {
  //   mass_dens = results[0];
  // } else {
  //   molar_dens = results[0];
  //   mass_dens = results[1];
  // }

  // if (wrt_key == pres_key_) {
  //   if (molar_dens != Teuchos::null) {
  //     // evaluate MolarDensity()
  //     for (CompositeVector::name_iterator comp=molar_dens->begin();
  //          comp!=molar_dens->end(); ++comp) {
  //       const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
  //       const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
  //       Epetra_MultiVector& dens_v = *(molar_dens->ViewComponent(*comp,false));

  //       int count = dens_v.MyLength();
  //       for (int i=0; i!=count; ++i) {
  //         dens_v[0][i] = eos_->DMolarDensityDp(temp_v[0][i], pres_v[0][i]);
  //       }
  //     }
  //   }

  //   if (mass_dens != Teuchos::null) {
  //     for (CompositeVector::name_iterator comp=mass_dens->begin();
  //          comp!=mass_dens->end(); ++comp) {
  //       if (mode_ == EOS_MODE_BOTH && eos_->IsConstantMolarMass() &&
  //           molar_dens->HasComponent(*comp)) {
  //         // calculate MassDensity from MolarDensity and molar mass.
  //         double M = eos_->MolarMass();

  //         mass_dens->ViewComponent(*comp,false)->Update(M,
  //                 *molar_dens->ViewComponent(*comp,false), 0.);
  //       } else {
  //         // evaluate MassDensity() directly
  //         const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
  //         const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
  //         Epetra_MultiVector& dens_v = *(mass_dens->ViewComponent(*comp,false));

  //         int count = dens_v.MyLength();
  //         for (int i=0; i!=count; ++i) {
  //           dens_v[0][i] = eos_->DMassDensityDp(temp_v[0][i], pres_v[0][i]);
  //         }
  //       }
  //     }
  //   }

  // } else if (wrt_key == temp_key_) {

  //   if (molar_dens != Teuchos::null) {
  //     // evaluate MolarDensity()
  //     for (CompositeVector::name_iterator comp=molar_dens->begin();
  //          comp!=molar_dens->end(); ++comp) {
  //       const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
  //       const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
  //       Epetra_MultiVector& dens_v = *(molar_dens->ViewComponent(*comp,false));

  //       int count = dens_v.MyLength();
  //       for (int i=0; i!=count; ++i) {
  //         dens_v[0][i] = eos_->DMolarDensityDT(temp_v[0][i], pres_v[0][i]);
  //       }
  //     }
  //   }

  //   if (mass_dens != Teuchos::null) {
  //     for (CompositeVector::name_iterator comp=mass_dens->begin();
  //          comp!=mass_dens->end(); ++comp) {
  //       if (mode_ == EOS_MODE_BOTH && eos_->IsConstantMolarMass() &&
  //           molar_dens->HasComponent(*comp)) {
  //         // calculate MassDensity from MolarDensity and molar mass.
  //         double M = eos_->MolarMass();

  //         mass_dens->ViewComponent(*comp,false)->Update(M,
  //                 *molar_dens->ViewComponent(*comp,false), 0.);
  //       } else {
  //         // evaluate MassDensity() directly
  //         const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
  //         const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
  //         Epetra_MultiVector& dens_v = *(mass_dens->ViewComponent(*comp,false));

  //         int count = dens_v.MyLength();
  //         for (int i=0; i!=count; ++i) {
  //           dens_v[0][i] = eos_->DMassDensityDT(temp_v[0][i], pres_v[0][i]);
  //         }
  //       }
  //     }
  //   }

  // } else {
  //   ASSERT(0);
  // }
}

} // namespace
} // namespace
