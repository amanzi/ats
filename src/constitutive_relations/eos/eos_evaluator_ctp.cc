/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EOSFieldEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Daniil Svyatsky(dasvyat@lanl.gov)
*/

#include "eos_factory.hh"
#include "eos_evaluator_ctp.hh"

namespace Amanzi {
namespace Relations {

EOSEvaluatorCTP::EOSEvaluatorCTP(Teuchos::ParameterList& plist) :
  EOSEvaluator(plist) {
  
  Key name = plist_.get<std::string>("evaluator name");
  // Set up my dependencies.
  Key domain_name = Keys::getDomain(name);
  
  // -- pressure
  pres_key_ = Keys::readKey(plist_, domain_name, "pressure", "pressure");
  dependencies_.insert(pres_key_);
  
  // -- temperature
  temp_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(temp_key_);

  // -- concentration
  conc_key_ = Keys::readKey(plist_, domain_name, "concentration", "total_component_concentration");  
  dependencies_.insert(conc_key_);  

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
  eos_sw_ = Teuchos::rcp_dynamic_cast<EOS_SW>(eos_);
};


EOSEvaluatorCTP::EOSEvaluatorCTP(const EOSEvaluatorCTP& other) :
    EOSEvaluator(other),
    temp_key_(other.temp_key_),
    conc_key_(other.conc_key_),    
    pres_key_(other.pres_key_)
 {
   eos_sw_ = other.eos_sw_;
 }


Teuchos::RCP<FieldEvaluator> EOSEvaluatorCTP::Clone() const {
  return Teuchos::rcp(new EOSEvaluatorCTP(*this));
}

void EOSEvaluatorCTP::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  
  int num_dep = dependencies_.size();  
  std::vector<double> eos_params(num_dep);

  Teuchos::RCP<const CompositeVector> conc = S->GetPtr<CompositeVector>(conc_key_);  
  Teuchos::RCP<const CompositeVector> temp = S->GetPtr<CompositeVector>(temp_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetPtr<CompositeVector>(pres_key_);

  // Pull dependencies out of state.  

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

      const Epetra_MultiVector& conc_v = *(conc->ViewComponent(*comp,false));
      const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
      const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
      Epetra_MultiVector& dens_v = *(molar_dens->ViewComponent(*comp,false));

      int count = dens_v.MyLength();
      for (int id=0; id!=count; ++id) {

        eos_params[0] = conc_v[0][id];        
        eos_params[1] = temp_v[0][id];
        eos_params[2] = pres_v[0][id];
       
        dens_v[0][id] = eos_sw_->MolarDensity(eos_params);
        
        AMANZI_ASSERT(dens_v[0][id] > 0.);
      }
    }
  }

  if (mass_dens != Teuchos::null) {
    for (CompositeVector::name_iterator comp=mass_dens->begin();
         comp!=mass_dens->end(); ++comp) {
      if (mode_ == EOS_MODE_BOTH && eos_->IsConstantMolarMass() &&
          molar_dens->HasComponent(*comp)) {
        // calculate MassDensity from MolarDensity and molar mass.
        double M = eos_sw_->MolarMass();

        mass_dens->ViewComponent(*comp,false)->Update(M,
                *molar_dens->ViewComponent(*comp,false), 0.);
      } else {
        // evaluate MassDensity() directly
        const Epetra_MultiVector& conc_v = *(conc->ViewComponent(*comp,false));
        const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
        const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
        Epetra_MultiVector& dens_v = *(mass_dens->ViewComponent(*comp,false));

        int count = dens_v.MyLength();
        for (int id=0; id!=count; ++id) {

          eos_params[0] = conc_v[0][id];        
          eos_params[1] = temp_v[0][id];
          eos_params[2] = pres_v[0][id];
          
          dens_v[0][id] = eos_sw_->MassDensity(eos_params);
          AMANZI_ASSERT(dens_v[0][id] > 0.);
        }
      }
    }
  }
}

  

void EOSEvaluatorCTP::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
                                                   Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  
  int num_dep = dependencies_.size();  
  std::vector<double> eos_params(num_dep);
   
  // Pull dependencies out of state.  
  Teuchos::RCP<const CompositeVector> temp = S->GetPtr<CompositeVector>(temp_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetPtr<CompositeVector>(pres_key_);
  Teuchos::RCP<const CompositeVector> conc = S->GetPtr<CompositeVector>(conc_key_);    

  
  Teuchos::Ptr<CompositeVector> molar_dens, mass_dens;
  if (mode_ == EOS_MODE_MOLAR) {
    molar_dens = results[0];
  } else if (mode_ == EOS_MODE_MASS) {
    mass_dens = results[0];
  } else {
    molar_dens = results[0];
    mass_dens = results[1];
  }

  if (wrt_key == pres_key_) {
    if (molar_dens != Teuchos::null) {
      // evaluate MolarDensity()
      for (CompositeVector::name_iterator comp=molar_dens->begin();
           comp!=molar_dens->end(); ++comp) {
        const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
        const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
        const Epetra_MultiVector& conc_v = *(conc->ViewComponent(*comp,false));        
        Epetra_MultiVector& dens_v = *(molar_dens->ViewComponent(*comp,false));

        int count = dens_v.MyLength();
        for (int i=0; i!=count; ++i) {
          eos_params[0] = conc_v[0][i];          
          eos_params[1] = temp_v[0][i];
          eos_params[2] = pres_v[0][i];
            
          dens_v[0][i] = eos_sw_->DMolarDensityDp(eos_params);
        }
      }
    }

    if (mass_dens != Teuchos::null) {
      for (CompositeVector::name_iterator comp=mass_dens->begin();
           comp!=mass_dens->end(); ++comp) {
        if (mode_ == EOS_MODE_BOTH && eos_sw_->IsConstantMolarMass() &&
            molar_dens->HasComponent(*comp)) {
          // calculate MassDensity from MolarDensity and molar mass.
          double M = eos_sw_->MolarMass();

          mass_dens->ViewComponent(*comp,false)->Update(M,
                  *molar_dens->ViewComponent(*comp,false), 0.);
        } else {
          // evaluate MassDensity() directly
          const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
          const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
          const Epetra_MultiVector& conc_v = *(conc->ViewComponent(*comp,false));                  
          Epetra_MultiVector& dens_v = *(mass_dens->ViewComponent(*comp,false));

          int count = dens_v.MyLength();
          for (int i=0; i!=count; ++i) {
            eos_params[0] = conc_v[0][i];          
            eos_params[1] = temp_v[0][i];
            eos_params[2] = pres_v[0][i];            

            dens_v[0][i] = eos_sw_->DMassDensityDp(eos_params);
          }
        }
      }
    }

  } else if (wrt_key == temp_key_) {

    if (molar_dens != Teuchos::null) {
      // evaluate MolarDensity()
      for (CompositeVector::name_iterator comp=molar_dens->begin();
           comp!=molar_dens->end(); ++comp) {
        const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
        const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
        const Epetra_MultiVector& conc_v = *(conc->ViewComponent(*comp,false));        
        Epetra_MultiVector& dens_v = *(molar_dens->ViewComponent(*comp,false));

        int count = dens_v.MyLength();
        for (int i=0; i!=count; ++i) {
          eos_params[0] = conc_v[0][i];          
          eos_params[1] = temp_v[0][i];
          eos_params[2] = pres_v[0][i];
                                
          dens_v[0][i] = eos_sw_->DMolarDensityDT(eos_params);
        }
      }
    }

    if (mass_dens != Teuchos::null) {
      for (CompositeVector::name_iterator comp=mass_dens->begin();
           comp!=mass_dens->end(); ++comp) {
        if (mode_ == EOS_MODE_BOTH && eos_sw_->IsConstantMolarMass() &&
            molar_dens->HasComponent(*comp)) {
          // calculate MassDensity from MolarDensity and molar mass.
          double M = eos_sw_->MolarMass();

          mass_dens->ViewComponent(*comp,false)->Update(M,
                                                        *molar_dens->ViewComponent(*comp,false), 0.);
        } else {
          // evaluate MassDensity() directly
          const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
          const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
          const Epetra_MultiVector& conc_v = *(conc->ViewComponent(*comp,false));          
          Epetra_MultiVector& dens_v = *(mass_dens->ViewComponent(*comp,false));

          int count = dens_v.MyLength();
          for (int i=0; i!=count; ++i) {
            eos_params[0] = conc_v[0][i];          
            eos_params[1] = temp_v[0][i];
            eos_params[2] = pres_v[0][i];
            
            dens_v[0][i] = eos_sw_->DMassDensityDT(eos_params);
          }
        }
      }
    }

  } else if (wrt_key == conc_key_) {

    if (molar_dens != Teuchos::null) {
      // evaluate MolarDensity()
      for (CompositeVector::name_iterator comp=molar_dens->begin();
           comp!=molar_dens->end(); ++comp) {
        
        const Epetra_MultiVector& conc_v = *(conc->ViewComponent(*comp,false));
        const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
        const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
        Epetra_MultiVector& dens_v = *(molar_dens->ViewComponent(*comp,false));

        int count = dens_v.MyLength();
        for (int i=0; i!=count; ++i) {
          eos_params[0] = conc_v[0][i];          
          eos_params[1] = temp_v[0][i];
          eos_params[2] = pres_v[0][i];
            
          dens_v[0][i] = eos_sw_->DMolarDensityDC(eos_params);
        }
      }
    }

    if (mass_dens != Teuchos::null) {
      for (CompositeVector::name_iterator comp=mass_dens->begin();
           comp!=mass_dens->end(); ++comp) {
        if (mode_ == EOS_MODE_BOTH && eos_sw_->IsConstantMolarMass() &&
            molar_dens->HasComponent(*comp)) {
          // calculate MassDensity from MolarDensity and molar mass.
          double M = eos_sw_->MolarMass();

          mass_dens->ViewComponent(*comp,false)->Update(M,
                  *molar_dens->ViewComponent(*comp,false), 0.);
        } else {
          // evaluate MassDensity() directly
          const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
          const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
          const Epetra_MultiVector& conc_v = *(conc->ViewComponent(*comp,false));          
          Epetra_MultiVector& dens_v = *(mass_dens->ViewComponent(*comp,false));

          int count = dens_v.MyLength();
          for (int i=0; i!=count; ++i) {
            eos_params[0] = conc_v[0][i];          
            eos_params[1] = temp_v[0][i];
            eos_params[2] = pres_v[0][i];
            
            dens_v[0][i] = eos_sw_->DMassDensityDC(eos_params);
          }
        }
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }

} 

} // namespace
} // namespace
