/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EOSFieldEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_factory.hh"
#include "eos_evaluator_tp.hh"

namespace Amanzi {
namespace Relations {

EOSEvaluatorTP::EOSEvaluatorTP(Teuchos::ParameterList& plist) :
  EOSEvaluator(plist) {
  
  Key name = plist_.get<std::string>("evaluator name");
  // Set up my dependencies.
  Key domain_name = Keys::getDomain(name);

  // -- temperature
  temp_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(temp_key_);
  
  // -- pressure
  pres_key_ = Keys::readKey(plist_, domain_name, "pressure", "pressure");
  dependencies_.insert(pres_key_);
  


  // -- logging
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    for (KeySet::const_iterator dep=dependencies_.begin();
         dep!=dependencies_.end(); ++dep) {
      *vo_->os() << " dep: " << *dep << std::endl;
    }
  }

  // // Construct my EOS model
  // AMANZI_ASSERT(plist_.isSublist("EOS parameters"));
  // EOSFactory eos_fac;
  // eos_ = eos_fac.createEOS(plist_.sublist("EOS parameters"));
};


EOSEvaluatorTP::EOSEvaluatorTP(const EOSEvaluatorTP& other) :
    EOSEvaluator(other),
    temp_key_(other.temp_key_),
    pres_key_(other.pres_key_)
 {}


Teuchos::RCP<FieldEvaluator> EOSEvaluatorTP::Clone() const {
  return Teuchos::rcp(new EOSEvaluatorTP(*this));
}


void EOSEvaluatorTP::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  
  int num_dep = dependencies_.size();  
  std::vector<double> eos_params(num_dep);

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);

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


      const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
      const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
      Epetra_MultiVector& dens_v = *(molar_dens->ViewComponent(*comp,false));

      int count = dens_v.MyLength();
      for (int id=0; id!=count; ++id) {
        
        eos_params[0] = temp_v[0][id];
        eos_params[1] = pres_v[0][id];
       
        dens_v[0][id] = eos_->MolarDensity(eos_params);
                  
        if (dens_v[0][id] < 0.){
          Errors::Message msg;
          msg<<"Values of pressure and temperature result in negative density\n"<<
            "Pressure: "<< pres_key_ <<", value : "<<pres_v[0][id]<<"\n"<<
            "Temperature: "<< temp_key_ <<", value : "<<temp_v[0][id]<<"\n"<<
            "Density "<< dens_v[0][id]<<"\n";
          Exceptions::amanzi_throw(msg);
        }
        
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

        const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
        const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
        Epetra_MultiVector& dens_v = *(mass_dens->ViewComponent(*comp,false));

        int count = dens_v.MyLength();
        for (int id=0; id!=count; ++id) {
          
          eos_params[0] = temp_v[0][id];
          eos_params[1] = pres_v[0][id];          
          dens_v[0][id] = eos_->MassDensity(eos_params);
          AMANZI_ASSERT(dens_v[0][id] > 0.);
        }
      }
    }
  }
}

  
void EOSEvaluatorTP::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
                                                   Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  
  int num_dep = dependencies_.size();  
  std::vector<double> eos_params(num_dep);
   
  // Pull dependencies out of state.  
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);  

  
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
        Epetra_MultiVector& dens_v = *(molar_dens->ViewComponent(*comp,false));

        int count = dens_v.MyLength();
        for (int i=0; i!=count; ++i) {
          eos_params[0] = temp_v[0][i];
          eos_params[1] = pres_v[0][i];
            
          dens_v[0][i] = eos_->DMolarDensityDp(eos_params);
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
          const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
          const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
          Epetra_MultiVector& dens_v = *(mass_dens->ViewComponent(*comp,false));

          int count = dens_v.MyLength();
          for (int i=0; i!=count; ++i) {
            eos_params[0] = temp_v[0][i];
            eos_params[1] = pres_v[0][i];            
            dens_v[0][i] = eos_->DMassDensityDp(eos_params);
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
        Epetra_MultiVector& dens_v = *(molar_dens->ViewComponent(*comp,false));

        int count = dens_v.MyLength();
        for (int i=0; i!=count; ++i) {
          eos_params[0] = temp_v[0][i];
          eos_params[1] = pres_v[0][i];          
          dens_v[0][i] = eos_->DMolarDensityDT(eos_params);
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
          const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
          const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
          Epetra_MultiVector& dens_v = *(mass_dens->ViewComponent(*comp,false));

          int count = dens_v.MyLength();
          for (int i=0; i!=count; ++i) {
            eos_params[0] = temp_v[0][i];
            eos_params[1] = pres_v[0][i];            
            dens_v[0][i] = eos_->DMassDensityDT(eos_params);
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
