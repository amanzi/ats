/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "fractional_conductance_evaluator.hh"
#include "subgrid.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

FractionalConductanceEvaluator::FractionalConductanceEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  Key domain = Keys::getDomain(my_key_);

  vpd_key_ = Keys::readKey(plist_, domain, "volumetric ponded depth", "volumetric_ponded_depth");
  dependencies_.insert(vpd_key_); 

  pdd_key_ = Keys::readKey(plist_, domain, "mobile depth", "mobile_depth");
  dependencies_.insert(pdd_key_);
  
  delta_max_key_ = Keys::readKey(plist_, domain, "microtopographic relief", "microtopographic_relief");
  dependencies_.insert(delta_max_key_);

  delta_ex_key_ = plist_.get<std::string>("excluded volume key", Keys::getKey(domain,"excluded_volume"));
  dependencies_.insert(delta_ex_key_);

  depr_depth_key_ = plist_.get<std::string>("depression depth key", Keys::getKey(domain,"depression_depth"));
  dependencies_.insert(depr_depth_key_);
}


Teuchos::RCP<FieldEvaluator>
FractionalConductanceEvaluator::Clone() const {
  return Teuchos::rcp(new FractionalConductanceEvaluator(*this));
}


void FractionalConductanceEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  const Epetra_MultiVector& vpd = *S->GetFieldData(vpd_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& del_max = *S->GetFieldData(delta_max_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& del_ex = *S->GetFieldData(delta_ex_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& depr_depth = *S->GetFieldData(depr_depth_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& mobile_depth = *S->GetFieldData(pdd_key_)->ViewComponent("cell",false);
  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  
  int ncells = res.MyLength();
  for (int c=0; c!=ncells; ++c) {
    if (mobile_depth[0][c] <= 0.0) {
      res[0][c] = 0;
    } else {
      double vol_depr_depth = subgrid_VolumetricDepth(depr_depth[0][c], del_max[0][c], del_ex[0][c]);
      res[0][c] = (vpd[0][c] - vol_depr_depth) / (mobile_depth[0][c]);
    }
  }
}


void
FractionalConductanceEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  const Epetra_MultiVector& vpd = *S->GetFieldData(vpd_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& del_max = *S->GetFieldData(delta_max_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& del_ex = *S->GetFieldData(delta_ex_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& depr_depth = *S->GetFieldData(depr_depth_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& mobile_depth = *S->GetFieldData(pdd_key_)->ViewComponent("cell",false);
  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  
  if (wrt_key == pdd_key_) {
    int ncells = res.MyLength();
    for (int c=0; c!=ncells; ++c) {
      if (mobile_depth[0][c] <= 0.0) {
        res[0][c] = 0;
      } else {
        double vol_depr_depth = subgrid_VolumetricDepth(depr_depth[0][c], del_max[0][c], del_ex[0][c]);
        res[0][c] = - (vpd[0][c] - vol_depr_depth) * std::pow(mobile_depth[0][c],-2.);
      }
    }

  } else if (wrt_key == vpd_key_) {
    int ncells = res.MyLength();
    for (int c=0; c!=ncells; ++c) {
      if (mobile_depth[0][c] <= 0.0) {
        res[0][c] = 0;
      } else {
        res[0][c] =  1.0 / mobile_depth[0][c];
      }
    }
    
  } else {
    Errors::Message msg("VolumetricHeightSubgridEvaluator: Not Implemented: no derivatives implemented other than ponded depth.");
    Exceptions::amanzi_throw(msg);
  }
}

} //namespace
} //namespace
} //namespace
