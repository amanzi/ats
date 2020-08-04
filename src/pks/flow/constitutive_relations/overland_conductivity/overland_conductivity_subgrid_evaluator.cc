/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow subgrid model.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "independent_variable_field_evaluator.hh"
#include "overland_conductivity_subgrid_evaluator.hh"
#include "manning_conductivity_model.hh"

namespace Amanzi {
namespace Flow {

OverlandConductivitySubgridEvaluator::OverlandConductivitySubgridEvaluator(Teuchos::ParameterList& plist)
    : SecondaryVariableFieldEvaluator(plist)
{
  Key domain = Keys::getDomain(my_key_);

  if (plist_.isParameter("height key") ||
      plist_.isParameter("ponded depth key") ||
      plist_.isParameter("height key suffix") ||
      plist_.isParameter("ponded depth key suffix")) {
    Errors::Message message("OverlandConductivitySubgrid: only use \"depth key\" or \"depth key suffix\", not \"height key\" or \"ponded depth key\".");
    Exceptions::amanzi_throw(message);
  }

  depth_key_ = Keys::readKey(plist_, domain, "mobile depth", "mobile_depth");
  dependencies_.insert(depth_key_);

  slope_key_ = Keys::readKey(plist_, domain, "slope", "slope_magnitude");
  dependencies_.insert(slope_key_);

  coef_key_ = Keys::readKey(plist_, domain, "coefficient", "manning_coefficient");
  dependencies_.insert(coef_key_);

  dens_key_ = Keys::readKey(plist_, domain, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(dens_key_);

  frac_cond_key_ = Keys::readKey(plist_, domain, "fractional conductance", "fractional_conductance");
  dependencies_.insert(frac_cond_key_); 
  
  drag_exp_key_ = Keys::readKey(plist_, domain, "drag exponent", "drag_exponent");
  dependencies_.insert(drag_exp_key_); 
  
  // create the model
  Teuchos::ParameterList& sublist = plist_.sublist("overland conductivity model");
  model_ = Teuchos::rcp(new ManningConductivityModel(sublist));
}


Teuchos::RCP<FieldEvaluator>
OverlandConductivitySubgridEvaluator::Clone() const {
  return Teuchos::rcp(new OverlandConductivitySubgridEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void OverlandConductivitySubgridEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> depth = S->GetFieldData(depth_key_);
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> coef = S->GetFieldData(coef_key_);

  Teuchos::RCP<const CompositeVector> frac_cond = S->GetFieldData(frac_cond_key_);
  Teuchos::RCP<const CompositeVector> drag = S->GetFieldData(drag_exp_key_);

  Teuchos::RCP<const CompositeVector> dens = S->GetFieldData(dens_key_);
  
  for (const auto& comp : *result) {
    const Epetra_MultiVector& depth_v = *depth->ViewComponent(comp,false);
    const Epetra_MultiVector& slope_v = *slope->ViewComponent(comp,false);
    const Epetra_MultiVector& coef_v = *coef->ViewComponent(comp,false);

    const Epetra_MultiVector& frac_cond_v = *frac_cond->ViewComponent(comp,false);
    const Epetra_MultiVector& drag_v = *drag->ViewComponent(comp,false);
    const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent(comp,false);    
    Epetra_MultiVector& result_v = *result->ViewComponent(comp,false);
    
    int ncomp = result->size(comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][i], coef_v[0][i]);
      result_v[0][i] *= dens_v[0][i] * std::pow(frac_cond_v[0][i], drag_v[0][i] + 1);
    }
  }
}
  

void OverlandConductivitySubgridEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> depth = S->GetFieldData(depth_key_);
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> coef = S->GetFieldData(coef_key_);

  Teuchos::RCP<const CompositeVector> frac_cond = S->GetFieldData(frac_cond_key_);
  Teuchos::RCP<const CompositeVector> drag = S->GetFieldData(drag_exp_key_);

  Teuchos::RCP<const CompositeVector> dens = S->GetFieldData(dens_key_);
  
  if (wrt_key == depth_key_) {
    for (const auto& comp : *result) {
      const Epetra_MultiVector& depth_v = *depth->ViewComponent(comp,false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(comp,false);

      const Epetra_MultiVector& frac_cond_v = *frac_cond->ViewComponent(comp,false);
      const Epetra_MultiVector& drag_v = *drag->ViewComponent(comp,false);
      const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent(comp,false);    
      Epetra_MultiVector& result_v = *result->ViewComponent(comp,false);
    
      int ncomp = result->size(comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DConductivityDDepth(depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        result_v[0][i] *= dens_v[0][i] * std::pow(frac_cond_v[0][i], drag_v[0][i] + 1);
      }
    }

  } else if (wrt_key == dens_key_) {
    for (const auto& comp : *result) {
      const Epetra_MultiVector& depth_v = *depth->ViewComponent(comp,false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(comp,false);

      const Epetra_MultiVector& frac_cond_v = *frac_cond->ViewComponent(comp,false);
      const Epetra_MultiVector& drag_v = *drag->ViewComponent(comp,false);
      const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent(comp,false);    
      Epetra_MultiVector& result_v = *result->ViewComponent(comp,false);
    
      int ncomp = result->size(comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        result_v[0][i] *= std::pow(frac_cond_v[0][i], drag_v[0][i] + 1);
      }
    }

  } else if (wrt_key == frac_cond_key_) {
    for (const auto& comp : *result) {
      const Epetra_MultiVector& depth_v = *depth->ViewComponent(comp,false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(comp,false);

      const Epetra_MultiVector& frac_cond_v = *frac_cond->ViewComponent(comp,false);
      const Epetra_MultiVector& drag_v = *drag->ViewComponent(comp,false);
      const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent(comp,false);    
      Epetra_MultiVector& result_v = *result->ViewComponent(comp,false);
    
      int ncomp = result->size(comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        result_v[0][i] *= dens_v[0][i] * (drag_v[0][i] + 1) *
                          std::pow(frac_cond_v[0][i], drag_v[0][i]);
      }
    }
  } else {
    result->PutScalar(0.);
  }
}  

} // namespace Flow
} // namespace Amanzi

