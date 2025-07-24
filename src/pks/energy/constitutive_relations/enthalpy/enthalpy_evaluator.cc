/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/* -----------------------------------------------------------------------------
ATS

Evaluator for enthalpy.
----------------------------------------------------------------------------- */


#include "enthalpy_evaluator.hh"

namespace Amanzi {
namespace Energy {

EnthalpyEvaluator::EnthalpyEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;
  include_work_ = plist_.get<bool>("include work term", false);

  // -- pressure
  if (include_work_) {
    pres_key_ = Keys::readKey(plist_, domain_name, "pressure", "pressure");
    dependencies_.insert(KeyTag{ pres_key_, tag });

    dens_key_ = Keys::readKey(plist_, domain_name, "molar density liquid", "molar_density_liquid");
    dependencies_.insert(KeyTag{ dens_key_, tag });
  }

  ie_key_ = Keys::readKey(plist_, domain_name, "internal energy liquid", "internal_energy_liquid");
  dependencies_.insert(KeyTag{ ie_key_, tag });
};


Teuchos::RCP<Evaluator>
EnthalpyEvaluator::Clone() const
{
  return Teuchos::rcp(new EnthalpyEvaluator(*this));
};


void
EnthalpyEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Teuchos::OSTab tab = vo_.getOSTab();
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> u_l = S.GetPtr<CompositeVector>(ie_key_, tag);
  *result[0] = *u_l;

  if (include_work_) {
    Teuchos::RCP<const CompositeVector> pres = S.GetPtr<CompositeVector>(pres_key_, tag);
    Teuchos::RCP<const CompositeVector> n_l = S.GetPtr<CompositeVector>(dens_key_, tag);

    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        // 1.e-6 converts to MJoules
        result_v[0][i] += 1.e-6 * pres_v[0][i] / nl_v[0][i];
      }
    }
  }
};


void
EnthalpyEvaluator::EvaluatePartialDerivative_(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag,
                                              const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  // not implemented
  if (wrt_key == ie_key_) {
    result[0]->PutScalar(1.);
  } else if (wrt_key == pres_key_) {
    AMANZI_ASSERT(include_work_);

    Teuchos::RCP<const CompositeVector> n_l = S.GetPtr<CompositeVector>(dens_key_, tag);

    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        // 1.e-6 converts to MJoules
        result_v[0][i] = 1.e-6 / nl_v[0][i];
      }
    }

  } else if (wrt_key == dens_key_) {
    AMANZI_ASSERT(include_work_);

    Teuchos::RCP<const CompositeVector> pres = S.GetPtr<CompositeVector>(pres_key_, tag);
    Teuchos::RCP<const CompositeVector> n_l = S.GetPtr<CompositeVector>(dens_key_, tag);

    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp, false);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        // 1.e-6 converts to MJoules
        result_v[0][i] = -1.e-6 * pres_v[0][i] / std::pow(nl_v[0][i], 2);
      }
    }
  }
};

} // namespace Energy
} // namespace Amanzi
