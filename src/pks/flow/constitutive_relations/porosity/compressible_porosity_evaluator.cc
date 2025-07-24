/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates the porosity, given a small compressibility of rock.

*/

#include "compressible_porosity_evaluator.hh"
#include "compressible_porosity_model.hh"

namespace Amanzi {
namespace Flow {

CompressiblePorosityEvaluator::CompressiblePorosityEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  std::string domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;
  pres_key_ = Keys::readKey(plist_, domain_name, "pressure", "pressure");
  dependencies_.insert(KeyTag{ pres_key_, tag });

  poro_key_ = Keys::readKey(plist_, domain_name, "base porosity", "base_porosity");
  dependencies_.insert(KeyTag{ poro_key_, tag });

  AMANZI_ASSERT(plist_.isSublist("compressible porosity model parameters"));
  models_ = createCompressiblePorosityModelPartition(
    plist_.sublist("compressible porosity model parameters"));
}

Teuchos::RCP<Evaluator>
CompressiblePorosityEvaluator::Clone() const
{
  return Teuchos::rcp(new CompressiblePorosityEvaluator(*this));
}


// Required methods from EvaluatorSecondaryMonotypeCV
void
CompressiblePorosityEvaluator::Evaluate_(const State& S,
                                         const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  // Initialize the MeshPartition
  if (!models_->first->initialized()) {
    models_->first->Initialize(result[0]->Mesh(), -1);
    models_->first->Verify();
  }

  Teuchos::RCP<const CompositeVector> pres = S.GetPtr<CompositeVector>(pres_key_, tag);
  Teuchos::RCP<const CompositeVector> poro = S.GetPtr<CompositeVector>(poro_key_, tag);
  const double& patm = S.Get<double>("atmospheric_pressure", Tags::DEFAULT);

  // evaluate the model
  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    AMANZI_ASSERT(
      *comp ==
      "cell"); // partition on cell only, could add boundary_face if needed (but not currently)
    const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp, false));
    const Epetra_MultiVector& poro_v = *(poro->ViewComponent(*comp, false));
    Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

    int count = result[0]->size(*comp);
    for (int id = 0; id != count; ++id) {
      result_v[0][id] =
        models_->second[(*models_->first)[id]]->Porosity(poro_v[0][id], pres_v[0][id], patm);
    }
  }
}


void
CompressiblePorosityEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  // Initialize the MeshPartition
  if (!models_->first->initialized()) {
    models_->first->Initialize(result[0]->Mesh(), -1);
    models_->first->Verify();
  }

  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> pres = S.GetPtr<CompositeVector>(pres_key_, tag);
  Teuchos::RCP<const CompositeVector> poro = S.GetPtr<CompositeVector>(poro_key_, tag);
  const double& patm = S.Get<double>("atmospheric_pressure", Tags::DEFAULT);

  if (wrt_key == pres_key_) {
    // evaluate the model
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      AMANZI_ASSERT(
        *comp ==
        "cell"); // partition on cell only, could add boundary_face if needed (but not currently)
      const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp, false));
      const Epetra_MultiVector& poro_v = *(poro->ViewComponent(*comp, false));
      Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

      int count = result[0]->size(*comp);
      for (int id = 0; id != count; ++id) {
        result_v[0][id] = models_->second[(*models_->first)[id]]->DPorosityDPressure(
          poro_v[0][id], pres_v[0][id], patm);
      }
    }

  } else if (wrt_key == poro_key_) {
    // evaluate the model
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      AMANZI_ASSERT(
        *comp ==
        "cell"); // partition on cell only, could add boundary_face if needed (but not currently)
      const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp, false));
      const Epetra_MultiVector& poro_v = *(poro->ViewComponent(*comp, false));
      Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

      int count = result[0]->size(*comp);
      for (int id = 0; id != count; ++id) {
        result_v[0][id] = models_->second[(*models_->first)[id]]->DPorosityDBasePorosity(
          poro_v[0][id], pres_v[0][id], patm);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} // namespace Flow
} // namespace Amanzi
