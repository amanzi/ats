/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  The manning coefficient with variable litter evaluator is an algebraic evaluator of a given model.
  Manning's coefficient that varies based on litter thickness and ponded depth.
  Generated via evaluator_generator.
*/

#include "boost/algorithm/string/predicate.hpp"

#include "manning_coefficient_litter_evaluator.hh"
#include "manning_coefficient_litter_model.hh"
#include "manning_coefficient_litter_constant_model.hh"
#include "manning_coefficient_litter_variable_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
ManningCoefficientLitterEvaluator::ManningCoefficientLitterEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("manning coefficient parameters");
  models_ = createManningCoefPartition(sublist);
  InitializeFromPlist_();
}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
ManningCoefficientLitterEvaluator::Clone() const
{
  return Teuchos::rcp(new ManningCoefficientLitterEvaluator(*this));
}


// Initialize by setting up dependencies
void
ManningCoefficientLitterEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // - pull Keys from plist
  // dependency: litter_thickness
  Key litter_domain = Keys::readDomainHint(plist_, domain, "surface", "litter");

  ld_key_ = Keys::readKey(plist_, litter_domain, "litter thickness", "thickness");
  dependencies_.insert(KeyTag{ ld_key_, tag });

  // dependency: ponded_depth
  pd_key_ = Keys::readKey(plist_, domain, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ pd_key_, tag });
}


void
ManningCoefficientLitterEvaluator::Evaluate_(const State& S,
                                             const std::vector<CompositeVector*>& result)
{
  // Initialize the MeshPartition
  if (!models_->first->initialized()) {
    models_->first->Initialize(result[0]->Mesh(), -1);
    models_->first->Verify();
  }

  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> ld = S.GetPtr<CompositeVector>(ld_key_, tag);
  Teuchos::RCP<const CompositeVector> pd = S.GetPtr<CompositeVector>(pd_key_, tag);

  // cell values
  {
    const Epetra_MultiVector& ld_v = *ld->ViewComponent("cell", false);
    const Epetra_MultiVector& pd_v = *pd->ViewComponent("cell", false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent("cell", false);

    int ncomp = result[0]->size("cell", false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] =
        models_->second[(*models_->first)[i]]->ManningCoefficient(ld_v[0][i], pd_v[0][i]);
    }
  }

  // potential boundary face values
  if (result[0]->HasComponent("boundary_face")) {
    const Epetra_MultiVector& ld_v = *ld->ViewComponent("boundary_face", false);
    const Epetra_MultiVector& pd_v = *pd->ViewComponent("boundary_face", false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent("boundary_face", false);

    // Need to get boundary face's inner cell to specify the WRM.
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = result[0]->Mesh();
    const Epetra_Map& vandelay_map = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false);
    const Epetra_Map& face_map = mesh->getMap(AmanziMesh::Entity_kind::FACE,false);

    int ncomp = result[0]->size("boundary_face", false);
    for (int bf = 0; bf != ncomp; ++bf) {
      // given a boundary face, we need the internal cell to choose the right model
      AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
      auto cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      AMANZI_ASSERT(cells.size() == 1);

      int index = (*models_->first)[cells[0]];
      result_v[0][bf] = models_->second[index]->ManningCoefficient(ld_v[0][bf], pd_v[0][bf]);
    }
  }
}


void
ManningCoefficientLitterEvaluator::EvaluatePartialDerivative_(
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
  Teuchos::RCP<const CompositeVector> ld = S.GetPtr<CompositeVector>(ld_key_, tag);
  Teuchos::RCP<const CompositeVector> pd = S.GetPtr<CompositeVector>(pd_key_, tag);
  {
    // cell values
    const Epetra_MultiVector& ld_v = *ld->ViewComponent("cell", false);
    const Epetra_MultiVector& pd_v = *pd->ViewComponent("cell", false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent("cell", false);

    int ncomp = result[0]->size("cell", false);
    if (wrt_key == ld_key_) {
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = models_->second[(*models_->first)[i]]->DManningCoefficientDLitterThickness(
          ld_v[0][i], pd_v[0][i]);
      }

    } else if (wrt_key == pd_key_) {
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = models_->second[(*models_->first)[i]]->DManningCoefficientDPondedDepth(
          ld_v[0][i], pd_v[0][i]);
      }
    } else {
      AMANZI_ASSERT(0);
    }
  }

  // potential boundary face values
  if (result[0]->HasComponent("boundary_face")) {
    const Epetra_MultiVector& ld_v = *ld->ViewComponent("boundary_face", false);
    const Epetra_MultiVector& pd_v = *pd->ViewComponent("boundary_face", false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent("boundary_face", false);

    // Need to get boundary face's inner cell to specify the WRM.
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = result[0]->Mesh();
    const Epetra_Map& vandelay_map = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false);
    const Epetra_Map& face_map = mesh->getMap(AmanziMesh::Entity_kind::FACE,false);

    int ncomp = result[0]->size("boundary_face", false);
    if (wrt_key == ld_key_) {
      for (int bf = 0; bf != ncomp; ++bf) {
        // given a boundary face, we need the internal cell to choose the right model
        AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
        auto cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
        AMANZI_ASSERT(cells.size() == 1);

        int index = (*models_->first)[cells[0]];
        result_v[0][bf] =
          models_->second[index]->DManningCoefficientDLitterThickness(ld_v[0][bf], pd_v[0][bf]);
      }

    } else if (wrt_key == pd_key_) {
      for (int bf = 0; bf != ncomp; ++bf) {
        // given a boundary face, we need the internal cell to choose the right model
        AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
        auto cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
        AMANZI_ASSERT(cells.size() == 1);

        int index = (*models_->first)[cells[0]];
        result_v[0][bf] =
          models_->second[index]->DManningCoefficientDPondedDepth(ld_v[0][bf], pd_v[0][bf]);
      }

    } else {
      AMANZI_ASSERT(0);
    }
  }
}


// Non-member factory
Teuchos::RCP<ManningCoefPartition>
createManningCoefPartition(Teuchos::ParameterList& plist)
{
  std::vector<Teuchos::RCP<ManningCoefficientLitterModel>> models;
  std::vector<std::string> region_list;

  for (Teuchos::ParameterList::ConstIterator lcv = plist.begin(); lcv != plist.end(); ++lcv) {
    std::string name = lcv->first;
    if (plist.isSublist(name)) {
      Teuchos::ParameterList sublist = plist.sublist(name);
      region_list.push_back(sublist.get<std::string>("region"));

      std::string coef_type = sublist.get<std::string>("manning coefficient model type");
      if (coef_type == "constant") {
        models.push_back(Teuchos::rcp(new ManningCoefficientLitterConstantModel(sublist)));
      } else if (coef_type == "variable") {
        models.push_back(Teuchos::rcp(new ManningCoefficientLitterVariableModel(sublist)));
      } else {
        Errors::Message message("ManningCoefficient: unknown model type");
        Exceptions::amanzi_throw(message);
      }

    } else {
      Errors::Message message("ManningCoefficient: incorrectly formed input parameter list");
      Exceptions::amanzi_throw(message);
    }
  }

  Teuchos::RCP<Functions::MeshPartition> part =
    Teuchos::rcp(new Functions::MeshPartition(AmanziMesh::Entity_kind::CELL, region_list));

  return Teuchos::rcp(new ManningCoefPartition(part, models));
}


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
