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
ManningCoefficientLitterEvaluator::ManningCoefficientLitterEvaluator(Teuchos::ParameterList& plist) :
  EvaluatorSecondaryMonotypeCV(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("manning coefficient parameters");
  models_ = createManningCoefPartition(sublist);
  InitializeFromPlist_();
}


// Copy constructor
ManningCoefficientLitterEvaluator::ManningCoefficientLitterEvaluator(const ManningCoefficientLitterEvaluator& other) :
  EvaluatorSecondaryMonotypeCV(other),
  ld_key_(other.ld_key_),
  pd_key_(other.pd_key_),
  models_(other.models_) {}


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
  Key domain_name = Keys::getDomain(my_key_);

  // - pull Keys from plist
  // dependency: litter_thickness
  Key litter_domain_name;
  if (domain_name == "surface") {
    litter_domain_name = "litter";
    litter_domain_name = plist_.get<std::string>("litter domain name", litter_domain_name);
  } else if (boost::starts_with(domain_name, "surface")) {
    litter_domain_name = Key("litter")+domain_name.substr(7,domain_name.size());
    litter_domain_name = plist_.get<std::string>("litter domain name", litter_domain_name);
  }


  ld_key_ = plist_.get<std::string>("litter thickness key",
          Keys::getKey(litter_domain_name, "thickness"));
  dependencies_.insert(ld_key_);

  // dependency: ponded_depth
  pd_key_ = plist_.get<std::string>("ponded depth key",
          Keys::getKey(domain_name,"ponded_depth"));
  dependencies_.insert(pd_key_);
}


void
ManningCoefficientLitterEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  // Initialize the MeshPartition
  if (!models_->first->initialized()) {
    models_->first->Initialize(result->Mesh(), -1);
    models_->first->Verify();
  }

  Teuchos::RCP<const CompositeVector> ld = S->GetPtr<CompositeVector>(ld_key_);
  Teuchos::RCP<const CompositeVector> pd = S->GetPtr<CompositeVector>(pd_key_);

  // cell values
  {
    const Epetra_MultiVector& ld_v = *ld->ViewComponent("cell", false);
    const Epetra_MultiVector& pd_v = *pd->ViewComponent("cell", false);
    Epetra_MultiVector& result_v = *result->ViewComponent("cell", false);

    int ncomp = result->size("cell", false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = models_->second[(*models_->first)[i]]->ManningCoefficient(ld_v[0][i], pd_v[0][i]);
    }
  }

  // potential boundary face values
  if (result->HasComponent("boundary_face")) {
    const Epetra_MultiVector& ld_v = *ld->ViewComponent("boundary_face", false);
    const Epetra_MultiVector& pd_v = *pd->ViewComponent("boundary_face", false);
    Epetra_MultiVector& result_v = *result->ViewComponent("boundary_face", false);

    // Need to get boundary face's inner cell to specify the WRM.
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();
    const Epetra_Map& vandelay_map = mesh->exterior_face_map(false);
    const Epetra_Map& face_map = mesh->face_map(false);
    AmanziMesh::Entity_ID_List cells;

    int ncomp = result->size("boundary_face", false);
    for (int bf=0; bf!=ncomp; ++bf) {
      // given a boundary face, we need the internal cell to choose the right model
      AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
      mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      AMANZI_ASSERT(cells.size() == 1);

      int index = (*models_->first)[cells[0]];
      result_v[0][bf] = models_->second[index]->ManningCoefficient(ld_v[0][bf], pd_v[0][bf]);
    }
  }
}


void
ManningCoefficientLitterEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  // Initialize the MeshPartition
  if (!models_->first->initialized()) {
    models_->first->Initialize(result->Mesh(), -1);
    models_->first->Verify();
  }

  Teuchos::RCP<const CompositeVector> ld = S->GetPtr<CompositeVector>(ld_key_);
  Teuchos::RCP<const CompositeVector> pd = S->GetPtr<CompositeVector>(pd_key_);

  {
    // cell values
    const Epetra_MultiVector& ld_v = *ld->ViewComponent("cell", false);
    const Epetra_MultiVector& pd_v = *pd->ViewComponent("cell", false);
    Epetra_MultiVector& result_v = *result->ViewComponent("cell", false);

    int ncomp = result->size("cell", false);
    if (wrt_key == ld_key_) {
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = models_->second[(*models_->first)[i]]
          ->DManningCoefficientDLitterThickness(ld_v[0][i], pd_v[0][i]);
      }

    } else if (wrt_key == pd_key_) {
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = models_->second[(*models_->first)[i]]
          ->DManningCoefficientDPondedDepth(ld_v[0][i], pd_v[0][i]);
      }
    } else {
      AMANZI_ASSERT(0);
    }
  }

  // potential boundary face values
  if (result->HasComponent("boundary_face")) {
    const Epetra_MultiVector& ld_v = *ld->ViewComponent("boundary_face", false);
    const Epetra_MultiVector& pd_v = *pd->ViewComponent("boundary_face", false);
    Epetra_MultiVector& result_v = *result->ViewComponent("boundary_face", false);

    // Need to get boundary face's inner cell to specify the WRM.
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();
    const Epetra_Map& vandelay_map = mesh->exterior_face_map(false);
    const Epetra_Map& face_map = mesh->face_map(false);
    AmanziMesh::Entity_ID_List cells;

    int ncomp = result->size("boundary_face", false);
    if (wrt_key == ld_key_) {
      for (int bf=0; bf!=ncomp; ++bf) {
        // given a boundary face, we need the internal cell to choose the right model
        AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
        mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        AMANZI_ASSERT(cells.size() == 1);

        int index = (*models_->first)[cells[0]];
        result_v[0][bf] = models_->second[index]->DManningCoefficientDLitterThickness(ld_v[0][bf], pd_v[0][bf]);
      }

    } else if (wrt_key == pd_key_) {
      for (int bf=0; bf!=ncomp; ++bf) {
        // given a boundary face, we need the internal cell to choose the right model
        AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
        mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        AMANZI_ASSERT(cells.size() == 1);

        int index = (*models_->first)[cells[0]];
        result_v[0][bf] = models_->second[index]->DManningCoefficientDPondedDepth(ld_v[0][bf], pd_v[0][bf]);
      }

    } else {
      AMANZI_ASSERT(0);
    }
  }
}



// Non-member factory
Teuchos::RCP<ManningCoefPartition>
createManningCoefPartition(Teuchos::ParameterList& plist) {

  std::vector<Teuchos::RCP<ManningCoefficientLitterModel> > models;
  std::vector<std::string> region_list;

  for (Teuchos::ParameterList::ConstIterator lcv=plist.begin();
       lcv!=plist.end(); ++lcv) {
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
    Teuchos::rcp(new Functions::MeshPartition(AmanziMesh::CELL,region_list));

  return Teuchos::rcp(new ManningCoefPartition(part, models));
}


} //namespace
} //namespace
} //namespace
