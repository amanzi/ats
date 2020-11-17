/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "ColumnSumEvaluator.hh"

namespace Amanzi {
namespace Relations {



ColumnSumEvaluator::ColumnSumEvaluator(Teuchos::ParameterList& plist)
    : SecondaryVariableFieldEvaluator(plist)
{
  surf_domain_ = Keys::getDomain(my_key_);
  if (surf_domain_ == "surface") {
    domain_ = "";
    domain_ = plist.get<std::string>("column domain name", domain_);
  } else if (Keys::starts_with(surf_domain_, "surface_")) {
    domain_ = surf_domain_.substr(8, surf_domain_.size());
    domain_ = plist.get<std::string>("column domain name", domain_);
  } else {
    domain_ = plist.get<std::string>("column domain name");
  }

  dep_key_ = Keys::readKey(plist_, domain_, "evaluator dependency", Keys::getKey(domain_, Keys::getVarName(my_key_)));
  dependencies_.insert(dep_key_);

  Key pname = dep_key_ + " coefficient";
  coef_ = plist.get<double>(pname, 1.0);

  // dependency: cell volume, surface cell volume
  if (plist_.get<bool>("include volume factor", true)) {
    cv_key_ = Keys::readKey(plist_, domain_, "cell volume", "cell_volume");
    dependencies_.insert(cv_key_);

    surf_cv_key_ = Keys::readKey(plist_, surf_domain_, "surface cell volume", "cell_volume");
    dependencies_.insert(surf_cv_key_);
  }

  if (plist_.get<bool>("divide by density", true)) {
    molar_dens_key_ = Keys::readKey(plist_, domain_, "molar density", "molar_density_liquid");
    dependencies_.insert(molar_dens_key_);
  }
}


ColumnSumEvaluator::ColumnSumEvaluator(const ColumnSumEvaluator& other) :
  SecondaryVariableFieldEvaluator(other),
  coef_(other.coef_),
  dep_key_(other.dep_key_),
  cv_key_(other.cv_key_),
  surf_cv_key_(other.surf_cv_key_),
  molar_dens_key_(other.molar_dens_key_) {}


Teuchos::RCP<FieldEvaluator>
ColumnSumEvaluator::Clone() const
{
  return Teuchos::rcp(new ColumnSumEvaluator(*this));
}


void
ColumnSumEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& dep_c = *S->GetFieldData(dep_key_)->ViewComponent("cell", false);

  if (cv_key_ != "") {
    const Epetra_MultiVector& cv = *S->GetFieldData(cv_key_)->ViewComponent("cell", false);
    const Epetra_MultiVector& surf_cv = *S->GetFieldData(surf_cv_key_)->ViewComponent("cell", false);

    if (molar_dens_key_ != "") {
      const Epetra_MultiVector& dens = *S->GetFieldData(molar_dens_key_)->ViewComponent("cell",false);

      Teuchos::RCP<const AmanziMesh::Mesh> subsurf_mesh = S->GetMesh(domain_);
      for (int c=0; c!=res_c.MyLength(); ++c) {
        double sum = 0;
        for (auto i : subsurf_mesh->cells_of_column(c)) {
          sum += dep_c[0][i] * cv[0][i] / dens[0][i];
        }
        res_c[0][c] = coef_ * sum / surf_cv[0][c];
      }
    } else {

      Teuchos::RCP<const AmanziMesh::Mesh> subsurf_mesh = S->GetMesh(domain_);
      for (int c=0; c!=res_c.MyLength(); ++c) {
        double sum = 0;
        for (auto i : subsurf_mesh->cells_of_column(c)) {
          sum += dep_c[0][i] * cv[0][i];
        }
        res_c[0][c] = coef_ * sum / surf_cv[0][c];
      }
    }

  } else {
    if (molar_dens_key_ != "") {
      const Epetra_MultiVector& dens = *S->GetFieldData(molar_dens_key_)->ViewComponent("cell",false);

      Teuchos::RCP<const AmanziMesh::Mesh> subsurf_mesh = S->GetMesh(domain_);
      for (int c=0; c!=res_c.MyLength(); ++c) {
        double sum = 0;
        for (auto i : subsurf_mesh->cells_of_column(c)) {
          sum += dep_c[0][i] / dens[0][i];
        }
        res_c[0][c] = coef_ * sum;
      }
    } else {

      Teuchos::RCP<const AmanziMesh::Mesh> subsurf_mesh = S->GetMesh(domain_);
      for (int c=0; c!=res_c.MyLength(); ++c) {
        double sum = 0;
        for (auto i : subsurf_mesh->cells_of_column(c)) {
          sum += dep_c[0][i];
        }
        res_c[0][c] = coef_ * sum;
      }
    }
  }
}


void
ColumnSumEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  AMANZI_ASSERT(false);
}


// Custom EnsureCompatibility forces this to be updated once.
bool
ColumnSumEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
        Key request)
{
  bool changed = SecondaryVariableFieldEvaluator::HasFieldChanged(S,request);

  if (!updated_once_) {
    UpdateField_(S);
    updated_once_ = true;
    return true;
  }
  return changed;
}

void
ColumnSumEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  if (my_fac->Mesh() != Teuchos::null) {
    // Recurse into the tree to propagate info to leaves.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      S->RequireFieldEvaluator(*key)->EnsureCompatibility(S);
    }
  }
}


} //namespace
} //namespace
