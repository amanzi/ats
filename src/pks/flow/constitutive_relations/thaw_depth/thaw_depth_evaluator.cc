/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! An evaluator for calculating the depth to frozen soil/permafrost, relative to the surface.

#include "thaw_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {

ThawDepthEvaluator::ThawDepthEvaluator(Teuchos::ParameterList& plist)
    : SecondaryVariableFieldEvaluator(plist)
{
  domain_ = Keys::getDomain(my_key_);
  domain_ss_ = Keys::readDomainHint(plist, domain_, "surface", "domain");
  temp_key_ = Keys::readKey(plist, domain_ss_, "temperature", "temperature");
  dependencies_.insert(temp_key_);

  trans_width_ =  plist.get<double>("transition width [K]", 0.2);
}


void
ThawDepthEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  const auto& temp_c = *S->GetFieldData(temp_key_)->ViewComponent("cell", false);

  // search through the column and find the first frozen cell
  const auto& surf_mesh = S->GetMesh(domain_);
  const auto& subsurf_mesh = S->GetMesh(domain_ss_);
  int z_dim = subsurf_mesh->space_dimension() - 1;

  double trans_temp = 273.15 + 0.5*trans_width_;

  for (AmanziMesh::Entity_ID sc=0; sc!=res_c.MyLength(); ++sc) {
    AmanziMesh::Entity_ID top_f = surf_mesh->entity_get_parent(AmanziMesh::CELL, sc);
    double top_z = subsurf_mesh->face_centroid(top_f)[z_dim];

    double thaw_z = std::numeric_limits<double>::quiet_NaN();
    for (int i=0; i!=subsurf_mesh->cells_of_column(sc).size(); ++i) {
      AmanziMesh::Entity_ID c = subsurf_mesh->cells_of_column(sc)[i];
      if (temp_c[0][c] < trans_temp) {
        // find the face above
        AmanziMesh::Entity_ID f_up = subsurf_mesh->faces_of_column(sc)[i];
        thaw_z = subsurf_mesh->face_centroid(f_up)[z_dim];
        break;
      }
    }
    res_c[0][sc] = top_z - thaw_z;
  }
}


// Custom EnsureCompatibility forces this to be updated once.
bool
ThawDepthEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
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
ThawDepthEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  AMANZI_ASSERT(my_key_ != std::string(""));

  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, true);
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
