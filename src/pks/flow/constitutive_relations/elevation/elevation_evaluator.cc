/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The elevation evaluator gets the surface elevation, slope, and updates pres + elev.

*/

#include "elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {

ElevationEvaluator::ElevationEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist), updated_once_(false), dynamic_mesh_(false)
{
  auto domain = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;
  my_keys_.clear(); // clear and re-insert to ensure proper order

  Key elev_key = Keys::readKey(plist_, domain, "elevation", "elevation");
  my_keys_.emplace_back(KeyTag{ elev_key, tag });
  Key slope_key = Keys::readKey(plist_, domain, "slope magnitude", "slope_magnitude");
  my_keys_.emplace_back(KeyTag{ slope_key, tag });
  Key aspect_key = Keys::readKey(plist_, domain, "aspect", "aspect");
  my_keys_.emplace_back(KeyTag{ aspect_key, tag });

  // If the mesh changes dynamically (e.g. due to the presence of a deformation
  // pk, then we must make sure that elevation is recomputed every time the
  // mesh has been deformed. The indicator for the mesh deformation event is the
  // the deformation field.
  dynamic_mesh_ = plist_.get<bool>("dynamic mesh", false);
  if (dynamic_mesh_) {
    deformation_key_ = Keys::readKey(plist_, domain, "deformation indicator", "deformation");
    dependencies_.insert(KeyTag{ deformation_key_, tag });
  }
}


void
ElevationEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  EvaluateElevationAndSlope_(S, results);

  // If boundary faces are requested, grab the slopes on the internal cell
  CompositeVector* slope = results[1];

  if (slope->HasComponent("boundary_face")) {
    const Epetra_Map& vandelay_map = slope->Mesh()->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false);
    const Epetra_Map& face_map = slope->Mesh()->getMap(AmanziMesh::Entity_kind::FACE,false);
    Epetra_MultiVector& slope_bf = *slope->ViewComponent("boundary_face", false);
    const Epetra_MultiVector& slope_c = *slope->ViewComponent("cell", false);

    // calculate boundary face values
    int nbfaces = slope_bf.MyLength();
    for (int bf = 0; bf != nbfaces; ++bf) {
      // given a boundary face, we need the internal cell to choose the right WRM
      AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
      auto cells = slope->Mesh()->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      AMANZI_ASSERT(cells.size() == 1);

      slope_bf[0][bf] = slope_c[0][cells[0]];
    }
  }
}

// This is hopefully never called?
void
ElevationEvaluator::EvaluatePartialDerivative_(const State& S,
                                               const Key& wrt_key,
                                               const Tag& wrt_tag,
                                               const std::vector<CompositeVector*>& results)
{
  AMANZI_ASSERT(0);
}

// Custom EnsureCompatibility forces this to be updated once.
bool
ElevationEvaluator::Update(State& S, const Key& request)
{
  bool changed = EvaluatorSecondaryMonotypeCV::Update(S, request);
  if (!updated_once_) {
    Update_(S);
    updated_once_ = true;
    return true;
  }
  return changed;
}

} // namespace Flow
} // namespace Amanzi
