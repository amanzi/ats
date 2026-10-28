/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Mesh.hh"
#include "CompositeVector.hh"
#include "State.hh"

#include "evaluator_upwind_arithmetic_mean.hh"

namespace Amanzi {
namespace Relations {

const std::string EvaluatorUpwindArithmeticMean::eval_type = "upwind: arithmetic mean";

EvaluatorUpwindArithmeticMean::EvaluatorUpwindArithmeticMean(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorUpwinding(plist)
{}


Teuchos::RCP<Evaluator>
EvaluatorUpwindArithmeticMean::Clone() const
{
  return Teuchos::rcp(new EvaluatorUpwindArithmeticMean(*this));
}


void
EvaluatorUpwindArithmeticMean::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  CompositeVector& face_coef = *result[0];

  const auto& field = S.Get<CompositeVector>(field_key_, tag);

  face_coef.getComponent("face", true)->putScalar(0.0);
  if (face_coef.hasComponent("cell")) { face_coef.getComponent("cell", true)->putScalar(1.0); }

  field.scatterMasterToGhosted("cell");

  // A boundary face has no second neighbor cell; its "other side" of the
  // mean is field_key_'s own "boundary_face" component, seeded by
  // SeedBoundaryFaces_. When field_key_ has no "boundary_face" component at
  // all, the boundary face's value is just the single interior cell's
  // value.
  bool have_boundary_candidate = field.hasComponent("boundary_face");
  SeedBoundaryFaces_(S, face_coef, /*is_derivative=*/false);

  const AmanziMesh::MeshCache& m = face_coef.getMesh()->getCache();
  auto coef_faces = face_coef.viewComponent("face", false);
  const auto coef_cells = field.viewComponent("cell", true);

  Kokkos::parallel_for(
    "EvaluatorUpwindArithmeticMean::Evaluate_", coef_faces.extent(0), KOKKOS_LAMBDA(const int& f) {
      auto cells = m.getFaceCells(f);
      if (cells.size() == 1) {
        coef_faces(f, 0) = have_boundary_candidate ?
          0.5 * (coef_cells(cells(0), 0) + coef_faces(f, 0)) :
          coef_cells(cells(0), 0);
      } else {
        double sum = 0.0;
        for (const auto& c : cells) sum += coef_cells(c, 0);
        coef_faces(f, 0) = sum / cells.size();
      }
    });
  face_coef.scatterMasterToGhosted("face");

  ApplyClobber_(S, face_coef, /*is_derivative=*/false);
}


void
EvaluatorUpwindArithmeticMean::EvaluatePartialDerivative_(const State& S,
                                                          const Key& wrt_key,
                                                          const Tag& wrt_tag,
                                                          const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(0);
}


void
EvaluatorUpwindArithmeticMean::UpwindDerivative_(const State& S,
                                                 const CompositeVector& dfield_dp,
                                                 CompositeVector& result) const
{
  dfield_dp.scatterMasterToGhosted("cell");

  // As in Evaluate_: a boundary face's derivative has a second term only
  // when field_key_'s own derivative has a "boundary_face" component.
  bool have_boundary_candidate = dfield_dp.hasComponent("boundary_face");

  const AmanziMesh::MeshCache& m = result.getMesh()->getCache();
  auto dcoef_faces = result.viewComponent("face", false);
  const auto dfield_cells = dfield_dp.viewComponent("cell", true);

  Kokkos::parallel_for(
    "EvaluatorUpwindArithmeticMean::UpwindDerivative_",
    dcoef_faces.extent(0),
    KOKKOS_LAMBDA(const int& f) {
      auto cells = m.getFaceCells(f);
      if (cells.size() == 1) {
        // boundary face: slot 0 is d(value)/d(p_cell), slot 1 is
        // d(value)/d(p_boundary) -- the boundary candidate's own
        // derivative, pre-loaded into dcoef_faces(f,0) by
        // SeedBoundaryFaces_ when have_boundary_candidate is true.
        double dcoef_boundary = have_boundary_candidate ? dcoef_faces(f, 0) : 0.0;
        dcoef_faces(f, 0) = 0.5 * dfield_cells(cells(0), 0);
        dcoef_faces(f, 1) = have_boundary_candidate ? 0.5 * dcoef_boundary : 0.0;
      } else {
        dcoef_faces(f, 0) = 0.5 * dfield_cells(cells(0), 0);
        dcoef_faces(f, 1) = 0.5 * dfield_cells(cells(1), 0);
      }
    });
  result.scatterMasterToGhosted("face");
}

} // namespace Relations
} // namespace Amanzi
