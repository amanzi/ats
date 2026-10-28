/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Mesh.hh"
#include "MeshHelpers_decl.hh"
#include "CompositeVector.hh"
#include "State.hh"

#include "evaluator_upwind_potential_difference.hh"

namespace Amanzi {
namespace Relations {

const std::string EvaluatorUpwindPotentialDifference::eval_type = "upwind: potential difference";

EvaluatorUpwindPotentialDifference::EvaluatorUpwindPotentialDifference(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorUpwinding(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  potential_key_ = Keys::readKey(*plist, domain, "potential", "pressure");
  dependencies_.insert(KeyTag{ potential_key_, tag });
}


Teuchos::RCP<Evaluator>
EvaluatorUpwindPotentialDifference::Clone() const
{
  return Teuchos::rcp(new EvaluatorUpwindPotentialDifference(*this));
}


bool
EvaluatorUpwindPotentialDifference::IsDifferentiableWRT(const State& S,
                                                        const Key& wrt_key,
                                                        const Tag& wrt_tag) const
{
  if (wrt_key == potential_key_) return false;
  return EvaluatorUpwinding::IsDifferentiableWRT(S, wrt_key, wrt_tag);
}


void
EvaluatorUpwindPotentialDifference::EnsureCompatibility_ToDeps_(State& S)
{
  EvaluatorUpwinding::EnsureCompatibility_ToDeps_(S);

  Tag tag = my_keys_.front().second;
  CompositeVectorSpace pot_fac;
  pot_fac.SetMesh(S.GetMesh(Keys::getDomain(potential_key_)))
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);
  S.Require<CompositeVector, CompositeVectorSpace>(potential_key_, tag).Update(pot_fac);
}


void
EvaluatorUpwindPotentialDifference::Evaluate_(const State& S,
                                              const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  CompositeVector& face_coef = *result[0];

  const auto& field = S.Get<CompositeVector>(field_key_, tag);
  const auto& potential = S.Get<CompositeVector>(potential_key_, tag);

  if (face_coef.hasComponent("cell")) { face_coef.getComponent("cell", true)->putScalar(1.0); }

  double eps = 1.e-16;
  field.scatterMasterToGhosted("cell");
  potential.scatterMasterToGhosted("cell");

  // A boundary face's "other side" coefficient candidate is field_key_'s
  // own "boundary_face" component, seeded into face_coef_f(f,0) and read
  // back inside the loop below.
  bool have_boundary_candidate = field.hasComponent("boundary_face");
  SeedBoundaryFaces_(S, face_coef, /*is_derivative=*/false);

  // The boundary-side POTENTIAL is potential_key_'s own "boundary_face"
  // component (the imposed Dirichlet pressure, same mesh) -- required
  // unconditionally in EnsureCompatibility_ToDeps_ above.
  const auto potential_bf = potential.viewComponent("boundary_face", false);

  const AmanziMesh::MeshCache& m = face_coef.getMesh()->getCache();
  auto face_coef_f = face_coef.viewComponent("face", false);
  const auto potential_c = potential.viewComponent("cell", true);
  const auto cell_coef_c = field.viewComponent("cell", true);

  int nfaces_local = face_coef_f.extent(0);
  Kokkos::parallel_for(
    "EvaluatorUpwindPotentialDifference::Evaluate_", nfaces_local, KOKKOS_LAMBDA(const int& f) {
      auto cells = m.getFaceCells(f);

      if (cells.size() == 1) {
        if (!have_boundary_candidate) {
          // no boundary-side coefficient candidate (field_key_ has no
          // "boundary_face" component) -- the boundary potential alone has
          // nothing to compare it against, so the face value is just the
          // interior cell's value.
          face_coef_f(f, 0) = cell_coef_c(cells[0], 0);
          return;
        }

        // boundary-side potential: the imposed Dirichlet pressure at this
        // same face, read via the mesh's own face<->boundary_face map.
        AmanziMesh::Entity_ID bf = AmanziMesh::getFaceOnBoundaryBoundaryFace(m, f);
        double p_bnd = potential_bf(bf, 0);
        double coef_bnd = face_coef_f(f, 0); // seeded above

        double ol0 = Kokkos::fmax(0., potential_c(cells[0], 0));
        double ol1 = Kokkos::fmax(0., p_bnd);
        double flow_eps = 0.0;
        if ((ol0 > 0) || (ol1 > 0)) { flow_eps = (ol0 * ol1) / (ol0 + ol1); }
        flow_eps = Kokkos::fmax(flow_eps, eps);

        if (potential_c(cells[0], 0) - p_bnd > flow_eps) {
          face_coef_f(f, 0) = cell_coef_c(cells[0], 0);
        } else if (p_bnd - potential_c(cells[0], 0) > flow_eps) {
          face_coef_f(f, 0) = coef_bnd;
        } else {
          double param;
          if (flow_eps < 2 * eps) {
            param = 0.5;
          } else {
            param = (p_bnd - potential_c(cells[0], 0)) / (2 * flow_eps) + 0.5;
          }
          face_coef_f(f, 0) = coef_bnd * param + cell_coef_c(cells[0], 0) * (1. - param);
        }
      } else {
        // overlap region sized by the potential field itself.
        double ol0 = Kokkos::fmax(0., potential_c(cells[0], 0));
        double ol1 = Kokkos::fmax(0., potential_c(cells[1], 0));

        double flow_eps = 0.0;
        if ((ol0 > 0) || (ol1 > 0)) { flow_eps = (ol0 * ol1) / (ol0 + ol1); }
        flow_eps = Kokkos::fmax(flow_eps, eps);

        if (potential_c(cells[0], 0) - potential_c(cells[1], 0) > flow_eps) {
          face_coef_f(f, 0) = cell_coef_c(cells[0], 0);
        } else if (potential_c(cells[1], 0) - potential_c(cells[0], 0) > flow_eps) {
          face_coef_f(f, 0) = cell_coef_c(cells[1], 0);
        } else {
          double param;
          if (flow_eps < 2 * eps) {
            param = 0.5;
          } else {
            param = (potential_c(cells[1], 0) - potential_c(cells[0], 0)) / (2 * flow_eps) + 0.5;
          }
          face_coef_f(f, 0) = cell_coef_c(cells[1], 0) * param + cell_coef_c(cells[0], 0) * (1. - param);
        }
      }
    });
  face_coef.scatterMasterToGhosted("face");

  ApplyClobber_(S, face_coef, /*is_derivative=*/false);
}


void
EvaluatorUpwindPotentialDifference::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(0);
}


void
EvaluatorUpwindPotentialDifference::UpwindDerivative_(const State& S,
                                                      const CompositeVector& dfield_dp,
                                                      CompositeVector& result) const
{
  Tag tag = my_keys_.front().second;
  const auto& potential = S.Get<CompositeVector>(potential_key_, tag);
  double eps = 1.e-16;

  dfield_dp.scatterMasterToGhosted("cell");
  potential.scatterMasterToGhosted("cell");

  bool have_boundary_candidate = dfield_dp.hasComponent("boundary_face");
  const auto potential_bf = potential.viewComponent("boundary_face", false);

  const AmanziMesh::MeshCache& m = result.getMesh()->getCache();
  auto dcoef_faces = result.viewComponent("face", false);
  const auto dfield_cells = dfield_dp.viewComponent("cell", true);
  const auto potential_c = potential.viewComponent("cell", true);

  Kokkos::parallel_for(
    "EvaluatorUpwindPotentialDifference::UpwindDerivative_",
    dcoef_faces.extent(0),
    KOKKOS_LAMBDA(const int& f) {
      auto cells = m.getFaceCells(f);

      if (cells.size() == 1) {
        if (!have_boundary_candidate) {
          dcoef_faces(f, 0) = dfield_cells(cells[0], 0);
          dcoef_faces(f, 1) = 0.0;
          return;
        }

        // slot 0 is d(value)/d(p_cell), slot 1 is d(value)/d(p_boundary) --
        // the boundary candidate's own derivative, pre-loaded into
        // dcoef_faces(f,0) by SeedBoundaryFaces_.
        double dcoef_bnd = dcoef_faces(f, 0);
        AmanziMesh::Entity_ID bf = AmanziMesh::getFaceOnBoundaryBoundaryFace(m, f);
        double p_bnd = potential_bf(bf, 0);

        double ol0 = Kokkos::fmax(0., potential_c(cells[0], 0));
        double ol1 = Kokkos::fmax(0., p_bnd);
        double flow_eps = 0.0;
        if ((ol0 > 0) || (ol1 > 0)) { flow_eps = (ol0 * ol1) / (ol0 + ol1); }
        flow_eps = Kokkos::fmax(flow_eps, eps);

        if (potential_c(cells[0], 0) - p_bnd > flow_eps) {
          dcoef_faces(f, 0) = dfield_cells(cells[0], 0);
          dcoef_faces(f, 1) = 0.0;
        } else if (p_bnd - potential_c(cells[0], 0) > flow_eps) {
          dcoef_faces(f, 0) = 0.0;
          dcoef_faces(f, 1) = dcoef_bnd;
        } else {
          double param;
          if (flow_eps < 2 * eps) {
            param = 0.5;
          } else {
            param = (p_bnd - potential_c(cells[0], 0)) / (2 * flow_eps) + 0.5;
          }
          dcoef_faces(f, 0) = (1. - param) * dfield_cells(cells[0], 0);
          dcoef_faces(f, 1) = param * dcoef_bnd;
        }
        return;
      }

      double ol0 = Kokkos::fmax(0., potential_c(cells[0], 0));
      double ol1 = Kokkos::fmax(0., potential_c(cells[1], 0));
      double flow_eps = 0.0;
      if ((ol0 > 0) || (ol1 > 0)) { flow_eps = (ol0 * ol1) / (ol0 + ol1); }
      flow_eps = Kokkos::fmax(flow_eps, eps);

      double p0 = potential_c(cells[0], 0);
      double p1 = potential_c(cells[1], 0);

      if (p0 - p1 > flow_eps) {
        dcoef_faces(f, 0) = dfield_cells(cells[0], 0);
        dcoef_faces(f, 1) = 0.0;
      } else if (p1 - p0 > flow_eps) {
        dcoef_faces(f, 0) = 0.0;
        dcoef_faces(f, 1) = dfield_cells(cells[1], 0);
      } else {
        double param;
        if (flow_eps < 2 * eps) {
          param = 0.5;
        } else {
          param = (p1 - p0) / (2 * flow_eps) + 0.5;
        }
        dcoef_faces(f, 0) = (1. - param) * dfield_cells(cells[0], 0);
        dcoef_faces(f, 1) = param * dfield_cells(cells[1], 0);
      }
    });
  result.scatterMasterToGhosted("face");
}

} // namespace Relations
} // namespace Amanzi
