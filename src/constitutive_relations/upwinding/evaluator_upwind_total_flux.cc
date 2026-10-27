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

#include "evaluator_upwind_total_flux.hh"

namespace Amanzi {
namespace Relations {

const std::string EvaluatorUpwindTotalFlux::eval_type = "upwind: total flux";

EvaluatorUpwindTotalFlux::EvaluatorUpwindTotalFlux(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorUpwinding(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  flux_dir_key_ = Keys::readKey(*plist, domain, "flux direction", "water_flux_direction");
  dependencies_.insert(KeyTag{ flux_dir_key_, tag });

  flux_eps_ = plist->get<double>("flux epsilon", 1.e-5);
}


Teuchos::RCP<Evaluator>
EvaluatorUpwindTotalFlux::Clone() const
{
  return Teuchos::rcp(new EvaluatorUpwindTotalFlux(*this));
}


bool
EvaluatorUpwindTotalFlux::IsDifferentiableWRT(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag) const
{
  // never differentiate through the flux direction -- it is a direction
  // indicator only.
  if (wrt_key == flux_dir_key_) return false;
  return EvaluatorUpwinding::IsDifferentiableWRT(S, wrt_key, wrt_tag);
}


void
EvaluatorUpwindTotalFlux::EnsureCompatibility_ToDeps_(State& S)
{
  EvaluatorUpwinding::EnsureCompatibility_ToDeps_(S);

  Tag tag = my_keys_.front().second;
  CompositeVectorSpace flux_fac;
  flux_fac.SetMesh(S.GetMesh(Keys::getDomain(flux_dir_key_)))
    ->SetGhosted()
    ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  S.Require<CompositeVector, CompositeVectorSpace>(flux_dir_key_, tag).Update(flux_fac);
}


void
EvaluatorUpwindTotalFlux::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  CompositeVector& face_coef = *result[0];

  const auto& field = S.Get<CompositeVector>(field_key_, tag);
  const auto& flux_dir = S.Get<CompositeVector>(flux_dir_key_, tag);

  if (face_coef.hasComponent("cell")) { face_coef.getComponent("cell", true)->putScalar(1.0); }

  field.scatterMasterToGhosted("cell");

  // Unconditionally seed boundary faces from field_key_'s own boundary_face
  // component BEFORE the upwind loop: a boundary face has no second
  // neighbor cell, so the fallback read of coef_faces(f,0) below (for the
  // missing side, in the hard-upwind and blended cases) needs this primed.
  SeedBoundaryFaces_(S, face_coef, /*is_derivative=*/false);

  const AmanziMesh::MeshCache& m = face_coef.getMesh()->getCache();
  const auto flux_v = flux_dir.viewComponent("face", false);
  auto coef_faces = face_coef.viewComponent("face", false);
  const auto coef_cells = field.viewComponent("cell", true);
  double flow_eps = flux_eps_;

  int nfaces_local = coef_faces.extent(0);
  Kokkos::parallel_for(
    "EvaluatorUpwindTotalFlux::Evaluate_", nfaces_local, KOKKOS_LAMBDA(const int& f) {
      auto fcells = m.getFaceCells(f);

      int uw = -1, dw = -1;
      int c0 = fcells(0);
      int orientation = 0;
      m.getFaceNormal(f, c0, &orientation);
      if (flux_v(f, 0) * orientation > 0) {
        uw = c0;
        if (fcells.size() == 2) dw = fcells(1);
      } else {
        dw = c0;
        if (fcells.size() == 2) uw = fcells(1);
      }
      assert(!((uw == -1) && (dw == -1)));

      double coefs[2];
      coefs[0] = (uw == -1) ? coef_faces(f, 0) : coef_cells(uw, 0);
      coefs[1] = (dw == -1) ? coef_faces(f, 0) : coef_cells(dw, 0);

      if (Kokkos::abs(flux_v(f, 0)) >= flow_eps) {
        coef_faces(f, 0) = coefs[0];
      } else {
        double param = Kokkos::abs(flux_v(f, 0)) / (2 * flow_eps) + 0.5;
        coef_faces(f, 0) = coefs[0] * param + coefs[1] * (1. - param);
      }
    });
  face_coef.scatterMasterToGhosted("face");

  ApplyClobber_(S, face_coef, /*is_derivative=*/false);
}


void
EvaluatorUpwindTotalFlux::EvaluatePartialDerivative_(const State& S,
                                                     const Key& wrt_key,
                                                     const Tag& wrt_tag,
                                                     const std::vector<CompositeVector*>& result)
{
  // Not used: EvaluatorUpwinding::UpdateDerivative_ is overridden to call
  // UpwindDerivative_ directly rather than going through the generic
  // elementWiseMultiply chain-rule loop that would call this.
  AMANZI_ASSERT(0);
}


void
EvaluatorUpwindTotalFlux::UpwindDerivative_(const State& S,
                                            const CompositeVector& dfield_dp,
                                            CompositeVector& result) const
{
  Tag tag = my_keys_.front().second;
  const auto& flux_dir = S.Get<CompositeVector>(flux_dir_key_, tag);

  dfield_dp.scatterMasterToGhosted("cell");

  const AmanziMesh::MeshCache& m = result.getMesh()->getCache();
  const auto flux_v = flux_dir.viewComponent("face", false);
  auto dcoef_faces = result.viewComponent("face", false);
  const auto dfield_cells = dfield_dp.viewComponent("cell", true);
  double flow_eps = flux_eps_;

  int nfaces_local = dcoef_faces.extent(0);
  Kokkos::parallel_for(
    "EvaluatorUpwindTotalFlux::UpwindDerivative_", nfaces_local, KOKKOS_LAMBDA(const int& f) {
      auto fcells = m.getFaceCells(f);

      if (fcells.size() == 1) {
        // boundary face: slot 0 is d(value)/d(p_cell), slot 1 is
        // d(value)/d(p_boundary) -- the boundary candidate's own
        // derivative, pre-loaded into dcoef_faces(f,0) by
        // SeedBoundaryFaces_.
        double dcoef_boundary = dcoef_faces(f, 0);
        int c0 = fcells(0);
        int orientation = 0;
        m.getFaceNormal(f, c0, &orientation);
        if (Kokkos::abs(flux_v(f, 0)) >= flow_eps) {
          if (flux_v(f, 0) * orientation > 0) {
            dcoef_faces(f, 0) = dfield_cells(c0, 0);
            dcoef_faces(f, 1) = 0.0;
          } else {
            dcoef_faces(f, 0) = 0.0;
            dcoef_faces(f, 1) = dcoef_boundary;
          }
        } else {
          double param = Kokkos::abs(flux_v(f, 0)) / (2 * flow_eps) + 0.5;
          if (flux_v(f, 0) * orientation > 0) {
            dcoef_faces(f, 0) = param * dfield_cells(c0, 0);
            dcoef_faces(f, 1) = (1. - param) * dcoef_boundary;
          } else {
            dcoef_faces(f, 0) = (1. - param) * dfield_cells(c0, 0);
            dcoef_faces(f, 1) = param * dcoef_boundary;
          }
        }
        return;
      }

      // interior face, 2 neighbor cells, ordered by getFaceCells(f) --
      // cells[0]/cells[1] is a fixed topological convention.
      int c0 = fcells(0);
      int orientation = 0;
      m.getFaceNormal(f, c0, &orientation);
      bool c0_is_upwind = (flux_v(f, 0) * orientation > 0);

      double dcell[2] = { dfield_cells(fcells(0), 0), dfield_cells(fcells(1), 0) };

      if (Kokkos::abs(flux_v(f, 0)) >= flow_eps) {
        // hard upwind: derivative goes entirely to the upwind cell's slot.
        if (c0_is_upwind) {
          dcoef_faces(f, 0) = dcell[0];
          dcoef_faces(f, 1) = 0.0;
        } else {
          dcoef_faces(f, 0) = 0.0;
          dcoef_faces(f, 1) = dcell[1];
        }
      } else {
        // blended zone: differentiate the linear blend, ignoring the
        // (small, flux-eps-scale) dependence of param itself on p.
        double param = Kokkos::abs(flux_v(f, 0)) / (2 * flow_eps) + 0.5;
        if (c0_is_upwind) {
          dcoef_faces(f, 0) = param * dcell[0];
          dcoef_faces(f, 1) = (1. - param) * dcell[1];
        } else {
          dcoef_faces(f, 0) = (1. - param) * dcell[0];
          dcoef_faces(f, 1) = param * dcell[1];
        }
      }
    });
  result.scatterMasterToGhosted("face");
}

} // namespace Relations
} // namespace Amanzi
