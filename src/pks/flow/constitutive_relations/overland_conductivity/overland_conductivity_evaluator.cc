/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "MeshHelpers.hh"
#include "manning_conductivity_model.hh"
#include "overland_conductivity_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

const std::string OverlandConductivityEvaluator::eval_type = "overland conductivity";

OverlandConductivityEvaluator::OverlandConductivityEvaluator(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  if (plist_->isParameter("height key") || plist_->isParameter("ponded depth key") ||
      plist_->isParameter("depth key") || plist_->isParameter("height key suffix") ||
      plist_->isParameter("ponded depth key suffix") || plist_->isParameter("depth key suffix")) {
    Errors::Message message(
      "OverlandConductivity: only use \"mobile depth key\" or \"mobile depth key suffix\", not "
      "\"height key\" or \"ponded depth key\" or \"depth key\".");
    Exceptions::amanzi_throw(message);
  }

  mobile_depth_key_ = Keys::readKey(*plist_, domain, "mobile depth", "ponded_depth");
  dependencies_.insert(KeyTag{ mobile_depth_key_, tag });

  slope_key_ = Keys::readKey(*plist_, domain, "slope", "slope_magnitude");
  dependencies_.insert(KeyTag{ slope_key_, tag });

  coef_key_ = Keys::readKey(*plist_, domain, "coefficient", "manning_coefficient");
  dependencies_.insert(KeyTag{ coef_key_, tag });

  dens_key_ = Keys::readKey(*plist_, domain, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ dens_key_, tag });

  // create the model
  auto sublist = Teuchos::sublist(plist_, "overland conductivity model");
  model_ = Teuchos::rcp(new ManningConductivityModel(sublist));
}


Teuchos::RCP<Evaluator>
OverlandConductivityEvaluator::Clone() const
{
  return Teuchos::rcp(new OverlandConductivityEvaluator(*this));
}


void
OverlandConductivityEvaluator::Evaluate_(const State& S,
                                         const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> depth = S.GetPtr<CompositeVector>(mobile_depth_key_, tag);
  Teuchos::RCP<const CompositeVector> slope = S.GetPtr<CompositeVector>(slope_key_, tag);
  Teuchos::RCP<const CompositeVector> coef = S.GetPtr<CompositeVector>(coef_key_, tag);
  Teuchos::RCP<const CompositeVector> dens = S.GetPtr<CompositeVector>(dens_key_, tag);

  // Note, the logic here splits vectors into two groups -- those that will be
  // defined on all components and those that do not make sense to define on
  // not-cells.  These vectors are used in a wierd way on boundary faces,
  // getting the internal cell and using that.  This currently cannot be used
  // on any other components than "cell" and "boundary_face".
  const AmanziMesh::MeshCache& mesh = result[0]->getMesh()->getCache();
  const ManningConductivityModel& model = *model_;

  for (const auto& comp : *result[0]) {
    AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
    bool is_internal_comp = comp == "boundary_face";
    Key internal_comp = is_internal_comp ? "cell" : comp;

    auto slope_v = slope->viewComponent(internal_comp, false);
    auto coef_v = coef->viewComponent(internal_comp, false);
    auto depth_v = depth->viewComponent(comp, false);
    auto dens_v = dens->viewComponent(comp, false);

    auto result_v = result[0]->viewComponent(comp, false);

    Kokkos::parallel_for(
      "OverlandConductivityEvaluator::Evaluate", result_v.extent(0), KOKKOS_LAMBDA(const int& i) {
        int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
        result_v(i, 0) = dens_v(i, 0) * depth_v(i, 0) *
                         model.Conductivity(depth_v(i, 0), slope_v(ii, 0), coef_v(ii, 0));
      });
  }
}


void
OverlandConductivityEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  // NOTE, we can only differentiate with respect to quantities that exist on
  // all entities, not just cell entities.
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> depth = S.GetPtr<CompositeVector>(mobile_depth_key_, tag);
  Teuchos::RCP<const CompositeVector> slope = S.GetPtr<CompositeVector>(slope_key_, tag);
  Teuchos::RCP<const CompositeVector> coef = S.GetPtr<CompositeVector>(coef_key_, tag);
  Teuchos::RCP<const CompositeVector> dens = S.GetPtr<CompositeVector>(dens_key_, tag);

  const AmanziMesh::MeshCache& mesh = result[0]->getMesh()->getCache();
  const ManningConductivityModel& model = *model_;

  if (wrt_key == mobile_depth_key_) {
    for (const auto& comp : *result[0]) {
      AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
      bool is_internal_comp = comp == "boundary_face";
      Key internal_comp = is_internal_comp ? "cell" : comp;

      auto slope_v = slope->viewComponent(internal_comp, false);
      auto coef_v = coef->viewComponent(internal_comp, false);
      auto depth_v = depth->viewComponent(comp, false);
      auto dens_v = dens->viewComponent(comp, false);

      auto result_v = result[0]->viewComponent(comp, false);

      Kokkos::parallel_for(
        "OverlandConductivityEvaluator::EvaluatePartialDerivative",
        result_v.extent(0),
        KOKKOS_LAMBDA(const int& i) {
          int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
          result_v(i, 0) =
            dens_v(i, 0) * (model.Conductivity(depth_v(i, 0), slope_v(ii, 0), coef_v(ii, 0)) +
                            depth_v(i, 0) * model.DConductivityDDepth(
                                              depth_v(i, 0), slope_v(ii, 0), coef_v(ii, 0)));
        });
    }

  } else if (wrt_key == dens_key_) {
    for (const auto& comp : *result[0]) {
      AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
      bool is_internal_comp = comp == "boundary_face";
      Key internal_comp = is_internal_comp ? "cell" : comp;

      auto slope_v = slope->viewComponent(internal_comp, false);
      auto coef_v = coef->viewComponent(internal_comp, false);
      auto depth_v = depth->viewComponent(comp, false);
      auto dens_v = dens->viewComponent(comp, false);

      auto result_v = result[0]->viewComponent(comp, false);

      Kokkos::parallel_for(
        "OverlandConductivityEvaluator::EvaluatePartialDerivative",
        result_v.extent(0),
        KOKKOS_LAMBDA(const int& i) {
          int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
          result_v(i, 0) =
            depth_v(i, 0) * model.Conductivity(depth_v(i, 0), slope_v(ii, 0), coef_v(ii, 0));
        });
    }

  } else {
    // FIX ME -- need to add derivatives of conductivity model wrt slope, coef
    // -- probably never called, but maybe eventually Manning's n?
    result[0]->putScalar(0.);
  }
}


void
OverlandConductivityEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  // Ensure my field exists.  Requirements should be already set.
  const auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first,
                                                                        my_keys_.front().second);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (my_fac.getMesh() != Teuchos::null) {
    // Create an unowned factory to check my dependencies.
    Teuchos::RCP<CompositeVectorSpace> dep_fac = Teuchos::rcp(new CompositeVectorSpace(my_fac));
    dep_fac->SetOwned(false);

    Teuchos::RCP<CompositeVectorSpace> no_bf_dep_fac;
    if (dep_fac->hasComponent("boundary_face")) {
      no_bf_dep_fac = Teuchos::rcp(new CompositeVectorSpace());
      no_bf_dep_fac->SetMesh(dep_fac->getMesh())
        ->SetGhosted(true)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    } else {
      no_bf_dep_fac = dep_fac;
    }

    // Loop over my dependencies, ensuring they meet the requirements.
    for (const auto& key_tag : dependencies_) {
      if (key_tag.first == coef_key_ || key_tag.first == slope_key_) {
        S.Require<CompositeVector, CompositeVectorSpace>(key_tag.first, key_tag.second)
          .Update(*no_bf_dep_fac);
      } else {
        S.Require<CompositeVector, CompositeVectorSpace>(key_tag.first, key_tag.second)
          .Update(*dep_fac);
      }
    }
  }
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
