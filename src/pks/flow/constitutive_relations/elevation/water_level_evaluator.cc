/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Bo Gao (gaob@ornl.gov)
*/

/*
  Evaluator for determining water level, a combined variable of
  water table and ponded depth.

*/

#include "height_model.hh"
#include "water_level_evaluator.hh"

namespace Amanzi {
namespace Flow {


WaterLevelEvaluator::WaterLevelEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  std::string domain_surf = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;
  Key domain_ss = Keys::readDomainHint(plist, domain_surf, "surface", "subsurface");

  // my dependencies
  dens_key_ = Keys::readKey(plist_, domain_surf, "mass density", "mass_density_liquid");
  dependencies_.insert(KeyTag{ dens_key_, tag });

  pres_key_ = Keys::readKey(plist_, domain_surf, "pressure", "pressure");
  dependencies_.insert(KeyTag{ pres_key_, tag });

  sat_key_ = Keys::readKey(plist_, domain_ss, "gas saturation", "saturation_gas");
  dependencies_.insert(KeyTag{ sat_key_, tag });

  cv_key_ = Keys::readKey(plist, domain_ss, "subsurface cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });

  cv_surf_key_ = Keys::readKey(plist, domain_surf, "surface cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_surf_key_, tag });

  // model
  Teuchos::ParameterList model_plist = plist_.sublist("height model parameters");
  model_ = Teuchos::rcp(new HeightModel(model_plist));
}


Teuchos::RCP<Evaluator>
WaterLevelEvaluator::Clone() const
{
  return Teuchos::rcp(new WaterLevelEvaluator(*this));
}


void
WaterLevelEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> pres = S.GetPtr<CompositeVector>(pres_key_, tag);

  // this is rather hacky.  surface_pressure is a mixed field vector -- it has
  // pressure on cells and ponded depth on faces.
  // -- copy the faces over directly
  if (result[0]->HasComponent("face") && pres->HasComponent("face"))
    *result[0]->ViewComponent("face", false) = *pres->ViewComponent("face", false);

  // -- cells need the function eval for ponded depth
  const Epetra_MultiVector& pres_c = *pres->ViewComponent("cell", false);
  const Epetra_MultiVector& rho =
    *S.GetPtr<CompositeVector>(dens_key_, tag)->ViewComponent("cell", false);
  double p_atm = S.Get<double>("atmospheric_pressure", Tags::DEFAULT);
  const AmanziGeometry::Point& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double gz = -gravity[2];

  // -- cells need the function eval for water table depth
  Teuchos::RCP<const CompositeVector> sat = S.GetPtr<CompositeVector>(sat_key_, tag);
  Teuchos::RCP<const CompositeVector> cv = S.GetPtr<CompositeVector>(cv_key_, tag);
  Teuchos::RCP<const CompositeVector> cv_surf = S.GetPtr<CompositeVector>(cv_surf_key_, tag);
  const Epetra_MultiVector& sat_c = *sat->ViewComponent("cell", false);
  const Epetra_MultiVector& cv_c = *cv->ViewComponent("cell", false);
  const Epetra_MultiVector& cv_surf_c = *cv_surf->ViewComponent("cell", false);
  auto mesh = result[0]->Mesh()->getParentMesh();

  const Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell", false);
  int ncells = res_c.MyLength();
  for (int sc = 0; sc != ncells; ++sc) {
    if (pres_c[0][sc] < p_atm) {
      // for each column, loop over cells until stop is requested or column end
      AmanziGeometry::Point val(0., 0.);
      auto col_cell = mesh->columns.getCells(sc);
      int i = 0;
      while (sat_c[0][col_cell[i]] > 0.0 && i != col_cell.size()) {
        val[0] += cv_c[0][col_cell[i]];
        i++;
      }
      // val[1] is typically e.g. cell volume, but can be 0 to indicate no
      // denominator.  Coefficient provides a hook for column-wide multiples
      // (e.g. 1/surface area).
      if (val[1] > 0.)
        res_c[0][sc] = -val[1] / cv_surf_c[0][sc] * val[0];
      else
        res_c[0][sc] = -1. / cv_surf_c[0][sc] * val[0];
    } else {
      // ponded depth
      res_c[0][sc] = model_->Height(pres_c[0][sc], rho[0][sc], p_atm, gz);
    }
  }
}


void
WaterLevelEvaluator::EvaluatePartialDerivative_(const State& S,
                                                const Key& wrt_key,
                                                const Tag& wrt_tag,
                                                const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  // -- cells need the function eval for ponded depth
  const Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell", false);
  const Epetra_MultiVector& pres_c =
    *S.GetPtr<CompositeVector>(pres_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& rho =
    *S.GetPtr<CompositeVector>(dens_key_, tag)->ViewComponent("cell", false);
  double p_atm = S.Get<double>("atmospheric_pressure", Tags::DEFAULT);
  const AmanziGeometry::Point& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double gz = -gravity[2];

  int ncells = res_c.MyLength();
  for (int sc = 0; sc != ncells; ++sc) {
    if (pres_c[0][sc] >= p_atm) {
      if (wrt_key == pres_key_) {
        res_c[0][sc] = model_->DHeightDPressure(pres_c[0][sc], rho[0][sc], p_atm, gz);
      } else if (wrt_key == dens_key_) {
        res_c[0][sc] = model_->DHeightDRho(pres_c[0][sc], rho[0][sc], p_atm, gz);
      } else {
        AMANZI_ASSERT(0);
      }
    } else {
      AMANZI_ASSERT(0);
    }
  }
}


void
WaterLevelEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  const auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first,
                                                                     my_keys_.front().second);
  if (fac.Mesh() != Teuchos::null) {
    CompositeVectorSpace dep_fac;
    dep_fac.SetMesh(fac.Mesh()->getParentMesh())
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    for (const auto& dep : dependencies_) {
      if (Keys::getDomain(dep.first) == Keys::getDomain(my_keys_.front().first)) {
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(fac);
      } else {
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(dep_fac);
      }
    }
  }
}

} // namespace Flow
} // namespace Amanzi
