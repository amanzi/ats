/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Distributes and downregulates potential transpiration to the rooting zone.
#include "Brent.hh"
#include "transpiration_distribution_relperm_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

SoilPlantFluxFunctor::SoilPlantFluxFunctor(AmanziMesh::Entity_ID sc_,
        const AmanziMesh::Entity_ID_View& cells_of_col_,
        const LandCover& lc_,
        const Epetra_MultiVector& soil_wp_,
        const Epetra_MultiVector& soil_kr_,
        const Epetra_MultiVector& f_root_,
        const Epetra_MultiVector& pet_,
        double c0_,
        double krp_)
  : sc(sc_),
    cells_of_col(cells_of_col_),
    lc(lc_),
    soil_kr(soil_kr_),
    soil_wp(soil_wp_),
    f_root(f_root_),
    pet(pet_),
    c0(c0_),
    krp(krp_)
{}


// Functor used in the solve for plant water potential
//
// Note the sign here -- pet should not change, as plant water potential
// increases, soil->plant flux should increase, so total_trans should increase,
// so this function is hopefully monotonically decreasing with plant_wp
double
SoilPlantFluxFunctor::operator()(double plant_wp) const
{
  double beta = transpirationReductionFunction(plant_wp);
  double total_trans = 0.;
  for (auto c : cells_of_col) {
    total_trans += soilPlantFlux(plant_wp, c);
  }
  return beta * pet[0][sc] - total_trans;
}


// right hand side
double
SoilPlantFluxFunctor::soilPlantFlux(double plant_wp, AmanziMesh::Entity_ID c) const
{
  double kr = (plant_wp > soil_wp[0][c]) ? soil_kr[0][c] : krp;
  return c0 * f_root[0][c] * kr * (plant_wp - soil_wp[0][c]);
}


// beta in left hand side
double
SoilPlantFluxFunctor::transpirationReductionFunction(double plant_wp) const
{
  if (plant_wp <= lc.stomata_open_water_potential) {
    return 1.0;
  } else if (plant_wp >= lc.stomata_closed_water_potential) {
    return 0.0;
  } else {
    return (plant_wp - lc.stomata_closed_water_potential)
      / (lc.stomata_open_water_potential - lc.stomata_closed_water_potential);
  }
}


// Constructor from ParameterList
TranspirationDistributionRelPermEvaluator::TranspirationDistributionRelPermEvaluator(
  Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  InitializeFromPlist_();
}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
TranspirationDistributionRelPermEvaluator::Clone() const
{
  return Teuchos::rcp(new TranspirationDistributionRelPermEvaluator(*this));
}


// Initialize by setting up dependencies
void
TranspirationDistributionRelPermEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  domain_sub_ = Keys::getDomain(my_keys_.front().first);
  domain_surf_ = Keys::readDomainHint(plist_, domain_sub_, "domain", "surface");
  Tag tag = my_keys_.front().second;

  // my_keys_.front() will always (for now?) be transpiration...
  Key plant_wp_key =
    Keys::readKey(plist_, domain_surf_, "plant water potential", "plant_water_potential");
  my_keys_.emplace_back(std::make_pair(plant_wp_key, tag));

  // - pull Keys from plist
  // dependency: soil p
  soil_wp_key_ =
    Keys::readKey(plist_, domain_sub_, "soil water potential", "capillary_pressure_gas_liq");
  dependencies_.insert(KeyTag{ soil_wp_key_, tag });

  // dependency: soil rel perm
  soil_kr_key_ =
    Keys::readKey(plist_, domain_sub_, "soil relative permeability", "relative_permeability");
  dependencies_.insert(KeyTag{ soil_kr_key_, tag });

  // dependency: rooting_depth_fraction
  f_root_key_ =
    Keys::readKey(plist_, domain_sub_, "rooting depth fraction", "rooting_depth_fraction");
  dependencies_.insert(KeyTag{ f_root_key_, tag });

  // dependency: transpiration
  potential_trans_key_ =
    Keys::readKey(plist_, domain_surf_, "potential transpiration", "potential_transpiration_mols");
  dependencies_.insert(KeyTag{ potential_trans_key_, tag });

  // dependency: cell volume, surface cell volume
  cv_key_ = Keys::readKey(plist_, domain_sub_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
  surf_cv_key_ = Keys::readKey(plist_, domain_surf_, "surface cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ surf_cv_key_, tag });

  // c0, maximal trans rate
  c0_ = plist_.get<double>("total maximal conductance [mol m^-2 s^-1 MPa^-1]", 10.) * 1.e-6; // per MPa --> per Pa
  krp_ = plist_.get<double>("plant relative conductance [-]", 0.);
}


void
TranspirationDistributionRelPermEvaluator::Evaluate_(const State& S,
                                                     const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  // on the subsurface
  const Epetra_MultiVector& soil_wp =
    *S.Get<CompositeVector>(soil_wp_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& soil_kr =
    *S.Get<CompositeVector>(soil_kr_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& f_root =
    *S.Get<CompositeVector>(f_root_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);

  // on the surface
  const Epetra_MultiVector& potential_trans =
    *S.Get<CompositeVector>(potential_trans_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& surf_cv =
    *S.Get<CompositeVector>(surf_cv_key_, tag).ViewComponent("cell", false);

  Epetra_MultiVector& trans_v = *result[0]->ViewComponent("cell", false);
  Epetra_MultiVector& plant_wp_v = *result[1]->ViewComponent("cell", false);

  auto& subsurf_mesh = *S.GetMesh(domain_sub_);
  auto& surf_mesh = *S.GetMesh(domain_surf_);

  for (const auto& region_lc : land_cover_) {
    auto lc_ids = surf_mesh.getSetEntities(
      region_lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    if (lc_ids.size() > 0) {
      for (int sc : lc_ids) {
        if (potential_trans[0][sc] > 0. || krp_ > 0.) {
          SoilPlantFluxFunctor func(sc, subsurf_mesh.columns.getCells(sc),
                  region_lc.second, soil_wp, soil_kr, f_root, potential_trans, c0_, krp_);

          // bracket the root -- linear to the left of 1, log to the right of 1
          int itrs = 1.e4;
          std::pair<double,double> ab;
          if (func(1) > 0) {
            // f(1) is positive, and function is decreasing, so look to the right of 1 in log-space
            auto func2 = [&](double x) { return func(std::pow(10, x)); };
            ab = Amanzi::Utils::bracketRoot(func2, 1, 1, &itrs);
            ab.first = std::pow(10, ab.first);
            ab.second = std::pow(10, ab.second);
          } else {
            // f(1) is negative, and function is decreasing, so look to the left of 1 in linear space
            ab = Amanzi::Utils::bracketRoot(func, 1, 1.e4, &itrs);
          }
          if (itrs >= 1.e4) {
            std::cout << "failing to bracket root:" << std::endl
                      << "  start = " << 0. << " f_start = " << func(0.) << std::endl
                      << "  a = " << ab.first << " fa = " << func(ab.first) << std::endl
                      << "  b = " << ab.second << " fb = " << func(ab.second) << std::endl;
            AMANZI_ASSERT(itrs < 1.e4);
          }

          // compute the value
          itrs = 100;
          plant_wp_v[0][sc] = Amanzi::Utils::computeRootBrent(func, ab.first, ab.second, 1.e-6, &itrs);
          AMANZI_ASSERT(itrs > 0 && itrs < 100);

          int i = 0;
          for (auto c : subsurf_mesh.columns.getCells(sc)) {
            trans_v[0][c] = surf_cv[0][sc] * func.soilPlantFlux(plant_wp_v[0][sc], i++) / cv[0][c];
          }
        } else {
          for (auto c : subsurf_mesh.columns.getCells(sc)) {
            trans_v[0][c] = 0.;
          }
        }
      }
    }
  }
}


void
TranspirationDistributionRelPermEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(
    0.); // this would be a nontrivial calculation, as it is technically nonlocal due to rescaling issues?
}


void
TranspirationDistributionRelPermEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  Tag tag = my_keys_.front().second;

  // new state!
  if (land_cover_.size() == 0)
    land_cover_ =
      getLandCover(S.ICList().sublist("land cover types"),
                   { "stomata_closed_water_potential", "stomata_open_water_potential" });

  Key domain = Keys::getDomain(my_keys_.front().first);

  // Create an unowned factory to check my dependencies.
  // -- first those on the subsurface mesh
  CompositeVectorSpace dep_fac;
  dep_fac.SetMesh(S.GetMesh(domain))->AddComponent("cell", AmanziMesh::CELL, 1);
  S.Require<CompositeVector, CompositeVectorSpace>(f_root_key_, tag).Update(dep_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(soil_wp_key_, tag).Update(dep_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(soil_kr_key_, tag).Update(dep_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(cv_key_, tag).Update(dep_fac);

  // -- next those on the surface mesh
  CompositeVectorSpace surf_fac;
  surf_fac.SetMesh(S.GetMesh(Keys::getDomain(surf_cv_key_)))
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S.Require<CompositeVector, CompositeVectorSpace>(potential_trans_key_, tag).Update(surf_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(surf_cv_key_, tag).Update(surf_fac);
}


void
TranspirationDistributionRelPermEvaluator::EnsureCompatibility_Structure_(State& S)
{
  Tag tag = my_keys_.front().second;

  // my_keys_[0] is transpiration
  S.Require<CompositeVector,CompositeVectorSpace>(my_keys_.front().first, tag)
    .SetMesh(S.GetMesh(domain_sub_))
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // my_keys_[1] is plant water potential
  S.Require<CompositeVector,CompositeVectorSpace>(my_keys_.back().first, tag)
    .SetMesh(S.GetMesh(domain_surf_))
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
