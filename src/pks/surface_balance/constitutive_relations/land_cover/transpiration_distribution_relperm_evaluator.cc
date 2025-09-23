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
namespace ATS_Physics {
namespace SurfaceBalance {
namespace Relations {

SoilPlantFluxFunctor::SoilPlantFluxFunctor(AmanziMesh::Entity_ID sc_,
                                           const AmanziMesh::Entity_ID_View& cells_of_col_,
                                           const LandCover& lc_,
                                           const Epetra_MultiVector& soil_pc_,
                                           const Epetra_MultiVector& soil_kr_,
                                           const Epetra_MultiVector& f_root_,
                                           const Epetra_MultiVector& pet_,
                                           const Epetra_MultiVector& rho_,
                                           const Epetra_MultiVector& nliq_,
                                           const Epetra_MultiVector& visc_,
                                           const Epetra_MultiVector& cv_,
                                           const Epetra_MultiVector& sa_,
                                           double K_,
                                           double krp_, // includes permeability_rescaling factor
                                           double g_)
  : sc(sc_),
    cells_of_col(cells_of_col_),
    lc(lc_),
    soil_kr(soil_kr_),
    soil_pc(soil_pc_),
    f_root(f_root_),
    pet(pet_),
    rho(rho_),
    nliq(nliq_),
    visc(visc_),
    cv(cv_),
    sa(sa_),
    K(K_),
    krp(krp_),
    g(g_)
{}


// Functor used in the solve for plant water potential
//
// Note the sign here -- pet should not change, as plant water potential
// increases, soil->plant flux should increase, so total_trans should increase,
// so this function is hopefully monotonically decreasing with plant_pc
double
SoilPlantFluxFunctor::operator()(double plant_pc) const
{
  double total_trans = computeSoilPlantFluxes(plant_pc);
  return pet[0][sc] - total_trans;
}


// right hand side
double
SoilPlantFluxFunctor::computeSoilPlantFlux(double root_pc, AmanziMesh::Entity_ID c) const
{
  double kr = (root_pc > soil_pc[0][c]) ? soil_kr[0][c] : (krp * nliq[0][c] / visc[0][c]);
  return -K * f_root[0][c] * kr * (soil_pc[0][c] - root_pc);
}


// all right hand sides
double
SoilPlantFluxFunctor::computeSoilPlantFluxes(double plant_pc, Epetra_MultiVector* res) const
{
  double total_trans = 0.;
  double root_pc = plant_pc;
  for (auto c : cells_of_col) {
    double Mg_dz_on_2 = rho[0][c] * g * cv[0][c] / (sa[0][sc] * 2);

    // top half-cell
    root_pc += Mg_dz_on_2;

    // compute flux
    double local_trans = computeSoilPlantFlux(root_pc, c);
    total_trans += local_trans;
    if (res) (*res)[0][c] = local_trans * sa[0][sc] / cv[0][c];

    // bottom half-cell
    root_pc += Mg_dz_on_2;
  }
  return total_trans;
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

  // tolerances for root finding
  tol_ = plist_.get<double>("tolerance", 1.e-12);
  nits_ = plist_.get<int>("maximum number of iterations", 100);

  // my_keys_.front() will always (for now?) be transpiration...
  Key plant_pc_key =
    Keys::readKey(plist_, domain_surf_, "plant capillary pressure", "capillary_pressure_plant");
  my_keys_.emplace_back(std::make_pair(plant_pc_key, tag));

  // - pull Keys from plist
  // dependency: soil p
  soil_pc_key_ =
    Keys::readKey(plist_, domain_sub_, "soil capillary pressure", "capillary_pressure_gas_liq");
  dependencies_.insert(KeyTag{ soil_pc_key_, tag });

  // dependency: soil rel perm
  soil_kr_key_ =
    Keys::readKey(plist_, domain_sub_, "soil relative permeability", "relative_permeability");
  dependencies_.insert(KeyTag{ soil_kr_key_, tag });

  // dependency: rooting_depth_fraction
  f_root_key_ = Keys::readKey(plist_, domain_sub_, "plant rooting fraction", "root_fraction");
  dependencies_.insert(KeyTag{ f_root_key_, tag });

  // dependency: liquid molar density
  rho_key_ = Keys::readKey(plist_, domain_sub_, "mass density", "mass_density_liquid");
  dependencies_.insert(KeyTag{ rho_key_, tag });

  // dependency: liquid molar density
  nliq_key_ = Keys::readKey(plist_, domain_sub_, "molar density", "molar_density_liquid");
  dependencies_.insert(KeyTag{ nliq_key_, tag });

  // dependency: viscosity
  visc_key_ = Keys::readKey(plist_, domain_sub_, "viscosity", "viscosity_liquid");
  dependencies_.insert(KeyTag{ visc_key_, tag });

  // dependency: cv, sa
  cv_key_ = Keys::readKey(plist_, domain_sub_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
  sa_key_ = Keys::readKey(plist_, domain_surf_, "surface area", "cell_volume");
  dependencies_.insert(KeyTag{ sa_key_, tag });

  // dependency: transpiration
  potential_trans_key_ =
    Keys::readKey(plist_, domain_surf_, "potential transpiration", "potential_transpiration_mols");
  dependencies_.insert(KeyTag{ potential_trans_key_, tag });

  // c0, maximal trans rate
  K_ = plist_.get<double>("plant permeability per m [m]", 1.e-12); // ? what is the equivalent?
  krp_ = plist_.get<double>("plant relative conductance [-]", 0.); // 1 = redistribution, 0 = none
}


void
TranspirationDistributionRelPermEvaluator::Evaluate_(const State& S,
                                                     const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  // on the subsurface
  const Epetra_MultiVector& soil_pc =
    *S.Get<CompositeVector>(soil_pc_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& soil_kr =
    *S.Get<CompositeVector>(soil_kr_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& f_root =
    *S.Get<CompositeVector>(f_root_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& rho =
    *S.Get<CompositeVector>(rho_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& nliq =
    *S.Get<CompositeVector>(nliq_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& visc =
    *S.Get<CompositeVector>(visc_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& sa = *S.Get<CompositeVector>(sa_key_, tag).ViewComponent("cell", false);

  const auto& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double g = gravity[gravity.dim() - 1];

  double perm_scale = S.Get<double>("permeability_rescaling", Tags::DEFAULT);

  // on the surface
  const Epetra_MultiVector& potential_trans =
    *S.Get<CompositeVector>(potential_trans_key_, tag).ViewComponent("cell", false);

  Epetra_MultiVector& trans_v = *result[0]->ViewComponent("cell", false);
  Epetra_MultiVector& plant_pc_v = *result[1]->ViewComponent("cell", false);

  auto& subsurf_mesh = *S.GetMesh(domain_sub_);
  auto& surf_mesh = *S.GetMesh(domain_surf_);

  for (const auto& region_lc : land_cover_) {
    auto lc_ids = surf_mesh.getSetEntities(
      region_lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    if (lc_ids.size() > 0) {
      for (int sc : lc_ids) {
        if (potential_trans[0][sc] > 0. || krp_ > 0.) {
          SoilPlantFluxFunctor func(sc,
                                    subsurf_mesh.columns.getCells(sc),
                                    region_lc.second,
                                    soil_pc,
                                    soil_kr,
                                    f_root,
                                    potential_trans,
                                    rho,
                                    nliq,
                                    visc,
                                    cv,
                                    sa,
                                    K_ * perm_scale,
                                    krp_ / perm_scale,
                                    g);

          // compute the flux at max pc
          double pc_plant_max = region_lc.second.maximum_xylem_capillary_pressure;
          double total_trans = func.computeSoilPlantFluxes(pc_plant_max, &trans_v);
          if (total_trans <= potential_trans[0][sc]) {
            // insufficient water
            plant_pc_v[0][sc] = pc_plant_max;
          } else {
            // sufficient water -- solve for plant_pc

            // bracket the root -- linear to the left of 1, log to the right of 1
            std::pair<double, double> ab;
            if (func(0.) > 0.) {
              ab.first = 0.;
              ab.second = 1.e4;
              while (func(ab.second) > 0) {
                ab.first = ab.second;
                ab.second *= 10;
                if (ab.second > 1.e14) {
                  Errors::Message msg;
                  msg << "TranspirationDistributionRelPermEvaluator:: Failing to bracket root: "
                      << "start = " << 0. << " f_start = " << func(0.) << "  stop = " << ab.second
                      << " f_stop = " << func(ab.second);
                  Exceptions::amanzi_throw(msg);
                }
              }
            } else {
              ab.second = 0.;
              ab.first = -1.e4;
              while (func(ab.first) < 0) {
                ab.second = ab.first;
                ab.first *= 10;
                if (ab.first < -1.e14) {
                  Errors::Message msg;
                  msg << "TranspirationDistributionRelPermEvaluator:: Failing to bracket root: "
                      << "start = " << 0. << " f_start = " << func(0.) << "  stop = " << ab.first
                      << " f_stop = " << func(ab.first);
                  Exceptions::amanzi_throw(msg);
                }
              }
            }

            // compute the plant capillary pressure using a root-finder
            int itrs = nits_;
            plant_pc_v[0][sc] =
              Amanzi::Utils::findRootBrent(func, ab.first, ab.second, tol_, &itrs);
            AMANZI_ASSERT(itrs > 0 && itrs <= nits_);

            // compute the distributed transpiration fluxes for each grid cell
            func.computeSoilPlantFluxes(plant_pc_v[0][sc], &trans_v);
          }
        } else {
          for (auto c : subsurf_mesh.columns.getCells(sc)) {
            trans_v[0][c] = 0.;
          }
        }
      }
    }
  }
  AMANZI_ASSERT(!std::isnan(trans_v[0][0]));
}


void
TranspirationDistributionRelPermEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  // this would be a nontrivial calculation, as it is technically nonlocal due
  // to rescaling issues?
  result[0]->PutScalar(0.);
}


void
TranspirationDistributionRelPermEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  Tag tag = my_keys_.front().second;

  // new state!
  if (land_cover_.size() == 0) {
    land_cover_ = getLandCover(S.GetModelParameters("land cover types"),
                               { "maximum_xylem_capillary_pressure" });
  }

  // Create an unowned factory to check my dependencies.
  // -- first those on the subsurface mesh
  CompositeVectorSpace dep_fac;
  dep_fac.SetMesh(S.GetMesh(domain_sub_))->AddComponent("cell", AmanziMesh::CELL, 1);
  S.Require<CompositeVector, CompositeVectorSpace>(f_root_key_, tag).Update(dep_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(soil_pc_key_, tag).Update(dep_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(soil_kr_key_, tag).Update(dep_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(rho_key_, tag).Update(dep_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(nliq_key_, tag).Update(dep_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(visc_key_, tag).Update(dep_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(cv_key_, tag).Update(dep_fac);

  // -- next those on the surface mesh
  CompositeVectorSpace surf_fac;
  surf_fac.SetMesh(S.GetMesh(domain_surf_))->AddComponent("cell", AmanziMesh::CELL, 1);
  S.Require<CompositeVector, CompositeVectorSpace>(potential_trans_key_, tag).Update(surf_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(sa_key_, tag).Update(surf_fac);

  // gravity
  S.Require<AmanziGeometry::Point>("gravity", Tags::DEFAULT);

  // perm rescaling
  S.Require<double>("permeability_rescaling", Tags::DEFAULT);
}


void
TranspirationDistributionRelPermEvaluator::EnsureCompatibility_Structure_(State& S)
{
  Tag tag = my_keys_.front().second;

  // my_keys_[0] is transpiration
  S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first, tag)
    .SetMesh(S.GetMesh(domain_sub_))
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // my_keys_[1] is plant water potential
  S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.back().first, tag)
    .SetMesh(S.GetMesh(domain_surf_))
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace ATS_Physics
} // namespace Amanzi
