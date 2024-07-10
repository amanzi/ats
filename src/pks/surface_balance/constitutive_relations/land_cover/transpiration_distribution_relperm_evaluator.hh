/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Distributes potential transpiration to the rooting zone using a rel perm weighting method.
/*!

The transpiration distribution evaluator looks to take a potential
evapotranspiration and distribute it across the vertical column based on water
availability and rooting depths.  It also potentially limits the transpiration
to avoid taking water where it is not available (thereby crashing the code).

This is a new model developed by Painter and Coon that uses a competing rates
approach.  It assumes that the water potential in the plant is constant, and
that the flux of water between any given grid cell and the plant is given by a
potential difference multiplied by a conductance:

.. math::
   Q_{T} = Q_{pot T} \beta(\Psi_p) = \sum_{i} Q_{i} = c_i k_r(\Psi_i, \Psi_p + \rho g z_i)[\Psi_i - \Psi_p]

where :math:`Q_{T}` is the total actual transpiration, :math:`Q_{pot T}` the
potential transpiration, :math:`\beta` a transpiration reduction function that
is only a function of the plant water potential, :math:`\Psi_p` and
:math:`\Psi_i` the water potential (capillary pressure) in the plant and the
soil cell :math:`i`, respectively, :math:`\rho` is the density of water,
:math:`g` the gravitational force (so that :math:`\Psi_p + \rho g z_i` is the
pressure in the root in soil cell :math:`i`), :math:`Q_i` is the transpiration
taken from each grid cell, :math:`c_i` is a maximum conductance of cell i,
which is proportional to the root fraction in cell :math:`i`, and :math:`k_r`
is the relative permeability (which is often upwinded between soil relative
permeability and a plant relative permeability that can be used to control
hydraulic redistribution).

The second and fourth terms above are set equal, and solved using a scalar
root-finding algorithm for :math:`\Psi_p`.  Doing so requires the functional
form of both :math:`\beta` and :math:`c_i`:

.. math::
   \beta(\Psi_p) = max(0, min(1, \frac{\Psi_p - \Psi_{ft}}{\Psi_{wp} - \Psi_{ft}} ))

where :math:`\Psi_{ft}` is the water potential at full turgor, e.g. the point
at which stomates are completely open, and :math:`\Psi_{wp}` is the water
potential at the wilting point, e.g. the point at which stomates are completely
closed.  Both are properties of the plant.

.. math::
   c_i = c0 \rho_i

where :math:`c0` is a maximal conductance, here assumed to be constant as
prescribed by Verhoef and Egea (2014), and :math:`\rho_i` is the rooting
fraction in cell :math:`i`.

Additionally, since :math:`k_r` is upwinded, we need a plant relative
permeability describing what value is used when the potential in the plant is
higher than that in the soil.  If this is set to 0, no hydraulic
redistribution, e.g. flow from plant to soil, is allowed.  If it is set to 1,
hydraulic redistribution is the maximal flow rate, e.g. no regulation by the
plant is assumed.


.. _transpiration-distribution-relperm-evaluator-spec:
.. admonition:: transpiration-distribution-relperm-evaluator-spec

   * `"total maximal conductance [mol m^-2 s^-1 MPa^-1]`" ``[double]`` **10**
     :math:`c0` above, the total maximal conductance
   * `"plant relative conductance [-]`" ``[double]`` **0** Relative conductance
     of the plant -- 0 indicates no hydraulic redistribution, 1 indicates
     maximal redistribution.
   * `"tolerance`" ``[double]`` **1.e-12** Tolerance of the root-finding
     algorithm, which is a mixed absolute and relative tolerance.  Note that
     the default is likely fine for most problems.
   * `"maximum number of iterations`" ``[int]`` **100** Maximum number of
     iterations allowed for root finding.  Typically this is solved in < 10, so
     100 is quite safe.

   KEYS:

   - `"soil water potential`" **DOMAIN-capillary_pressure_gas_liq**
   - `"soil relative permeability`" **DOMAIN-relative_permeability**
   - `"rooting depth fraction`" **DOMAIN-rooting_depth_fraction**
   - `"potential transpiration`" **DOMAIN_SURF-potential_transpiration**
   - `"cell volume`" **DOMAIN-cell_volume**
   - `"surface cell volume`" **DOMAIN_SURF-cell_volume**
   - `"depth`" **DOMAIN-depth**

   Note this also uses `"water potential at fully closed stomata [Pa]`" and
   `"water potential at fully open stomata [Pa]`" from the land cover.

*/

#pragma once

#include "Factory.hh"
#include "BlockVector_decl.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {

namespace SurfaceBalance {
namespace Relations {

struct SoilPlantFluxFunctor {
  using cView_type = BlockVector<double>::cView_type;
  using View_type = BlockVector<double>::View_type;

  KOKKOS_INLINE_FUNCTION
  SoilPlantFluxFunctor(AmanziMesh::Entity_ID sc,
                       const AmanziMesh::MeshCache::cEntity_ID_View& cells_of_col,
                       const LandCover& lc,
                       const cView_type& soil_pc,
                       const cView_type& soil_kr,
                       const cView_type& f_root,
                       const cView_type& pet,
                       const cView_type& rho,
                       const cView_type& nliq,
                       const cView_type& visc,
                       const cView_type& cv,
                       const cView_type& sa,
                       double K,
                       double krp,
                       double g);

  // error function used for rootfinder
  KOKKOS_INLINE_FUNCTION
  double operator()(double plant_pc) const;

  // right hand side
  KOKKOS_INLINE_FUNCTION
  double computeSoilPlantFlux(double root_pc, AmanziMesh::Entity_ID c) const;

  KOKKOS_INLINE_FUNCTION
  double computeSoilPlantFluxes(double root_pc, const View_type* trans = nullptr) const;

  // beta in left hand side
  KOKKOS_INLINE_FUNCTION
  double computeTranspirationReductionFunction(double plant_pc) const;

  LandCover lc;
  cView_type soil_pc;
  cView_type soil_kr;
  cView_type f_root;
  cView_type pet;
  cView_type rho;
  cView_type nliq;
  cView_type visc;
  cView_type cv;
  cView_type sa;
  double K, krp, g;

  AmanziMesh::Entity_ID sc;
  AmanziMesh::MeshCache::cEntity_ID_View cells_of_col;
};


class TranspirationDistributionRelPermEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit TranspirationDistributionRelPermEvaluator(
    const Teuchos::RCP<Teuchos::ParameterList>& plist);
  TranspirationDistributionRelPermEvaluator(
    const TranspirationDistributionRelPermEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    // calculate of derivatives of this is a tricky thing to do, with
    // non-cell-local terms due to rescaling.  Just turn off derivatives
    // instead.
    return false;
  }

  std::string getType() const override { return eval_type; }

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  static const std::string eval_type;

  // need a custom EnsureCompatibility as some vectors cross meshes.
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  //  need a custom EnsureCompatibility as my_keys have different structure
  virtual void EnsureCompatibility_Structure_(State& S) override;

  void InitializeFromPlist_();
  bool TranspirationPeriod_(double time, double leaf_on_doy, double leaf_off_doy);

 protected:
  Key domain_surf_;
  Key domain_sub_;

  Key soil_pc_key_;
  Key soil_kr_key_;
  Key f_root_key_;
  Key potential_trans_key_;
  Key rho_key_;
  Key nliq_key_;
  Key visc_key_;
  Key sa_key_;
  Key cv_key_;

  double K_;
  double krp_;
  double rho_;

  double tol_;
  int nits_;

  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<Evaluator, TranspirationDistributionRelPermEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
