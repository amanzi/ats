/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Distributes potential transpiration to the rooting zone using a rel perm weighting method.
/*!

The transpiration distribution evaluator looks to take a potential
evapotranspiration and distribute it across the vertical column based on water
availability and rooting depths.  It also potentially limits the transpiration
to avoid taking water where it is not available (thereby crashing the code).

This is a new model developed by Painter and Coon that uses a competing rates
approach.  Given a capillary pressure at the plant collar, :math:`\Psi_p`,

.. math::
   Q_{T}(\Psi_p) = \sum_{i} Q_{i} = \frac{n}{\nu}} K r_i k_r(\Psi_i, \Psi_p + \rho g z_i)[\Psi_i - (\Psi_p + \rho g z_i)]

where:

- :math:`Q_{T}` is the total actual transpiration formed as a sum over grid
  cells :math:`i`,
- :math:`n` and :math:`\nu` are the liquid water density and
  viscosity, respectively,
- :math:`K` is the permeability of the root-soil interface,
- :math:`r_i` is the fraction of the plant's roots in cell :math:`i`,
- k_r is the upwinded relative permeability, which may be either the soil
  relative peremability or the root relative permeability (see below),
- :math:`\Psi_p` and :math:`\Psi_i` the water potential (capillary pressure) at
  the plant collar and the soil cell :math:`i`, respectively,
- :math:`\rho` is the density of water and :math:`g` the gravitational force
  (so that :math:`\Psi_p + \rho g z_i` is the pressure in the root in soil cell
  :math:`i`),
- and :math:`Q_i` is the transpiration taken from each grid cell.

Given a potential transpiration :math:`Q_{pot T}`, which is provided by the
user through another model (e.g. Priestley-Taylor or comparable), we first
compare to a maximum plant capillary pressure, :math:`\Psi_{max}`, which is a
parameter of the plant (e.g. the wilting point).  If :math:`Q_{T}(\Psi_p =
\Psi_{max})` results in a value less than the potential transpiration, this is
used directly and :math:`Q_i` are computed directly.  If this results in a
tranpiration greater than the potential, then we set

.. math::
   Q_{T} = Q_{pot T}

and the above equations are solved using a scalar root-finding algorithm for
:math:`\Psi_p`.

Additionally, :math:`k_r` is upwinded (e.g. :math:`k_r` of the soil is used if
capillary pressure in the root is smaller than that of the soil; otherwise the
plant's relative permeability is used).  The plant relative permeability is
chosen by the user; if it is set to 0, no hydraulic redistribution, e.g. flow
from plant to soil, is allowed.  If it is set to 1, hydraulic redistribution is
the maximal flow rate, e.g. no regulation by the plant is assumed.


.. _transpiration_distribution_relperm_evaluator-spec:
.. admonition:: transpiration_distribution_relperm_evaluator-spec

   * `"plant permeability per m [m]`" ``[double]`` **1.e-12**
     :math:`K` above, the total plant permeability.
   * `"plant relative permeability [-]`" ``[double]`` **0** Relative
     permeability of the plant -- 0 indicates no hydraulic redistribution, 1
     indicates maximal redistribution.
   * `"tolerance`" ``[double]`` **1.e-12** Tolerance of the root-finding
     algorithm, which is a mixed absolute and relative tolerance.  Note that
     the default is likely fine for most problems.
   * `"maximum number of iterations`" ``[int]`` **100** Maximum number of
     iterations allowed for root finding.  Typically this is solved in < 10, so
     100 is quite safe.

   KEYS:

   - `"soil water potential`" **DOMAIN-capillary_pressure_gas_liq**
   - `"soil relative conductance`" **DOMAIN-relative_permeability**
   - `"rooting depth fraction`" **DOMAIN-rooting_depth_fraction**
   - `"density`" **DOMAIN-molar_density_liquid**
   - `"viscosity`" **DOMAIN-viscosity_liquid**
   - `"rooting depth fraction`" **DOMAIN-rooting_depth_fraction**
   - `"potential transpiration`" **DOMAIN_SURF-potential_transpiration**
   - `"cell volume`" **DOMAIN-cell_volume**
   - `"surface cell volume`" **DOMAIN_SURF-cell_volume**
   - `"depth`" **DOMAIN-depth**

   Note this also uses `"maximum xylem capillary pressure [Pa]`"from the land
   cover parameters.

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {

namespace SurfaceBalance {
namespace Relations {

struct SoilPlantFluxFunctor {
  SoilPlantFluxFunctor(AmanziMesh::Entity_ID sc,
                       const AmanziMesh::Entity_ID_View& cells_of_col,
                       const LandCover& lc,
                       const Epetra_MultiVector& soil_pc,
                       const Epetra_MultiVector& soil_kr,
                       const Epetra_MultiVector& f_root,
                       const Epetra_MultiVector& pet,
                       const Epetra_MultiVector& rho,
                       const Epetra_MultiVector& nliq,
                       const Epetra_MultiVector& visc,
                       const Epetra_MultiVector& cv,
                       const Epetra_MultiVector& sa,
                       double K,
                       double krp,
                       double g);

  // error function used for rootfinder
  double operator()(double plant_pc) const;

  // right hand side
  double computeSoilPlantFlux(double root_pc, AmanziMesh::Entity_ID c) const;
  double computeSoilPlantFluxes(double root_pc, Epetra_MultiVector* trans=nullptr) const;

  const LandCover& lc;
  const Epetra_MultiVector& soil_pc;
  const Epetra_MultiVector& soil_kr;
  const Epetra_MultiVector& f_root;
  const Epetra_MultiVector& pet;
  const Epetra_MultiVector& rho;
  const Epetra_MultiVector& nliq;
  const Epetra_MultiVector& visc;
  const Epetra_MultiVector& cv;
  const Epetra_MultiVector& sa;
  double K, krp, g;

  AmanziMesh::Entity_ID sc;
  AmanziMesh::Entity_ID_View cells_of_col;
};


class TranspirationDistributionRelPermEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit TranspirationDistributionRelPermEvaluator(Teuchos::ParameterList& plist);
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

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  // need a custom EnsureCompatibility as some vectors cross meshes.
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  //  need a custom EnsureCompatibility as my_keys have different structure
  virtual void EnsureCompatibility_Structure_(State& S) override;

  void InitializeFromPlist_();
  bool TranspirationPeriod_(double time, double leaf_on_doy, double leaf_off_doy);

 protected:
  Key domain_surf_;
  Key domain_sub_;

  Key rho_key_;
  Key nliq_key_;
  Key visc_key_;
  Key soil_pc_key_;
  Key soil_kr_key_;
  Key f_root_key_;
  Key potential_trans_key_;
  Key sa_key_;
  Key cv_key_;

  double K_;
  double krp_;

  double tol_;
  int nits_;

  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<Evaluator, TranspirationDistributionRelPermEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
