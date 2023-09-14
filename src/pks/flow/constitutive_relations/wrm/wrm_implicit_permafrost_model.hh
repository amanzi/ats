/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Painter's original, implicitly defined permafrost model.
/*!

.. _wrm-implicit-permafrost-spec:
.. admonition:: wrm-implicit-permafrost-spec

    * `"converged tolerance`" ``[double]`` **1.e-12** Convergence tolerance of the implicit solve.
    * `"max iterations`" ``[int]`` **100** Maximum allowable iterations of the implicit solve.
    * `"solver algorithm [bisection/toms]`" ``[string]`` **bisection** Use bisection or the TOMS algorithm from boost.

*/

#ifndef AMANZI_FLOWRELATIONS_WRM_IMPLICIT_PERMAFROST_MODEL_
#define AMANZI_FLOWRELATIONS_WRM_IMPLICIT_PERMAFROST_MODEL_

#include "boost/cstdint.hpp"
#include "boost/math/tools/roots.hpp"
#include "boost/cstdint.hpp"

#include "wrm_permafrost_model.hh"
#include "wrm_permafrost_factory.hh"

namespace Amanzi {
namespace Flow {

class WRM;

class WRMImplicitPermafrostModel : public WRMPermafrostModel {
 public:
  explicit WRMImplicitPermafrostModel(Teuchos::ParameterList& plist);

  // required methods from the base class
  // sats[0] = sg, sats[1] = sl, sats[2] = si
  virtual bool freezing(double T, double pc_liq, double pc_ice) { return T < 273.15; }
  virtual void saturations(double pc_liq, double pc_ice, double (&sats)[3]);
  virtual void dsaturations_dpc_liq(double pc_liq, double pc_ice, double (&dsats)[3]);
  virtual void dsaturations_dpc_ice(double pc_liq, double pc_ice, double (&dsats)[3]);

 protected:
  // calculation if unfrozen
  bool sats_unfrozen_(double pc_liq, double pc_ice, double (&sats)[3]);
  bool dsats_dpc_liq_unfrozen_(double pc_liq, double pc_ice, double (&dsats)[3]);
  bool dsats_dpc_ice_unfrozen_(double pc_liq, double pc_ice, double (&dsats)[3]);

  // calculation if saturated
  bool sats_saturated_(double pc_liq, double pc_ice, double (&sats)[3]);
  bool dsats_dpc_liq_saturated_(double pc_liq, double pc_ice, double (&dsats)[3]);
  bool dsats_dpc_ice_saturated_(double pc_liq, double pc_ice, double (&dsats)[3]);

  // calculation if unfrozen and saturated
  bool sats_frozen_unsaturated_(double pc_liq, double pc_ice, double (&sats)[3]);
  bool dsats_dpc_liq_frozen_unsaturated_(double pc_liq, double pc_ice, double (&dsats)[3]);
  bool dsats_dpc_ice_frozen_unsaturated_(double pc_liq, double pc_ice, double (&dsats)[3]);

  double si_frozen_unsaturated_(double pc_liq, double pc_ice);
  double dsi_dpc_liq_frozen_unsaturated_(double pc_liq, double pc_ice, double si);
  double dsi_dpc_ice_frozen_unsaturated_(double pc_liq, double pc_ice, double si);

  double si_frozen_unsaturated_nospline_(double pc_liq, double pc_ice, bool throw_ok = false);
  double dsi_dpc_liq_frozen_unsaturated_nospline_(double pc_liq, double pc_ice, double si);
  double dsi_dpc_ice_frozen_unsaturated_nospline_(double pc_liq, double pc_ice, double si);

  bool DetermineSplineCutoff_(double pc_liq, double pc_ice, double& cutoff, double& si);
  bool FitSpline_(double pc_ice, double cutoff, double si_cutoff, double (&coefs)[4]);


 protected:
  double eps_;
  boost::uintmax_t max_it_;
  double deriv_regularization_;
  std::string solver_;

 private:
  // Functor for ice saturation, gets used within a root-finding algorithm
  class SatIceFunctor_ {
   public:
    SatIceFunctor_(double pc_liq, double pc_ice, const Teuchos::RCP<WRM>& wrm)
      : pc_liq_(pc_liq), pc_ice_(pc_ice), wrm_(wrm)
    {}

    double operator()(double si)
    {
      double tmp = (1.0 - si) * wrm_->saturation(pc_liq_);
      double result = tmp - wrm_->saturation(pc_ice_ + wrm_->capillaryPressure(tmp + si));
      return result;
    }

   private:
    double pc_liq_;
    double pc_ice_;
    Teuchos::RCP<WRM> wrm_;
  };

  // Convergence criteria for root-finding
  struct Tol_ {
    Tol_(double eps) : eps_(eps) {}
    bool operator()(const double& a, const double& b) const { return std::abs(a - b) <= eps_; }
    double eps_;
  };

 private:
  // factory registration
  static Utils::RegisteredFactory<WRMPermafrostModel, WRMImplicitPermafrostModel> factory_;
};


} // namespace Flow
} // namespace Amanzi

#endif
