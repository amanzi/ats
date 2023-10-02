/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Plant wilting factor provides a moisture availability-based limiter on transpiration.
/*!

Also known as Beta, or the water availability factor, or the plant wilting
factor, or the transpiration reduction function.

.. math::
   Beta =  (p_closed - p) / (p_closed - p_open)

where p is the capillary pressure or water potential, and closed
and open indicate the values at which stomates are fully open or fully
closed (the wilting point).

Note this makes use of LandCover objects for water potential of fully open and
fully closed stomata.

Note the challenges of using this model with arbitrary van Genuchten WRMs.  See
Verhoef & Egea, Ag. & Forest Meteorology, 2014
https://doi.org/10.1016/j.agrformet.2014.02.009


.. _plant-wilting-factor-evaluator-spec:
.. admonition:: plant-wilting-factor-evaluator-spec

   KEYS:

   - `"capillary pressure`" **DOMAIN-capillary_pressure_gas_liq**


*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class PlantWiltingFactorModel;

class PlantWiltingFactorEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit PlantWiltingFactorEvaluator(Teuchos::ParameterList& plist);
  PlantWiltingFactorEvaluator(const PlantWiltingFactorEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

 protected:
  Key pc_key_;
  Key domain_surf_;
  Key domain_sub_;

  LandCoverMap land_cover_;
  std::map<std::string, Teuchos::RCP<PlantWiltingFactorModel>> models_;

 private:
  static Utils::RegisteredFactory<Evaluator, PlantWiltingFactorEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
