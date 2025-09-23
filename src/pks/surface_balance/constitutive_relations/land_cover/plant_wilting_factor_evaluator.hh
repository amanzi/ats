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

   \beta =  \frac{p_{closed} - p}{p_{closed} - p_{open}}

where p is the capillary pressure or water potential, and closed
and open indicate the values at which stomates are fully open or fully
closed (the wilting point).

Note the challenges of using this model with arbitrary van Genuchten WRMs.  See
`Verhoef & Egea, Ag. & Forest Meteorology, 2014
<https://doi.org/10.1016/j.agrformet.2014.02.009>`_


.. _evaluator-plant-wilting-factor-spec:
.. admonition:: evaluator-plant-wilting-factor-spec

   KEYS:

   - `"capillary pressure`" **DOMAIN-capillary_pressure_gas_liq**

.. note::

   This evaluator also uses the :ref:`Land Cover` types.  From that struct, it
   requires the value of the following parameters:

   - `"capillary pressure at fully closed stomata [Pa]`"
   - `"capillary pressure at fully open stomata [Pa]`"


*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace ATS_Physics {
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
} // namespace ATS_Physics
} // namespace Amanzi
