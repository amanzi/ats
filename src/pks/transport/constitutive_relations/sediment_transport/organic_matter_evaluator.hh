/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*!

The organic evaluator gets the organic matter deposition rates.

..math::

   Q_{db} = Q_{db_0} \frac{B}{B_{max}}

where :math:`Q_{db_0}` is a typical deposition rate, empirical coefficient with
the dimensions of [L/T], which is derived empirically from field measurements,
:math:`B` is the current biomass, :math:`B_{max}` is the maximum value of the
biomass.

`"evaluator type`" = `"organic matter rate`"

.. _evaluator-organic-matter-rate-spec:
.. admonition:: evaluator-organic-matter-rate-spec

   * `"empirical coefficient`" ``[double]`` **9.5129e-11**
   * `"maximum biomass`" ``[double]`` **2000**

   DEPENDENCIES:

   - `"biomass`" **DOMAIN-biomass**

*/

#pragma once

// TPLs
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

class OrganicMatterRateEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit OrganicMatterRateEvaluator(Teuchos::ParameterList& plist);
  OrganicMatterRateEvaluator(const OrganicMatterRateEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;


  double Bmax_;
  double Q_db0_;
  double Q_on_Bmax_;
  Key biomass_key_;

  static Utils::RegisteredFactory<Evaluator, OrganicMatterRateEvaluator> factory_;
};

} // namespace Amanzi
