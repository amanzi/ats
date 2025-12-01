/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Interfrost water content: liquid with compressibility, and ice.
/*!

.. math::
  \Theta = (n_l s_l (1 + \beta (pressure - 101325) ) + n_i s_i) \phi V

`"evaluator type`" = `"interfrost water content`"

.. _evaluator-interfrost-water-content-spec:
.. admonition:: evaluator-interfrost-water-content-spec

   DEPENDENCIES:

   - `"porosity`"
   - `"molar density liquid`"
   - `"saturation liquid`"
   - `"molar density ice`"
   - `"saturation ice`"
   - `"cell volume`"
   - `"pressure`"

*/


#ifndef AMANZI_INTERFROST_WATER_CONTENT_HH_
#define AMANZI_INTERFROST_WATER_CONTENT_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class InterfrostWaterContent : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit InterfrostWaterContent(Teuchos::ParameterList& wc_plist);
  InterfrostWaterContent(const InterfrostWaterContent& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key phi_key_;
  Key sl_key_;
  Key nl_key_;
  Key si_key_;
  Key ni_key_;
  Key cv_key_;
  Key pres_key_;

  double beta_;

 private:
  static Utils::RegisteredFactory<Evaluator, InterfrostWaterContent> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

#endif
