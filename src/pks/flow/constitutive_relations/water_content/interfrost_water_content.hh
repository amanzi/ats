/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/* -----------------------------------------------------------------------------
ATS

Evaluator for water content.

INTERFROST's comparison uses a very odd compressibility term that doesn't
quite fit into either compressible porosity or into a compressible density, so
it needs a special evaluator.

----------------------------------------------------------------------------- */


#ifndef AMANZI_INTERFROST_WATER_CONTENT_HH_
#define AMANZI_INTERFROST_WATER_CONTENT_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace ATS_Physics {
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
  double beta_;

 private:
  static Utils::RegisteredFactory<Evaluator, InterfrostWaterContent> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
