/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the unfrozen fraction of water.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_MODEL_
#define AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_MODEL_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {

class UnfrozenFractionModel {
 public:
  UnfrozenFractionModel(Teuchos::ParameterList& list);

  double UnfrozenFraction(double temp) const;
  double DUnfrozenFractionDT(double temp) const;

 protected:
  double halfwidth_;
  double T0_;
  double min_uf_;
  const double pi_;

  Teuchos::ParameterList plist_;
};

} // namespace Flow
} // namespace Amanzi

#endif
