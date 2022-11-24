/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow as a function of ponded
  depth using Manning's model. The denominator in the model is evaluated separately.

*/

#ifndef AMANZI_FLOWRELATIONS_SPLIT_DENOMINATOR_CONDUCTIVITY_MODEL_
#define AMANZI_FLOWRELATIONS_SPLIT_DENOMINATOR_CONDUCTIVITY_MODEL_

#include "Teuchos_ParameterList.hpp"
#include "overland_conductivity_model.hh"

namespace Amanzi {
namespace Flow {

class SplitDenominatorConductivityModel : public OverlandConductivityModel {
 public:
  explicit SplitDenominatorConductivityModel(Teuchos::ParameterList& plist);

  virtual double Conductivity(double depth, double slope, double coef);

  virtual double DConductivityDDepth(double depth, double slope, double coef);

 protected:
  Teuchos::ParameterList plist_;

  double manning_exp_;
};

} // namespace Flow
} // namespace Amanzi

#endif
