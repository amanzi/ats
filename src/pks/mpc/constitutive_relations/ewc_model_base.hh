/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

EWCModelBase provides some of the functionality of EWCModel for inverse
evaluating.

------------------------------------------------------------------------- */

#ifndef AMANZI_EWC_MODEL_BASE_HH_
#define AMANZI_EWC_MODEL_BASE_HH_

#include "Tensor.hh"
#include "Point.hh"

#include "ewc_model.hh"

namespace Amanzi {

class EWCModelBase : public EWCModel {
 public:
  EWCModelBase() {}
  virtual ~EWCModelBase() = default;

  virtual int Evaluate(double T, double p, double& energy, double& wc) override;
  virtual int InverseEvaluate(double energy,
                              double wc,
                              double& T,
                              double& p,
                              bool verbose = false) override;
  virtual int InverseEvaluateEnergy(double energy, double p, double& T) override;

 protected:
  virtual int EvaluateEnergyAndWaterContent_(double T, double p, AmanziGeometry::Point& result) = 0;

  int EvaluateEnergyAndWaterContentAndJacobian_(double T,
                                                double p,
                                                AmanziGeometry::Point& result,
                                                WhetStone::Tensor& jac);

  int EvaluateEnergyAndWaterContentAndJacobian_FD_(double T,
                                                   double p,
                                                   AmanziGeometry::Point& result,
                                                   WhetStone::Tensor& jac);
};

} // namespace Amanzi

#endif
