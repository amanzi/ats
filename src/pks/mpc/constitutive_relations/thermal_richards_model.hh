/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Ugly hackjob to enable direct evaluation of the full model, on a single
  WRM/region.  This is bypassing much of the "niceness" of the framework, but
  seems necessary for solving a cell-wise correction equation.

*/

#ifndef AMANZI_PERMAFROST_MODEL_
#define AMANZI_PERMAFROST_MODEL_

#include "Teuchos_ParameterList.hpp"

#include "Tensor.hh"
#include "Point.hh"

#include "ewc_model.hh"

namespace Amanzi {
namespace ATS_Physics {

namespace Flow {
namespace Flow {
class WRM;
}
}
namespace Energy {
namespace Energy {
class IEM;
class IEMWaterVapor;
}
}
namespace Relations {
class EOS;
class VaporPressureRelation;
}


class ThermalRichardsModel : public EWCModel {
 public:
  ThermalRichardsModel() {}
  virtual void InitializeModel(const Teuchos::Ptr<State>& S);
  virtual void UpdateModel(const Teuchos::Ptr<State>& S, int c);
  virtual int Evaluate(double T,
                       double p,
                       double base_poro,
                       double& energy,
                       double& wc,
                       bool verbose = false);
  virtual int InverseEvaluate(double energy, double wc, double base_poro, double& T, double& p);

 protected:
  bool IsSetUp_();

  int EvaluateEnergyAndWaterContent_(double T,
                                     double p,
                                     double base_poro,
                                     AmanziGeometry::Point& result);

  int EvaluateEnergyAndWaterContentAndJacobian_(double T,
                                                double p,
                                                double base_poro,
                                                AmanziGeometry::Point& result,
                                                WhetStone::Tensor& jac);

  int EvaluateEnergyAndWaterContentAndJacobian_FD_(double T,
                                                   double p,
                                                   double base_poro,
                                                   AmanziGeometry::Point& result,
                                                   WhetStone::Tensor& jac);

 protected:
  Teuchos::RCP<Flow::WRM> wrm_;
  Teuchos::RCP<Relations::EOS> liquid_eos_;
  Teuchos::RCP<Relations::EOS> gas_eos_;
  Teuchos::RCP<Flow::PCLiqAtm> pc_l_;
  Teuchos::RCP<Relations::VaporPressureRelation> vpr_;
  Teuchos::RCP<Energy::IEM> liquid_iem_;
  Teuchos::RCP<Energy::IEMWaterVapor> gas_iem_;
  Teuchos::RCP<Energy::IEM> rock_iem_;
  double rho_rock_;
  double p_atm_;
  double poro_;
};


#endif
