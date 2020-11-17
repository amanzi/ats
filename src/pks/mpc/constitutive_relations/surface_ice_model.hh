/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Ugly hackjob to enable direct evaluation of the full model, on a single
  WRM/region.  This is bypassing much of the "niceness" of the framework, but
  seems necessary for solving a cell-wise correction equation.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_SURFACE_ICE_MODEL_
#define AMANZI_SURFACE_ICE_MODEL_

#include "Teuchos_ParameterList.hpp"

#include "Tensor.hh"
#include "Point.hh"

#include "ewc_model_base.hh"

namespace Amanzi {

namespace Flow {
class IcyHeightModel;
class UnfrozenFractionModel;
}

namespace Energy {
class IEM;
}

namespace Relations {
class EOS;
}

class SurfaceIceModel : public EWCModelBase {

 public:
  SurfaceIceModel() {}
  virtual void InitializeModel(const Teuchos::Ptr<State>& S, Teuchos::ParameterList& plist);
  virtual void UpdateModel(const Teuchos::Ptr<State>& S, int c);

  virtual bool Freezing(double T, double p) { return T < 273.15; }
  virtual int EvaluateSaturations(double T, double p, double& s_gas, double& s_liq, double& s_ice) {
    AMANZI_ASSERT(0);
    return 1;
  }

 protected:
  bool IsSetUp_();

  int EvaluateEnergyAndWaterContent_(double T, double p,
          AmanziGeometry::Point& result);

 protected:
  Teuchos::RCP<Flow::IcyHeightModel> pd_;
  Teuchos::RCP<Flow::UnfrozenFractionModel> uf_;
  Teuchos::RCP<Relations::EOS> liquid_eos_;
  Teuchos::RCP<Relations::EOS> ice_eos_;
  Teuchos::RCP<Energy::IEM> liquid_iem_;
  Teuchos::RCP<Energy::IEM> ice_iem_;

  double p_atm_;
  double gz_;
  double M_;

  Key ice_dens_key_;
};


}


#endif
