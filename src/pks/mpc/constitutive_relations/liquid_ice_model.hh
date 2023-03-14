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

#ifndef AMANZI_LIQUID_ICE_MODEL_
#define AMANZI_LIQUID_ICE_MODEL_

#include "Teuchos_ParameterList.hpp"

#include "Tensor.hh"
#include "Point.hh"

#include "wrm_partition.hh"
#include "compressible_porosity_model_partition.hh"
#include "compressible_porosity_leijnse_model_partition.hh"
#include "ewc_model_base.hh"

namespace Amanzi {

namespace Flow {
class PCIceWater;
class PCLiqAtm;
} // namespace Flow

namespace Energy {
class IEM;
}
namespace Relations {
class EOS;
}

class LiquidIceModel : public EWCModelBase {
 public:
  LiquidIceModel() {}
  virtual ~LiquidIceModel(){};

  virtual void InitializeModel(const Teuchos::Ptr<State>& S,
                               const Tag& tag,
                               Teuchos::ParameterList& plist) override;
  virtual void UpdateModel(const Teuchos::Ptr<State>& S, int c) override;
  virtual bool Freezing(double T, double p) override;
  virtual int
  EvaluateSaturations(double T, double p, double& s_gas, double& s_liq, double& s_ice) override;

 protected:
  bool IsSetUp_();

  int EvaluateEnergyAndWaterContent_(double T, double p, AmanziGeometry::Point& result) override;

 protected:
  Teuchos::RCP<Flow::WRMPermafrostModelPartition> wrms_;
  Teuchos::RCP<Flow::WRMPermafrostModel> wrm_;
  Teuchos::RCP<Relations::EOS> liquid_eos_;
  Teuchos::RCP<Relations::EOS> gas_eos_;
  Teuchos::RCP<Relations::EOS> ice_eos_;
  Teuchos::RCP<Flow::PCIceWater> pc_i_;
  Teuchos::RCP<Flow::PCLiqAtm> pc_l_;
  Teuchos::RCP<Energy::IEM> liquid_iem_;
  Teuchos::RCP<Energy::IEM> ice_iem_;
  Teuchos::RCP<Energy::IEM> rock_iem_;
  Teuchos::RCP<Flow::CompressiblePorosityModelPartition> poro_models_;
  Teuchos::RCP<Flow::CompressiblePorosityModel> poro_model_;

  Teuchos::RCP<Flow::CompressiblePorosityLeijnseModelPartition> poro_leij_models_;
  Teuchos::RCP<Flow::CompressiblePorosityLeijnseModel> poro_leij_model_;

  double p_atm_;
  double poro_;
  double rho_rock_;
  bool poro_leij_;
  Key domain;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Tag tag_;

  bool use_pc_ice_;
};

} // namespace Amanzi


#endif
