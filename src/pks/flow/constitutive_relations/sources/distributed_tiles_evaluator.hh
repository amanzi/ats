/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

//! Evaluates water/solute source which represent effect of distributed subsurface tiles.
/*!

Distributed, subsurface sources due to tile drains.

.. _distributed-tiles-spec:
.. admonition:: distributed-tiles-spec

   * `"number of ditches`" ``[int]`` Number of ditches, corresponding to the number of unique IDs.
   * `"tile permeability [m^2]`" ``[double]`` Permeability of the tile/pipe connecting soil to ditch.
   * `"number of components`" ``[int]`` **1** Number of components in the source/sink pair.
   * `"entering pressure [Pa]`" ``[double]`` **101325** Pressure required to enter the tile drain.

   KEYS:
   - `"accumulated source`" **DOMAIN-accumulated_source** Source to the ditch from the tile.
   - `"distributed source`" **DOMAIN-distributed_source** Source to the soil from the tile (likely negative, sink).
   - `"catchment ID`" **DOMAIN-catchments_id** ID indicating which ditch a given cell drains to.
   - `"pressure`"
   - `"molar density liquid`"
   - `"factor field`" **NO FACTOR** Scalar factor?

*/

#ifndef AMANZI_FLOW_RELATIONS_DISTTILES_EVALUATOR_HH_
#define AMANZI_FLOW_RELATIONS_DISTTILES_EVALUATOR_HH_

// TPLs
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"


#include "Factory.hh"
#include "EvaluatorSecondary.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace Relations {

class DistributedTilesRateEvaluator : public EvaluatorSecondary {
 public:
  explicit DistributedTilesRateEvaluator(Teuchos::ParameterList& plist);
  DistributedTilesRateEvaluator(const DistributedTilesRateEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new DistributedTilesRateEvaluator(*this));
  }

  virtual void EnsureCompatibility(State& S) override;

  // derivatives aren't implemented here
  bool IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    return false;
  }

 protected:
  // Required methods from EvaluatorSecondary
  virtual void Update_(State& S) override;
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override {};

 protected:
  Key domain_;

  // my keys
  Key dist_sources_key_;
  Key acc_sources_key_;

  // my dependencies
  Key catch_id_key_;
  Key pres_key_;
  Key mol_dens_key_, mass_dens_key_, visc_key_;
  Key factor_key_;

  double p_enter_;
  double ka_, kb_, d_, L_, th_;
  int num_ditches_;
  int num_components_;

 private:
  static Utils::RegisteredFactory<Evaluator, DistributedTilesRateEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
#endif
