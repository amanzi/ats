/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Bo Gao (gaob@ornl.gov)
*/

/*!
  Sakagucki-Zeng soil resistance model refered to Sakaguchi and Zeng (2009).
  Note that `"dessicated_zone_thickness`" is given by soil types.
  If it is not declared in `"WRM paramters`" through `"model parameters`"
  under state, default value 0.1 m is used for all soil types.
*/

#pragma once

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "StateDefs.hh"
#include "State.hh"

#include "EvaluatorModelCVByMaterial.hh"
#include "EvaluatorModelCV.hh"
#include "soil_resistance_sakagucki_zeng_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

template <class Evaluator_type>
class SoilResistanceSakaguckiZengEvaluator_ : public Evaluator_type {
 public:
  SoilResistanceSakaguckiZengEvaluator_(const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : Evaluator_type(plist)
  {}
  SoilResistanceSakaguckiZengEvaluator_(const SoilResistanceSakaguckiZengEvaluator_& other) =
    default;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new SoilResistanceSakaguckiZengEvaluator_(*this));
  }

 protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override
  {
    Teuchos::RCP<const AmanziMesh::Mesh> parent =
      S.GetMesh(Keys::getDomain(my_keys_.front().first))->getParentMesh();
    AMANZI_ASSERT(parent != Teuchos::null);

    for (const auto& dep : dependencies_) {
      S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second)
        .SetMesh(parent)
        ->SetGhosted()
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    }
  }

  using Evaluator_type::eval_type;
  using Evaluator_type::my_keys_;
  using Evaluator_type::dependencies_;

 private:
  static Utils::RegisteredFactory<Evaluator, SoilResistanceSakaguckiZengEvaluator_> reg_;
};


using SoilResistanceSakaguckiZengEvaluator =
  SoilResistanceSakaguckiZengEvaluator_<EvaluatorModelCV<SoilResistanceSakaguckiZengModel>>;

using SoilResistanceSakaguckiZengEvaluatorByMaterial = SoilResistanceSakaguckiZengEvaluator_<
  EvaluatorModelCVByMaterial<SoilResistanceSakaguckiZengModel>>;

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
