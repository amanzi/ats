/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Bo Gao (gaob@ornl.gov)
*/

/*
  Evaluator for determining water level, a combined variable of
  water table and ponded depth.
*/

/*!
When surface pressure is over 101325 Pa, water level = surface ponded depth;
otherwise water level = subsurface water table depth * (-1)

Computes ponded depth from surface water pressure.

.. math::
   h = \frac{H(p - p_{atm})}{\rho g}

where :math:`H` is the Heaviside function.

Computes water table by looping through cells in the column until saturation of gas is 0.

`"evaluator type`" = `"ponded depth`"

.. _evaluator-height-spec:
.. admonition:: evaluator-height-spec

   KEYS:

   - `"mass density`"
   - `"pressure`"
   - `"gas saturation`"
   - `"subsurface cell volume`"
   - `"surface cell volume`"
*/

#ifndef AMANZI_FLOW_RELATIONS_WATER_LEVEL_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_WATER_LEVEL_EVALUATOR_

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class HeightModel;

class WaterLevelEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit WaterLevelEvaluator(Teuchos::ParameterList& plist);
  WaterLevelEvaluator(const WaterLevelEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<HeightModel> get_Model() { return model_; }


 protected:
  // Needs a special EnsureCompatibility to get around trying to find face
  // values and derivatives of face values.
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key dens_key_;
  Key pres_key_;
  Key gravity_key_;
  Key patm_key_;
  Key sat_key_;
  Key cv_key_;
  Key cv_surf_key_;

  Teuchos::RCP<HeightModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, WaterLevelEvaluator> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
