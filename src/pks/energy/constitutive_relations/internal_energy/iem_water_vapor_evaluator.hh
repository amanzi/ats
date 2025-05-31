/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The IEM Evaluator simply calls the IEM with the correct arguments.

*/

/*!

Computes (specific) internal energy of as a function of temperature and molar
fraction of water vapor in the gaseous phase.

`"evaluator type`" = `"iem water vapor`"

.. _evaluator-iem-water-vapor-spec:
.. admonition:: evaluator-iem-water-vapor-spec

   * `"IEM parameters`" ``[iem-water-vapor-spec]``

   KEYS:

   - `"temperature`"
   - `"vapor molar fraction`"

*/

#ifndef AMANZI_ENERGY_RELATIONS_IEM_WATER_VAPOR_EVALUATOR_
#define AMANZI_ENERGY_RELATIONS_IEM_WATER_VAPOR_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "iem_water_vapor.hh"

namespace Amanzi {
namespace Energy {

class IEMWaterVaporEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit IEMWaterVaporEvaluator(Teuchos::ParameterList& plist);
  IEMWaterVaporEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<IEMWaterVapor>& iem);
  IEMWaterVaporEvaluator(const IEMWaterVaporEvaluator& other) = default;

  Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<IEMWaterVapor> get_IEM() { return iem_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  void InitializeFromPlist_();

  Key temp_key_;
  Key mol_frac_key_;
  Teuchos::RCP<IEMWaterVapor> iem_;

 private:
  static Utils::RegisteredFactory<Evaluator, IEMWaterVaporEvaluator> factory_;
};

} // namespace Energy
} // namespace Amanzi

#endif
