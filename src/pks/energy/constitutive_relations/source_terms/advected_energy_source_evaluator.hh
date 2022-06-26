/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Source term evaluator for enthalpy of mass source.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ENERGY_RELATIONS_ADVECTED_ENERGY_SOURCE_EVALUATOR_
#define AMANZI_ENERGY_RELATIONS_ADVECTED_ENERGY_SOURCE_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Energy {

class AdvectedEnergySourceEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  // constructor format for all derived classes
  explicit
  AdvectedEnergySourceEvaluator(Teuchos::ParameterList& plist);
  AdvectedEnergySourceEvaluator(const AdvectedEnergySourceEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& results) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& results) override;

 protected:
  void InitializeFromPlist_();

  Key internal_enthalpy_key_;
  Key external_enthalpy_key_;
  Key water_source_key_;
  Key internal_density_key_;
  Key external_density_key_;
  Key conducted_source_key_;
  Key cell_vol_key_;

  bool include_conduction_;
  enum SourceUnits {
    SOURCE_UNITS_METERS_PER_SECOND,
    SOURCE_UNITS_MOLS_PER_SECOND,
    SOURCE_UNITS_MOLS_PER_SECOND_PER_METERSD
  };

  SourceUnits source_units_;

 private:
  static Utils::RegisteredFactory<Evaluator,AdvectedEnergySourceEvaluator> factory_;

};

} //namespace
} //namespace

#endif
