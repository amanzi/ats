/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

This evaluates an EOS as a function of temperature, pressure, and/or mole
ratio/concentration.

.. _evaluator-eos-spec:
.. admonition:: evaluator-eos-spec

   * `"EOS parameters`" ``[eos-typedinline-spec]``

   KEYS:

   - `"molar density`" Note that molar density is [mol H2O m^-3].  This is a
     _component_ quantity.
   - `"mass density`" Note that mass density is [kg m^-3].  This is a _phase_
     quantity.

   DEPENDENCIES:

   - `"temperature`"
   - `"pressure`"
   - `"mole fraction`"

*/


#ifndef AMANZI_RELATIONS_EOS_EVALUATOR_HH_
#define AMANZI_RELATIONS_EOS_EVALUATOR_HH_

#include "eos.hh"
#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

class EOSEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  enum EOSMode { EOS_MODE_MASS, EOS_MODE_MOLAR, EOS_MODE_BOTH };

  // constructor format for all derived classes
  explicit EOSEvaluator(Teuchos::ParameterList& plist);
  EOSEvaluator(const EOSEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<EOS> get_EOS() { return eos_; }

 protected:
  // ensures a given structure of all of my_keys
  virtual void EnsureCompatibility_Structure_(State& S) override
  {
    EnsureCompatibility_StructureSame_(S);
  }


  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  void ParsePlistKeys_();
  void ParsePlistTemp_();
  void ParsePlistPres_();
  void ParsePlistConc_();

 protected:
  // the actual model
  Teuchos::RCP<EOS> eos_;

  KeyTag temp_key_, pres_key_, conc_key_;
  EOSMode mode_;
  bool updated_once_;

 private:
  static Utils::RegisteredFactory<Evaluator, EOSEvaluator> fac_;
};

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi

#endif
