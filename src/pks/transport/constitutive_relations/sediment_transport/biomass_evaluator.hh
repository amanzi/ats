/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  This biomass model evaluates
*/

#ifndef AMANZI_BIOMASS_EVALUATOR
#define AMANZI_BIOMASS_EVALUATOR

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {

class BiomassEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit BiomassEvaluator(Teuchos::ParameterList& plist);
  BiomassEvaluator(const BiomassEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
                            const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override {
    AMANZI_ASSERT(0);
  }

  virtual void EnsureCompatibility_Structure_(State& S) {
    EnsureCompatibility_StructureSame_(S);
  }

  void InitializeFromPlist_();

 protected:

  int nspecies_, type_;
  std::vector<double> alpha_n, alpha_h, alpha_d, alpha_a;
  std::vector<double> beta_n, beta_h, beta_d, beta_a;
  std::vector<double> Bmax, zmax, zmin;

  Key biomass_key_, stem_density_key_, stem_height_key_, stem_diameter_key_, plant_area_key_;
  Key elev_key_, msl_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, BiomassEvaluator> factory_;
};

} // namespace Amanzi

#endif
