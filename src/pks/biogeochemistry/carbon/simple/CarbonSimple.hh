/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Process kernel for energy equation for Richard's flow.
------------------------------------------------------------------------- */

#ifndef PKS_CARBON_SIMPLE_HH_
#define PKS_CARBON_SIMPLE_HH_

#include "PK_Factory.hh"
#include "pk_physical_explicit_default.hh"
#include "PK.hh"

namespace Amanzi {
namespace BGC {

class CarbonSimple : public PK_Physical_Explicit_Default {
 public:
  CarbonSimple(Teuchos::ParameterList& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& glist,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& solution);

  // EnergyBase is a PK
  // -- Setup data
  virtual void Setup() override;

  // -- Initialize owned (dependent) variables.
  // default ok?
  // virtual void initialize(const State& S);

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override {};

  // -- Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Tag& tag) override;

  // EnergyBase is a BDFFnBase
  // computes the non-linear functional f = f(t,u,udot)
  virtual void
  FunctionalTimeDerivative(const double t, const TreeVector& u, TreeVector& f) override;

 protected:
  virtual void ApplyDiffusion_(const Teuchos::Ptr<CompositeVector>& g);
  virtual void AddSources_(const Teuchos::Ptr<CompositeVector>& g);
  virtual void AddDecomposition_(const Teuchos::Ptr<CompositeVector>& g);


 protected:
  int npools_;

  Key cell_vol_key_;

  bool is_diffusion_;
  Key div_diff_flux_key_;

  bool is_source_;
  Key source_key_;

  bool is_decomp_;
  Key decomp_key_;

 private:
  // factory registration
  static RegisteredPKFactory<CarbonSimple> reg_;
};


} // namespace BGC
} // namespace Amanzi


#endif
