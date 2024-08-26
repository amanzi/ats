/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Weakly coupled N PKs, potentially subcycling any of them.

  * `"subcycle`" ``[Array(bool)]`` Array of the same length as sub_pks.
  * `"minimum subcycled relative dt`" ``[double]`` **1.e-5** Sets the minimum
    time step size of the subcycled PKs, as a multiple of the minimum of the
    non-subcycled PKs' timestep sizes.

  INCLUDES:
  - ``[mpc-spec]``

*/

#ifndef ATS_AMANZI_SUBCYCLED_MPC_HH_
#define ATS_AMANZI_SUBCYCLED_MPC_HH_

#include "PK.hh"
#include "mpc.hh"

namespace Amanzi {

class MPCSubcycled : public MPC<PK> {
 public:
  MPCSubcycled(Teuchos::ParameterList& pk_tree_or_fe_list,
               const Teuchos::RCP<Teuchos::ParameterList>& global_list,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  void parseParameterList() override;

  // -- dt is the minimum of the sub pks
  double get_dt() override;
  void set_dt(double dt) override;
  void set_tags(const Tag& current, const Tag& next) override;

  void Setup() override;
  void Initialize() override;

  // -- advance each sub pk dt.
  bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;
  void CommitStep(double t_old, double t_new, const Tag& tag) override;

 protected:
  // advance the ith sub_pk
  bool AdvanceStep_i_(std::size_t i, double t_old, double t_new, bool reinit);

  Teuchos::Array<int> subcycling_;
  double dt_, target_dt_;
  std::vector<double> dts_;
  std::vector<std::pair<Tag, Tag>> tags_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCSubcycled> reg_;
};

} // namespace Amanzi

#endif
