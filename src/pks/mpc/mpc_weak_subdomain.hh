/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Weak MPC for subdomain model MPCs.
/*!

  A weak MPC that couples the same PK across many subdomains.  Note that this
  means that the number of PKs is not known a priori -- it depends on a domain
  set.

 */

#pragma once
#include "mpc.hh"

namespace Amanzi {

class MPCWeakSubdomain : public MPC<PK> {
 public:
  MPCWeakSubdomain(Teuchos::ParameterList& FElist,
                   const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution);

  // PK methods
  void parseParameterList() override;

  // -- dt is the minimum of the sub pks
  double get_dt() override;
  void set_dt(double dt) override;
  void set_tags(const Tag& current, const Tag& next) override;

  void Setup() override;
  void Initialize() override;

  // -- advance each sub pk dt.
  bool AdvanceStep(double t_old, double t_new, bool reinit) override;
  void CommitStep(double t_old, double t_new, const Tag& tag_next) override;

 protected:
  void init_();

  bool AdvanceStep_Standard_(double t_old, double t_new, bool reinit);
  bool AdvanceStep_InternalSubcycling_(double t_old, double t_new, bool reinit);

  Tag get_ds_tag_next_(const std::string& subdomain);
  Tag get_ds_tag_current_(const std::string& subdomain);

  Comm_ptr_type comm_;
  bool internal_subcycling_;
  double internal_subcycling_target_dt_;
  double cycle_dt_;
  Key ds_name_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCWeakSubdomain> reg_;
};

} // namespace Amanzi
