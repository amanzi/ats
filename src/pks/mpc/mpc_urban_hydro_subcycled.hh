/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

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

#ifndef ATS_AMANZI_URBAN_HYDRO_SUBCYCLED_MPC_HH_
#define ATS_AMANZI_URBAN_HYDRO_SUBCYCLED_MPC_HH_

#include "PK.hh"
#include "mpc.hh"
#include "mpc_subcycled.hh"

namespace Amanzi {

class MPCUrbanHydroSub : public MPCSubcycled {

public:
  MPCUrbanHydroSub(Teuchos::ParameterList& pk_tree_or_fe_list,
                      const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt() override;
  virtual void set_dt(double dt) override;
  virtual void set_tags(const Tag& current, const Tag& next) override;

  virtual void Setup() override;
  virtual void Initialize() override;

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;

 protected:

  // advance the ith sub_pk
  // bool AdvanceStep_i_(std::size_t i, double t_old, double t_new, bool reinit);

  // Teuchos::Array<int> subcycling_;
  // double dt_, target_dt_;
  // std::vector<double> dts_;
  // std::vector<std::pair<Tag,Tag>> tags_;
  Key pipe_flow_src_key_, pipe_domain_;
  Key master_src_key_, master_domain_;
  int master_id_, subcycled_id_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCUrbanHydroSub> reg_;
};

}  // namespace Amanzi

#endif
