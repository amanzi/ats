/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

//! Implements a predictor-corrector strategy.

/*!
See additional documentation in the base class src/pks/mpc/MPC.hh
*/

#ifndef PKS_MPC_PREDICTOR_CORRECTOR_HH_
#define PKS_MPC_PREDICTOR_CORRECTOR_HH_

#include <vector>

#include "strong_mpc.hh"
#include "pk_bdf_default.hh"


namespace Amanzi {

template <class PK_t>
class MPCPredictorCorrector :  public StrongMPC<PK_t> {

public:

  MPCPredictorCorrector(Teuchos::ParameterList& FElist,
            const Teuchos::RCP<Teuchos::ParameterList>& plist,
            const Teuchos::RCP<State>& S,
            const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~MPCPredictorCorrector() {}

  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<State>& S_inter,
                  const Teuchos::RCP<State>& S_next);

  virtual double get_dt() {
    return std::min(StrongMPC<PK_t>::get_dt(), predictor_pk_->get_dt());
  }
  virtual void set_dt(double dt) {
    StrongMPC<PK_t>::set_dt(dt);
    predictor_pk_->set_dt(dt);
  }
  
  
  virtual void CommitStep(double t_old, double t_new,
                          const Teuchos::RCP<State>& S);
  
  // -- Modify the predictor.
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u);

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);
  
protected:
  Teuchos::RCP<PK> predictor_pk_;
  Teuchos::RCP<TreeVector> predictor_soln_;
  bool predictor_fail_;

  int pred_fail_gi_fail_;
  int pred_fail_gi_good_;
  int pred_good_gi_fail_;
  int pred_good_gi_itr0_;
  int pred_good_gi_itrN_;
  
  
private:
  // factory registration
  static RegisteredPKFactory<MPCPredictorCorrector> reg_;

};

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
template<class PK_t>
MPCPredictorCorrector<PK_t>::MPCPredictorCorrector(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln) :
    PK(pk_tree, global_list, S, soln),
    StrongMPC<PK_t>(pk_tree, global_list, S, soln),
    pred_fail_gi_fail_(0),
    pred_fail_gi_good_(0),
    pred_good_gi_fail_(0),
    pred_good_gi_itr0_(0),
    pred_good_gi_itrN_(0)
{
  // create the predictor PK
  Teuchos::RCP<Teuchos::ParameterList> pk_tree_list = Teuchos::sublist(this->plist_, "predictor PK tree");
  Teuchos::ParameterList::ConstIterator pk_item = pk_tree_list->begin();
  const std::string &pk_name = pk_tree_list->name(pk_item);

  predictor_soln_ = Teuchos::rcp(new Amanzi::TreeVector());
  
  Amanzi::PKFactory pk_factory;
  predictor_pk_ = pk_factory.CreatePK(pk_name, *pk_tree_list, global_list, S, predictor_soln_);

}


// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
template<class PK_t>
void MPCPredictorCorrector<PK_t>::Setup(const Teuchos::Ptr<State>& S) {
  StrongMPC<PK_t>::Setup(S);
  predictor_pk_->Setup(S);
};


// -----------------------------------------------------------------------------
// Required unique initialize(), just calls both of its base class
// initialize() methods.
// -----------------------------------------------------------------------------
template<class PK_t>
void MPCPredictorCorrector<PK_t>::Initialize(const Teuchos::Ptr<State>& S) {
  StrongMPC<PK_t>::Initialize(S);
  predictor_pk_->Initialize(S);
};

template<class PK_t>
void MPCPredictorCorrector<PK_t>::set_states(const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<State>& S_inter,
                                 const Teuchos::RCP<State>& S_next){
  StrongMPC<PK_t>::set_states(S, S_inter, S_next);
  predictor_pk_->set_states(S, S_inter, S_next);
} 

template<class PK_t>
void MPCPredictorCorrector<PK_t>::CommitStep(double t_old, double t_new,
        const Teuchos::RCP<State>& S) {
  // commit the BDF step
  StrongMPC<PK_t>::CommitStep(t_old, t_new, S);

  // copy the BDF solution to the explicit solution
  // NOTE special purpose
  *predictor_soln_->SubVector(1)->SubVector(1) = *this->solution_->SubVector(0)->SubVector(1);
  *predictor_soln_->SubVector(1)->SubVector(0) = *this->solution_->SubVector(0)->SubVector(0);
  predictor_pk_->ChangedSolutionPK(S.ptr());

  // now commit the predictor, storing this solution
  predictor_pk_->CommitStep(t_old, t_new, S);
}


// -----------------------------------------------------------------------------
// Modify predictor from each sub pk.
// -----------------------------------------------------------------------------
template<class PK_t>
bool MPCPredictorCorrector<PK_t>::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
        Teuchos::RCP<TreeVector> u) {
  predictor_fail_ = predictor_pk_->AdvanceStep(this->S_inter_->time(), this->S_next_->time(), false);
  
  // copy explicit prediction into BDF
  // NOTE special purpose
  *u->SubVector(0)->SubVector(0) = *predictor_soln_->SubVector(1)->SubVector(0);
  *u->SubVector(0)->SubVector(1) = *predictor_soln_->SubVector(1)->SubVector(1);

  return true;
};


template<class PK_t>
bool MPCPredictorCorrector<PK_t>::AdvanceStep(double t_old, double t_new, bool reinit) {
  bool fail = StrongMPC<PK_t>::AdvanceStep(t_old, t_new, reinit);
  if (predictor_fail_) {
    if (fail) pred_fail_gi_fail_++;
    else pred_fail_gi_good_++;
  } else {
    if (fail) pred_good_gi_fail_++;
    else if (this->time_stepper_->number_nonlinear_steps() == 0) pred_good_gi_itr0_++;
    else pred_good_gi_itrN_++;
  }

  if (this->vo_->os_OK(Teuchos::VERB_HIGH)) {
    *this->vo_->os() << "Tallies: [P_FAIL,GI_FAIL=" << pred_fail_gi_fail_ << "][P_FAIL,GI_GOOD=" << pred_fail_gi_good_
                     << "][P_GOOD,GI_FAIL=" << pred_good_gi_fail_ << "][P_GOOD,GI_0=" << pred_good_gi_itr0_ << "][P_GOOD,GI_N=" << pred_good_gi_itrN_ << "]" << std::endl;
  }
  return fail;
}


} // close namespace Amanzi

#endif
