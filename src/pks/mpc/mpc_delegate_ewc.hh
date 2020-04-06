/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Globalization for nonlinearity associated with phase change and latent heat.

/*!

The EWC delegate works to deal with strong nonlinearities associated with
latent heat and phase change.  Provided a change in primary variables pressure
and temperature, it works by first multiplying those changes by the local
Jacobian matrix, :math:`\frac{\partial \left\{ \Theta, E \right\} }{ \partial
\left\{ p, T \right\} }` to calculate changes in water content and energy, then
calculating the new water content and energy and inverting the functions
:math:`\Theta(p,T), E(p,T)` to determine what pressure and temperature would
have resulted in those values.  This provides a corrected change in the primary
variables.

Conceptually, this is a "more robust" choice in nonlinearities associated with
phase change, where the derivatives go from small to large to small again, and
small changes in pressure and temperature result in large changes in water
content and energy.

This delegate manages these globalization strategies, which can be used both in
modifying the correction supplied by a nonlinear iterate, and in modifying a
predictor, the extrapolated projection (from previous timesteps) that
provides the initial guess to the nonlinear solve.

.. _mpc-delegate-ewc-spec:
.. admonition:: mpc-delegate-ewc-spec

    * `"verbose object`" ``[verbose-object-spec]`` See `Verbose Object`_.
    
    * `"PK name`" ``[string]`` Name of the owning PK -- simply for logging and
      debugging.
    * `"domain name`" ``[string]`` **"domain"** The mesh.

    * `"preconditioner type`" ``[string]`` When to use EWC on the nonlinear
      iterate's correction.  One of:

      - `"none`" Never do EWC
      - `"ewc`" Always do EWC
      - `"smart ewc`" Attempt EWC when it seems likely it will be useful and
        take the EWC correction if it is smaller than the standard correction.

    * `"predictor type`" ``[string]`` When to use EWC on the predictor.  One
      of:
      
      - `"none`" Never do EWC
      - `"ewc`" Always do EWC
      - `"smart ewc`" Attempt EWC when it seems likely it will be useful and
        take the EWC correction if it is smaller than the standard correction.

    * `"freeze-thaw cusp width [K]`" ``[double]`` Controls a width over which
      to assume we are close to the latent heat cliff, and begins applying the
      EWC algorithm in `"ewc smarter`".

    * `"freeze-thaw cusp width (freezing) [K]`" ``[double]`` Controls a width
      over which to assume we are close to the latent heat cliff as we get
      colder, and begins applying the EWC algorithm in `"ewc smarter`".

    * `"freeze-thaw cusp width (thawing) [K]`" ``[double]`` Controls a width
      over which to assume we are close to the latent heat cliff as we get
      warmer, and begins applying the EWC algorithm in `"ewc smarter`".
        
    * `"pressure key`" ``[string]`` **DOMAIN-pressure**
    * `"temperature key`" ``[string]`` **DOMAIN-temperature**
    * `"water content key`" ``[string]`` **DOMAIN-water_content**
    * `"energy key`" ``[string]`` **DOMAIN-energy**
    * `"cell volume key`" ``[string]`` **DOMAIN-cell_volume**

    INCLUDES

    - ``[debugger-spec]`` Uses a Debugger_
    
*/

#ifndef MPC_DELEGATE_EWC_HH_
#define MPC_DELEGATE_EWC_HH_

#include "VerboseObject.hh"
#include "Debugger.hh"
#include "Tensor.hh"
#include "State.hh"
#include "TreeVector.hh"

namespace Amanzi {

class EWCModel;

class MPCDelegateEWC {

 public:

  MPCDelegateEWC(Teuchos::ParameterList& plist);
  virtual ~MPCDelegateEWC() = default;

  virtual void setup(const Teuchos::Ptr<State>& S);
  virtual void initialize(const Teuchos::Ptr<State>& S);
  virtual void set_states(const Teuchos::RCP<State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);
  virtual bool ModifyPredictor(double h, Teuchos::RCP<TreeVector> up);
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  void set_model(const Teuchos::RCP<EWCModel>& model) { model_ = model; }

 protected:
  virtual bool modify_predictor_smart_ewc_(double h, Teuchos::RCP<TreeVector> up) = 0;
  virtual void precon_ewc_(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<TreeVector> Pu) = 0;

  virtual void update_precon_ewc_(double t, Teuchos::RCP<const TreeVector> up, double h);



 protected:
  Teuchos::RCP<Teuchos::ParameterList> plist_;
  Teuchos::RCP<VerboseObject> vo_;
  Teuchos::RCP<Debugger> db_;

  // model
  Teuchos::RCP<EWCModel> model_;

  enum PredictorType {
    PREDICTOR_NONE = 0,
    PREDICTOR_EWC,
    PREDICTOR_SMART_EWC
  };

  enum PreconditionerType {
    PRECON_NONE = 0,
    PRECON_EWC,
    PRECON_SMART_EWC,
  };

  // control flags
  PreconditionerType precon_type_;
  PredictorType predictor_type_;

  // extra data
  std::vector<WhetStone::Tensor> jac_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<Epetra_MultiVector> wc_prev2_;
  Teuchos::RCP<Epetra_MultiVector> e_prev2_;
  double time_prev2_;

  // parameters for heuristic
  double cusp_size_T_freezing_;
  double cusp_size_T_thawing_;

  // states
  Teuchos::RCP<State> S_next_;
  Teuchos::RCP<State> S_inter_;
  Teuchos::RCP<const State> S_;

  // keys
  Key pres_key_;
  Key temp_key_;
  Key e_key_;
  Key wc_key_;
  Key cv_key_;

};

} // namespace

#endif
