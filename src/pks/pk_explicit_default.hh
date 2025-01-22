/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   ATS

   Default base with default implementations of methods for a PK integrated using
   Explicit.
   ------------------------------------------------------------------------- */

#ifndef ATS_PK_EXPLICIT_DEFAULT_HH_
#define ATS_PK_EXPLICIT_DEFAULT_HH_

#include "Teuchos_TimeMonitor.hpp"
#include "Explicit_TI_RK.hh"
#include "PK_Explicit.hh"
#include "TreeVector.hh"
#include "Epetra_Vector.h"

namespace Amanzi {

class TreeVector;
class PK_Explicit_Default : public PK_Explicit<TreeVector> {
 public:
  PK_Explicit_Default(Teuchos::ParameterList& pk_tree,
                      const Teuchos::RCP<Teuchos::ParameterList>& glist,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& solution)
    : PK(pk_tree, glist, S, solution), PK_Explicit<TreeVector>(pk_tree, glist, S, solution)
  {}

  // Virtual destructor
  virtual ~PK_Explicit_Default(){};

  // Default implementations of PK methods.
  // -- setup
  virtual void Setup() override;

  // -- initialize
  virtual void Initialize() override;

  // -- Choose a timestep compatible with physics.
  virtual double get_dt() override { return dt_; }
  virtual void set_dt(double dt) override { dt_ = dt; }

  // -- Advance from state S0 to state S1 at time S0.time + dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override;

 protected: //data  timestep control
  double dt_;
  Teuchos::RCP<Explicit_TI::RK<TreeVector>> time_stepper_;

  // timing
  Teuchos::RCP<Teuchos::Time> step_walltime_;

  // solution at the old timestep
  Teuchos::RCP<TreeVector> solution_old_;
};

} // namespace Amanzi

#endif
