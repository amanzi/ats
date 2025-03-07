/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Default implementation of both BDF and Physical PKs.
/*!

A base class for all PKs that are both physical, in the sense that they
implement an equation and are not couplers, and BDF, in the sense that they
support the implicit integration interface.  This largely just supplies a
default error norm based on error in conservation relative to the extent of the
conserved quantity.

By default, the error norm used by solvers is given by:

:math:`ENORM(u, du) = |du| / ( a_tol + r_tol * |u| )`

The defaults here are typically good, or else good defaults are set in the
code, so usually are not supplied by the user.


.. _pk-physical-bdf-default-spec:
.. admonition:: pk-physical-bdf-default-spec

    * `"conserved quantity key`" ``[string]`` Name of the conserved quantity.
      Usually a sane default is set by the PK.

    * `"max valid change`" ``[double]`` **-1** Sets a limiter on what is a
      valid change in a single timestep.  Changes larger than this are declared
      invalid and the timestep shrinks.  By default, any change is valid.
      Units are the same as the primary variable.

    * `"absolute error tolerance`" ``[double]`` **1.0** Absolute tolerance,
      :math:`a_tol` in the equation above.  Unit are the same as the conserved
      quantity.  Note that this default is often overridden by PKs with more
      physical values, and very rarely are these set by the user.

    * `"relative error tolerance`" ``[double]`` **1.0** Relative tolerance,
      :math:`r_tol` in the equation above.  ``[-]`` Note that this default can
      be overridden by PKs with more physical values, and very rarely are these
      set by the user.

    * `"flux error tolerance`" ``[double]`` **1.0** Relative tolerance on the
      flux.  Note that this default is often overridden by PKs with more physical
      values, and very rarely are these set by the user.

    INCLUDES:

    - ``[pk-bdf-default-spec]`` *Is a* `PK: BDF`_
    - ``[pk-physical-default-spec]`` *Is a* `PK: Physical`_


*/

#ifndef ATS_PK_PHYSICAL_BDF_BASE_HH_
#define ATS_PK_PHYSICAL_BDF_BASE_HH_

#include "errors.hh"
#include "pk_bdf_default.hh"
#include "PK_Physical.hh"

#include "BCs.hh"
#include "Operator.hh"

namespace Amanzi {

class PK_PhysicalBDF_Default : public PK_BDF_Default, public PK_Physical {
 public:
  PK_PhysicalBDF_Default(Teuchos::ParameterList& pk_tree,
                         const Teuchos::RCP<Teuchos::ParameterList>& glist,
                         const Teuchos::RCP<State>& S,
                         const Teuchos::RCP<TreeVector>& solution)
    : PK(pk_tree, glist, S, solution),
      PK_BDF_Default(pk_tree, glist, S, solution),
      PK_Physical(pk_tree, glist, S, solution),
      max_valid_change_(-1.0)
  {}

  virtual void parseParameterList() override;

  virtual void Setup() override;

  // initialize.  Note both BDFBase and PhysicalBase have initialize()
  // methods, so we need a unique overrider.
  virtual void Initialize() override;

  // Default preconditioner is Picard
  virtual int
  ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) override;

  // updates the preconditioner, default does nothing
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override
  {}

  // Default implementations of BDFFnBase methods.
  // -- Compute a norm on u-du and return the result.
  virtual double
  ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) override;

  bool IsValid(const Teuchos::RCP<const TreeVector>& up) override;

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;
  virtual void FailStep(double t_old, double t_new, const Tag& tag) override;

  // -- Experimental approach -- calling this indicates that the time
  //    integration scheme is changing the value of the solution in
  //    state.
  virtual void ChangedSolution(const Tag& tag) override;
  virtual void ChangedSolution() override { return ChangedSolution(tag_next_); }

  // PC operator access
  Teuchos::RCP<Operators::Operator> preconditioner() { return preconditioner_; }

  // BC access
  std::vector<int>& bc_markers() { return bc_->bc_model(); }
  std::vector<double>& bc_values() { return bc_->bc_value(); }
  Teuchos::RCP<Operators::BCs> BCs() { return bc_; }

 protected:
  // PC
  Teuchos::RCP<Operators::Operator> preconditioner_;

  // BCs
  Teuchos::RCP<Operators::BCs> bc_;

  // step validity
  double max_valid_change_;

  // error criteria
  Key conserved_key_;
  Key cell_vol_key_;
  double atol_, rtol_, fluxtol_;
};


} // namespace Amanzi

#endif
