/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

This is a :ref:`Strong MPC` which uses a preconditioner in which the
block-diagonal cell-local matrix is dense.  If the system looks something
like:

.. math::

   A( y_1, y_2, x, t ) = 0
   B( y_1, y_2, x, t ) = 0

where :math:`y_1, y_2` are spatially varying unknowns that are discretized
using the MFD method (and therefore have both cell and face unknowns), an
approximation to the Jacobian is written as

.. math::

   \begin{bmatrix}
     \frac{dA^c}{dy_1^c}  & \frac{dA^c}{dy_1^f} & \frac{dA^c}{dy_2^c} & 0                   \\
     \frac{dA^f}{dy_1^c}  & \frac{dA^f}{dy_1^f} & 0                   & 0                   \\
     \frac{dB^c}{dy_1^c}  & 0                   & \frac{dB^c}{dy_2^c} & \frac{dB^c}{dy_2^f} \\
     0                    & 0                   & \frac{dB^f}{dy_2^c} & \frac{dB^f}{dy_2^f}
   \end{bmatrix}


Note that the upper left block is the standard preconditioner for the A system,
and the lower right block is the standard precon for the B system, and we have
simply added cell-based couplings, :math:`\frac{dA^c}{dy_2^c}` and
:math:`\frac{dB^c}{dy_1^c}`.

Most commonly this is used to couple flow and energy equations on the same
mesh.  In the temperature/pressure system, these extra blocks correspond to

.. math::
    \frac{\partial \Theta}{\partial T} \; , \; \frac{\partial E}{\partial p}

`"PK type`" = `"mpc coupled cells`"

.. _pk-mpc-coupled-cells-spec:
.. admonition:: pk-mpc-coupled-cells-spec

   * `"domain name`" ``[string]`` Domain of simulation
   * `"conserved quantity A`" ``[string]`` Key of the first sub-PK's conserved quantity.
   * `"conserved quantity B`" ``[string]`` Key of the second sub-PK's conserved quantity.
   * `"primary variable A`" ``[string]`` Key of the first sub-PK's primary variable.
   * `"primary variable B`" ``[string]`` Key of the second sub-PK's primary variable.
   * `"no dA/dy2 block`" ``[bool]`` **false** Excludes the dA^c/dy_2^c block above.
   * `"no dB/dy1 block`" ``[bool]`` **false** Excludes the dB^c/dy1^c block above.

   INCLUDES:

   - ``[strong-mpc-spec]`` *Is a* :ref:`Strong MPC`.

*/

#ifndef MPC_COUPLED_CELLS_HH_
#define MPC_COUPLED_CELLS_HH_

#include "pk_physical_bdf_default.hh"
#include "strong_mpc.hh"

namespace Amanzi {

namespace Operators {
class TreeOperator;
class PDE_Accumulation;
} // namespace Operators

class MPCCoupledCells : public StrongMPC<ATS_Physics::PK_PhysicalBDF_Default> {
 public:
  MPCCoupledCells(Teuchos::ParameterList& FElist,
                  const Teuchos::RCP<Teuchos::ParameterList>& plist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& solution)
    : PK(FElist, plist, S, solution), StrongMPC<ATS_Physics::PK_PhysicalBDF_Default>(FElist, plist, S, solution)
  {}

  void parseParameterList() override;
  virtual void Setup() override;

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<TreeVector> Pu) override;

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

 protected:
  Key A_key_;
  Key y2_key_;
  Key B_key_;
  Key y1_key_;
  Teuchos::RCP<Operators::TreeOperator> preconditioner_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  Teuchos::RCP<Operators::PDE_Accumulation> dA_dy2_;
  Teuchos::RCP<Operators::PDE_Accumulation> dB_dy1_;

  // cruft for easier global debugging
  Teuchos::RCP<Debugger> db_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCCoupledCells> reg_;
};


} // namespace Amanzi
#endif
