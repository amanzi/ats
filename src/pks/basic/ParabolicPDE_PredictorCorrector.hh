/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! A parabolic PDE solved using a predictor-corrector scheme.

/*!

This is the canonical nonlinear parabolic PDE, in mixed form for
predictor-corrector methods.

.. math::
    \frac{\partial \Psi(u) }{\partial t} - \nabla \cdot K(u) \nabla \Phi(u) = Q(u,x,t)

where:

- :math:`u` the primary variable, key_
- :math:`\Psi` the conserved quantity, conserved_quantity_key_
- :math:`\Phi` the potential field, potential_key_
- :math:`K` the diffusion coefficient
- :math:`Q` any source term, source_key_

Any of these may be a function of the primary variable :math:`u`.

Note that frequently \Phi(u) = u (e.g. the potential field is the primary
variable) but not always -- specifically for surface water where the potential
field is :math:`\delta(p) + z` for primary variable pressure.

This simply uses the Explicit scheme as a predictor for the implicit scheme.

.._pk-parabolic-pde-predictor-corrector-spec
.. admonition:: pk-parabolic-pde-predictor-corrector-spec

   INCLUDES:
   - mixin-parabolic-pde-mixed-form-implicit-spec
   - mixin-parabolic-pde-primary-form-explicit-spec
   - mixin-conservation-equation-spec
   - mixin-explicit-spec
   - mixin-leaf-spec
   - pk-spec

*/

#pragma once

#include "TreeVector.hh"
#include "Operator.hh"
#include "Operator_Factory.hh"
#include "Inverse.hh"
#include "SolverDefs.hh"

#include "PK_Adaptors.hh"
#include "ParabolicPDE_PrimaryFormExplicit.hh"
#include "ParabolicPDE_MixedFormImplicit.hh"
#include "PK_MixinExplicit.hh"
#include "PK_MixinImplicit.hh"
#include "PK_MixinPredictorCorrector.hh"
#include "PK_MixinLeaf.hh"
#include "PK_Default.hh"

namespace ATS {
namespace Basic {

using PK_ParabolicPDE_PredictorCorrector =
    PK_ImplicitExplicit_Adaptor<ParabolicPDE_MixedFormImplicit<
                          ParabolicPDE_PrimaryFormExplicit<
                            PK_MixinConservationEquation<
                              PK_MixinPredictorCorrector<
                                PK_MixinLeafCompositeVector<
                                  PK_Default>>>>>>;

} // namespace Basic
} // namespace ATS
