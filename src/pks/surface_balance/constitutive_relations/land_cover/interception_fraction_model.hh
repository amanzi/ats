/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Fraction of incoming water that is intercepted.
/*!

Based on CLM 4.5 and Lawrence et al 2007:

.. math::
  f_I = \alpha * (1 - exp(-.5 LAI)))

The interception fraction is everything here after the precip.

.. _interception-fraction-model-spec:
.. admonition:: interception-fraction-model-spec

   * `"leaf area interception fraction [-]`" ``[double]`` **0.25** The alpha
      term, this describes the fraction of leaf area that intercepts water.

*/

#pragma once

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class InterceptionFractionModel {
 public:
  explicit InterceptionFractionModel(Teuchos::ParameterList& plist);

  double InterceptionFraction(double ai) const;

  double DInterceptionFractionDAreaIndex(double ai) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:
  double alpha_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
