/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! A model for providing the functional form of the exchange flux.

/*!

.. _micropore-macropore-flux-model-spec
.. admonition:: micropore-macropore-flux-model-spec

   * `"gamma [-]`" ``[double]`` The exchange coefficient.
   * `"delta [m]`" ``[double]`` Typical length scale between pores.

*/

#pragma once

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class MicroporeMacroporeFluxModel {

 public:
  explicit
  MicroporeMacroporeFluxModel(Teuchos::ParameterList& plist);

  double MicroporeMacroporeFlux(double pm, double pM, double krM, double krm, double K) const;

  double DMicroporeMacroporeFluxDMicroporePressure(double pm, double pM, double krM, double krm, double K) const;
  double DMicroporeMacroporeFluxDPressure(double pm, double pM, double krM, double krm, double K) const;
  double DMicroporeMacroporeFluxDRelativePermeability(double pm, double pM, double krM, double krm, double K) const;
  double DMicroporeMacroporeFluxDMicroporeRelativePermeability(double pm, double pM, double krM, double krm, double K) const;
  double DMicroporeMacroporeFluxDMicroporeAbsolutePermeability(double pm, double pM, double krM, double krm, double K) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:

  double gamma_;
  double delta_;

};

} //namespace
} //namespace
} //namespace


