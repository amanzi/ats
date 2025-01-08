/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The interfrost dtheta_dpressure model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Interfrost water content portion sl.

*/

#ifndef AMANZI_FLOW_INTERFROST_DTHETA_DPRESSURE_MODEL_HH_
#define AMANZI_FLOW_INTERFROST_DTHETA_DPRESSURE_MODEL_HH_

namespace Amanzi {
namespace Flow {
namespace Relations {

class InterfrostDthetaDpressureModel {
 public:
  explicit InterfrostDthetaDpressureModel(Teuchos::ParameterList& plist);

  double DThetaDpCoef(double nl, double sl, double phi) const;

  double DDThetaDpCoefDMolarDensityLiquid(double nl, double sl, double phi) const;
  double DDThetaDpCoefDSaturationLiquid(double nl, double sl, double phi) const;
  double DDThetaDpCoefDPorosity(double nl, double sl, double phi) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:
  double beta_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

#endif
