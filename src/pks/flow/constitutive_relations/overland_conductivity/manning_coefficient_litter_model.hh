/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The manning coefficient with variable litter model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Manning's coefficient that varies based on litter thickness and ponded depth.

*/

#ifndef AMANZI_FLOW_MANNING_COEFFICIENT_LITTER_MODEL_HH_
#define AMANZI_FLOW_MANNING_COEFFICIENT_LITTER_MODEL_HH_

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace Relations {

class ManningCoefficientLitterModel {
 public:
  virtual ~ManningCoefficientLitterModel() {}

  virtual double ManningCoefficient(double litter_depth, double ponded_depth) const = 0;

  virtual double DManningCoefficientDLitterThickness(double litter_depth,
                                                     double ponded_depth) const = 0;

  virtual double DManningCoefficientDPondedDepth(double litter_depth,
                                                 double ponded_depth) const = 0;
};

} // namespace Relations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
