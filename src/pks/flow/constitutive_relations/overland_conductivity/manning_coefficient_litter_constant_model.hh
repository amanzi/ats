/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

A constant Manning's n

`"Manning coefficient model type`" = `"constant`"

.. manning-coefficient-constant-spec:
.. admonition:: manning-coefficient-constant-spec

   * `"Manning coefficient [s m^-1/3]`" ``[double]``

*/

#ifndef AMANZI_FLOW_MANNING_LITTER_COEFFICIENT_CONSTANT_MODEL_HH_
#define AMANZI_FLOW_MANNING_LITTER_COEFFICIENT_CONSTANT_MODEL_HH_

#include "manning_coefficient_litter_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class ManningCoefficientLitterConstantModel : public ManningCoefficientLitterModel {
 public:
  ManningCoefficientLitterConstantModel(Teuchos::ParameterList& plist)
  {
    n_ = plist.get<double>("Manning coefficient [s m^-1/3]");
  }


  double ManningCoefficient(double litter_depth, double ponded_depth) const { return n_; }

  double DManningCoefficientDLitterThickness(double litter_depth, double ponded_depth) const
  {
    return 0;
  }

  double DManningCoefficientDPondedDepth(double litter_depth, double ponded_depth) const
  {
    return 0;
  }

 protected:
  double n_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

#endif
