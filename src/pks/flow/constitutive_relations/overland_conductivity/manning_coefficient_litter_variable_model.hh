/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

A Manning coefficient with variable litter thickness.  Manning's n is taken to
vary with litter depth.  If ponded depth is less than the litter depth, then n
is given by litter n.  If ponded depth is greater than litter depth, it is
approaches bare ground n for ponded depth >> litter depth.

`"Manning coefficient model type`" = `"variable`"

.. manning-coefficient-variable-spec:
.. admonition:: manning-coefficient-variable-spec

   * `"Manning coefficient bare ground [s m^-1/3]`" ``[double]`` **0.02**
   * `"Manning coefficient litter [s m^-1/3]`" ``[double]`` **0.1**

*/

#ifndef AMANZI_FLOW_MANNING_LITTER_COEFFICIENT_VARIABLE_MODEL_HH_
#define AMANZI_FLOW_MANNING_LITTER_COEFFICIENT_VARIABLE_MODEL_HH_

#include "manning_coefficient_litter_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class ManningCoefficientLitterVariableModel : public ManningCoefficientLitterModel {
 public:
  ManningCoefficientLitterVariableModel(Teuchos::ParameterList& plist)
  {
    n_bg_ = plist.get<double>("Manning coefficient bare ground [s m^-1/3]", 0.02);
    n_l_ = plist.get<double>("Manning coefficient litter [s m^-1/3]", 0.1);
  }


  double ManningCoefficient(double ld, double pd) const
  {
    double n;
    if (ld > 0.) {
      n = n_l_;
    } else {
      n = n_bg_;
    }

    if (pd > 0) {
      if (pd <= ld) {
        n = n_l_;
      } else {
        n = (ld * n_l_ + n_bg_ * (-ld + pd)) / pd;
      }
    }

    return n;
  }


  double DManningCoefficientDLitterThickness(double ld, double pd) const
  {
    double n = 0.;

    if (pd > 0 && pd > ld) {
      n = n_l_ / pd - n_bg_ / pd;
    }
    return n;
  }


  double DManningCoefficientDPondedDepth(double ld, double pd) const
  {
    double n = 0.;
    if (pd > 0 && pd > ld) {
      n = n_bg_ / pd - (ld * n_l_ + n_bg_ * (-ld + pd)) / std::pow(pd, 2);
    }
    return n;
  }

 protected:
  double n_bg_;
  double n_l_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

#endif
