/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

.. _overland-conductivity-manning-spec:
.. admonition:: overland-conductivity-manning-spec

   * `"Manning exponent`" ``[double]`` **2/3**
   * `"slope regularization epsilon`" ``[double]`` **1.e-8** In ATS's
     implementation of the diffusion wave equation, it is expected that |S| >
     0.  This may be arbitrarily small, but it keeps slope from being exactly
     0, which crashes the code.
   * `"maximum ponded depth [m]`" ``[double]`` **1.e8** Arbitrarily large
     ponded depth creates arbitrarily large flowing velocities -- sometimes we
     wish to use this model (even though it is incorrect) for larger rivers.
     This limits the velocity from being unbounded.

 */

#ifndef AMANZI_FLOWRELATIONS_MANNING_CONDUCTIVITY_MODEL_
#define AMANZI_FLOWRELATIONS_MANNING_CONDUCTIVITY_MODEL_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class ManningConductivityModel {
 public:
  explicit ManningConductivityModel(Teuchos::ParameterList& plist);

  double Conductivity(double depth, double slope, double coef);
  double DConductivityDDepth(double depth, double slope, double coef);

 protected:
  double slope_regularization_;
  double manning_exp_;
  double depth_max_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
