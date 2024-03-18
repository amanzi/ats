/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates the unfrozen fraction of water.

*/
/*!

.. _unfrozen-fraction-model-spec:
.. admonition:: unfrozen-fraction-model-spec

   * `"transition width [K]`" ``[double]`` **0.2** Degrees over which to
     transition from no ice to all ice.

   * `"freezing point [K]`" ``[double]`` **273.15** Center of the transition,
     at this point unfrozen fraction is 0.5.

   * `"minimum unfrozen fraction [-]`" ``[double]`` **0** Sets a minimum value.

*/

#ifndef AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_MODEL_
#define AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_MODEL_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {

class UnfrozenFractionModel {
 public:
  UnfrozenFractionModel(Teuchos::ParameterList& list);

  double UnfrozenFraction(double temp) const;
  double DUnfrozenFractionDT(double temp) const;

 protected:
  double halfwidth_;
  double T0_;
  double min_uf_;
  const double pi_;

  Teuchos::ParameterList plist_;
};

} // namespace Flow
} // namespace Amanzi

#endif
