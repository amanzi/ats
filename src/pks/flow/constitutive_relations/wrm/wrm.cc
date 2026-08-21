/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "errors.hh"
#include "exceptions.hh"
#include "Brent.hh"

#include "wrm.hh"

namespace Amanzi {
namespace Flow {

double
computeFieldCapacity(WRM& wrm, double k_rel_fc)
{
  double sr = wrm.residualSaturation();
  if (!(k_rel_fc > 0. && k_rel_fc < 1.)) {
    Errors::Message msg;
    msg << "computeFieldCapacity: \"field capacity relative permeability [-]\" must be in "
           "(0,1), got "
        << k_rel_fc;
    Exceptions::amanzi_throw(msg);
  }

  auto func = [&wrm, k_rel_fc](double s) { return wrm.k_relative(s) - k_rel_fc; };

  int itr = 100;
  double tol = 1.e-8;
  double result = Utils::findRootBrent(func, sr, 1.0, tol, &itr);
  if (itr < 0) {
    Errors::Message msg;
    msg << "computeFieldCapacity: failed to bracket a root of k_relative(s) - " << k_rel_fc
        << " on [" << sr << ", 1]";
    Exceptions::amanzi_throw(msg);
  }
  return result;
}

} // namespace Flow
} // namespace Amanzi
