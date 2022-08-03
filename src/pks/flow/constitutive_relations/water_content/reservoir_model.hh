/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
/*

Base class and factory for a reservoir model.

Also supplies the dumbest reservoir model ever for testing.

*/

#pragma once

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {
namespace Relations {


class ReservoirModel {
 public:
  virtual ~ReservoirModel() = default;
  virtual double computeDischarge(double wc) = 0;
};

Teuchos::RCP<ReservoirModel>
createReservoirModel(Teuchos::ParameterList& plist);


// and a really dumb one
class DumbReservoirModel : public ReservoirModel {
  virtual double computeDischarge(double wc) override {
    return 0.1 * wc;
  }
};


} //namespace
} //namespace
} //namespace

