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

#include "exceptions.hh"
#include "errors.hh"
#include "reservoir_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


Teuchos::RCP<ReservoirModel>
createReservoirModel(Teuchos::ParameterList& plist) {
  std::string reservoir_model_type = plist.get<std::string>("reservoir model type");
  Teuchos::RCP<ReservoirModel> model;
  if (reservoir_model_type == "dumb") {
    model = Teuchos::rcp(new DumbReservoirModel());
  } else {
    Errors::Message msg;
    msg << "ReservoirModel: unknown model type \"" << reservoir_model_type << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return model;
}


} //namespace
} //namespace
} //namespace

