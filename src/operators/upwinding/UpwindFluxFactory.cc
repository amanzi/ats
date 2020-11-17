// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
//
// Factory for taking coefficients for div-grad operators from cells to
// faces.
// -----------------------------------------------------------------------------

#include "errors.hh"

#include "Teuchos_ParameterList.hpp"

#include "UpwindFluxFactory.hh"
#include "upwind_total_flux.hh"
#include "upwind_flux_harmonic_mean.hh"
#include "upwind_flux_split_denominator.hh"
#include "upwind_flux_fo_cont.hh"
#include "upwind_cell_centered.hh"

namespace Amanzi {
namespace Operators {

Teuchos::RCP<Upwinding>
UpwindFluxFactory::Create(Teuchos::ParameterList& oplist,
                          std::string pkname,
                          std::string cell_coef,
                          std::string face_coef,
                          std::string flux)
{
  std::string model_type = oplist.get<std::string>("upwind type", "manning upwind");
  double flux_eps = oplist.get<double>("upwind flux epsilon", 1.e-8);

  if (model_type == "manning upwind") {
    return Teuchos::rcp(new UpwindTotalFlux(pkname, cell_coef, face_coef, flux, flux_eps));

  } else if (model_type == "manning harmonic mean") {
    // // this is dangerous because it can result in 0 flux when there is a
    // // neighboring dry cell, which almost is never what is really wanted.  This
    // // option has shown up in a few input files and broken user code, so for
    // // now I am removing it as an option.
    // std::stringstream message;
    // message << "============== WARNING WARNING WARNING WARNING WARNING =============" << std::endl
    //           << "Are you certain you intend to use the option:" << std::endl
    //           << "  \"overland conductivity model\" = \"manning harmonic mean\"?"  << std::endl
    //           << "As in, really really certain?  You might be better off removing this" << std::endl
    //           << "option and using the default value of \"manning\" unless you know" << std::endl
    //           << "what you are doing!" << std::endl
    //           << "============== WARNING WARNING WARNING WARNING WARNING =============" << std::endl;
    // Errors::Message msg(message.str());
    // Exceptions::amanzi_throw(msg);
    return Teuchos::rcp(new UpwindFluxHarmonicMean(pkname, cell_coef, face_coef, flux, flux_eps));

  } else if (model_type == "manning split denominator") {
    Key domain = Keys::getDomain(cell_coef);
    std::string slope = Keys::readKey(oplist, domain, "slope", "slope_magnitude");
    std::string manning_coef = Keys::readKey(oplist, domain, "coefficient", "manning_coefficient");
    std::string ponded_depth = Keys::readKey(oplist, domain, "height", "ponded_depth");
    double slope_regularization = oplist.get<double>("slope regularization epsilon", 1.e-8);
    return Teuchos::rcp(new UpwindFluxSplitDenominator(pkname, cell_coef, face_coef, flux, flux_eps, slope, manning_coef, slope_regularization, ponded_depth));

  } else if (model_type == "manning ponded depth passthrough") {
    Key domain = Keys::getDomain(cell_coef);
    std::string slope = Keys::readKey(oplist, domain, "slope", "slope_magnitude");
    std::string manning_coef = Keys::readKey(oplist, domain, "coefficient", "manning_coefficient");
    std::string elev = Keys::readKey(oplist, domain, "elevation", "elevation");
    double slope_regularization = oplist.get<double>("slope regularization epsilon", 1.e-8);
    double manning_exp = oplist.get<double>("Manning exponent");
    return Teuchos::rcp(new UpwindFluxFOCont(pkname, cell_coef, face_coef, flux, slope, manning_coef, elev, slope_regularization, manning_exp));

  } else if (model_type == "manning cell centered") {
    return Teuchos::rcp(new UpwindCellCentered(pkname, cell_coef, face_coef));
  } else {
    AMANZI_ASSERT(0);
    return Teuchos::null;
  }
}
  
}  // namespace Operators
}  // namespace Amanzi


