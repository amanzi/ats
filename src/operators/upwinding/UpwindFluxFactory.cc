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

#include "State.hh"
#include "UpwindFluxFactory.hh"
#include "upwind_total_flux.hh"
#include "upwind_flux_harmonic_mean.hh"
#include "upwind_flux_split_denominator.hh"
#include "upwind_elevation_stabilized.hh"
#include "upwind_flux_fo_cont.hh"
#include "upwind_cell_centered.hh"

namespace Amanzi {
namespace Operators {

Teuchos::RCP<Upwinding>
UpwindFluxFactory::Create(Teuchos::ParameterList& oplist,
                          const Teuchos::Ptr<State>& S,
                          std::string pkname,
                          std::string cell_coef,
                          std::string face_coef,
                          std::string flux)
{
  std::string model_type = oplist.get<std::string>("upwind type", "manning upwind");
  double flux_eps = oplist.get<double>("upwind flux epsilon", 1.e-8);
  auto domain = Keys::getDomain(face_coef);

  if (model_type == "manning upwind") {
    S->RequireField(flux)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    S->RequireField(cell_coef)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    return Teuchos::rcp(new UpwindTotalFlux(pkname, cell_coef, face_coef, flux, flux_eps));

  } else if (model_type == "manning harmonic mean") {
    // this is dangerous because it can result in 0 flux when there is a
    // neighboring dry cell, which almost is never what is really wanted.  This
    // option has shown up in a few input files and broken user code, so for
    // now I am removing it as an option.
    // std::stringstream message;
    // message << "============== WARNING WARNING WARNING WARNING WARNING =============" << std::endl
    //           << "Are you certain you intend to use the option:" << std::endl
    //           << "  \"overland conductivity model\" = \"manning harmonic mean\"?"  << std::endl
    //           << "As in, really really certain?  You might be better off removing this" << std::endl
    //           << "option and using the default value of \"manning upwind\" unless you know" << std::endl
    //           << "what you are doing!" << std::endl
    //           << "============== WARNING WARNING WARNING WARNING WARNING =============" << std::endl;
    // Errors::Message msg(message.str());
    // Exceptions::amanzi_throw(msg);
    S->RequireField(flux)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    S->RequireField(cell_coef)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    return Teuchos::rcp(new UpwindFluxHarmonicMean(pkname, cell_coef, face_coef, flux, flux_eps));

  } else if (model_type == "manning split denominator") {
    Key domain = Keys::getDomain(cell_coef);
    std::string slope = Keys::readKey(oplist, domain, "slope", "slope_magnitude");
    std::string manning_coef = Keys::readKey(oplist, domain, "coefficient", "manning_coefficient");
    double slope_regularization = oplist.get<double>("slope regularization epsilon", 1.e-2);

    S->RequireField(flux)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    S->RequireField(cell_coef)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireField(slope)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireField(manning_coef)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    return Teuchos::rcp(new UpwindFluxSplitDenominator(pkname, cell_coef, face_coef, flux, flux_eps, slope, manning_coef, slope_regularization));

  } else if (model_type == "manning elevation stabilized") {
    Key domain = Keys::getDomain(cell_coef);
    std::string slope = Keys::readKey(oplist, domain, "slope", "slope_magnitude");
    std::string manning_coef = Keys::readKey(oplist, domain, "coefficient", "manning_coefficient");
    std::string ponded_depth = Keys::readKey(oplist, domain, "ponded depth", "ponded_depth");
    std::string elev = Keys::readKey(oplist, domain, "elevation", "elevation");
    std::string dens = Keys::readKey(oplist, domain, "molar density liquid", "molar_density_liquid");
    double slope_regularization = oplist.get<double>("slope regularization epsilon", 1.e-2);
    double manning_exp = oplist.get<double>("Manning exponent", 2.0/3);

    S->RequireField(slope)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireField(elev)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireField(manning_coef)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    // add boundary face components?

    return Teuchos::rcp(new UpwindElevationStabilized(pkname, face_coef, slope, manning_coef, ponded_depth, elev, dens, slope_regularization, manning_exp));

  } else if (model_type == "manning ponded depth passthrough") {
    Key domain = Keys::getDomain(cell_coef);
    std::string slope = Keys::readKey(oplist, domain, "slope", "slope_magnitude");
    std::string manning_coef = Keys::readKey(oplist, domain, "coefficient", "manning_coefficient");
    std::string elev = Keys::readKey(oplist, domain, "elevation", "elevation");
    double slope_regularization = oplist.get<double>("slope regularization epsilon", 1.e-8);
    double manning_exp = oplist.get<double>("Manning exponent");

    S->RequireField(flux)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    S->RequireField(cell_coef)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireField(slope)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireField(elev)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireField(manning_coef)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    return Teuchos::rcp(new UpwindFluxFOCont(pkname, cell_coef, face_coef, flux, slope, manning_coef, elev, slope_regularization, manning_exp));

  } else if (model_type == "manning cell centered") {
    S->RequireField(cell_coef)->SetGhosted()->SetMesh(S->GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    return Teuchos::rcp(new UpwindCellCentered(pkname, cell_coef, face_coef));
  } else {
    Errors::Message msg;
    msg << "Unknown \"upwind type\" value \"" << model_type << ",\" must be one of \"manning upwind\", \"manning harmonic mean\", or \"manning ponded depth passthrough.\"";
    Exceptions::amanzi_throw(msg);
  }
  return Teuchos::null;
}

}  // namespace Operators
}  // namespace Amanzi


