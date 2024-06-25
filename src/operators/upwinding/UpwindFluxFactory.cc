/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Factory for taking coefficients for div-grad operators from cells to faces.
#include "errors.hh"

#include "Teuchos_ParameterList.hpp"

#include "State.hh"
#include "UpwindFluxFactory.hh"
#include "upwind_total_flux.hh"
#include "upwind_flux_harmonic_mean.hh"
#include "upwind_flux_split_denominator.hh"
//#include "upwind_elevation_stabilized.hh"
#include "upwind_flux_fo_cont.hh"
#include "upwind_cell_centered.hh"

namespace Amanzi {
namespace Operators {
namespace UpwindFactory {

Teuchos::RCP<Upwinding>
Create(Teuchos::ParameterList& oplist,
       State& S,
       const std::string& pkname,
       const Tag& tag,
       const Key& flux)
{
  std::string model_type = oplist.get<std::string>("upwind type", "manning upwind");
  double flux_eps = oplist.get<double>("upwind flux epsilon", 1.e-8);
  auto domain = Keys::getDomain(flux);

  if (model_type == "manning upwind") {
    S.Require<CompositeVector, CompositeVectorSpace>(flux, tag)
      .SetGhosted()
      ->SetMesh(S.GetMesh(domain))
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    return Teuchos::rcp(new UpwindTotalFlux(pkname, tag, flux, flux_eps));

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
    S.Require<CompositeVector, CompositeVectorSpace>(flux, tag)
      .SetGhosted()
      ->SetMesh(S.GetMesh(domain))
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    return Teuchos::rcp(new UpwindFluxHarmonicMean(pkname, tag, flux, flux_eps));

  } else if (model_type == "manning split denominator") {
    Key slope = Keys::readKey(oplist, domain, "slope", "slope_magnitude");
    Key manning_coef = Keys::readKey(oplist, domain, "coefficient", "manning_coefficient");
    double slope_regularization = oplist.get<double>("slope regularization epsilon", 1.e-2);

    S.Require<CompositeVector, CompositeVectorSpace>(flux, tag)
      .SetGhosted()
      ->SetMesh(S.GetMesh(domain))
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    S.Require<CompositeVector, CompositeVectorSpace>(slope, tag)
      .SetGhosted()
      ->SetMesh(S.GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S.Require<CompositeVector, CompositeVectorSpace>(manning_coef, tag)
      .SetGhosted()
      ->SetMesh(S.GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    return Teuchos::rcp(new UpwindFluxSplitDenominator(
      pkname, tag, flux, slope, manning_coef, flux_eps, slope_regularization));

  // } else if (model_type == "manning elevation stabilized") {
  //   Key slope = Keys::readKey(oplist, domain, "slope", "slope_magnitude");
  //   Key manning_coef = Keys::readKey(oplist, domain, "coefficient", "manning_coefficient");
  //   Key ponded_depth = Keys::readKey(oplist, domain, "ponded depth", "ponded_depth");
  //   Key elev = Keys::readKey(oplist, domain, "elevation", "elevation");
  //   Key dens = Keys::readKey(oplist, domain, "molar density liquid", "molar_density_liquid");
  //   double slope_regularization = oplist.get<double>("slope regularization epsilon", 1.e-2);
  //   double manning_exp = oplist.get<double>("Manning exponent", 2.0 / 3);

  //   S.Require<CompositeVector, CompositeVectorSpace>(slope, tag)
  //     .SetGhosted()
  //     ->SetMesh(S.GetMesh(domain))
  //     ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  //   S.Require<CompositeVector, CompositeVectorSpace>(elev, tag)
  //     .SetGhosted()
  //     ->SetMesh(S.GetMesh(domain))
  //     ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  //   S.Require<CompositeVector, CompositeVectorSpace>(manning_coef, tag)
  //     .SetGhosted()
  //     ->SetMesh(S.GetMesh(domain))
  //     ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  //   // add boundary face components?

  //   return Teuchos::rcp(new UpwindElevationStabilized(pkname,
  //                                                     tag,
  //                                                     slope,
  //                                                     manning_coef,
  //                                                     ponded_depth,
  //                                                     elev,
  //                                                     dens,
  //                                                     slope_regularization,
  //                                                     manning_exp));

  } else if (model_type == "manning ponded depth passthrough") {
    Key slope = Keys::readKey(oplist, domain, "slope", "slope_magnitude");
    Key manning_coef = Keys::readKey(oplist, domain, "coefficient", "manning_coefficient");
    Key elev = Keys::readKey(oplist, domain, "elevation", "elevation");
    double slope_regularization = oplist.get<double>("slope regularization epsilon", 1.e-8);
    double manning_exp = oplist.get<double>("Manning exponent");

    S.Require<CompositeVector, CompositeVectorSpace>(flux, tag)
      .SetGhosted()
      ->SetMesh(S.GetMesh(domain))
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    S.Require<CompositeVector, CompositeVectorSpace>(slope, tag)
      .SetGhosted()
      ->SetMesh(S.GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S.Require<CompositeVector, CompositeVectorSpace>(elev, tag)
      .SetGhosted()
      ->SetMesh(S.GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S.Require<CompositeVector, CompositeVectorSpace>(manning_coef, tag)
      .SetGhosted()
      ->SetMesh(S.GetMesh(domain))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    return Teuchos::rcp(new UpwindFluxFOCont(
      pkname, tag, flux, slope, manning_coef, elev, slope_regularization, manning_exp));

  } else if (model_type == "manning cell centered") {
    return Teuchos::rcp(new UpwindCellCentered(pkname, tag));
  } else {
    Errors::Message msg;
    msg << "Unknown \"upwind type\" value \"" << model_type
        << ",\" must be one of \"manning upwind\", \"manning harmonic mean\", or \"manning ponded "
           "depth passthrough.\"";
    Exceptions::amanzi_throw(msg);
  }
  return Teuchos::null;
}

} // namespace UpwindFactory
} // namespace Operators
} // namespace Amanzi
