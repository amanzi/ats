/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (berndt@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include "Kokkos_MathematicalFunctions.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "Spline.hh"

#include "wrm_van_genuchten.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


/* ******************************************************************
 * Setup fundamental parameters for this model.
 ****************************************************************** */
WRMVanGenuchten::WRMVanGenuchten(Teuchos::ParameterList& plist)
{
  std::string fname = plist.get<std::string>("Krel function name", "Mualem");
  if (fname == std::string("Mualem")) {
    function_ = RelPermFunction_kind::MUALEM;
  } else if (fname == std::string("Burdine")) {
    function_ = RelPermFunction_kind::BURDINE;
  } else {
    AMANZI_ASSERT(0);
  }

  // DEPRECATION ERROR
  if (plist.isParameter("van Genuchten alpha")) {
    Errors::Message message(
      "WRM: DEPRECATION: parameter \"van Genuchten alpha\" is now \"van Genuchten alpha [Pa^-1]\"");
    Exceptions::amanzi_throw(message);
  }
  if (plist.isParameter("van Genuchten n")) {
    Errors::Message message(
      "WRM: DEPRECATION: parameter \"van Genuchten n\" is now \"van Genuchten n [-]\"");
    Exceptions::amanzi_throw(message);
  }
  if (plist.isParameter("van Genuchten m")) {
    Errors::Message message(
      "WRM: DEPRECATION: parameter \"van Genuchten m\" is now \"van Genuchten m [-]\"");
    Exceptions::amanzi_throw(message);
  }
  if (plist.isParameter("residual saturation")) {
    Errors::Message message(
      "WRM: DEPRECATION: parameter \"residual saturation\" is now \"residual saturation [-]\"");
    Exceptions::amanzi_throw(message);
  }
  if (plist.isParameter("Mualem exponent l")) {
    Errors::Message message(
      "WRM: DEPRECATION: parameter \"Mualem exponent l\" is now \"Mualem exponent l [-]\"");
    Exceptions::amanzi_throw(message);
  }
  if (plist.isParameter("smoothing interval width")) {
    Errors::Message message("WRM: DEPRECATION: option \"smoothing interval width\" is now "
                            "\"smoothing interval width [saturation]\"");
    Exceptions::amanzi_throw(message);
  }
  if (plist.isParameter("van Genuchten residual saturation")) {
    Errors::Message message("WRM: DEPRECATION: option \"van Genuchten residual saturation\" is now "
                            "\"residual saturation [-]\"");
    Exceptions::amanzi_throw(message);
  }

  alpha_ = plist.get<double>("van Genuchten alpha [Pa^-1]");
  sr_ = plist.get<double>("residual saturation [-]");
  l_ = plist.get<double>("Mualem exponent l [-]", 0.5);

  // map to n,m
  if (plist.isParameter("van Genuchten m [-]")) {
    m_ = plist.get<double>("van Genuchten m [-]");
    if (function_ == RelPermFunction_kind::MUALEM) {
      n_ = 1.0 / (1.0 - m_);
    } else {
      n_ = 2.0 / (1.0 - m_);
    }
  } else {
    n_ = plist.get<double>("van Genuchten n [-]");
    if (function_ == RelPermFunction_kind::MUALEM) {
      m_ = 1.0 - 1.0 / n_;
    } else {
      m_ = 1.0 - 2.0 / n_;
    }
  }

  s0_ = 1.0 - plist.get<double>("smoothing interval width [saturation]", 0.0);
  if (s0_ < 1.) { fit_kr_.Setup(s0_, k_relative(s0_), d_k_relative(s0_), 1.0, 1.0, 0.0); }

  pc0_ = plist.get<double>("saturation smoothing interval [Pa]", 0.0);
  if (pc0_ > 0.) { fit_s_.Setup(0.0, 1.0, 0.0, pc0_, saturation(pc0_), d_saturation(pc0_)); }
};


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
