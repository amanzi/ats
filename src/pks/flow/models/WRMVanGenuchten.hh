/*
  Copyright 2010-201x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

//! Water Retention Model, van Genuchten curve.

/*!

van Genuchten's water retention curve.

* `"alpha [Pa^-1]`" ``[double]`` van Genuchten's alpha

ONE OF:
* `"n`" ``[double]`` van Genuchten's n
OR
* `"m`" ``[double]`` van Genuchten's m, m = 1 - 1/n
END

* `"residual saturation [-]`" ``[double]`` **0.0**

* `"smoothing interval width [saturation]`" ``[double]`` **0.0**

* `"Mualem exponent l`" ``[double]`` **0.5**

* `"Krel function name`" ``[string]`` **Mualem**  `"Mualem`" or `"Burdine`"

*/

#pragma once

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Key.hh"
#include "Factory.hh"
#include "StateDefs.hh"
#include "EvaluatorModelByMaterial.hh"
#include "EvaluatorModel_CompositeVector.hh"

namespace ATS {
namespace Flow {
namespace Relations {

const double FLOW_WRM_TOLERANCE = 1e-10;


using namespace Amanzi;


template <class InView_type, class OutView_type>
class WRMVanGenuchten {
 public:
  static const int n_dependencies = 1;
  static const std::string name;

  WRMVanGenuchten(Teuchos::ParameterList& plist)
      : alpha_(plist.sublist("model parameters").get<double>("van Genuchten alpha [Pa^-1]")),
        sr_(plist.sublist("model parameters").get<double>("residual saturation [-]", 0.0))
  {
    if (plist.sublist("model parameters").isParameter("van Genuchten m [-]")) {
      m_ = plist.sublist("model parameters").get<double>("van Genuchten m [-]");
      n_ = 1.0 / (1.0 - m_);
    } else {
      n_ = plist.sublist("model parameters").get<double>("van Genuchten n [-]");
      m_ = 1.0 - 1.0/n_;
    }

    sat_key_ = Keys::cleanPListName(plist.name());
    std::string domain = Keys::getDomain(sat_key_);
    pc_key_ = Keys::readKey(plist, domain, "capillary pressure", "capillary_pressure_gas_liq");
  }

  void SetViews(const std::vector<InView_type>& dependency_views,
                const std::vector<OutView_type>& result_views)
  {
    AMANZI_ASSERT(dependency_views.size() == 1);
    AMANZI_ASSERT(result_views.size() == 1);

    sat_ = result_views[0];
    pc_ = dependency_views[0];
  }

  KeyVector my_keys() const { return { sat_key_ }; }
  KeyVector dependencies() const { return { pc_key_ }; }

  // the model
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    sat_(i,0) =  pc_(i,0) <= 0.0
                 ? 1.0
                 : pow(1.0 + pow(alpha_*pc_(i,0), n_), -m_) * (1.0 - sr_) + sr_;
  }

  // derivatives
  //
  // NOTE: the order of these function tags (i.e. Deriv<0>, ...) is set by the
  // above call to dependencies().  Deriv<I> must correspond to the derivative
  // with respect to dependencies()[I];

  // d/dB
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const
  {
    sat_(i,0) =  pc_(i,0) <= 0.0
                 ? 0.0
                 : -m_*n_ * pow(1.0 + pow(alpha_*pc_(i,0), n_), -m_-1.0) * pow(alpha_*pc_(i,0), n_-1) * alpha_ * (1.0 - sr_);
  }

 private:
  OutView_type sat_;
  InView_type pc_;

  Key sat_key_;
  Key pc_key_;

  double alpha_, sr_, n_, m_;

  static Utils::RegisteredFactory<Evaluator, EvaluatorModelByMaterial<WRMVanGenuchten>> by_material_reg_;
  static Utils::RegisteredFactory<Evaluator, EvaluatorModel_CompositeVector<WRMVanGenuchten>> global_reg_;

};


template <class InView_type, class OutView_type>
class WRMVanGenuchten_Kr {
 public:
  static const int n_dependencies = 1;
  static const std::string name;

  WRMVanGenuchten_Kr(Teuchos::ParameterList& plist)
      : alpha_(plist.sublist("model parameters").get<double>("van Genuchten alpha [Pa^-1]")),
        sr_(plist.sublist("model parameters").get<double>("residual saturation [-]", 0.0))
  {
    if (plist.sublist("model parameters").isParameter("van Genuchten m [-]")) {
      m_ = plist.sublist("model parameters").get<double>("van Genuchten m [-]");
      n_ = 1.0 / (1.0 - m_);
    } else {
      n_ = plist.sublist("model parameters").get<double>("van Genuchten n [-]");
      m_ = 1.0 - 1.0/n_;
    }

    kr_key_ = Keys::cleanPListName(plist.name());
    std::string domain = Keys::getDomain(kr_key_);
    sat_key_ = Keys::readKey(plist, domain, "saturation", "saturation_liquid");
  }

  void SetViews(const std::vector<InView_type>& dependency_views,
                const std::vector<OutView_type>& result_views)
  {
    AMANZI_ASSERT(dependency_views.size() == 1);
    AMANZI_ASSERT(result_views.size() == 1);

    kr = result_views[0];
    sat = dependency_views[0];
  }

  KeyVector my_keys() const { return { kr_key_ }; }
  KeyVector dependencies() const { return { sat_key_ }; }

  // the model
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    double se = (sat(i,0) - sr_)/(1-sr_);
    kr(i,0) =  sat(i,0) < 1.0
                          ? sqrt(se) * pow(1.0 - pow(1.0 - pow(se, 1.0/m_), m_), 2.0)
                          : 1.0;
  }

  // derivatives
  //
  // NOTE: the order of these function tags (i.e. Deriv<0>, ...) is set by the
  // above call to dependencies().  Deriv<I> must correspond to the derivative
  // with respect to dependencies()[I];

  // d/dB
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const
  {
    double se = (sat(i,0) - sr_)/(1-sr_);
    if (sat(i,0) < 1.0) {
      kr(i,0) =  0.0;
    } else {
      double x = pow(se, 1.0 / m_);
      if (fabs(1.0 - x) < FLOW_WRM_TOLERANCE) {
        kr(i,0) = 0.;
      }  else {
        double y = pow(1.0 - x, m_);
        double dkdse = (1.0 - y) * (0.5 * (1.0 - y) + 2 * x * y / (1.0 - x)) * pow(se, -0.5);
        kr(i,0) = dkdse / (1 - sr_);
      }
    }
  }

 private:
  OutView_type kr;
  InView_type sat;

  Key sat_key_;
  Key kr_key_;

  double alpha_, sr_, n_, m_;

  static Utils::RegisteredFactory<Evaluator, EvaluatorModelByMaterial<WRMVanGenuchten_Kr>> by_material_reg_;
  static Utils::RegisteredFactory<Evaluator, EvaluatorModel_CompositeVector<WRMVanGenuchten_Kr>> global_reg_;

};




} // namespace Relations
} // namespace Flow
} // namespace ATS
