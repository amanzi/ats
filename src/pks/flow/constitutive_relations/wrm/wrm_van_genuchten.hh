/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! WRMVanGenuchten : water retention model using van Genuchten's parameterization
/*!

van Genuchten's water retention curve.

.. _WRM-van-Genuchten-spec
.. admonition:: WRM-van-Genuchten-spec

    * `"van Genuchten alpha [Pa^-1]`" ``[double]`` van Genuchten's alpha

    ONE OF:

    * `"van Genuchten n [-]`" ``[double]`` van Genuchten's n

    OR

    * `"van Genuchten m [-]`" ``[double]`` van Genuchten's m, m = 1 - 1/n

    END

    * `"residual saturation [-]`" ``[double]`` **0.0**
    * `"smoothing interval width [saturation]`" ``[double]`` **0.0**
    * `"Mualem exponent l [-]`" ``[double]`` **0.5**
    * `"Krel function name`" ``[string]`` **Mualem**  `"Mualem`" or `"Burdine`"

Example:

.. code-block:: xml

    <ParameterList name="moss" type="ParameterList">
      <Parameter name="van Genuchten alpha [Pa^-1]" type="double" value="0.002" />
      <Parameter name="van Genuchten m [-]" type="double" value="0.2" />
      <Parameter name="residual saturation [-]" type="double" value="0.0" />
      <Parameter name="smoothing interval width [saturation]" type="double" value=".05" />
    </ParameterList>

*/

#ifndef ATS_FLOWRELATIONS_WRM_VAN_GENUCHTEN_
#define ATS_FLOWRELATIONS_WRM_VAN_GENUCHTEN_

#include "Teuchos_ParameterList.hpp"
#include "Spline.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

enum class RelPermFunction_kind { MUALEM, BURDINE };
const static double FLOW_WRM_TOLERANCE = 1e-10;


class WRMVanGenuchten {
 public:
  static const std::string eval_type;

  explicit WRMVanGenuchten(Teuchos::ParameterList& plist);

  // required methods from the base class
  KOKKOS_INLINE_FUNCTION double k_relative(double s) const {
    if (s <= s0_) {
      double se = (s - sr_) / (1 - sr_);
      if (function_ == RelPermFunction_kind::MUALEM) {
        return Kokkos::pow(se, l_) * Kokkos::pow(1.0 - Kokkos::pow(1.0 - Kokkos::pow(se, 1.0 / m_), m_), 2.0);
      } else {
        return se * se * (1.0 - Kokkos::pow(1.0 - Kokkos::pow(se, 1.0 / m_), m_));
      }
    } else if (s == 1.0) {
      return 1.0;
    } else {
      return fit_kr_(s);
    }
  }

  KOKKOS_INLINE_FUNCTION double d_k_relative(double s) const {
    if (s <= s0_) {
      double se = (s - sr_) / (1 - sr_);

      double x = Kokkos::pow(se, 1.0 / m_);
      if (fabs(1.0 - x) < FLOW_WRM_TOLERANCE) return 0.0;
      if (fabs(x) < FLOW_WRM_TOLERANCE) return 0.0;

      double y = Kokkos::pow(1.0 - x, m_);
      double dkdse;
      if (function_ == RelPermFunction_kind::MUALEM)
        dkdse = (1.0 - y) * (l_ * (1.0 - y) + 2 * x * y / (1.0 - x)) * Kokkos::pow(se, l_ - 1.0);
      else
        dkdse = (2 * (1.0 - y) + x / (1.0 - x)) * se;

      return dkdse / (1 - sr_);

    } else if (s == 1.0) {
      return 0.0;
    } else {
      return fit_kr_.Derivative(s);
    }
  }

  KOKKOS_INLINE_FUNCTION double saturation(double pc) const {
    if (pc > pc0_) {
      return Kokkos::pow(1.0 + Kokkos::pow(alpha_ * pc, n_), -m_) * (1.0 - sr_) + sr_;
    } else if (pc <= 0.) {
      return 1.0;
    } else {
      return fit_s_(pc);
    }
  }

  KOKKOS_INLINE_FUNCTION double d_saturation(double pc) const {
    if (pc > pc0_) {
      return -m_ * n_ * Kokkos::pow(1.0 + Kokkos::pow(alpha_ * pc, n_), -m_ - 1.0) *
        Kokkos::pow(alpha_ * pc, n_ - 1) * alpha_ * (1.0 - sr_);
    } else if (pc <= 0.) {
      return 0.0;
    } else {
      return fit_s_.Derivative(pc);
    }
  }

  KOKKOS_INLINE_FUNCTION double capillaryPressure(double s) const {
    double se = (s - sr_) / (1.0 - sr_);
    se = Kokkos::min(se, 1.0);
    se = Kokkos::max(se, 1.e-40);
    if (se < 1.e-8) {
      return Kokkos::pow(se, -1.0 / (m_ * n_)) / alpha_;
    } else {
      return (Kokkos::pow(Kokkos::pow(se, -1.0 / m_) - 1.0, 1 / n_)) / alpha_;
    }
  }

  KOKKOS_INLINE_FUNCTION double d_capillaryPressure(double s) const {
    double se = (s - sr_) / (1.0 - sr_);
    se = Kokkos::min(se, 1.0);
    se = Kokkos::max(se, 1.e-40);
    if (se < 1.e-8) {
      return -1.0 / (m_ * n_ * alpha_) * Kokkos::pow(se, -1.0 / (m_ * n_) - 1.) / (1.0 - sr_);
    } else {
      return -1.0 / (m_ * n_ * alpha_) * Kokkos::pow(Kokkos::pow(se, -1.0 / m_) - 1.0, 1 / n_ - 1.0) *
        Kokkos::pow(se, -1.0 / m_ - 1.0) / (1.0 - sr_);
    }
  }

  KOKKOS_INLINE_FUNCTION double residualSaturation() const { return sr_; }

 private:
  double m_; // van Genuchten parameters: m, n, alpha
  double n_;
  double l_;
  double alpha_;
  double sr_; // van Genuchten residual saturation

  RelPermFunction_kind function_;
  double s0_; // regularization threshold in saturation
  Amanzi::Utils::Spline fit_kr_;

  double pc0_;
  Amanzi::Utils::Spline fit_s_;
};

inline const std::string WRMVanGenuchten::eval_type = "van Genuchten";



} // namespace Relations
} // namespace Flow
} // namespace Amanzi

#endif
