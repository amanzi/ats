/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

The Leinjnse model is an exponential model of porosity as a function of
pressure, based on (insert citation!):

.. math:
   p_\lambda = p - p_\text{ref}
   \phi = 1 - (1-\phi_\text{base}) * (exp(- \alpha (p_\lambda - \delta)) - 0.5 \alpha \delta), p_\lambda > \delta
   \phi = 1 - (1-\phi_\text{base}) * (1 - 0.5 \alpha / \delta p_\lambda^2)

where :math:`\alpha` is the provided
compressibility, and :math:`\delta` is the cutoff (inflection point).

If the inflection point is set to zero, the above function is exact.  However,
then the porosity function is not smooth (has discontinuous derivatives).

.. _compressible-porosity-leijnse-model-spec:
.. admonition:: compressible-porosity-leijnse-model-spec

   * `"pore compressibility [Pa^-1]`" ``[double]`` :math:`\alpha` as described above
   * `"pore compressibility inflection point [Pa]`" ``[double]`` **1000** The
     inflection point above which the function is linear.

NOTE: additionally the user should provide a parameter in the `EWC
Globalization Delegate`_ to turn Leijnse model ON in the EWC calculations.

.. code-block:: xml

   <Parameter name="porosity leijnse model" type="bool" value="true"/>

*/

#ifndef AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_LEIJNSE_MODEL_HH_
#define AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_LEIJNSE_MODEL_HH_

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"

namespace Amanzi {
namespace Flow {

class CompressiblePorosityLeijnseModel {
 public:
  explicit CompressiblePorosityLeijnseModel(Teuchos::ParameterList& plist) : plist_(plist)
  {
    InitializeFromPlist_();
  }

  double Porosity(double base_poro, double pres, double patm)
  {
    double poro = base_poro;
    double p_over = pres - patm;
    if (p_over > cutoff_) {
      double a1 = -compressibility_ * (p_over - cutoff_);
      double a2 = 0.5 * compressibility_ * cutoff_;
      poro = 1.0 - (1.0 - base_poro) * (exp(a1) - a2);
    } else if (p_over > 0.) {
      double b1 = 1 - compressibility_ / (2 * cutoff_) * std::pow(p_over, 2.0);
      poro = 1.0 - (1.0 - base_poro) * b1;
    }

    //    return poro > 1. ? 1. : poro;
    return poro;
  }

  double DPorosityDPressure(double base_poro, double pres, double patm)
  {
    double p_over = pres - patm;
    double poro = Porosity(base_poro, pres, patm);
    if (poro == 1.0) {
      return 0.;
    } else if (p_over > cutoff_) {
      double a1 = -compressibility_ * (p_over - cutoff_);
      return (compressibility_ * (1 - base_poro) * exp(a1));
    } else if (p_over > 0.) {
      return (1 - base_poro) * compressibility_ * p_over / cutoff_;
    }

    return 0.;
  }

  double DPorosityDBasePorosity(double base_poro, double pres, double patm)
  {
    std::cout
      << "Derivative of Porosity w.r.t base porosity not implemented for Leijnse model \n!!";
    abort();
    return pres > patm ? (Porosity(base_poro, pres, patm) > 1.0 ? 0. : 1.) : 1.;
  }

 protected:
  void InitializeFromPlist_()
  {
    compressibility_ = plist_.get<double>("pore compressibility [Pa^-1]");
    cutoff_ = plist_.get<double>("pore compressibility inflection point [Pa]", 1000.);
  }

 protected:
  Teuchos::ParameterList plist_;
  double compressibility_;
  double cutoff_;
};

} // namespace Flow
} // namespace Amanzi

#endif
