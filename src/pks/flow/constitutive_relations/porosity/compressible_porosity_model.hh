/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A simple model for allowing porosity to vary with pressure.
/*!

Based on a linear increase, i.e.

.. math::

   \phi = \phi_{base} + H(p - p_{atm}) * \alpha

where :math:`H` is the heaviside function and :math:`\alpha` is the provided
compressibility.  If the inflection point is set to zero, the above function is
exact.  However, then the porosity function is not smooth (has discontinuous
derivatives), so the inflection point smooths this with a quadratic that
matches the value and derivative at the inflection point and is 0 with 0 slope
at atmospheric pressure.

.. _compressible-porosity-standard-model-spec:
.. admonition:: compressible-porosity-standard-model-spec

   * `"region`" ``[string]`` Region on which this is applied.
   * `"pore compressibility [Pa^-1]`" ``[double]``  :math:`\alpha` as described above
   * `"pore compressibility inflection point [Pa]`" ``[double]`` **1000**

  The inflection point above which the function is linear.

Example:

.. code-block:: xml

  <ParameterList name="soil" type="ParameterList">
    <Parameter name="region" type="string" value="soil" />
    <Parameter name="pore compressibility [Pa^-1]" type="double" value="1.e-9" />
    <Parameter name="pore compressibility inflection point [Pa]" type="double" value="1000." />
  </ParameterList>

*/

#ifndef AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_MODEL_HH_
#define AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_MODEL_HH_

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"

namespace Amanzi {
namespace Flow {

class CompressiblePorosityModel {
 public:
  explicit CompressiblePorosityModel(Teuchos::ParameterList& plist) : plist_(plist)
  {
    InitializeFromPlist_();
  }

  double Porosity(double base_poro, double pres, double patm)
  {
    double poro = base_poro;
    double p_over = pres - patm;
    if (p_over > cutoff_) {
      poro = base_poro + compressibility_ * (cutoff_ / 2. + (p_over - cutoff_));
    } else if (p_over > 0.) {
      poro = base_poro + compressibility_ * (std::pow(p_over, 2.) / 2. / cutoff_);
    }

    if (max_is_one_)
      return poro > 1. ? 1. : poro;
    else
      return poro;
  }

  double DPorosityDPressure(double base_poro, double pres, double patm)
  {
    double p_over = pres - patm;
    double poro = Porosity(base_poro, pres, patm);
    if (max_is_one_ && poro == 1.) {
      return 0.;
    } else if (p_over > cutoff_) {
      return compressibility_;
    } else if (p_over > 0.) {
      return compressibility_ * p_over / cutoff_;
    }

    return 0.;
  }

  double DPorosityDBasePorosity(double base_poro, double pres, double patm)
  {
    if (max_is_one_)
      return pres > patm ? (Porosity(base_poro, pres, patm) > 1.0 ? 0. : 1.) : 1.;
    else
      return 1;
  }

 protected:
  void InitializeFromPlist_()
  {
    compressibility_ = plist_.get<double>("pore compressibility [Pa^-1]");
    cutoff_ = plist_.get<double>("pore compressibility inflection point [Pa]", 1000.);
    max_is_one_ = plist_.get<bool>("cap porosity at 1", true);
  }

 protected:
  Teuchos::ParameterList plist_;
  double compressibility_;
  double cutoff_;
  bool max_is_one_;
};

} // namespace Flow
} // namespace Amanzi

#endif
