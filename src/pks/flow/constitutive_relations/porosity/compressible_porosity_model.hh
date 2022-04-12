/* -*-  mode: c++; indent-tabs-mode: nil -*- */

//! A simple model for allowing porosity to vary with pressure.

/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.


  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!

Based on a linear increase, i.e.

$$ \phi = \phi_{base} + H(p - p_{atm}) * \alpha $$

where $H$ is the heaviside function and $\alpha$ is the provided
compressibility.  If the inflection point is set to zero, the above function
is exact.  However, then the porosity function is not smooth (has
discontinuous derivatives).

* `"pore compressibility [Pa^-1]`" ``[double]``  $\alpha$ as described above
* `"pore compressibility inflection point [Pa]`" ``[double]`` **1000**
* `"region`" ``[string]`` Region on which this is applied.

  The inflection point above which the function is linear.

Example:

.. code-block::xml

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
  explicit
  CompressiblePorosityModel(Teuchos::ParameterList& plist) :
      plist_(plist) {
    InitializeFromPlist_();
  }

  double Porosity(double base_poro, double pres, double patm) {
    double poro = base_poro;
    double p_over = pres - patm;
    if (p_over > cutoff_) {
      poro = base_poro + compressibility_ * ( cutoff_ / 2. + (p_over - cutoff_));
    } else if (p_over > 0.) {
      poro = base_poro + compressibility_ * (std::pow(p_over,2.) / 2. / cutoff_);
    }

    return poro > 1. ? 1. : poro;
  }

  double DPorosityDPressure(double base_poro, double pres, double patm) {
    double p_over = pres - patm;
    double poro = Porosity(base_poro, pres, patm);
    if (poro == 1.) {
      return 0.;
    } else if (p_over > cutoff_) {
      return compressibility_;
    } else if (p_over > 0.) {
      return compressibility_ * p_over / cutoff_;
    }

    return 0.;
  }

  double DPorosityDBasePorosity(double base_poro, double pres, double patm) {
    return pres > patm ? (Porosity(base_poro, pres, patm) > 1.0 ? 0. : 1.) : 1.;
  }

 protected:
  void InitializeFromPlist_() {
    compressibility_ = plist_.get<double>("pore compressibility [Pa^-1]");
    cutoff_ = plist_.get<double>("pore compressibility inflection point [Pa]", 1000.);
  }

 protected:

  Teuchos::ParameterList plist_;
  double compressibility_;
  double cutoff_;

};

} // namespace
} // namespace

#endif
