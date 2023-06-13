/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A simple model for allowing porosity to vary with pressure.
/*!

Compressibility based on a linear increase, i.e.

.. math::

   \phi = \phi_{base} + H(p - p_{atm}) * \alpha

where :math:`H` is the heaviside function and :math:`\alpha` is the provided
compressibility.  If the inflection point is set to zero, the above function is
exact.  However, then the porosity function is not differentiable, so the
inflection point can be used to smooth this with a quadratic that matches the
value and derivative at the inflection point and matches the value and slope at
atmospheric pressure.

type : `"compressible porosity linear`"

.. _compressible-porosity-linear-model-spec
.. admonition:: compressible-porosity-linear-model-spec

   * `"region`" ``[string]`` Region on which this is applied.
   * `"pore compressibility [Pa^-1]`" ``[double]``  :math:`\alpha` as described above
   * `"pore compressibility inflection point [Pa]`" ``[double]`` **1000**

   KEYS:
   - `"pressure`"
   - `"base porosity`"

  The inflection point above which the function is linear.

Example:

.. code-block:: xml

  <ParameterList name="soil" type="ParameterList">
    <Parameter name="region" type="string" value="soil" />
    <Parameter name="pore compressibility [Pa^-1]" type="double" value="1.e-9" />
    <Parameter name="pore compressibility inflection point [Pa]" type="double" value="1000." />
  </ParameterList>

*/

#pragma once

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "Key.hh"
#include "StateDefs.hh"
#include "State.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

template <class cView_type, class View_type>
class CompressiblePorosityLinearModel {
 public:
  static const int n_dependencies = 2;
  static const std::string name;

  CompressiblePorosityLinearModel(Teuchos::ParameterList& plist) {
    my_key_ = Keys::cleanPListName(plist);
    auto domain = Keys::getDomain(my_key_);
    phi_key_ = Keys::readKey(plist, domain, "porosity", "porosity");
    p_key_ = Keys::readKey(plist, domain, "pressure", "pressure");

    compressibility_ = plist.get<double>("pore compressibility [Pa^-1]");
    cutoff_ = plist.get<double>("pore compressibility inflection point [Pa]", 1000.);
    max_is_one_ = plist.get<bool>("cap porosity at 1", true);
  }

  void setViews(const std::vector<cView_type>& deps,
                const std::vector<View_type>& res,
                const State& s) {
    res_ = res[0];
    phi_ = deps[0];
    p_ = deps[1];
    patm_ = s.Get<double>("atmospheric_pressure", Tags::DEFAULT);
  }

  KeyVector getMyKeys() const { return { my_key_ }; }
  KeyVector getDependencies() const { return { phi_key_, p_key_ }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
    double poro = phi_(i,0);
    double p_over = p_(i,0) - patm_;
    if (p_over > cutoff_) {
      poro = poro + compressibility_ * (cutoff_ / 2. + (p_over - cutoff_));
    } else if (p_over > 0.) {
      poro = poro + compressibility_ * (std::pow(p_over, 2.) / 2. / cutoff_);
    }

    if (max_is_one_) res_(i,0) = poro > 1. ? 1. : poro;
    else res_(i,0) = poro;
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const {
    operator()(i);
    double poro = res_(i,0);

    if (max_is_one_) res_(i,0) = p_(i,0) > patm_ ? (poro > 1.0 ? 0. : 1.) : 1.;
    else res_(i,0) = 1.0;
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const {
    operator()(i);
    double poro = res_(i,0);
    double p_over = p_(i,0) - patm_;
    if (max_is_one_ && poro == 1.) {
      res_(i,0) = 0.;
    } else if (p_over > cutoff_) {
      res_(i,0) = compressibility_;
    } else if (p_over > 0.) {
      res_(i,0) = compressibility_ * p_over / cutoff_;
    } else {
      res_(i,0) = 0.;
    }
  }

 private:
  View_type res_;
  cView_type phi_, p_;
  Key my_key_, phi_key_, p_key_;

  double patm_;

  double compressibility_;
  double cutoff_;
  bool max_is_one_;
};


template <class cView_type, class View_type>
const std::string CompressiblePorosityLinearModel<cView_type, View_type>::name = "compressible porosity linear";


} // namespace Relations
} // namespace Flow
} // namespace Amanzi


