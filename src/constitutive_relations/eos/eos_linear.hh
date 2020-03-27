/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)

*/

//! EOS that is linear in pressure.

/*!

An Equation of State that is linear in pressure, i.e.

.. math::
    \rho = \rho_0 + \beta \cdot max(p - p_atm, 0)

Parameters are:

* `"atmospheric pressure [Pa]`" ``[double]`` **101325.0**
* `"density [kg m^-3]`" ``[double]`` **1000.0**
* `"compressibility [Pa^-1]`" ``[double]``
  
*/

#pragma once

#include "Teuchos_ParameterList.hpp"
#include "Keys.hh"

namespace ATS {
namespace Relations {

using namespace Amanzi;


template <class cView_type, class View_type>
class EOSLinearMass {
 public:
  static const int n_dependencies = 1;
  static const std::string name;
  static const Utils::RegisteredFactory<Evaluator,EvaluatorModel_CompositeVector<EOSLinear> > fac_cv_;
  // static Utils::RegisteredFactory<Evaluator,EvaluatorModel_MultiPatch<EOSLinear> > fac_patch_;

  
  EOSLinear(OutputVector_type<DeviceType> rho,
            InputVector_type<DeviceType> p,
            Teuchos::ParameterList& plist)
      : plist_(plist)
  {
    rho_key_ = Keys::cleanPListName(plist_.name());
    molar_ = rho_key_.find("molar", 0) != std::string::npos;
    molar_ = plist_.sublist("model parameters").get<bool>("is molar density", molar_);

    pres_key_ = Keys::readKey(plist_, Keys::getDomain(rho_key_), "pressure");
    
    if (molar_) {
      rho0_ = plist_.sublist("model parameters").get<double>("density [mol m^-3]", 55000.);
      patm_ = plist_.sublist("model parameters").get<double>("atmospheric pressure [Pa]", 101325.);
      beta_ = plist_.sublist("model parameters").get<double>("compressibility [Pa^-1]");
    } else {
      rho0_ = plist_.sublist("model parameters").get<double>("density [kg m^-3]", 1000.);
      patm_ = plist_.sublist("model parameters").get<double>("atmospheric pressure [Pa]", 101325.);
      beta_ = plist_.sublist("model parameters").get<double>("compressibility [Pa^-1]");
    }            

    
  }

  void SetViews(const std::vector<cView_type>& dependency_views,
                const std::vector<View_type>& result_views)
  {
    AMANZI_ASSERT(dependency_views.size() == 1);
    AMANZI_ASSERT(result_views.size() == 1);

    rho_ = result_views[0];
    p_ = dependency_views[0];
  }

  KeyVector my_keys() const { return { rho_key_, }; }
  KeyVector dependencies() const { return { pres_key_, }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int& i) const
  {
    rho_(i,0) = rho0_ + beta_ * max(p_(i,0) - p_atm_, 0.);
  }

  class dd1 {};
  KOKKOS_INLINE_FUNCTION void operator()(dAdB, const int& i) const
  {
    rho_(i,0) = p_(i,0) > p_atm_ ? beta_ : 0.;
  }

 private:
  OutputVector_type<DeviceType> rho_;
  InputVector_type<DeviceType> p_;

  Teuchos::ParameterList plist_;
  Key rho_key_;
  Key pres_key_;

  double rho0_;
  double beta_;
  double p_atm_;

  
};

template <class cView_type, class View_type>
const std::string EOSLinear<cView_type, View_type>::name("eos linear");


} // namespace Relations
} // namespace ATS
