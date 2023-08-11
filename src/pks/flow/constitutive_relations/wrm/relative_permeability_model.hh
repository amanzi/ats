/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Uses WRM class instances to compute the scalar portion of the hydraulic conductance.
/*!

This computes the product,

.. math::
   k = \frac{n}{\mu} k_r

which is the scalar portion of the hydraulic conductance (excludes the absolute
permeability).  Relative permeability, k_r, is computed using a water retention
model, which provides generic methods for this quantity based on internal
parameters (e.g. van Genuchten/Mualem, Brooks & Corey, etc.).

While we call this "relative permeability", it is actually the product above,
and might better be called the scalar hydraulic conductance?

type : `"relative permeability`"

.. _relative-permeability-model-spec
.. admonition:: relative-permeability-model-spec

   * `"model parameters`" ``[wrm-spec-list]``

   KEYS
   - `"saturation`" **DOMAIN-saturation_liquid**
   - `"density`" **DOMAIN-molar_density_liquid**
   - `"viscosity`" **DOMAIN-viscosity**  
                
*/

#pragma once

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "StateDefs.hh"
#include "wrm_van_genuchten.hh"

namespace Amanzi {

class State;

namespace Flow {
namespace Relations {

template <class cView_type, class View_type, class WRM_type>
class RelativePermeabilityModel {
 public:
  static const int n_dependencies = 1;
  static const std::string name;

  RelativePermeabilityModel(Teuchos::ParameterList& plist)
    : model_(plist.sublist("model parameters"))
  {
    my_key_ = Keys::cleanPListName(plist);
    auto domain = Keys::getDomain(my_key_);
    s_key_ = Keys::readKey(plist, domain, "saturation", "saturation_liquid");
    dens_key_ = Keys::readKey(plist, domain, "density", "molar_density_liquid");
    visc_key_ = Keys::readKey(plist, domain, "viscosity", "viscosity");

    rescaling_ = 1.0 / plist.get<double>("permeability rescaling", 1.0);
  }

  void setViews(const std::vector<cView_type>& deps,
                const std::vector<View_type>& res,
                const State& s) {
    res_ = res[0];
    s_ = deps[0];
    dens_ = deps[1];
    visc_ = deps[2];
  }

  KeyVector getMyKeys() const { return { my_key_ }; }
  KeyVector getDependencies() const { return { s_key_, dens_key_, visc_key_ }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
    res_(i,0) = rescaling_ * model_.k_relative(s_(i,0)) * dens_(i,0) / visc_(i,0);
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const {
    res_(i,0) = rescaling_ * model_.d_k_relative(s_(i,0)) * dens_(i,0) / visc_(i,0);
  }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const {
    res_(i,0) = rescaling_ * model_.k_relative(s_(i,0)) / visc_(i,0);
  }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<2>, const int i) const {
    res_(i,0) = - rescaling_ * dens_(i,0) * model_.k_relative(s_(i,0)) / pow(visc_(i,0), 2);
  }

 private:
  View_type res_;
  cView_type s_, dens_, visc_;
  Key my_key_, s_key_, dens_key_, visc_key_;
  double rescaling_;

  WRM_type model_;
};


template <class cView_type, class View_type, class WRM_type>
const std::string RelativePermeabilityModel<cView_type, View_type, WRM_type>::name = "relative permeability " + WRM_type::name;


template<class cView_type, class View_type>
using RelativePermeabilityVanGenuchtenModel = RelativePermeabilityModel<cView_type, View_type, WRMVanGenuchten>;


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
