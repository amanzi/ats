/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Bo Gao (gaob@ornl.gov)
*/

/*!
  Sakagucki-Zeng soil resistance model refered to Sakaguchi and Zeng (2009).
  Note that `"dessicated_zone_thickness`" is given by soil types.
  If it is not declared in `"WRM paramters`" through `"model parameters`"
  under state, default value 0.1 m is used for all soil types.
*/

#pragma once

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "StateDefs.hh"
#include "State.hh"
#include "Mesh.hh"
#include "MeshHelpers.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

template <class cView_type, class View_type>
class SoilResistanceSakaguckiZengModel {
 public:
  static const int n_results = 1;
  static const int n_dependencies = 2;
  static const bool provides_derivatives = false;
  static const std::string eval_type;

  explicit SoilResistanceSakaguckiZengModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    my_key_ = { Keys::cleanPListName(*plist), Tag{ plist->get<std::string>("tag") } };

    Key domain = Keys::getDomain(my_key_.first);
    Key domain_sub = Keys::readDomainHint(*plist, domain, "surface", "domain");

    // dependencies
    sat_gas_key_ = Keys::readKeyTag(*plist, domain_sub, "gas saturation", "saturation_gas", my_key_.second);
    poro_key_ = Keys::readKeyTag(*plist, domain_sub, "porosity", "porosity", my_key_.second);

    Teuchos::ParameterList& model_list = plist->sublist("model parameters");
    d_ = model_list.get<double>("dessicated zone thickness [m]", 0.1);
    sr_ = model_list.get<double>("residual saturation [-]", 0.0);
    if (model_list.get<std::string>("wrm type") == "van Genuchten") {
      if (model_list.isParameter("van Genuchten m [-]")) {
        double m = model_list.get<double>("van Genuchten m [-]");
        double n = 1.0 / (1.0 - m);
        double lambda = (n - 1) * (1 - std::pow(0.5, n / (n - 1)));
        b_ = 1. / lambda;
      } else {
        double n = model_list.get<double>("van Genuchten n [-]");
        double lambda = (n - 1) * (1 - std::pow(0.5, n / (n - 1)));
        b_ = 1. / lambda;
      }
    } else if (model_list.get<std::string>("wrm type") == "Brooks-Corey") {
      double lambda = model_list.get<double>("Brooks-Corey lambda [-]");
      b_ = 1. / lambda;
    } else {
      b_ = model_list.get<double>("Clapp-Hornberger b [-]");
    }
  }


  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    AMANZI_ASSERT(deps.size() == n_dependencies);
    AMANZI_ASSERT(res.size() == n_results);

    res_ = res[0];
    sat_gas_ = deps[0];
    poro_ = deps[1];

    auto mesh = s.GetMesh(Keys::getDomain(my_key_.first));
    mesh_surf_ = mesh->getCache();
    mesh_sub_ = mesh->getParentMesh()->getCache();
  }

  void freeViews()
  {
    res_ = View_type();
    sat_gas_ = cView_type();
    poro_ = cView_type();
  }

  KeyTagVector getMyKeys() const
  {
    return {
      my_key_,
    };
  }

  KeyTagVector getDependencies() const
  {
    return {
      sat_gas_key_, poro_key_
    };
  }

  KOKKOS_INLINE_FUNCTION void operator()(const int sc) const
  {
    AmanziMesh::Entity_ID f = mesh_surf_.getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
    AmanziMesh::Entity_ID c = AmanziMesh::getFaceOnBoundaryInternalCell(mesh_sub_, f);

    double r_soil;
    if (sat_gas_(c, 0) == 0.) {
      r_soil = 0.; // ponded water
    } else {
      double vp_diffusion = 2.2e-5 * Kokkos::pow(poro_(c, 0), 2) * Kokkos::pow(1 - sr_, 2 + 3 * b_);
      double L_Rsoil = d_ * (Kokkos::exp(Kokkos::pow(sat_gas_(c, 0), 5)) - 1) / (Kokkos::exp(1) - 1);
      r_soil = L_Rsoil / vp_diffusion;
    }
    res_(sc, 0) = r_soil;
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const { assert(false); }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const { assert(false); }

 private:
  View_type res_;
  cView_type poro_, sat_gas_;

  double sr_, b_, d_;

  AmanziMesh::MeshCache mesh_surf_;
  AmanziMesh::MeshCache mesh_sub_;

  KeyTag my_key_;
  KeyTag sat_gas_key_, poro_key_;
};


template <class cView_type, class View_type>
const std::string SoilResistanceSakaguckiZengModel<cView_type, View_type>::eval_type =
  "soil resistance, Sakagucki-Zeng";

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
