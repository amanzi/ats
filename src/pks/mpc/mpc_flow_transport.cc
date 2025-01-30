/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "pk_helpers.hh"
#include "mpc_flow_transport.hh"

namespace Amanzi {

void
MPCFlowTransport::parseParameterList()
{
  // are we doing surface, subsurface, or integrated transport?
  auto flow_plist = getSubPKPlist_(0);
  bool surface(false), subsurface(false);
  bool chemistry(false);
  Key surf_lwc_key, sub_lwc_key;

  if (flow_plist->isParameter("domain name")) {
    std::string domain = Keys::readDomain(*flow_plist);
    if (Keys::in(domain, "surface")) {
      surface = true;
      surf_lwc_key = Keys::readKey(*getSubPKPlist_(1), domain, "water content", "water_content");
    } else {
      subsurface = true;
      sub_lwc_key = Keys::readKey(*getSubPKPlist_(1), domain, "water content", "water_content");
    }
  } else {
    // no domain name means it is an MPC, likely coupled water, but maybe permafrost
    surface = true;
    subsurface = true;
    auto transport_names = getSubPKPlist_(1)->get<Teuchos::Array<std::string>>("PKs order");

    // if reactive_transport_names[0] has a domain, it is transport only.
    // otherwise it is reactive transport, and we have to go another level down
    AMANZI_ASSERT(transport_names.size() == 2);

    if (!pks_list_->sublist(transport_names[0]).isParameter("domain name")) {
      chemistry = true;

      // chemistry first
      transport_names = pks_list_->sublist(transport_names[1]).get<Teuchos::Array<std::string>>("PKs order");
    }

    // now we have the actual transport names!
    auto subsurf_domain = Keys::readDomain(pks_list_->sublist(transport_names[0]));
    sub_lwc_key = Keys::readKey(pks_list_->sublist(transport_names[0]), subsurf_domain, "liquid water content", "water_content");

    auto surf_domain = Keys::readDomain(pks_list_->sublist(transport_names[1]));
    surf_lwc_key = Keys::readKey(pks_list_->sublist(transport_names[1]), surf_domain, "liquid water content", "water_content");
  }

  auto [flow_current_tag, flow_next_tag] = tags_[0];
  auto [transport_current_tag, transport_next_tag] = tags_[1];

  if (subsurface) {
    if (transport_next_tag != flow_next_tag) {
      // set the flow field evaluator as the flow's NEXT tag
      //
      // Note, we could be more careful here and readKey() the flow field's name
      // from the flow PK's sublist (which may be nested two deep).  Instead we
      // hard-code this as the default.  If this breaks in the future it can be
      // fixed. --ETC
      Teuchos::ParameterList& flux_list = S_->GetEvaluatorList(Keys::getKey("water_flux", transport_next_tag));
      if (!flux_list.isParameter("evaluator type")) {
        flux_list.set<std::string>("evaluator type", "alias");
        flux_list.set<std::string>("target", Keys::getKey("water_flux", flow_next_tag, true));
      }

      // velocity for dispersivity
      Teuchos::ParameterList& velo_list = S_->GetEvaluatorList(Keys::getKey("darcy_velocity", transport_next_tag));
      if (!velo_list.isParameter("evaluator type")) {
        velo_list.set<std::string>("evaluator type", "alias");
        velo_list.set<std::string>("target", Keys::getKey("darcy_velocity", flow_next_tag, true));
      }

      // now set the liquid water content as an interpolated field at next
      // note that current copy is kept by transport PK
      Teuchos::ParameterList& lwc_list_next = S_->GetEvaluatorList(Keys::getKey(sub_lwc_key, transport_next_tag));
      if (!lwc_list_next.isParameter("evaluator type")) {
        lwc_list_next.set<std::string>("evaluator type", "temporal interpolation");
        lwc_list_next.set<std::string>("current tag", flow_current_tag.get());
        lwc_list_next.set<std::string>("next tag", flow_next_tag.get());
      }

      // porosity used with velocity to compute particle velocity when dispersion is on
      Teuchos::ParameterList& poro_list_next = S_->GetEvaluatorList(Keys::getKey("porosity", transport_next_tag));
      if (!poro_list_next.isParameter("evaluator type")) {
        poro_list_next.set<std::string>("evaluator type", "temporal interpolation");
        poro_list_next.set<std::string>("current tag", flow_current_tag.get());
        poro_list_next.set<std::string>("next tag", flow_next_tag.get());
      }

      // chemistry uses (independently) density, saturation, and porosity
      if (chemistry) {
        Teuchos::ParameterList& dens_list_current = S_->GetEvaluatorList(Keys::getKey("mass_density_liquid", transport_current_tag));
        if (!dens_list_current.isParameter("evaluator type")) {
          dens_list_current.set<std::string>("evaluator type", "temporal interpolation");
          dens_list_current.set<std::string>("current tag", flow_current_tag.get());
          dens_list_current.set<std::string>("next tag", flow_next_tag.get());
        }

        Teuchos::ParameterList& sat_list_current = S_->GetEvaluatorList(Keys::getKey("saturation_liquid", transport_current_tag));
        if (!sat_list_current.isParameter("evaluator type")) {
          sat_list_current.set<std::string>("evaluator type", "temporal interpolation");
          sat_list_current.set<std::string>("current tag", flow_current_tag.get());
          sat_list_current.set<std::string>("next tag", flow_next_tag.get());
        }

        Teuchos::ParameterList& poro_list_current = S_->GetEvaluatorList(Keys::getKey("porosity", transport_current_tag));
        if (!poro_list_current.isParameter("evaluator type")) {
          poro_list_current.set<std::string>("evaluator type", "temporal interpolation");
          poro_list_current.set<std::string>("current tag", flow_current_tag.get());
          poro_list_current.set<std::string>("next tag", flow_next_tag.get());
        }
      }
    }
  }

  if (surface) {
    if (transport_next_tag != flow_next_tag) {
      // set the flow field evaluator as the flow's NEXT tag
      //
      // Note, we could be more careful here and readKey() the flow field's name
      // from the flow PK's sublist (which may be nested two deep).  Instead we
      // hard-code this as the default.  If this breaks in the future it can be
      // fixed. --ETC
      Teuchos::ParameterList& flux_list = S_->GetEvaluatorList(Keys::getKey("surface-water_flux", transport_next_tag));
      if (!flux_list.isParameter("evaluator type")) {
        flux_list.set<std::string>("evaluator type", "alias");
        flux_list.set<std::string>("target", Keys::getKey("surface-water_flux", flow_next_tag, true));
      }

      // velocity for dispersivity
      Teuchos::ParameterList& velo_list = S_->GetEvaluatorList(Keys::getKey("surface-water_velocity", transport_next_tag));
      if (!velo_list.isParameter("evaluator type")) {
        velo_list.set<std::string>("evaluator type", "alias");
        velo_list.set<std::string>("target", Keys::getKey("surface-water_velocity", flow_next_tag, true));
      }

      // set the liquid water content as an interpolated field
      // note that current copy is kept by transport PK
      Teuchos::ParameterList& lwc_list_next = S_->GetEvaluatorList(Keys::getKey(surf_lwc_key, transport_next_tag));
      if (!lwc_list_next.isParameter("evaluator type")) {
        lwc_list_next.set<std::string>("evaluator type", "temporal interpolation");
        lwc_list_next.set<std::string>("current tag", flow_current_tag.get());
        lwc_list_next.set<std::string>("next tag", flow_next_tag.get());
      }

      // porosity used with velocity to compute particle velocity when dispersion is on
      Teuchos::ParameterList& poro_list_next = S_->GetEvaluatorList(Keys::getKey("surface-porosity", transport_next_tag));
      if (!poro_list_next.isParameter("evaluator type")) {
        poro_list_next.set<std::string>("evaluator type", "temporal interpolation");
        poro_list_next.set<std::string>("current tag", flow_current_tag.get());
        poro_list_next.set<std::string>("next tag", flow_next_tag.get());
      }

      // chemistry uses (independently) density, saturation, and porosity at the current tag
      if (chemistry) {
        Teuchos::ParameterList& dens_list_current = S_->GetEvaluatorList(Keys::getKey("surface-mass_density_liquid", transport_current_tag));
        if (!dens_list_current.isParameter("evaluator type")) {
          dens_list_current.set<std::string>("evaluator type", "temporal interpolation");
          dens_list_current.set<std::string>("current tag", flow_current_tag.get());
          dens_list_current.set<std::string>("next tag", flow_next_tag.get());
        }

        Teuchos::ParameterList& sat_list_current = S_->GetEvaluatorList(Keys::getKey("surface-ponded_depth", transport_current_tag));
        if (!sat_list_current.isParameter("evaluator type")) {
          sat_list_current.set<std::string>("evaluator type", "temporal interpolation");
          sat_list_current.set<std::string>("current tag", flow_current_tag.get());
          sat_list_current.set<std::string>("next tag", flow_next_tag.get());
        }

        Teuchos::ParameterList& poro_list_current = S_->GetEvaluatorList(Keys::getKey("surface-porosity", transport_current_tag));
        if (!poro_list_current.isParameter("evaluator type")) {
          poro_list_current.set<std::string>("evaluator type", "temporal interpolation");
          poro_list_current.set<std::string>("current tag", flow_current_tag.get());
          poro_list_current.set<std::string>("next tag", flow_next_tag.get());
        }
      }

    }
  }

  MPCSubcycled::parseParameterList();

  if (subsurface) {
    // for interpolated evaluators, we must first require key@NEXT -- otherwise
    // it will get called to be an alias key@TRANSPORT_NEXT, which is an
    // interpolator that depends upon this.
    requireAtNext(sub_lwc_key, Tags::NEXT, *S_);
    requireAtNext("porosity", Tags::NEXT, *S_);
    if (chemistry) {
      requireAtNext("mass_density_liquid", Tags::NEXT, *S_);
      requireAtNext("saturation_liquid", Tags::NEXT, *S_);
    }

    // now require key@transport_next, which will be the interpolant
    requireAtNext(sub_lwc_key, transport_next_tag, *S_);
    requireAtNext("porosity", transport_next_tag, *S_);
    if (chemistry) {
      requireAtCurrent("mass_density_liquid", transport_current_tag, *S_);
      requireAtCurrent("saturation_liquid", transport_current_tag, *S_);
      requireAtCurrent("porosity", transport_current_tag, *S_);
    }
  }

  if (surface) {
    // for interpolated evaluators, we must first require key@NEXT -- otherwise
    // it will get called to be an alias key@TRANSPORT_NEXT, which is an
    // interpolator that depends upon this.
    requireAtNext(surf_lwc_key, Tags::NEXT, *S_);
    requireAtNext("surface-porosity", Tags::NEXT, *S_);
    if (chemistry) {
      requireAtNext("surface-mass_density_liquid", Tags::NEXT, *S_);
      requireAtNext("surface-ponded_depth", Tags::NEXT, *S_);
    }

    // now require key@transport_next, which will be the interpolant
    requireAtNext(surf_lwc_key, transport_next_tag, *S_);
    requireAtNext("surface-porosity", transport_next_tag, *S_);
    if (chemistry) {
      requireAtCurrent("surface-mass_density_liquid", transport_current_tag, *S_);
      requireAtCurrent("surface-ponded_depth", transport_current_tag, *S_);
      requireAtCurrent("surface-porosity", transport_current_tag, *S_);
    }
  }
}


} // namespace Amanzi
