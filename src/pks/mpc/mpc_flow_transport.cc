/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "PK_Helpers.hh"
#include "mpc_flow_transport.hh"

namespace Amanzi {

void
MPCFlowTransport::parseParameterList()
{
  // are we doing surface, subsurface, or integrated transport?
  auto flow_plist = getSubPKPlist_(0);
  Key surf_lwc_key, sub_lwc_key;

  if (flow_plist->isParameter("domain name")) {
    std::string domain = Keys::readDomain(*flow_plist);
    if (Keys::in(domain, "surface")) {
      surface_ = true;
      surf_lwc_key = Keys::readKey(*getSubPKPlist_(1), domain, "water content", "water_content");
    } else {
      subsurface_ = true;
      sub_lwc_key = Keys::readKey(*getSubPKPlist_(1), domain, "water content", "water_content");
    }
  } else {
    // no domain name means it is an MPC, likely coupled water, but maybe permafrost
    surface_ = true;
    subsurface_ = true;
    auto transport_names = getSubPKPlist_(1)->get<Teuchos::Array<std::string>>("PKs order");

    // if reactive_transport_names[0] has a domain, it is transport only.
    // otherwise it is reactive transport, and we have to go another level down
    AMANZI_ASSERT(transport_names.size() == 2);

    if (!pks_list_->sublist(transport_names[0]).isParameter("domain name")) {
      chemistry_ = true;

      // chemistry first
      transport_names =
        pks_list_->sublist(transport_names[1]).get<Teuchos::Array<std::string>>("PKs order");
    }

    // now we have the actual transport names!
    auto subsurf_domain = Keys::readDomain(pks_list_->sublist(transport_names[0]));
    sub_lwc_key = Keys::readKey(pks_list_->sublist(transport_names[0]),
                                subsurf_domain,
                                "liquid water content",
                                "water_content");

    auto surf_domain = Keys::readDomain(pks_list_->sublist(transport_names[1]));
    surf_lwc_key = Keys::readKey(
      pks_list_->sublist(transport_names[1]), surf_domain, "liquid water content", "water_content");
  }

  auto [flow_current_tag, flow_next_tag] = tags_[0];
  auto [transport_current_tag, transport_next_tag] = tags_[1];

  if (subsurface_) {
    if (transport_next_tag != flow_next_tag) {
      // set the flow field evaluator as the flow's NEXT tag
      //
      // Note, we could be more careful here and readKey() the flow field's name
      // from the flow PK's sublist (which may be nested two deep).  Instead we
      // hard-code this as the default.  If this breaks in the future it can be
      // fixed. --ETC
      Teuchos::ParameterList& flux_list =
        S_->GetEvaluatorList(Keys::getKey("water_flux", transport_next_tag));
      if (!flux_list.isParameter("evaluator type")) {
        flux_list.set<std::string>("evaluator type", "alias");
        flux_list.set<std::string>("target", Keys::getKey("water_flux", flow_next_tag, true));
      }

      // velocity for dispersivity
      Teuchos::ParameterList& velo_list =
        S_->GetEvaluatorList(Keys::getKey("darcy_velocity", transport_next_tag));
      if (!velo_list.isParameter("evaluator type")) {
        velo_list.set<std::string>("evaluator type", "alias");
        velo_list.set<std::string>("target", Keys::getKey("darcy_velocity", flow_next_tag, true));
      }

      // now set the liquid water content as an interpolated field at next
      // note that flow_current copy is kept by flow PK, and transport_current copy is kept by transport PK
      Teuchos::ParameterList& lwc_list_next =
        S_->GetEvaluatorList(Keys::getKey(sub_lwc_key, transport_next_tag));
      if (!lwc_list_next.isParameter("evaluator type")) {
        lwc_list_next.set<std::string>("evaluator type", "temporal interpolation");
        lwc_list_next.set<std::string>("current tag", flow_current_tag.get());
        lwc_list_next.set<std::string>("next tag", flow_next_tag.get());
      }

      // porosity used with velocity to compute particle velocity when dispersion is on
      // -- and an interpolation at transport's next
      Teuchos::ParameterList& poro_list_next =
        S_->GetEvaluatorList(Keys::getKey("porosity", transport_next_tag));
      if (!poro_list_next.isParameter("evaluator type")) {
        poro_list_next.set<std::string>("evaluator type", "temporal interpolation");
        poro_list_next.set<std::string>("current tag", flow_current_tag.get());
        poro_list_next.set<std::string>("next tag", flow_next_tag.get());
      }

      // chemistry uses (independently) density, saturation, and porosity
      if (chemistry_) {
        // mass density used by Alquimia
        // -- transport next is an interpolation
        Teuchos::ParameterList& dens_list_current =
          S_->GetEvaluatorList(Keys::getKey("mass_density_liquid", transport_current_tag));
        if (!dens_list_current.isParameter("evaluator type")) {
          dens_list_current.set<std::string>("evaluator type", "temporal interpolation");
          dens_list_current.set<std::string>("current tag", flow_current_tag.get());
          dens_list_current.set<std::string>("next tag", flow_next_tag.get());
        }

        // saturation
        // -- flow's current is done by flow
        // -- transport's next is an interoplation
        Teuchos::ParameterList& sat_list_current =
          S_->GetEvaluatorList(Keys::getKey("saturation_liquid", transport_current_tag));
        if (!sat_list_current.isParameter("evaluator type")) {
          sat_list_current.set<std::string>("evaluator type", "temporal interpolation");
          sat_list_current.set<std::string>("current tag", flow_current_tag.get());
          sat_list_current.set<std::string>("next tag", flow_next_tag.get());
        }

        // -- porosity at transport's current is an interpolation (could also be a copy of transport's next?)
        Teuchos::ParameterList& poro_list_current =
          S_->GetEvaluatorList(Keys::getKey("porosity", transport_current_tag));
        if (!poro_list_current.isParameter("evaluator type")) {
          poro_list_current.set<std::string>("evaluator type", "temporal interpolation");
          poro_list_current.set<std::string>("current tag", flow_current_tag.get());
          poro_list_current.set<std::string>("next tag", flow_next_tag.get());
        }
      }
    }
  }

  if (surface_) {
    if (transport_next_tag != flow_next_tag) {
      // set the flow field evaluator as the flow's NEXT tag
      //
      // Note, we could be more careful here and readKey() the flow field's name
      // from the flow PK's sublist (which may be nested two deep).  Instead we
      // hard-code this as the default.  If this breaks in the future it can be
      // fixed. --ETC
      Teuchos::ParameterList& flux_list =
        S_->GetEvaluatorList(Keys::getKey("surface-water_flux", transport_next_tag));
      if (!flux_list.isParameter("evaluator type")) {
        flux_list.set<std::string>("evaluator type", "alias");
        flux_list.set<std::string>("target",
                                   Keys::getKey("surface-water_flux", flow_next_tag, true));
      }

      // velocity for evaluators
      Teuchos::ParameterList& velo_list =
        S_->GetEvaluatorList(Keys::getKey("surface-velocity", transport_next_tag));
      if (!velo_list.isParameter("evaluator type")) {
        velo_list.set<std::string>("evaluator type", "alias");
        velo_list.set<std::string>("target", Keys::getKey("surface-velocity", flow_next_tag, true));
      }

      // set the liquid water content as an interpolated field
      // -- flow's current is a kept by flow
      // -- transport's next is an interpolation
      Teuchos::ParameterList& lwc_list_next =
        S_->GetEvaluatorList(Keys::getKey(surf_lwc_key, transport_next_tag));
      if (!lwc_list_next.isParameter("evaluator type")) {
        lwc_list_next.set<std::string>("evaluator type", "temporal interpolation");
        lwc_list_next.set<std::string>("current tag", flow_current_tag.get());
        lwc_list_next.set<std::string>("next tag", flow_next_tag.get());
      }

      // chemistry uses (independently) density & ponded depth at the current tag
      if (chemistry_) {
        // density used by Alquimia
        // -- transport's next is an interpolation
        Teuchos::ParameterList& dens_list_current =
          S_->GetEvaluatorList(Keys::getKey("surface-mass_density_liquid", transport_current_tag));
        if (!dens_list_current.isParameter("evaluator type")) {
          dens_list_current.set<std::string>("evaluator type", "temporal interpolation");
          dens_list_current.set<std::string>("current tag", flow_current_tag.get());
          dens_list_current.set<std::string>("next tag", flow_next_tag.get());
        }

        // ponded depth
        // -- flow's current is kept by flow
        // -- transport's next is an interpolation
        Teuchos::ParameterList& pd_list_current =
          S_->GetEvaluatorList(Keys::getKey("surface-ponded_depth", transport_current_tag));
        if (!pd_list_current.isParameter("evaluator type")) {
          pd_list_current.set<std::string>("evaluator type", "temporal interpolation");
          pd_list_current.set<std::string>("current tag", flow_current_tag.get());
          pd_list_current.set<std::string>("next tag", flow_next_tag.get());
        }
      }
    }
  }

  MPCSubcycled::parseParameterList();
  if (subcycling_[0]) {
    Errors::Message msg;
    msg << "In MPCFlowTransport PK \"" << name_ << "\", cannot subcycle flow PK, only transport PK";
    Exceptions::amanzi_throw(msg);
  }

  if (subsurface_) {
    // for interpolated evaluators, we must first require key@NEXT -- otherwise
    // it will get called to be an alias key@TRANSPORT_NEXT, which is an
    // interpolator that depends upon this.
    requireEvaluatorAtNext(sub_lwc_key, flow_next_tag, *S_);
    requireEvaluatorAtNext("porosity", flow_next_tag, *S_);
    requireEvaluatorAtCurrent("porosity", flow_current_tag, *S_, name_);

    if (chemistry_) {
      requireEvaluatorAtNext("mass_density_liquid", flow_next_tag, *S_);
      requireEvaluatorAtCurrent("mass_density_liquid", flow_current_tag, *S_, name_);

      requireEvaluatorAtNext("saturation_liquid", flow_next_tag, *S_);
    }

    // now require key@transport_next, which will be the interpolant
    requireEvaluatorAtNext(sub_lwc_key, transport_next_tag, *S_);
    requireEvaluatorAtNext("porosity", transport_next_tag, *S_);
    if (chemistry_) {
      requireEvaluatorAtCurrent("mass_density_liquid", transport_current_tag, *S_);
      requireEvaluatorAtCurrent("saturation_liquid", transport_current_tag, *S_);
      requireEvaluatorAtCurrent("porosity", transport_current_tag, *S_);
    }
  }

  if (surface_) {
    // for interpolated evaluators, we must first require key@NEXT -- otherwise
    // it will get called to be an alias key@TRANSPORT_NEXT, which is an
    // interpolator that depends upon this.
    requireEvaluatorAtNext(surf_lwc_key, flow_next_tag, *S_);
    if (chemistry_) {
      requireEvaluatorAtNext("surface-mass_density_liquid", flow_next_tag, *S_);
      requireEvaluatorAtCurrent("surface-mass_density_liquid", flow_current_tag, *S_, name_);

      requireEvaluatorAtNext("surface-ponded_depth", flow_next_tag, *S_);
    }

    // now require key@transport_next, which will be the interpolant
    requireEvaluatorAtNext(surf_lwc_key, transport_next_tag, *S_);
    if (chemistry_) {
      requireEvaluatorAtCurrent("surface-mass_density_liquid", transport_current_tag, *S_);
      requireEvaluatorAtCurrent("surface-ponded_depth", transport_current_tag, *S_);
      requireEvaluatorAtCurrent("surface-porosity", transport_current_tag, *S_);
    }
  }
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
void
MPCFlowTransport::CommitStep(double t_old, double t_new, const Tag& tag)
{
  MPCSubcycled::CommitStep(t_old, t_new, tag);

  // also save our flow quantities, needed for interpolation
  auto [flow_current_tag, flow_next_tag] = tags_[0];

  if (subsurface_) {
    assign("porosity", flow_current_tag, flow_next_tag, *S_);
    if (chemistry_) {
      assign("mass_density_liquid", flow_current_tag, flow_next_tag, *S_);
    }
  }

  if (surface_ && chemistry_) {
    assign("surface-mass_density_liquid", flow_current_tag, flow_next_tag, *S_);
  }
}


} // namespace Amanzi
