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
  Key surf_lwc_key, sub_lwc_key;

  if (flow_plist->isParameter("domain name")) {
    std::string domain = Keys::readDomain(*flow_plist);
    if (Keys::in(domain, "surface")) {
      surface = true;
      sub_lwc_key = Keys::readKey(*getSubPKPlist_(1), domain, "water content", "water_content");
    } else {
      subsurface = true;
      surf_lwc_key = Keys::readKey(*getSubPKPlist_(1), domain, "water content", "water_content");
    }
  } else {
    // no domain name means it is an MPC, likely coupled water, but maybe permafrost
    surface = true;
    subsurface = true;
    auto transport_names = getSubPKPlist_(1)->get<Teuchos::Array<std::string>>("PKs order");

    AMANZI_ASSERT(transport_names.size() == 2);
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
      flux_list.set<std::string>("evaluator type", "alias");
      flux_list.set<std::string>("target", Keys::getKey("water_flux", flow_next_tag, true));

      // velocity for dispersivity
      Teuchos::ParameterList& velo_list = S_->GetEvaluatorList(Keys::getKey("darcy_velocity", transport_next_tag));
      velo_list.set<std::string>("evaluator type", "alias");
      velo_list.set<std::string>("target", Keys::getKey("darcy_velocity", flow_next_tag, true));

      // now set the liquid water content as an interpolated field at next
      // note that current copy is kept by transport PK
      Teuchos::ParameterList& lwc_list_next = S_->GetEvaluatorList(Keys::getKey(sub_lwc_key, transport_next_tag));
      lwc_list_next.set<std::string>("evaluator type", "temporal interpolation");
      lwc_list_next.set<std::string>("current tag", flow_current_tag.get());
      lwc_list_next.set<std::string>("next tag", flow_next_tag.get());

      // porosity used with velocity to compute particle velocity when dispersion is on
      Teuchos::ParameterList& poro_list_next = S_->GetEvaluatorList(Keys::getKey("porosity", transport_next_tag));
      poro_list_next.set<std::string>("evaluator type", "alias");
      poro_list_next.set<std::string>("target", Keys::getKey("porosity", flow_next_tag, true));
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
      flux_list.set<std::string>("evaluator type", "alias");
      flux_list.set<std::string>("target", Keys::getKey("surface-water_flux", flow_next_tag, true));

      // velocity for dispersivity
      Teuchos::ParameterList& velo_list = S_->GetEvaluatorList(Keys::getKey("surface-water_velocity", transport_next_tag));
      velo_list.set<std::string>("evaluator type", "alias");
      velo_list.set<std::string>("target", Keys::getKey("surface-water_velocity", flow_next_tag, true));

      // set the liquid water content as an interpolated field
      // note that current copy is kept by transport PK
      Teuchos::ParameterList& lwc_list_next = S_->GetEvaluatorList(Keys::getKey(surf_lwc_key, transport_next_tag));
      lwc_list_next.set<std::string>("evaluator type", "temporal interpolation");
      lwc_list_next.set<std::string>("current tag", flow_current_tag.get());
      lwc_list_next.set<std::string>("next tag", flow_next_tag.get());
    }
  }

  MPCSubcycled::parseParameterList();

  if (subsurface) {
    // first require lwc@NEXT -- otherwise it will get called to be an alias
    // lwc@TRANSPORT_NEXT, which is an interpolator that depends upon this.
    requireAtNext(sub_lwc_key, Tags::NEXT, *S_);

    // now call at lwc@transport_next, which will be the interpolant
    requireAtNext(sub_lwc_key, transport_next_tag, *S_);
  }

  if (surface) {
    // first require lwc@NEXT -- otherwise it will get called to be an alias
    // lwc@TRANSPORT_NEXT, which is an interpolator that depends upon this.
    requireAtNext(surf_lwc_key, Tags::NEXT, *S_);

    // now call at lwc@transport_next, which will be the interpolant
    requireAtNext(surf_lwc_key, transport_next_tag, *S_);
  }

}


} // namespace Amanzi
