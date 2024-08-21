/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "mpc_flow_transport.hh"

namespace Amanzi {

void
MPCFlowTransport::parseParameterList()
{
  // are we doing surface, subsurface, or integrated transport?
  auto flow_plist = getSubPKPlist_(0);
  bool surface(false), subsurface(false);

  if (flow_plist->isParameter("domain name")) {
    std::string domain = Keys::readDomain(*flow_plist);
    if (Keys::in(domain, "surface")) {
      surface = true;
    } else {
      subsurface = true;
    }
  } else {
    surface = true;
    subsurface = true;
  }

  if (subsurface) {
    auto [flow_current_tag, flow_next_tag] = tags_[0];
    auto [transport_current_tag, transport_next_tag] = tags_[1];

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

      // set the saturation_liquid as an interpolated field
      Teuchos::ParameterList& sat_list_current = S_->GetEvaluatorList(Keys::getKey("saturation_liquid", transport_current_tag));
      sat_list_current.set<std::string>("evaluator type", "temporal interpolation");
      sat_list_current.set<std::string>("current tag", flow_current_tag.get());
      sat_list_current.set<std::string>("next tag", flow_next_tag.get());

      Teuchos::ParameterList& sat_list_next = S_->GetEvaluatorList(Keys::getKey("saturation_liquid", transport_next_tag));
      sat_list_next.set<std::string>("evaluator type", "temporal interpolation");
      sat_list_next.set<std::string>("current tag", flow_current_tag.get());
      sat_list_next.set<std::string>("next tag", flow_next_tag.get());

      Teuchos::ParameterList& poro_list_next = S_->GetEvaluatorList(Keys::getKey("porosity", transport_next_tag));
      poro_list_next.set<std::string>("evaluator type", "alias");
      poro_list_next.set<std::string>("target", Keys::getKey("porosity", flow_next_tag, true));
    }
  }

  if (surface) {
    auto [flow_current_tag, flow_next_tag] = tags_[0];
    auto [transport_current_tag, transport_next_tag] = tags_[1];

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

      // set the surface-ponded_depth as an interpolated field
      Teuchos::ParameterList& pd_list_current = S_->GetEvaluatorList(Keys::getKey("surface-ponded_depth", transport_current_tag));
      pd_list_current.set<std::string>("evaluator type", "temporal interpolation");
      pd_list_current.set<std::string>("current tag", flow_current_tag.get());
      pd_list_current.set<std::string>("next tag", flow_next_tag.get());

      Teuchos::ParameterList& pd_list_next = S_->GetEvaluatorList(Keys::getKey("surface-ponded_depth", transport_next_tag));
      pd_list_next.set<std::string>("evaluator type", "temporal interpolation");
      pd_list_next.set<std::string>("current tag", flow_current_tag.get());
      pd_list_next.set<std::string>("next tag", flow_next_tag.get());
    }
  }

  MPCSubcycled::parseParameterList();
}


} // namespace Amanzi
