/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (versions 1 & 2)  (nnc@lanl.gov)
           Ethan Coon (ATS version)
*/

#include "FunctionConstant.hh"
#include "FunctionFactory.hh"
#include "MultiFunction.hh"
#include "errors.hh"

#include "bc_factory.hh"

namespace Amanzi {

/* ******************************************************************
* Process Dirichet BC (pressure), step 1.
****************************************************************** */
Teuchos::RCP<Functions::BoundaryFunction>
BCFactory::CreateWithFunction(const std::string& list_name, const std::string& function_name) const
{
  Teuchos::RCP<Functions::BoundaryFunction> bc =
    Teuchos::rcp(new Functions::BoundaryFunction(mesh_));

  if (plist_.isParameter(list_name)) {
    if (plist_.isSublist(list_name)) {
      ProcessListWithFunction_(plist_.sublist(list_name), function_name, bc);
    } else {
      std::stringstream mstream;
      mstream << "BCFactory: " << list_name << " sublist error: not a sublist.";
      Errors::Message m(mstream.str());
      Exceptions::amanzi_throw(m);
    }
  }
  return bc;
}

/* ******************************************************************
* Process Dirichet BC (pressure), step 1.
****************************************************************** */
Teuchos::RCP<Functions::BoundaryFunction>
BCFactory::CreateWithoutFunction(const std::string& list_name) const
{
  Teuchos::RCP<Functions::BoundaryFunction> bc =
    Teuchos::rcp(new Functions::BoundaryFunction(mesh_));

  if (plist_.isParameter(list_name)) {
    if (plist_.isSublist(list_name)) {
      ProcessListWithoutFunction_(plist_.sublist(list_name), bc);
    } else {
      std::stringstream mstream;
      mstream << "BCFactory: " << list_name << " sublist error: not a sublist.";
      Errors::Message m(mstream.str());
      Exceptions::amanzi_throw(m);
    }
  }
  return bc;
}


Teuchos::RCP<Functions::DynamicBoundaryFunction>
BCFactory::CreateDynamicFunction(const std::string& list_name) const
{
  std::vector<std::string> regions;
  std::vector<std::string> bc_types;
  std::vector<std::string> bc_functions;

  Teuchos::RCP<Functions::BoundaryFunction> bc;

  Teuchos::RCP<Functions::DynamicBoundaryFunction> bcs =
    Teuchos::rcp(new Functions::DynamicBoundaryFunction(mesh_));

  if (plist_.isParameter(list_name)) {
    if (plist_.isSublist(list_name)) {
      if (plist_.sublist(list_name).isSublist("bcs")) {
        const Teuchos::ParameterList& list = plist_.sublist(list_name).sublist("bcs");
        // if (list.isParameter("regions")){
        //   regions = list.get<Teuchos::Array<std::string> >("regions").toVector();
        // } else {  // Parameter "regions" is missing.
        //   m << "parameter \"regions\" is missing";
        //   Errors::Message msg(m.str());
        //   Exceptions::amanzi_throw(msg);
        // }
        if (list.isParameter("bc types")) {
          bc_types = list.get<Teuchos::Array<std::string>>("bc types").toVector();
        } else { // Parameter "regions" is missing.
          Errors::Message msg;
          msg << "parameter \"bc types\" is missing";
          Exceptions::amanzi_throw(msg);
        }
        if (list.isParameter("bc functions")) {
          bc_functions = list.get<Teuchos::Array<std::string>>("bc functions").toVector();
        } else { // Parameter "regions" is missing.
          Errors::Message msg("parameter \"bc functions\" is missing");
          Exceptions::amanzi_throw(msg);
        }

        AMANZI_ASSERT(bc_types.size() == bc_functions.size());

        for (int i = 0; i < bc_types.size(); i++) {
          bc = Teuchos::rcp(new Functions::BoundaryFunction(mesh_));
          if (list.isParameter(bc_types[i])) {
            if (list.isSublist(bc_types[i])) {
              if ((bc_types[i] == "pressure") || (bc_types[i] == "water flux") ||
                  (bc_types[i] == "seepage face pressure") ||
                  (bc_types[i] == "seepage face head") || (bc_types[i] == "head") ||
                  (bc_types[i] == "fixed level")) {
                ProcessListWithFunction_(list.sublist(bc_types[i]), bc_functions[i], bc);
              } else {
                ProcessListWithoutFunction_(list.sublist(bc_types[i]), bc);
              }
            } else {
              Errors::Message msg;
              msg << "BCFactory: " << list_name << " sublist error: not a sublist.";
              Exceptions::amanzi_throw(msg);
            }
          }

          bcs->AddFunction(bc);
        }

      } else { // Parameter "regions" is missing.
        Errors::Message msg;
        msg << "sublist \"bcs \" is missing";
        Exceptions::amanzi_throw(msg);
      }

      if (plist_.sublist(list_name).isSublist("switch function")) {
        try {
          ProcessSpecWithFunction_(plist_.sublist(list_name), "switch function", bcs);
        } catch (Errors::Message& msg) {
          Errors::Message m;
          m << "in sublist " << "switch function" << ": " << msg.what();
          Exceptions::amanzi_throw(m);
        }
      } else { // ERROR -- parameter is not a sublist
        Errors::Message msg;
        msg << "parameter \"switch function\" is not a sublist";
        Exceptions::amanzi_throw(msg);
      }
    }
  }

  return bcs;
}


bool
BCFactory::CheckExplicitFlag(const std::string& list_name)
{
  bool is_explicit = false;
  if (plist_.isSublist(list_name)) {
    is_explicit = plist_.sublist(list_name).get<bool>("explicit time index", false);
    plist_.sublist(list_name).remove("explicit time index");
  }
  return is_explicit;
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 2.
****************************************************************** */
void
BCFactory::ProcessListWithFunction_(const Teuchos::ParameterList& list,
                                    const std::string& function_name,
                                    const Teuchos::RCP<Functions::BoundaryFunction>& bc) const
{
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList spec = list.sublist(name);
      try {
        ProcessSpecWithFunction_(spec, function_name, bc);
      } catch (Errors::Message& msg) {
        Errors::Message m;
        m << "in sublist " << spec.name().c_str() << ": " << msg.what();
        Exceptions::amanzi_throw(m);
      }
    }
  }
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 3.
****************************************************************** */
void
BCFactory::ProcessSpecWithFunction_(const Teuchos::ParameterList& list,
                                    const std::string& function_name,
                                    const Teuchos::RCP<Functions::BoundaryFunction>& bc) const
{
  std::stringstream m;
  std::vector<std::string> regions;

  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string>>("regions")) {
      regions = list.get<Teuchos::Array<std::string>>("regions").toVector();
    } else {
      m << "parameter \"regions\" is not of type Array string";
      Errors::Message msg(m.str());
      Exceptions::amanzi_throw(msg);
    }
  } else { // Parameter "regions" is missing.
    m << "parameter \"regions\" is missing";
    Errors::Message msg(m.str());
    Exceptions::amanzi_throw(msg);
  }

  Teuchos::ParameterList f_list;
  if (list.isParameter(function_name)) {
    if (list.isSublist(function_name)) {
      // validate function-factory sublist
      f_list = list.sublist(function_name);
    } else { // function_name is not a sublist
      m << "parameter " << function_name << " is not a sublist";
      Errors::Message msg(m.str());
      Exceptions::amanzi_throw(msg);
    }
  } else { // function_name sublist is missing.
    m << "sublist " << function_name << " is missing";
    Errors::Message msg(m.str());
    Exceptions::amanzi_throw(msg);
  }

  // Make the boundary pressure function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist " << function_name << ": " << msg.what();
    Errors::Message lmsg(m.str());
    Exceptions::amanzi_throw(lmsg);
  }

  // A bit hacky -- this entire code needs to be revisited in light of the new
  // options for mesh_functions.
  Teuchos::RCP<MultiFunction> func = Teuchos::rcp(new MultiFunction(f));


  // Add this BC specification to the boundary function.
  bc->Define(regions, func);
  bc->Finalize();
}


void
BCFactory::ProcessSpecWithFunctionRegions_(
  const Teuchos::ParameterList& list,
  const std::string& function_name,
  std::vector<std::string>& regions,
  const Teuchos::RCP<Functions::BoundaryFunction>& bc) const
{
  std::stringstream m;

  Teuchos::ParameterList f_list;
  if (list.isParameter(function_name)) {
    if (list.isSublist(function_name)) {
      // validate function-factory sublist
      f_list = list.sublist(function_name);
    } else { // function_name is not a sublist
      m << "parameter " << function_name << " is not a sublist";
      Errors::Message msg(m.str());
      Exceptions::amanzi_throw(msg);
    }
  } else { // function_name sublist is missing.
    m << "sublist " << function_name << " is missing";
    Errors::Message msg(m.str());
    Exceptions::amanzi_throw(msg);
  }

  // Make the boundary pressure function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(f_list));
  } catch (Errors::Message& msg) {
    m << "error in sublist " << function_name << ": " << msg.what();
    Errors::Message lmsg(m.str());
    Exceptions::amanzi_throw(lmsg);
  }

  // A bit hacky -- this entire code needs to be revisited in light of the new
  // options for mesh_functions.
  Teuchos::RCP<MultiFunction> func = Teuchos::rcp(new MultiFunction(f));

  // Add this BC specification to the boundary function.
  bc->Define(regions, func);
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 2.
****************************************************************** */
void
BCFactory::ProcessListWithoutFunction_(const Teuchos::ParameterList& list,
                                       const Teuchos::RCP<Functions::BoundaryFunction>& bc) const
{
  // Iterate through the BC specification sublists in the list.
  // All are expected to be sublists of identical structure.
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i) {
    std::string name = i->first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList spec = list.sublist(name);
      try {
        ProcessSpecWithoutFunction_(spec, bc);
      } catch (Errors::Message& msg) {
        std::stringstream m;
        m << "in sublist " << spec.name().c_str() << ": " << msg.what();
        Errors::Message lmsg(m.str());
        Exceptions::amanzi_throw(lmsg);
      }
    }
  }
}


/* ******************************************************************
* Process Dirichet BC (pressure), step 3.
****************************************************************** */
void
BCFactory::ProcessSpecWithoutFunction_(const Teuchos::ParameterList& list,
                                       const Teuchos::RCP<Functions::BoundaryFunction>& bc) const
{
  Errors::Message m;
  std::vector<std::string> regions;

  if (list.isParameter("regions")) {
    if (list.isType<Teuchos::Array<std::string>>("regions")) {
      regions = list.get<Teuchos::Array<std::string>>("regions").toVector();
    } else {
      m << "parameter \"regions\" is not of type \"Array string\"";
      Exceptions::amanzi_throw(m);
    }
  } else { // Parameter "regions" is missing.
    m << "parameter \"regions\" is missing";
    Exceptions::amanzi_throw(m);
  }

  // Make the boundary pressure function.
  Teuchos::RCP<Function> f = Teuchos::rcp(new FunctionConstant(0.));

  // A bit hacky -- this entire code needs to be revisited in light of the new
  // options for mesh_functions.
  Teuchos::RCP<MultiFunction> func = Teuchos::rcp(new MultiFunction(f));

  // Add this BC specification to the boundary function.
  bc->Define(regions, func);
  bc->Finalize();
}

} // namespace Amanzi
