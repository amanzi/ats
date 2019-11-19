/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  A collection of MDMs along with a mesh partition.
*/

#include "dbc.hh"
#include "MDMFactory.hh"
#include "MDMPartition.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Non-member factory.
****************************************************************** */
Teuchos::RCP<MDMPartition> CreateMDMPartition(
    Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    Teuchos::RCP<Teuchos::ParameterList> plist, bool& flag)
{
  MDMFactory factory;
  Teuchos::RCP<MDM> mdm;
  std::vector<Teuchos::RCP<MDM> > mdm_list;
  std::vector<std::vector<std::string> > regions;

  flag = false;
  for (Teuchos::ParameterList::ConstIterator lcv = plist->begin(); lcv != plist->end(); ++lcv) {
    std::string name = lcv->first;
    if (plist->isSublist(name)) {
      Teuchos::ParameterList sublist = plist->sublist(name);
      regions.push_back(sublist.get<Teuchos::Array<std::string> >("regions").toVector());

      mdm = factory.Create(sublist);
      mdm->set_dim(mesh->space_dimension());
      flag |= mdm->is_valid();
      mdm_list.push_back(mdm);
    } else {
      AMANZI_ASSERT(0);
    }
  }

  Teuchos::RCP<Functions::MeshPartition> partition = Teuchos::rcp(new Functions::MeshPartition());
  partition->Initialize(mesh, AmanziMesh::CELL, regions, -1);
  partition->Verify();

  return Teuchos::rcp(new MDMPartition(partition, mdm_list));
}

}  // namespace Transport
}  // namespace Amanzi

