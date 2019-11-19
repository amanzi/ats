/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov

  A collection of multiscale porosity models along with a mesh partition.
*/

#include "dbc.hh"
#include "MultiscaleTransportPorosityFactory.hh"
#include "MultiscaleTransportPorosityPartition.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Non-member factory.
****************************************************************** */
Teuchos::RCP<MultiscaleTransportPorosityPartition> CreateMultiscaleTransportPorosityPartition(
    Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    Teuchos::RCP<Teuchos::ParameterList> tplist)
{
  auto msp_list = Teuchos::sublist(tplist, "multiscale models", true);
  auto diff_list = Teuchos::sublist(tplist, "molecular diffusion", true);

  auto mol_diff = diff_list->get<Teuchos::Array<double> >("aqueous values");

  MultiscaleTransportPorosityFactory factory;
  std::vector<Teuchos::RCP<MultiscaleTransportPorosity> > msp;
  std::vector<std::vector<std::string> > regions;

  for (auto it = msp_list->begin(); it != msp_list->end(); ++it) {
    if (msp_list->isSublist(it->first)) {
      Teuchos::ParameterList plist = msp_list->sublist(it->first);
      plist.set<Teuchos::Array<double> >("molecular diffusion", mol_diff);

      regions.push_back(plist.get<Teuchos::Array<std::string> >("regions").toVector());
      msp.push_back(factory.Create(plist));
    } else {
      AMANZI_ASSERT(0);
    }
  }

  auto partition = Teuchos::rcp(new Functions::MeshPartition());
  partition->Initialize(mesh, AmanziMesh::CELL, regions, -1);
  partition->Verify();

  return Teuchos::rcp(new MultiscaleTransportPorosityPartition(partition, msp));
}


/* ******************************************************************
* Non-member function quering partition.
****************************************************************** */
int NumberMatrixNodes(Teuchos::RCP<MultiscaleTransportPorosityPartition>& msp)
{
  int nnodes(0);
  const auto& list = msp->second;

  for (int i = 0; i < list.size(); ++i) {
    nnodes = std::max(nnodes, list[i]->NumberMatrixNodes());
  }

  return nnodes;
}

}  // namespace Transport
}  // namespace Amanzi

