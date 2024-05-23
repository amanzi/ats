/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for water drainage through a pipe network.
  This depends on:
  1) flow depth above surface elevation (surface_depth_key_)
  2) hydraulic head in the pipe network (pressure_head_key_)
  3) pipe drain length that is determined by the bathymetry difference

  Authors: Giacomo Capodaglio (gcapodaglio@lanl.gov)
           Naren Vohra (vohra@lanl.gov)
*/

#include "pipe_drain_evaluator.hh"
#include "Geometry.hh"

namespace Amanzi {
namespace Flow {


PipeDrainEvaluator::PipeDrainEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  manhole_radius_ = plist_.get<double>("manhole radius", 0.24);
  energ_loss_coeff_weir_ = plist.get<double>("energy losses coeff weir", 0.54);
  energ_loss_coeff_subweir_ = plist.get<double>("energy losses coeff submerged weir", 0.056);
  energ_loss_coeff_orifice_ = plist.get<double>("energy losses coeff orifice", 0.167);

  sw_domain_name_ = plist.get<std::string>("surface domain name", "surface");
  pipe_domain_name_ = plist.get<std::string>("pipe domain name", "pipe");

  Tag tag = my_keys_.front().second;

  // my dependencies
  surface_depth_key_ = Keys::readKey(plist_, sw_domain_name_, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ surface_depth_key_, tag });

  // bathymetry
  surface_bathymetry_key_ =
    Keys::readKey(plist_, sw_domain_name_, "surface bathymetry", "bathymetry");
  pipe_bathymetry_key_ = Keys::readKey(plist_, pipe_domain_name_, "pipe bathymetry", "bathymetry");

  if (!pipe_domain_name_.empty()) {
    pressure_head_key_ = Keys::readKey(plist_, pipe_domain_name_, "pressure head", "pressure_head");
    dependencies_.insert(KeyTag{ pressure_head_key_, tag });
  }

  // figure out if SW or pipe is calling the evaluator
  auto domain = Keys::getDomain(my_keys_.front().first);

  if (domain == pipe_domain_name_) {
    pipe_flag_ = true;
    sw_flag_ = false;

    mask_key_ = Keys::readKey(plist_, pipe_domain_name_, "manhole locations", "manhole_locations");
    dependencies_.insert(KeyTag{ mask_key_, tag });

    sink_source_coeff_ = -1.0;
  } else if (domain == sw_domain_name_) {
    pipe_flag_ = false;
    sw_flag_ = true;

    mask_key_ = Keys::readKey(plist_, sw_domain_name_, "manhole locations", "manhole_locations");
    dependencies_.insert(KeyTag{ mask_key_, tag });

    sink_source_coeff_ = 1.0;
  } else {
    std::cout << "Unknown domain in pipe drain evaluator" << std::endl;
    AMANZI_ASSERT(0);
  }
}


Teuchos::RCP<Evaluator>
PipeDrainEvaluator::Clone() const
{
  return Teuchos::rcp(new PipeDrainEvaluator(*this));
}


void
PipeDrainEvaluator::CreateCellMap(const State& S)
{
  // grab the meshes
  Teuchos::RCP<const AmanziMesh::Mesh> pipe_mesh = S.GetMesh(pipe_domain_name_);
  Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh = S.GetMesh(sw_domain_name_);

  // get the manhole locations/ cell maps
  Tag tag = my_keys_.front().second;

  const Epetra_MultiVector& mnhMask =
    *S.GetPtr<CompositeVector>(mask_key_, tag)->ViewComponent("cell", false);

  // loop over mesh cells and create map
  int ncells_pipe =
    pipe_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int ncells_sw =
    surface_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  // SW map from pipe -> SW domain
  if (pipe_flag_ == true) {
    pipe_map_.resize(ncells_pipe);
    for (int c_pipe = 0; c_pipe < ncells_pipe; ++c_pipe) {
      if (std::abs(mnhMask[0][c_pipe] - 1.0) < 1.e-12) {
        const Amanzi::AmanziGeometry::Point& xc_pipe = pipe_mesh->getCellCentroid(c_pipe);
        for (int c_sw = 0; c_sw < ncells_sw; ++c_sw) {
          std::vector<AmanziGeometry::Point> coords;
          auto cnodes = surface_mesh->getCellNodes(c_sw);
          for (int node_sw = 0; node_sw < cnodes.size(); node_sw++) {
            coords.push_back(surface_mesh->getNodeCoordinate(node_sw));
          }

          if (AmanziGeometry::point_in_polygon(xc_pipe, coords) == true) {
            pipe_map_[c_pipe] = c_sw;
            break;
          }
        }
      }
    }
    pipe_map_created_ = true;
  }

  // Pipe map from SW -> pipe domain
  if (sw_flag_ == true) {
    sw_map_.resize(ncells_sw);
    for (int c_sw = 0; c_sw < ncells_sw; ++c_sw) {
      if (std::abs(mnhMask[0][c_sw] - 1.0) < 1.e-12) {
        const Amanzi::AmanziGeometry::Point& xc_sw = surface_mesh->getCellCentroid(c_sw);
        for (int c_pipe = 0; c_pipe < ncells_pipe; ++c_pipe) {
          std::vector<AmanziGeometry::Point> coords;
          auto cnodes = pipe_mesh->getCellNodes(c_pipe);
          for (int node_pipe = 0; node_pipe < cnodes.size(); node_pipe++) {
            coords.push_back(pipe_mesh->getNodeCoordinate(node_pipe));
          }

          if (AmanziGeometry::point_in_polygon(xc_sw, coords) == true) {
            sw_map_[c_sw] = c_pipe;
            break;
          }
        }
      }
    }
    sw_map_created_ = true;
  }
}

void
PipeDrainEvaluator::EnsureCompatibility_ToDeps_(State& S, const CompositeVectorSpace& fac)
{
  auto domain1 = Keys::getDomain(my_keys_.front().first);
  if (domain1 == pipe_domain_name_) {
    for (const auto& dep : dependencies_) {
      auto domain = Keys::getDomain(dep.first);
      if (pipe_domain_name_ == domain) {
        auto& dep_fac = S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second);
        dep_fac.Update(fac);
      }
    }
  }

  else {
    for (const auto& dep : dependencies_) {
      auto domain = Keys::getDomain(dep.first);
      if (sw_domain_name_ == domain) {
        auto& dep_fac = S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second);
        dep_fac.Update(fac);
      }
    }
  }
}

void
PipeDrainEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Epetra_MultiVector& res = *result[0]->ViewComponent("cell", false);

  const Epetra_MultiVector& srfcDepth =
    *S.GetPtr<CompositeVector>(surface_depth_key_, tag)->ViewComponent("cell", false);

  const Epetra_MultiVector& mnhMask =
    *S.GetPtr<CompositeVector>(mask_key_, tag)->ViewComponent("cell", false);

  // surface and pipe bathymetry needed for pipe drain length
  const Epetra_MultiVector& srfcB =
    *S.GetPtr<CompositeVector>(surface_bathymetry_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& pipeB =
    *S.GetPtr<CompositeVector>(pipe_bathymetry_key_, tag)->ViewComponent("cell", false);

  double g = norm(S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT));
  double pi = 3.14159265359;
  double mnhArea = pi * manhole_radius_ * manhole_radius_;
  double mnhPerimeter = 2.0 * pi * manhole_radius_;
  double sqrtTwoG = sqrt(2.0 * g);

  int ncells = res.MyLength();

  // generate cell maps
  if (pipe_map_created_ == false) {
    CreateCellMap(S);
    pipe_map_created_ = true;
  }
  if (sw_map_created_ == false) {
    CreateCellMap(S);
    sw_map_created_ = true;
  }

  double drain_length;

  if (!pipe_domain_name_.empty()) {
    const Epetra_MultiVector& pressHead =
      *S.GetPtr<CompositeVector>(pressure_head_key_, tag)->ViewComponent("cell", false);
    int c_pipe, c_sw;

    for (int c = 0; c != ncells; ++c) {
      // use cell map
      if (pipe_flag_ == true) {
        c_pipe = c;
        c_sw = pipe_map_[c];
      } else if (sw_flag_ == true) {
        c_sw = c;
        c_pipe = sw_map_[c];
      }

      // calculate pipe drain length using bathymetry
      drain_length = srfcB[0][c_sw] - pipeB[0][c_pipe];

      // throw error if bathymetrys are not physical
      if (drain_length < 0.0) {
        std::cout << "Pipe drain length negative; bathymetrys not physical" << std::endl;
        AMANZI_ASSERT(0);
      }

      if (pressHead[0][c_pipe] < drain_length) {
        res[0][c] = -mnhMask[0][c] * 2.0 / 3.0 * energ_loss_coeff_weir_ * mnhPerimeter * sqrtTwoG *
                    pow(srfcDepth[0][c_sw], 3.0 / 2.0);
      } else if (drain_length < pressHead[0][c_pipe] &&
                 pressHead[0][c_pipe] < (drain_length + srfcDepth[0][c_sw])) {
        res[0][c] = -mnhMask[0][c] * energ_loss_coeff_subweir_ * mnhArea * sqrtTwoG *
                    sqrt(srfcDepth[0][c_sw] + drain_length - pressHead[0][c_pipe]);
      } else if (pressHead[0][c_pipe] > (drain_length + srfcDepth[0][c_sw])) {
        res[0][c] = mnhMask[0][c] * energ_loss_coeff_orifice_ * mnhArea * sqrtTwoG *
                    sqrt(pressHead[0][c_pipe] - drain_length - srfcDepth[0][c_sw]);
      }
      res[0][c] *= sink_source_coeff_;
    }
  } else {
    for (int c = 0; c != ncells; ++c) {
      res[0][c] = -mnhMask[0][c] * 2.0 / 3.0 * energ_loss_coeff_weir_ * mnhPerimeter * sqrtTwoG *
                  pow(srfcDepth[0][c], 3.0 / 2.0);
      res[0][c] *= sink_source_coeff_;
    }
  }
}

void
PipeDrainEvaluator::EvaluatePartialDerivative_(const State& S,
                                               const Key& wrt_key,
                                               const Tag& wrt_tag,
                                               const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  Epetra_MultiVector& res = *result[0]->ViewComponent("cell", false);

  const Epetra_MultiVector& srfcDepth =
    *S.GetPtr<CompositeVector>(surface_depth_key_, tag)->ViewComponent("cell", false);

  const Epetra_MultiVector& mnhMask =
    *S.GetPtr<CompositeVector>(mask_key_, tag)->ViewComponent("cell", false);

  double g = norm(S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT));
  double pi = 3.14159265359;
  double mnhArea = pi * manhole_radius_ * manhole_radius_;
  double mnhPerimeter = 2.0 * pi * manhole_radius_;
  double sqrtTwoG = sqrt(2.0 * g);

  int ncells = res.MyLength();

  // surface and pipe bathymetry needed for pipe drain length
  const Epetra_MultiVector& srfcB =
    *S.GetPtr<CompositeVector>(surface_bathymetry_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& pipeB =
    *S.GetPtr<CompositeVector>(pipe_bathymetry_key_, tag)->ViewComponent("cell", false);

  double drain_length;

  if (!pipe_domain_name_.empty()) {
    const Epetra_MultiVector& pressHead =
      *S.GetPtr<CompositeVector>(pressure_head_key_, tag)->ViewComponent("cell", false);

    int c_pipe, c_sw;

    if (wrt_key == surface_depth_key_) {
      for (int c = 0; c != ncells; ++c) {
        // use cell map
        if (pipe_flag_ == true) {
          c_pipe = c;
          c_sw = pipe_map_[c];
        } else if (sw_flag_ == true) {
          c_sw = c;
          c_pipe = sw_map_[c];
        }
        // calculate pipe drain length using bathymetry
        drain_length = srfcB[0][c_sw] - pipeB[0][c_pipe];

        if (pressHead[0][c_pipe] < drain_length) {
          res[0][c] = -mnhMask[0][c] * energ_loss_coeff_weir_ * mnhPerimeter * sqrtTwoG *
                      sqrt(srfcDepth[0][c_sw]);
        } else if (drain_length < pressHead[0][c_pipe] &&
                   pressHead[0][c_pipe] < (drain_length + srfcDepth[0][c_sw])) {
          res[0][c] = -0.5 * mnhMask[0][c] * energ_loss_coeff_subweir_ * mnhArea * sqrtTwoG /
                      sqrt(srfcDepth[0][c] + drain_length - pressHead[0][c]);
        } else if (pressHead[0][c_pipe] > (drain_length + srfcDepth[0][c_sw])) {
          res[0][c] = -0.5 * mnhMask[0][c] * energ_loss_coeff_orifice_ * mnhArea * sqrtTwoG /
                      sqrt(pressHead[0][c_pipe] - drain_length - srfcDepth[0][c_sw]);
        }
      }
    } else if (wrt_key == pressure_head_key_) {
      for (int c = 0; c != ncells; ++c) {
        // use cell map
        if (pipe_flag_ == true) {
          c_pipe = c;
          c_sw = pipe_map_[c];
        } else if (sw_flag_ == true) {
          c_sw = c;
          c_pipe = sw_map_[c];
        }
        // calculate pipe drain length using bathymetry
        drain_length = srfcB[0][c_sw] - pipeB[0][c_pipe];

        if (pressHead[0][c_pipe] < drain_length) {
          res[0][c] = 0.0;
        } else if (drain_length < pressHead[0][c_pipe] &&
                   pressHead[0][c_pipe] < (drain_length + srfcDepth[0][c_sw])) {
          res[0][c] = 0.5 * mnhMask[0][c] * energ_loss_coeff_subweir_ * mnhArea * sqrtTwoG /
                      sqrt(srfcDepth[0][c_sw] + drain_length - pressHead[0][c_pipe]);
        } else if (pressHead[0][c_pipe] > (drain_length + srfcDepth[0][c_sw])) {
          res[0][c] = 0.5 * mnhMask[0][c] * energ_loss_coeff_orifice_ * mnhArea * sqrtTwoG /
                      sqrt(pressHead[0][c_pipe] - drain_length - srfcDepth[0][c_sw]);
        }
      }
    } else {
      AMANZI_ASSERT(0);
    }
  } else {
    if (wrt_key == surface_depth_key_) {
      for (int c = 0; c != ncells; ++c) {
        res[0][c] =
          -mnhMask[0][c] * energ_loss_coeff_weir_ * mnhPerimeter * sqrtTwoG * sqrt(srfcDepth[0][c]);
      }
    } else {
      AMANZI_ASSERT(0);
    }
  }
}


} // namespace Flow
} // namespace Amanzi
