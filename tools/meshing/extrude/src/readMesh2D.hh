#ifndef MESH_READER_HH_
#define MESH_READER_HH_

#include <map>

#include "Mesh2D.hh"

namespace Amanzi {
namespace AmanziGeometry {


struct PointFactory {
  PointFactory() {}

  bool addPoint(const Point& p, int& id) {
    double key = norm(p);
    auto range = points_sorted.equal_range(key);
    for (auto sp=range.first; sp!=range.second; ++sp) {
      if (sp->second.second == p) {
        id = sp->second.first;
        return false;
      }
    }

    id = points.size();
    points.push_back(p);
    points_sorted.insert(std::make_pair(key,
            std::make_pair(id, p)));
    return true;
  }

  std::multimap<double,std::pair<int,Point> > points_sorted;
  std::vector<Point> points;
};

Mesh2D readMesh2D_text(const std::string& filename,
                       std::vector<int>& soil_type,
                       std::vector<int>& bedrock_type,
                       std::vector<double>& depth_to_bedrock,
                       double cut_x=1.e80,
                       double cut_y=1.e80);

}
}


#endif
