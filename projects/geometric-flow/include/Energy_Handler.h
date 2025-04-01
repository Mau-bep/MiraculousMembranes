#pragma once

#include <Eigen/Core>
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "Beads.h"
// #include "Mem-3dg.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


class E_Handler {
    public: 
            
        ManifoldSurfaceMesh* mesh;
        VertexPositionGeometry* geometry;
        std::vector<std::string> Energies;
        std::vector<std::vector<double>> Energy_constants;
        std::vector<Bead *> Beads;
        
        
    
}