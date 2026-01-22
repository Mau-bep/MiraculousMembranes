

#ifndef __SHAPENLP_HPP__
#define __SHAPENLP_HPP__


#include "coin-or/IpTNLP.hpp"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "Mem-3dg.h"

using namespace Ipopt;


class ShapeNLP : public TNLP
{
public:
    ShapeNLP();
    virtual ~ShapeNLP();

    geometrycentral::surface::ManifoldSurfaceMesh* mesh;
    geometrycentral::surface::VertexPositionGeometry* geometry;
    Mem3DG* M3DG;


    virtual void set_M3DG(Mem3DG* inputM3DG){
        M3DG = inputM3DG;
      }

    /**@name Overloaded from TNLP */
   virtual void update_positions(
       const Number* x
    );
    
    virtual bool get_nlp_info(
       Index&          n,
       Index&          m,
       Index&          nnz_jac_g,
       Index&          nnz_h_lag,
       IndexStyleEnum& index_style
    );

    virtual bool get_bounds_info(
       Index   n,
       Number* x_l,
       Number* x_u,
       Index   m,
       Number* g_l,
       Number* g_u
    );

    virtual bool get_starting_point(
       Index   n,
       bool    init_x,
       Number* x,
       bool    init_z,
       Number* z_L,
       Number* z_U,
       Index   m,
       bool    init_lambda,
       Number* lambda
    );

    virtual bool eval_f(
       Index         n,
       const Number* x,
       bool          new_x,
       Number&       obj_value
    );

    virtual bool eval_grad_f(
       Index         n,
       const Number* x,
       bool          new_x,
       Number*       grad_f
    );

    virtual bool eval_g(
       Index         n,
       const Number* x,
       bool          new_x,
       Index         m,
       Number*       g
    );

    virtual bool eval_jac_g(
       Index         n,
       const Number* x,
       bool          new_x,
       Index         m,
       Index         nele_jac,
       Index*        iRow,
       Index*        jCol,
       Number*       values
    );

    virtual bool eval_h(
       Index         n,
       const Number* x,
       bool          new_x,
       Number        obj_factor,
       Index         m,
       const Number* lambda,
       bool          new_lambda,
       Index         nele_hess,
       Index*        iRow,
       Index*        jCol,
       Number*       values
    );
    virtual void finalize_solution(
       SolverReturn               status,
       Index                      n,
       const Number*              x,
       const Number*              z_L,
       const Number*              z_U,
       Index                      m,
       const Number*              g,
       const Number*              lambda,
       Number                     obj_value,
       const IpoptData*           ip_data,
       IpoptCalculatedQuantities* ip_cq
    );
   virtual bool intermediate_callback(
      AlgorithmMode              mode,
      Index                      iter,
      Number                     obj_value,
      Number                     inf_pr,
      Number                     inf_du,
      Number                     mu,
      Number                     d_norm,
      Number                     regularization_size,
      Number                     alpha_du,
      Number                     alpha_pr,
      Index                      ls_trials,
      const IpoptData*           ip_data,
      IpoptCalculatedQuantities* ip_cq
   );
private:

    ShapeNLP(
       const ShapeNLP&
    );

    ShapeNLP& operator=(
       const ShapeNLP&
    );
};

#endif
