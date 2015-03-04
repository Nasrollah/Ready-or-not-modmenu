/****************************************************************************** 

  (c) 2004-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE-SCOREC file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/

/**
 * @file CrvTriBlending.h
 * @author Qiukai Lu
 */

#ifndef CURVEMESHADAPT_CRVTRIBLENDING_H_
#define CURVEMESHADAPT_CRVTRIBLENDING_H_
 
#include "CrvFace.h"
#include "CrvEdge.h"

/** @brief The class of curved triangular faces using bounding edges and blending functions 
 */
class CrvTriBlending : public CrvFace {
 public:
    /** \brief ctor 1
     */
    CrvTriBlending(pMeshEnt in_face, int order);
    /** \brief ctor 2
     */
    CrvTriBlending(VtxPtrVec in_vert_vec, int order);
    /** \brief dtor
     */
    ~CrvTriBlending() {
      clear_edges();
    }

    int eval(double u, double v, double * pos); 

    int v_data_size() const {
      return 0;
    }

    std::string v_tag_name() const {
      return std::string( "crv_tri_blending" );
    }

    int eval_deriv1(double in_xi_1,
                    double in_xi_2,
                    Point3d & dxyz_dxi1,
                    Point3d & dxyz_dxi2);


 protected:

 private:
    /// set up the corresponding vertex and edges
    int setup_verts();

    /// this functions allocates memory
    int setup_edges();

    /// this function frees memory
    int clear_edges();

    /// eval first derivatives by forward differencing
    void eval_diff_fd(double in_xi1,
                      double in_xi2,
                      Point3d & dxyz_dxi1,
                      Point3d & dxyz_dxi2);

    /// eval first derivatives by analytic expressions
    void eval_diff_analytic(double in_xi1,
                            double in_xi2,
                            Point3d & dxyz_dxi1,
                            Point3d & dxyz_dxi2);

    /// eval first derivatives by analytic expressions
    void eval_diff_analytic_c2(double in_xi1,
                            double in_xi2,
                            Point3d & dxyz_dxi1,
                            Point3d & dxyz_dxi2);

    /// pointers to the mesh entities
    Point3d  m_vertex_100;
    Point3d  m_vertex_010;
    Point3d  m_vertex_001;
    CrvEdgePtr m_edge_P;
    CrvEdgePtr m_edge_Q;
    CrvEdgePtr m_edge_R;
};

#endif  // CURVEMESHADAPT_CRVTRIBLENDING_H_
