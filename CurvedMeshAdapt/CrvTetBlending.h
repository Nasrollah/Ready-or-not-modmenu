/****************************************************************************** 

  (c) 2004-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE-SCOREC file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/

#ifndef MESHADAPT_CRVTETBLENDING_H
#define MESHADAPT_CRVTETBLENDING_H

#include "CrvTet.h"
#include <stdlib.h>
#include <vector>

typedef class mEntity * pMeshEnt;
class ParametricCurve;

/** \brief This class implements the blending based volume mapping
 * The blending mapping formulation is coded according to the Dey1997
 * paper Geometric Representation Issues of p-Version FE Computations \cite Dey1997a.
 * The order of the linear blending functions and symbol for face/edges
 * are named according Equation (30) of that paper.
 * \image html img/bounding-entity-labeling-curved-tet.png "Labeling for the curved edges and faces bounding the tetrahedron according to [Dey1997]" width=10cm
 */
class CrvTetBlending : public CrvTet {
 public:
    /** \brief the only supported ctor
     */
    CrvTetBlending(VtxPtrVec in_vert_vec, int order);

    // the destructor
    ~CrvTetBlending();

 private:
   
    /// container of a ordered list of mesh edges 
    CrvEdgePtrVec m_crv_edges;

    /// container of a ordered list of mesh faces
    CrvFacePtrVec m_crv_faces;

    /// shape handle for the 4 faces
//    std::vector<pTriGeom> m_face_shape_vector;

    // evaluate the mapping (by blending)
    virtual int v_eval(Point3d in_xi, Point3d & out_pos);
    // evaluate first derivative
    virtual int v_eval_deriv1(Point3d in_xi,
                              Mat3x3 & out_jac);
    // evaluate first derivative
    virtual int v_eval_deriv1_old(Point3d in_xi,
                                Mat3x3 & out_jac);

    /// set up internal edge storage and create Parametric Curve objects
    void setup_edges(int order);

    /// set up internal face storage and create Parametric Tri objects
    void setup_faces(int order);

    /// blending functions w.r.t. the 4 faces
    void face_blending_eval(int index,
                                   double in_xi_1,
                                   double in_xi_2,
                                   double in_xi_3,
                                   double & out_fb);

    void face_blending_eval_deriv1(int index,
                                   double in_xi_1,
                                   double in_xi_2,
                                   double in_xi_3,
                                   Point3d & dfb_dxi);

    void projected_face_coords(int face_index,
                                   double in_xi1,
                                   double in_xi2,
                                   double in_xi3,
                                   Point2d & out_xip);

    void projected_face_coords_deriv1(int face_index,
                                      double in_xi1,
                                      double in_xi2,
                                      double in_xi3,
                                      Point3d & dxip1_dxi,
                                      Point3d & dxip2_dxi);
                                   
  /// blending functions w.r.t. the 6 edges
    void edge_blending_eval(int index,
                                   double in_xi_1,
                                   double in_xi_2,
                                   double in_xi_3,
                                   double & out_eb);

    void edge_blending_eval_deriv1(int index,
                                   double in_xi_1,
                                   double in_xi_2,
                                   double in_xi_3,
                                   Point3d & deb_dxi);

    void projected_edge_coords(int edge_index,
                                   double in_xi1,
                                   double in_xi2,
                                   double in_xi3,
                                   Point1d & out_xip);

    /** \brief \f$ \frac{\xi'(\xi)}{\xi} \f$
        \f{eqnarray*}{
           \frac{\xi'_1}{\xi_1},  \frac{\xi'_1}{\xi_2},  \frac{\xi'_1}{\xi_3}
        \f}
     */
    void projected_edge_coords_deriv1(int edge_index,
                                      double in_xi1,
                                      double in_xi2,
                                      double in_xi3,
                                      Point3d & dxip1_dxi);
};
#endif//MESHADAPT_CRVTETBLENDING_H
