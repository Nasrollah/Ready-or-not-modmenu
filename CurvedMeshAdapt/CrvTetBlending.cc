/****************************************************************************** 

  (c) 2004-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE-SCOREC file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "CrvTetBlending.h"
#include "LagrangeCurve.h"
#include "CADCurve.h"
#include "ParametricFace.h"
#include "LagrangeTri.h"
#include "CrvEdge.h"
#include "CrvEdge1.h"
#include "CrvEdge2.h"
#include "CrvEdge3.h"
#include "CrvEdge4.h"
#include "CrvEdge5.h"
#include "CrvEdgeCAD.h"
#include "CrvTri2.h"
#include "CrvTriBlending.h"
#include "curveMesh.h"
#include <pumi_mesh.h>

CrvTetBlending::CrvTetBlending(VtxPtrVec in_vert_vec, int order)
  : CrvTet(in_vert_vec)
{
  setup_edges(order);
  setup_faces(order);
}

CrvTetBlending::~CrvTetBlending()
{
}

/// set up internal edge storage and create Parametric Curve objects
void CrvTetBlending::setup_edges(int order) {
  // order of the edge shape objects is consistent with Dey's paper:
  // 0 => Edge W;
  // 1 => Edge T;
  // 2 => Edge Q;
  // 3 => Edge S;
  // 4 => Edge R;
  // 5 => Edge P;
  
  pVertex vert[4];
  for (int i = 0; i < 4; ++i)
    vert[i] = m_vert_vec[i];

  int edge_verts[6][2] = {
                             {2, 3},  // edge W
                             {1, 3},  // edge T
                             {1, 2},  // edge Q
                             {0, 3},  // S
                             {0, 2},  // R
                             {0, 1}   // P
                         };
  for (int i = 0; i < 6; ++i) {
      CrvEdge * crv_edge;
      int vid[2] = {
                       edge_verts[i][0],
                       edge_verts[i][1]
                   };
      VtxPtrVec verts;
      verts.push_back(vert[vid[0]]);
      verts.push_back(vert[vid[1]]);
      if (edge_on_model_bdry(vert[vid[0]], vert[vid[1]])) {
          crv_edge = new CrvEdge5(verts);
      } else {
          crv_edge = new CrvEdge5(verts);
      }
      m_crv_edges.push_back(crv_edge);
  }
}

/// set up internal face storage and create Parametric Tri objects
void CrvTetBlending::setup_faces(int order) {
  // order of the face shape objects is consistent with Dey's paper:
  // 0 => Face G;
  // 1 => Face E;
  // 2 => Face F;
  // 3 => Face D;

  pVertex vert[4];
  for(int i = 0; i < 4; ++i) {
    vert[i] = m_vert_vec[i];
  }

  VtxPtrVec verts;
  verts.clear();
  
  // set up face G
  pTriGeom face_g;

  verts.push_back(vert[1]);
  verts.push_back(vert[2]);
  verts.push_back(vert[3]);
  CrvFace * g = new CrvTriBlending(verts,order);
  m_crv_faces.push_back(g);
  verts.clear();

  // set up face E
  pTriGeom face_e;

  verts.push_back(vert[0]);
  verts.push_back(vert[2]);
  verts.push_back(vert[3]);
  CrvFace * e = new CrvTriBlending(verts,order);
  m_crv_faces.push_back(e);
  verts.clear();

  // set up face F
  pTriGeom face_f;

  verts.push_back(vert[0]);
  verts.push_back(vert[1]);
  verts.push_back(vert[3]);
  CrvFace * f = new CrvTriBlending(verts,order);
  m_crv_faces.push_back(f);
  verts.clear();

  // set up face D
  pTriGeom face_d;

  verts.push_back(vert[0]);
  verts.push_back(vert[1]);
  verts.push_back(vert[2]);
  CrvFace * d = new CrvTriBlending(verts,order);
  m_crv_faces.push_back(d);
  verts.clear();
}

int CrvTetBlending::v_eval(Point3d in_xi, Point3d & pos) {
  double xi1 = in_xi[0];
  double xi2 = in_xi[1];
  double xi3 = in_xi[2];
  double xi4 = 1.0 - xi1 - xi2 - xi3;
  pos = Point3d(0.0, 0.0, 0.0);

  // evaluate faces
  for(int i = 0; i < 4; ++i) {
    // evaluate blending function first
    double blend_eval;
    face_blending_eval(i, xi1, xi2, xi3, blend_eval);
    if (fabs(0.0 - blend_eval) > 1.0e-6) {
      Point2d fxip;
      projected_face_coords(i, xi1, xi2, xi3, fxip);
      Point3d in_fxip = Point3d(fxip[0], fxip[1], 0.0);
      Point3d face_eval;
      m_crv_faces[i]->eval_by_xi(in_fxip, face_eval);
      pos = pos + face_eval * blend_eval;
    }
    // if blending eval is zero, then no need to calculate face contribution
  }


  // evaluate edges
  for(int i = 0; i < 6; ++i) {
    // evaluate blending function first
    double blend_eval;
    edge_blending_eval(i, xi1, xi2, xi3, blend_eval);
    if (fabs(0.0 - blend_eval) > 1.0e-6) {
      Point1d exip;
      projected_edge_coords(i, xi1, xi2, xi3, exip);
      Point3d edge_eval;
      Point3d in(exip, 0.0, 0.0);
      m_crv_edges[i]->eval_by_xi(in, edge_eval);
      pos = pos - edge_eval * blend_eval;
    }
    // if blending eval is zero, then no need to calculate edge contribution
  }

  // evaluate vertices
  pVertex vert[4];
  double xi[4] = {xi1, xi2, xi3, xi4};
  for(int i = 0; i < 4; ++i) {
    vert[i] = m_vert_vec[i];
    double vcoords[3];
    V_coord(vert[i], vcoords);
    Point3d v(vcoords);
    pos = pos + v * xi[i] * xi[i];
  }

  return 0; 
}

/// the most up-do-date impl of derivatives eval based on blending
int CrvTetBlending::v_eval_deriv1(Point3d in_xi,
                                    Mat3x3 & dxyz_dxi)
{
  double tol = 1.0e-6;

  double xi1 = in_xi[0];
  double xi2 = in_xi[1];
  double xi3 = in_xi[2];
  double xi4 = 1.0 - xi1 - xi2 - xi3;
  // reset output
  for(int i = 0; i < 3; ++i)
  for(int j = 0; j < 3; ++j)
      dxyz_dxi[i][j] = 0.0;

  Point2d fxi; // projected face coords
  Point3d fxyz; // face eval
  Point3d df_dfxi1, df_dfxi2;

  // contribution from face G
  if (fabs(1.0 - xi1) > tol) {  // xi1 != 1.0
      projected_face_coords(0, xi1, xi2, xi3, fxi);
      Point3d fxi_pt = Point3d(fxi[0], fxi[1], 0.0);
      m_crv_faces[0]->eval_by_xi(fxi_pt, fxyz);
      m_crv_faces[0]->diff_by_xi(fxi_pt, 0, df_dfxi1);
      m_crv_faces[0]->diff_by_xi(fxi_pt, 1, df_dfxi2);
      for (int i = 0; i < 3; ++i) {
          // dxyz_dxi1
          dxyz_dxi[i][0] = dxyz_dxi[i][0] + fxyz[i] * (-2.0) * (1.0 - xi1)
                         + df_dfxi1[i] * xi2 + df_dfxi2[i] * xi3;
          // dxyz_dxi2
          dxyz_dxi[i][1] = dxyz_dxi[i][1] + (1.0 - xi1) * df_dfxi1[i];
          // dxyz_dxi3
          dxyz_dxi[i][2] = dxyz_dxi[i][2] + (1.0 - xi1) * df_dfxi2[i];
      }
  }

  // contribution from face E
  if (fabs(1.0 - xi2) > tol) {  // xi2 != 1.0
      projected_face_coords(1, xi1, xi2, xi3, fxi);
      Point3d fxi_pt = Point3d(fxi[0], fxi[1], 0.0);
      m_crv_faces[1]->eval_by_xi(fxi_pt, fxyz);
      m_crv_faces[1]->diff_by_xi(fxi_pt, 0, df_dfxi1);
      m_crv_faces[1]->diff_by_xi(fxi_pt, 1, df_dfxi2);
      for (int i = 0; i < 3; ++i) {
          // dxyz_dxi1
          dxyz_dxi[i][0] = dxyz_dxi[i][0] + (1.0 - xi2) * df_dfxi1[i];
          // dxyz_dxi2
          dxyz_dxi[i][1] = dxyz_dxi[i][1] + (-2.0) * (1.0 - xi2) * fxyz[i]
                         + xi3 * df_dfxi2[i] + xi1 * df_dfxi1[i];
          // dxyz_dxi3
          dxyz_dxi[i][2] = dxyz_dxi[i][2] + (1.0 - xi2) * df_dfxi2[i];
      }
  }

  // contribution from face F
  if (fabs(1.0 - xi3) > tol) {  // xi3 != 1.0
      projected_face_coords(2, xi1, xi2, xi3, fxi);
      Point3d fxi_pt = Point3d(fxi[0], fxi[1], 0.0);
      m_crv_faces[2]->eval_by_xi(fxi_pt, fxyz);
      m_crv_faces[2]->diff_by_xi(fxi_pt, 0, df_dfxi1);
      m_crv_faces[2]->diff_by_xi(fxi_pt, 1, df_dfxi2);
      for (int i = 0; i < 3; ++i) {
          // dxyz_dxi1
          dxyz_dxi[i][0] = dxyz_dxi[i][0] + (1.0 - xi3) * df_dfxi1[i];
          // dxyz_dxi2
          dxyz_dxi[i][1] = dxyz_dxi[i][1] + (1.0 - xi3) * df_dfxi2[i];
          // dxyz_dxi3
          dxyz_dxi[i][2] = dxyz_dxi[i][2] + (-2.0) * (1.0 - xi3) * fxyz[i]
                         + xi2 * df_dfxi2[i] + xi1 * df_dfxi1[i];
      }
  }

  // contribution from face D
  if (fabs(1.0 - xi4) > tol) {  // xi4 != 1.0
      projected_face_coords(3, xi1, xi2, xi3, fxi);
      Point3d fxi_pt = Point3d(fxi[0], fxi[1], 0.0);
      m_crv_faces[3]->eval_by_xi(fxi_pt, fxyz);
      m_crv_faces[3]->diff_by_xi(fxi_pt, 0, df_dfxi1);
      m_crv_faces[3]->diff_by_xi(fxi_pt, 1, df_dfxi2);
      for (int i = 0; i < 3; ++i) {
          // dxyz_dxi1
          dxyz_dxi[i][0] = dxyz_dxi[i][0] + 2.0 * (xi1 + xi2 + xi3) * fxyz[i]
                         - xi2 * df_dfxi2[i] + (xi2 + xi3) * df_dfxi1[i];
          // dxyz_dxi2
          dxyz_dxi[i][1] = dxyz_dxi[i][1] + 2.0 * (xi1 + xi2 + xi3) * fxyz[i]
                         + (xi1 + xi3) * df_dfxi2[i] - xi1 * df_dfxi1[i];
          // dxyz_dxi3
          dxyz_dxi[i][2] = dxyz_dxi[i][2] + 2.0 * (xi1 + xi2 + xi3) * fxyz[i]
                         - xi2 * df_dfxi2[i] - xi1 * df_dfxi1[i];
      }
  }

  Point1d exi;
  Point3d exyz;
  Point3d de_dexi;
  // contribution from edge W
  if (fabs(1.0 - xi1 - xi2) > tol) {  // xi1 + xi2 != 1.0
      projected_edge_coords(0, xi1, xi2, xi3, exi);
      Point3d exi_pt = Point3d(exi, 0.0, 0.0);
      m_crv_edges[0]->eval_by_xi(exi_pt, exyz);
      m_crv_edges[0]->diff_by_xi(exi_pt, 0, de_dexi);
      for (int i = 0; i < 3; ++i) {
          // dxyz_dxi1
          dxyz_dxi[i][0] = dxyz_dxi[i][0] - ( (-2.0) * (1.0 - xi1 - xi2) * exyz[i] + xi3 * de_dexi[i] );
          // dxyz_dxi2
          dxyz_dxi[i][1] = dxyz_dxi[i][1] - ( (-2.0) * (1.0 - xi1 - xi2) * exyz[i] + xi3 * de_dexi[i]  );
          // dxyz_dxi3
          dxyz_dxi[i][2] = dxyz_dxi[i][2] - ( (1.0 - xi1 - xi2) * de_dexi[i] );
      }
  }

  // contribution from edge T
  if (fabs(1.0 - xi1 - xi3) > tol) {  // xi1 + xi3 != 1.0
      projected_edge_coords(1, xi1, xi2, xi3, exi);
      Point3d exi_pt = Point3d(exi, 0.0, 0.0);
      m_crv_edges[1]->eval_by_xi(exi_pt, exyz);
      m_crv_edges[1]->diff_by_xi(exi_pt, 0, de_dexi);
      for (int i = 0; i < 3; ++i) {
          // dxyz_dxi1
          dxyz_dxi[i][0] = dxyz_dxi[i][0] - ( (-2.0) * (1.0 - xi1 - xi3) * exyz[i] + xi2 * de_dexi[i] );
          // dxyz_dxi2
          dxyz_dxi[i][1] = dxyz_dxi[i][1] - ( (1.0 - xi1 - xi3) * de_dexi[i] );
          // dxyz_dxi3
          dxyz_dxi[i][2] = dxyz_dxi[i][2] - ( (-2.0) * (1.0 - xi1 - xi3) * exyz[i] + xi2 * de_dexi[i] );
      }
  }

  // contribution from edge Q
  if (fabs(xi2 + xi3) > tol) {  // xi2 + xi3 != 0.0
      projected_edge_coords(2, xi1, xi2, xi3, exi);
      Point3d exi_pt = Point3d(exi, 0.0, 0.0);
      m_crv_edges[2]->eval_by_xi(exi_pt, exyz);
      m_crv_edges[2]->diff_by_xi(exi_pt, 0, de_dexi);
      for (int i = 0; i < 3; ++i) {
          // dxyz_dxi1
          dxyz_dxi[i][0] = dxyz_dxi[i][0];
          // dxyz_dxi2
          dxyz_dxi[i][1] = dxyz_dxi[i][1] - ( 2.0 * (xi2 + xi3) * exyz[i] + xi3 * de_dexi[i] );
          // dxyz_dxi3
          dxyz_dxi[i][2] = dxyz_dxi[i][2] - ( 2.0 * (xi2 + xi3) * exyz[i] - xi2 * de_dexi[i] );
      }
  }

  // contribution from edge S
  if (fabs(1.0 - xi2 - xi3) > tol) {  // xi2 + xi3 != 1.0
      projected_edge_coords(3, xi1, xi2, xi3, exi);
      Point3d exi_pt = Point3d(exi, 0.0, 0.0);
      m_crv_edges[3]->eval_by_xi(exi_pt, exyz);
      m_crv_edges[3]->diff_by_xi(exi_pt, 0, de_dexi);
      for (int i = 0; i < 3; ++i) {
          // dxyz_dxi1
          dxyz_dxi[i][0] = dxyz_dxi[i][0] - ( (1.0 - xi2 - xi3) * de_dexi[i] );
          // dxyz_dxi2
          dxyz_dxi[i][1] = dxyz_dxi[i][1] - ( (-2.0) * (1.0 - xi2 - xi3) * exyz[i] + xi1 * de_dexi[i] );
          // dxyz_dxi3
          dxyz_dxi[i][2] = dxyz_dxi[i][2] - ( (-2.0) * (1.0 - xi2 - xi3) * exyz[i] + xi1 * de_dexi[i] );
      }
  }

  // contribution from edge R
  if (fabs(xi1 + xi3) > tol) {  // xi1 + xi3 != 0.0
      projected_edge_coords(4, xi1, xi2, xi3, exi);
      Point3d exi_pt = Point3d(exi, 0.0, 0.0);
      m_crv_edges[4]->eval_by_xi(exi_pt, exyz);
      m_crv_edges[4]->diff_by_xi(exi_pt, 0, de_dexi);
      for (int i = 0; i < 3; ++i) {
          // dxyz_dxi1
          dxyz_dxi[i][0] = dxyz_dxi[i][0] - ( 2.0 * (xi1 + xi3) * exyz[i] + xi3 * de_dexi[i] );
          // dxyz_dxi2
          dxyz_dxi[i][1] = dxyz_dxi[i][1];
          // dxyz_dxi3
          dxyz_dxi[i][2] = dxyz_dxi[i][2] - ( 2.0 * (xi1 + xi3) * exyz[i] - xi1 * de_dexi[i] );
      }
  }

  // contribution from edge P
  if (fabs(xi1 + xi2) > tol) {  // xi1 + xi2 != 0.0
      projected_edge_coords(5, xi1, xi2, xi3, exi);
      Point3d exi_pt = Point3d(exi, 0.0, 0.0);
      m_crv_edges[5]->eval_by_xi(exi_pt, exyz);
      m_crv_edges[5]->diff_by_xi(exi_pt, 0, de_dexi);
      for (int i = 0; i < 3; ++i) {
          // dxyz_dxi1
          dxyz_dxi[i][0] = dxyz_dxi[i][0] - ( 2.0 * (xi1 + xi2) * exyz[i] + xi2 * de_dexi[i] );
          // dxyz_dxi2
          dxyz_dxi[i][1] = dxyz_dxi[i][1] - ( 2.0 * (xi1 + xi2) * exyz[i] - xi1 * de_dexi[i] );
          // dxyz_dxi3
          dxyz_dxi[i][2] = dxyz_dxi[i][2];
      }
  }

  // contribution from vertices
  double vcoords[4][3];
  Point3d v[4];
  for(int i = 0; i < 4; ++i) {
    V_coord(m_vert_vec[i], vcoords[i]);
    v[i] = Point3d(vcoords[i]);
  }

  for (int i = 0; i < 3; ++i) {
      // dxyz_dxi1
      dxyz_dxi[i][0] = dxyz_dxi[i][0] + 2.0 * xi1 * v[0][i]
                     - 2.0 * (1.0 - xi1 - xi2 - xi3) * v[3][i];
      // dxyz_dxi2
      dxyz_dxi[i][1] = dxyz_dxi[i][1] + 2.0 * xi2 * v[1][i]
                     - 2.0 * (1.0 - xi1 - xi2 - xi3) * v[3][i];
      // dxyz_dxi3
      dxyz_dxi[i][2] = dxyz_dxi[i][2] + 2.0 * xi3 * v[2][i]
                     - 2.0 * (1.0 - xi1 - xi2 - xi3) * v[3][i];
  }
}

/// evaluate first derivative -- set to retire
int CrvTetBlending::v_eval_deriv1_old(Point3d in_xi,
                                  Mat3x3 & dxyz_dxi) {
  double xi1 = in_xi[0];
  double xi2 = in_xi[1];
  double xi3 = in_xi[2];
//  double xi4 = 1.0 - xi1 - xi2 - xi3;
  // reset output
  for(int i = 0; i < 3; ++i)
  for(int j = 0; j < 3; ++j)
    dxyz_dxi[i][j] = 0.0;

  // evaluate faces
  double face_blend;
//  double dfb_dxi1, dfb_dxi2, dfb_dxi3;
  Point3d dfb_dxi;
  Point2d fxip; 
  Point3d face_eval;
  Point3d df_dxip1, df_dxip2;
  Point3d dxip1_dxi;
  Point3d dxip2_dxi;
  for(int iface = 0 ; iface < 4; ++iface) {
    // get face blend
    face_blending_eval(iface, xi1, xi2, xi3, face_blend);
    if (fabs(0.0 - face_blend) > 1.0e-6) {
      // get face blend deriv;
      face_blending_eval_deriv1(iface, xi1, xi3, xi3, dfb_dxi);
      // get face coord projection
      projected_face_coords(iface, xi1, xi2, xi3, fxip);
      // get face eval
//      m_face_shape_vector[iface]->eval(fxip, face_eval);
      Point3d in_fxip = Point3d(fxip[0], fxip[1], 0.0);
      m_crv_faces[iface]->eval_by_xi(in_fxip, face_eval);
      // get face eval deriv
      Mat3x3 df_dxip;
      m_crv_faces[iface]->diff_by_xi(in_fxip, df_dxip);
      for (int idim = 0; idim < 3; ++idim) {
          df_dxip1[idim] = df_dxip[idim][0];
          df_dxip2[idim] = df_dxip[idim][1];
      }
//      m_face_shape_vector[iface]->deriv1(fxip, df_dxip1, df_dxip2);
      // get coord projection deriv
      projected_face_coords_deriv1(iface, xi1, xi2, xi3, dxip1_dxi, dxip2_dxi);
      for(int i = 0; i < 3; ++i) {
          for(int j = 0; j < 3; ++j) {
              // face blend deriv * face eval + face blend * face deriv * project coord deriv
              dxyz_dxi[i][j] = dxyz_dxi[i][j] 
                       + dfb_dxi[j] * face_eval[i] 
                       + (dxip1_dxi[j] * df_dxip1[i] + dxip2_dxi[j] * df_dxip2[i]) * face_blend;
          }
      }
    }
  }


  // evaluate edges
  double edge_blend;
  Point3d deb_dxi;
  Point1d exip;  // xi1 prime for edge
  Point3d edge_eval;
  Point3d de_dxip1;
  Point3d dexip1_dxi;

  for(int iedge = 0 ; iedge < 6; ++iedge) {
    // get edge blend
    edge_blending_eval(iedge, xi1, xi2, xi3, edge_blend);

    if (fabs(0.0 - edge_blend) > 1.0e-6) {
      // get edge blend deriv;
      edge_blending_eval_deriv1(iedge, xi1, xi3, xi3, deb_dxi);

      // get edge coord projection
      projected_edge_coords(iedge, xi1, xi2, xi3, exip);
    
      // get edge eval
      Point3d in(exip, 0.0, 0.0);
      m_crv_edges[iedge]->eval_by_xi(in, edge_eval);
      // get edge eval deriv
//      m_edge_shape_vector[iedge]->deriv1(exip, de_dxip1);
      m_crv_edges[iedge]->diff_by_xi(in, 0, de_dxip1);
      // get coord projection deriv
      projected_edge_coords_deriv1(iedge, xi1, xi2, xi3, dexip1_dxi);

      for(int xyz = 0 ; xyz < 3; ++xyz) {
      // edge blend deriv * edge eval + edge blend * edge deriv * project coord deriv
        for (int xii = 0; xii < 3; ++xii) {
          dxyz_dxi[xyz][xii] = dxyz_dxi[xyz][xii]
              - (deb_dxi[xii] * edge_eval[xyz] + dexip1_dxi[xii] * de_dxip1[xyz] * edge_blend);
        }
      }
    }
  }

  // evaluate vertex derivatives here
  double vcoords[4][3];
  Point3d v[4];
  for(int i = 0; i < 4; ++i) {
    V_coord(m_vert_vec[i], vcoords[i]);
    v[i] = Point3d(vcoords[i]);
  }

  for (int xyz = 0; xyz < 3; ++xyz)
  for (int xi = 0; xi < 3; ++xi)
    dxyz_dxi[xyz][xi] = dxyz_dxi[xyz][xi] + v[xi][xyz] - v[3][xyz];


}

void CrvTetBlending::face_blending_eval(int index,
                                        double xi1,
                                        double xi2,
                                        double xi3,
                                        double & out) {
  // face blend uses quadratic functions to ensure c1 intra continuity
  assert(index >=0 && index <= 3);
  double xi4 = 1.0 - xi1 - xi2 -xi3;
  switch (index) {
    case 0:  // blend for G
      out = (1.0 - xi1) * (1.0 - xi1);
      break;
    case 1:  // blend for E
      out = (1.0 - xi2) * (1.0 - xi2);
      break;
    case 2:  // blend for F
      out = (1.0 - xi3) * (1.0 - xi3);
      break;
    case 3:  // blend for D
      out = (1.0 - xi4) * (1.0 - xi4);
      break;
    default:
      break;
  }
}

void CrvTetBlending::face_blending_eval_deriv1(int index,
                      double xi1, double xi2, double xi3,
                      Point3d & df_dxi) {
  assert(index >= 0 && index <= 3);
//  double xi4 = 1.0 - xi1 - xi2 - xi3;
  switch (index) {
    case 0:  // blend for G
      df_dxi = Point3d(-1.0, 0.0, 0.0);
      break;
    case 1:  // blend for E
      df_dxi = Point3d(0.0, -1.0, 0.0);
      break;
    case 2:  // blend for F
      df_dxi = Point3d(0.0, 0.0, -1.0);
      break;
    case 3:  // blend for D
      df_dxi = Point3d(1.0, 1.0, 1.0);
      break;
    default:
      break;
  }
}

void CrvTetBlending::projected_face_coords(int face_index,
                                               double xi1,
                                               double xi2,
                                               double xi3,
                                               Point2d & out) {
  assert(face_index >=0 && face_index <= 3);
  double xi4 = 1.0 - xi1 - xi2 - xi3;
  switch (face_index) {
    case 0:  // face G
      out[0] = xi2 / (xi2 + xi3 + xi4);
      out[1] = xi3 / (xi2 + xi3 + xi4);
      break;
    case 1:  // face E
      out[0] = xi1 / (xi1 + xi3 + xi4);
      out[1] = xi3 / (xi1 + xi3 + xi4);
      break;
    case 2:  // face F
      out[0] = xi1 / (xi1 + xi2 + xi4);
      out[1] = xi2 / (xi1 + xi2 + xi4);
      break;
    case 3:  // face D
      out[0] = xi1 / (xi1 + xi2 + xi3);
      out[1] = xi2 / (xi1 + xi2 + xi3);
      break;
    default:
      break;
  }
}

void CrvTetBlending::projected_face_coords_deriv1(int face_index,
                   double xi1, double xi2, double xi3,
                   Point3d & dxip1_dxi, Point3d & dxip2_dxi) {

  assert(face_index >=0 && face_index <= 3);
//  double xi4 = 1.0 - xi1 - xi2 - xi3;
  switch (face_index) {
    case 0:  // face G
      dxip1_dxi[0] = xi2 / pow((xi1 - 1.0), 2);
      dxip1_dxi[1] = 1.0 / (1.0 - xi1);
      dxip1_dxi[2] = 0.0;

      dxip2_dxi[0] = xi3 / pow((xi1 - 1.0), 2);
      dxip2_dxi[1] = 0.0;
      dxip2_dxi[2] = 1.0 / (1.0 - xi1);

      break;
    case 1:  // face E
      dxip1_dxi[0] = 1.0 / (1.0 - xi2);
      dxip1_dxi[1] = xi1 / pow((xi2 - 1.0), 2);
      dxip1_dxi[2] = 0.0;

      dxip2_dxi[0] = 0.0;
      dxip2_dxi[1] = xi3 / pow((xi2 - 1.0), 2);
      dxip2_dxi[2] = 1.0 / (1.0 - xi2);

      break;
    case 2:  // face F
      dxip1_dxi[0] = 1.0 / (1.0 - xi3);
      dxip1_dxi[1] = 0.0;
      dxip1_dxi[2] = xi1 / pow((xi3 - 1.0), 2);

      dxip2_dxi[0] = 0.0;
      dxip2_dxi[1] = 1.0 / (1.0 - xi3);
      dxip2_dxi[2] = xi2 / pow((xi3 - 1.0), 2);

      break;
    case 3:  // face D
      dxip1_dxi[0] = (xi2 + xi3) / pow((xi1 + xi2 + xi3), 2);
      dxip1_dxi[1] = -xi1 / pow((xi1 + xi2 + xi3), 2);
      dxip1_dxi[2] = -xi1 / pow((xi1 + xi2 + xi3), 2);

      dxip2_dxi[0] = -xi2 / pow((xi1 + xi2 + xi3), 2);
      dxip2_dxi[1] = (xi1 + xi3) / pow((xi1 + xi2 + xi3), 2);
      dxip2_dxi[2] = -xi2 / pow((xi1 + xi2 + xi3), 2);

      break;
    default:
      break;
  }

}

void CrvTetBlending::edge_blending_eval(int index,
                                   double xi1,
                                   double xi2,
                                   double xi3,
                                   double & out) {
  // edge blend uses quadratic functions to ensure c1
  assert(index >= 0 && index <= 5);
  double xi4 = 1.0 - xi1 - xi2 - xi3;
  switch (index) {
    case 0:  // edge W
      out = (1.0 - xi1 - xi2) * (1.0 - xi1 - xi2);
      break;
    case 1:  // edge T
      out = (1.0 - xi1 - xi3) * (1.0 - xi1 - xi3);
      break;
    case 2:  // edge Q
      out = (1.0 - xi1 - xi4) * (1.0 - xi1 - xi4);
      break;
    case 3:  // edge S
      out = (1.0 - xi2 - xi3) * (1.0 - xi2 - xi3);
      break;
    case 4:  // edge R
      out = (1.0 - xi2 - xi4) * (1.0 - xi2 - xi4);
      break;
    case 5:  // edge P
      out = (1.0 - xi3 - xi4) * (1.0 - xi3 - xi4);
      break;
    default:
      break;
  }
}

void CrvTetBlending::edge_blending_eval_deriv1(int index,
                                   double xi1,
                                   double xi2,
                                   double xi3,
                                   Point3d & de_dxi) {
  assert(index >= 0 && index <= 5);
//  double xi4 = 1.0 - xi1 - xi2 - xi3;
  switch (index) {
    case 0:  // edge W
      de_dxi[0] = -1.0;
      de_dxi[1] = -1.0;
      de_dxi[2] = 0.0;
      break;
    case 1:  // edge T
      de_dxi[0] = -1.0;
      de_dxi[1] = 0.0;
      de_dxi[2] = -1.0;
      break;
    case 2:  // edge Q
      de_dxi[0] = 0.0;
      de_dxi[1] = 1.0;
      de_dxi[2] = 1.0;
      break;
    case 3:  // edge S
      de_dxi[0] = 0.0;
      de_dxi[1] = -1.0;
      de_dxi[2] = -1.0;
      break;
    case 4:  // edge R
      de_dxi[0] = 1.0;
      de_dxi[1] = 0.0;
      de_dxi[2] = 1.0;
      break;
    case 5:  // edge P
      de_dxi[0] = 1.0;
      de_dxi[1] = 1.0;
      de_dxi[2] = 0.0;
      break;
    default:
      break;
  }
}

void CrvTetBlending::projected_edge_coords(int index,
                                   double xi1,
                                   double xi2,
                                   double xi3,
                                   double & xi1p) {
  assert(index >= 0 && index <= 5);
  double xi4 = 1.0 - xi1 - xi2 - xi3;
  switch (index) {
    case 0:  // edge W
      xi1p = xi3 / (xi3 + xi4);
      break;
    case 1:  // edge T
      xi1p = xi2 / (xi2 + xi4);
      break;
    case 2:  // edge Q
      xi1p = xi2 / (xi2 + xi3);
      break;
    case 3:  // edge S
      xi1p = xi1 / (xi1 + xi4);
      break;
    case 4:  // edge R
      xi1p = xi1 / (xi1 + xi3);
      break;
    case 5:  // edge P
      xi1p = xi1 / (xi1 + xi2);
      break;
    default:
      break;
  }

}

void CrvTetBlending::projected_edge_coords_deriv1(int index,
                                      double xi1,
                                      double xi2,
                                      double xi3,
                                      Point3d & dxip_dxi) {
  assert(index >= 0 && index <= 5);
//  double xi4 = 1.0 - xi1 - xi2 - xi3;
  switch (index) {
    case 0:  // edge W
      dxip_dxi[0] = xi3 / pow((xi1 + xi2 - 1.0), 2);
      dxip_dxi[1] = xi3 / pow((xi1 + xi2 - 1.0), 2);
      dxip_dxi[2] = 1.0 / (1.0 - xi1 - xi2);
      break;
    case 1:  // edge T
      dxip_dxi[0] = xi2 / pow((xi1 + xi3 - 1.0), 2);
      dxip_dxi[1] = 1.0 / (1.0 - xi1 - xi3);
      dxip_dxi[2] = xi2 / pow((xi1 + xi3 - 1.0), 2);
      break;
    case 2:  // edge Q
      dxip_dxi[0] = 0.0;
      dxip_dxi[1] = xi3 / pow((xi2 + xi3), 2);
      dxip_dxi[2] = -xi2 / pow((xi2 + xi3), 2);
      break;
    case 3:  // edge S
      dxip_dxi[0] = 1.0 / (1.0 - xi2 - xi3);
      dxip_dxi[1] = xi1 / pow((xi2 + xi3 - 1.0), 2);
      dxip_dxi[2] = xi1 / pow((xi2 + xi3 - 1.0), 2);
      break;
    case 4:  // edge R
      dxip_dxi[0] = xi3 / pow((xi1 + xi3), 2);
      dxip_dxi[1] = 0.0;
      dxip_dxi[2] = -xi1 / pow((xi1 + xi3), 2);
      break;
    case 5:  // edge P
      dxip_dxi[0] = xi2 / pow((xi1 + xi2), 2);
      dxip_dxi[1] = -xi1 / pow((xi1 + xi2), 2);
      dxip_dxi[2] = 0.0;
      break;
    default:
      break;
  }

}
