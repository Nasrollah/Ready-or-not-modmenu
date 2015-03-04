/****************************************************************************** 

  (c) 2004-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE-SCOREC file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/

/**
 * @file CrvTriBlending.cc
 * @author Qiukai Lu
 */

#include "curveUtil.h"
#include "CrvTriBlending.h"
#include "CrvEdgeCAD.h"
#include "CrvEdge1.h"
#include "CrvEdge2.h"
#include "CrvEdge3.h"
#include "CrvEdge4.h"
#include "CrvEdge5.h"
CrvTriBlending::CrvTriBlending(pMeshEnt in_face, int order)
  : CrvFace(in_face) {
  setup_verts();
  setup_edges();
}

CrvTriBlending::CrvTriBlending(VtxPtrVec in_vert_vec, int order)
  : CrvFace(in_vert_vec)
{
  double coords[3];

  std::cout << "building tri blending of order " << order << std::endl;

  V_coord(in_vert_vec[0], coords);
  m_vertex_100 = Point3d(coords);
  V_coord(in_vert_vec[1], coords);
  m_vertex_010 = Point3d(coords);
  V_coord(in_vert_vec[2], coords);
  m_vertex_001 = Point3d(coords);

  VtxPtrVec edge_verts;
  edge_verts.push_back(in_vert_vec[0]);
  edge_verts.push_back(in_vert_vec[1]);

  switch (order) {
    case 1:
    {
      m_edge_P = new CrvEdge1(edge_verts);
      break;
    }
    case 2:
    {
      m_edge_P = new CrvEdge2(edge_verts);
      break;
    }
    case 3:
    {
      m_edge_P = new CrvEdge3(edge_verts);
      break;
    }
    case 4:
    {
      m_edge_P = new CrvEdge4(edge_verts);
      break;
    }
    case 5:
    {
      m_edge_P = new CrvEdge5(edge_verts);
      break;
    }
    default:
      assert(0);
      break;
  }
  edge_verts.clear();

  edge_verts.push_back(in_vert_vec[1]);
  edge_verts.push_back(in_vert_vec[2]);
  switch (order) {
    case 1:
    {
      m_edge_Q = new CrvEdge1(edge_verts);
      break;
    }
    case 2:
    {
      m_edge_Q = new CrvEdge2(edge_verts);
      break;
    }
    case 3:
    {
      m_edge_Q = new CrvEdge3(edge_verts);
      break;
    }
    case 4:
    {
      m_edge_Q = new CrvEdge4(edge_verts);
      break;
    }
    case 5:
    {
      m_edge_Q = new CrvEdge5(edge_verts);
      break;
    }
    default:
      assert(0);
      break;
  }
  edge_verts.clear();

  edge_verts.push_back(in_vert_vec[2]);
  edge_verts.push_back(in_vert_vec[0]);
  switch (order) {
    case 1:
    {
      m_edge_R = new CrvEdge1(edge_verts);
      break;
    }
    case 2:
    {
      m_edge_R = new CrvEdge2(edge_verts);
      break;
    }
    case 3:
    {
      m_edge_R = new CrvEdge3(edge_verts);
      break;
    }
    case 4:
    {
      m_edge_R = new CrvEdge4(edge_verts);
      break;
    }
    case 5:
    {
      m_edge_R = new CrvEdge5(edge_verts);
      break;
    }
    default:
      assert(0);
      break;
  }
  edge_verts.clear();
}

int CrvTriBlending::setup_verts() {
  // set up vertex coordinates
  pMeshVtx vert;
  double coords[3];

  vert = F_vertex(get_mesh_ent(), 0);
  V_coord(vert, coords);
  m_vertex_100 = Point3d(coords);
 
  vert = F_vertex(get_mesh_ent(), 1);
  V_coord(vert, coords);
  m_vertex_010 = Point3d(coords);

  vert = F_vertex(get_mesh_ent(), 2);
  V_coord(vert, coords);
  m_vertex_001 = Point3d(coords);
 
  return 0; 
}

int CrvTriBlending::setup_edges() {
  FIXME_LATER();
}

int CrvTriBlending::clear_edges() {
  delete m_edge_P;
  delete m_edge_Q;
  delete m_edge_R;
  return 0;
}

int CrvTriBlending::eval(double u, double v, double * pos) {
  // refer to Eq (28) of [Dey et al 1997]
  // u -- xi1
  // v -- xi2
  // w -- xi3
  double w = 1.0 - u - v; 
  double xi1 = u;
  double xi2 = v;
  double xi3 = w;

  double tol = 1.0e-6;

  if(fabs(1.0 - xi1) < tol) {
      m_vertex_100.toArray(pos);
  } else if(fabs(1.0 - xi2) < tol) {
      m_vertex_010.toArray(pos);
  } else if(fabs(1.0 - xi3) < tol) {
      m_vertex_001.toArray(pos);
  } else {
    Point3d in, tmp[3];
    in = Point3d(xi1 / (xi1 + xi2), 0.0, 0.0);
    m_edge_P->eval_by_xi(in, tmp[0]);
    in = Point3d(xi2 / (xi2 + xi3), 0.0, 0.0);
    m_edge_Q->eval_by_xi(in, tmp[1]);
    in = Point3d(xi3 / (xi3 + xi1), 0.0, 0.0);
    m_edge_R->eval_by_xi(in, tmp[2]);
  
    Point3d pt_pos = tmp[0] * (1.0 - xi3) * (1.0 - xi3)
                   + tmp[1] * (1.0 - xi1) * (1.0 - xi1)
                   + tmp[2] * (1.0 - xi2) * (1.0 - xi2)
                   - m_vertex_100 * xi1 * xi1
                   - m_vertex_010 * xi2 * xi2
                   - m_vertex_001 * xi3 * xi3;
    pt_pos.toArray(pos);
  }
  return 0;
}

int CrvTriBlending::eval_deriv1(double in_xi1,
                                double in_xi2,
                                Point3d & dxyz_dxi1,
                                Point3d & dxyz_dxi2)
{
  CrvTriBlending::eval_diff_analytic(in_xi1, in_xi2, dxyz_dxi1, dxyz_dxi2);
  return 0;
}

void CrvTriBlending::eval_diff_fd(double in_xi_1,
                          double in_xi_2,
                          Point3d & dxyz_dxi1,
                          Point3d & dxyz_dxi2)
{
  // approximate by forward differencing for now
  double xi1 = in_xi_1;
  double xi2 = in_xi_2;
  double xi3 = 1.0 - xi1 - xi2;
  double delta = 1.0e-3;
  double out1[3], out0[3];

  if(fabs(0.0 - xi3) > delta)// xi3 != 0.0 -> forward difference
  {
    // eval in xi1 direction
    CrvTriBlending::eval(xi1 + delta, xi2, out1);
    CrvTriBlending::eval(xi1,         xi2, out0);
    dxyz_dxi1 = (Point3d(out1) - Point3d(out0)) / delta;

    // eval in xi2 direction
    CrvTriBlending::eval(xi1, xi2 + delta, out1);
    CrvTriBlending::eval(xi1, xi2,         out0);
    dxyz_dxi2 = (Point3d(out1) - Point3d(out0)) / delta;
  }
  else { // xi3 == 0.0
    if(fabs(0.0 - xi1) < delta) // xi1 == 0.0
    {
      xi1 = 0.0 + delta;
      xi2 = 1.0 - delta;
    }
    if(fabs(0.0 - xi2) < delta) // xi2 == 0.0
    {
      xi1 = 1.0 - delta;
      xi2 = 0.0 + delta;
    }
    // use backward difference
    // eval in xi1 direction
    CrvTriBlending::eval(xi1,         xi2, out1);
    CrvTriBlending::eval(xi1 - delta, xi2, out0);
    dxyz_dxi1 = (Point3d(out1) - Point3d(out0)) / delta;

    // eval in xi2 direction
    CrvTriBlending::eval(xi1,              xi2, out1);
    CrvTriBlending::eval(xi1, xi2 - delta, out0);
    dxyz_dxi2 = (Point3d(out1) - Point3d(out0)) / delta;
  }
}

void CrvTriBlending::eval_diff_analytic(double in_xi1,
                                        double in_xi2,
                                        Point3d & dxyz_dxi1,
                                        Point3d & dxyz_dxi2)
{
  double xi1 = in_xi1;
  if(fabs(xi1 - 0.0) < 1.0e-4)
    xi1 = 0.0;
  double xi2 = in_xi2;
  if(fabs(xi2 - 0.0) < 1.0e-4)
    xi2 = 0.0;
  double xi3 = 1.0 - xi1 - xi2;
  if(fabs(xi3 - 0.0) < 1.0e-4)
    xi3 = 0.0;

  Point3d in, tmp, tmp_diff;
  dxyz_dxi1 = Point3d(0.0, 0.0, 0.0);
  dxyz_dxi2 = Point3d(0.0, 0.0, 0.0);
  // eval edge P
  if(fabs(xi3 - 1.0) > 1.0e-4) { // xi3 != 1.0
    in = Point3d(xi1/(xi1+xi2), 0.0, 0.0);
    m_edge_P->eval_by_xi(in, tmp);
    m_edge_P->diff_by_xi(in, 0, tmp_diff);

    dxyz_dxi1 = dxyz_dxi1 + tmp * 2.0 * (xi1 + xi2) + tmp_diff * xi2;
    dxyz_dxi2 = dxyz_dxi2 + tmp * 2.0 * (xi1 + xi2) - tmp_diff * xi1;
  }

  // eval edge Q
  if(fabs(xi1 - 1.0) > 1.0e-4) { // xi1 != 1.0
    in = Point3d(xi2/(xi2+xi3), 0.0, 0.0);
    m_edge_Q->eval_by_xi(in, tmp);
    m_edge_Q->diff_by_xi(in, 0, tmp_diff);

    dxyz_dxi1 = dxyz_dxi1 - tmp * 2.0 * (xi2 + xi3) + tmp_diff * xi2;
    dxyz_dxi2 = dxyz_dxi2 + tmp_diff * (xi2 + xi3);
  }

  // eval edge R
  if(fabs(xi2 - 1.0) > 1.0e-4) { // xi2 != 1.0
    in = Point3d(xi3/(xi3+xi1), 0.0, 0.0);
    m_edge_R->eval_by_xi(in, tmp);
    m_edge_R->diff_by_xi(in, 0, tmp_diff);

    dxyz_dxi1 = dxyz_dxi1 - tmp_diff * (xi1 + xi3);
    dxyz_dxi2 = dxyz_dxi2 - tmp * 2.0 * (xi1 + xi3) - tmp_diff * xi1;
  }

  // eval at vertices

  dxyz_dxi1 = dxyz_dxi1 - m_vertex_100 * 2.0 * xi1 + m_vertex_001 * 2.0 * xi3;
  dxyz_dxi2 = dxyz_dxi2 - m_vertex_010 * 2.0 * xi2 + m_vertex_001 * 2.0 * xi3;

}

void CrvTriBlending::eval_diff_analytic_c2(double in_xi1,
                                        double in_xi2,
                                        Point3d & dxyz_dxi1,
                                        Point3d & dxyz_dxi2)
{
  double xi1 = in_xi1;
  if(fabs(xi1 - 0.0) < 1.0e-4)
    xi1 = 0.0;
  double xi2 = in_xi2;
  if(fabs(xi2 - 0.0) < 1.0e-4)
    xi2 = 0.0;
  double xi3 = 1.0 - xi1 - xi2;
  if(fabs(xi3 - 0.0) < 1.0e-4)
    xi3 = 0.0;

  Point3d in, tmp, tmp_diff;
  dxyz_dxi1 = Point3d(0.0, 0.0, 0.0);
  dxyz_dxi2 = Point3d(0.0, 0.0, 0.0);
  // eval edge P
  if(fabs(xi3 - 1.0) > 1.0e-4) { // xi3 != 1.0
    in = Point3d(xi1/(xi1+xi2), 0.0, 0.0);
    m_edge_P->eval_by_xi(in, tmp);
    m_edge_P->diff_by_xi(in, 0, tmp_diff);

    dxyz_dxi1 = dxyz_dxi1 + tmp * 3.0 * (xi1 + xi2) * (xi1 + xi2) + tmp_diff * xi2 * (xi1 + xi2);
    dxyz_dxi2 = dxyz_dxi2 + tmp * 3.0 * (xi1 + xi2) * (xi1 + xi2) - tmp_diff * xi1 * (xi1 + xi2);
  }

  // eval edge Q
  if(fabs(xi1 - 1.0) > 1.0e-4) { // xi1 != 1.0
    in = Point3d(xi2/(xi2+xi3), 0.0, 0.0);
    m_edge_Q->eval_by_xi(in, tmp);
    m_edge_Q->diff_by_xi(in, 0, tmp_diff);

    dxyz_dxi1 = dxyz_dxi1 - tmp * 3.0 * (xi2 + xi3) * (xi2 + xi3) + tmp_diff * xi2 * (xi2 + xi3);
    dxyz_dxi2 = dxyz_dxi2 + tmp_diff * (xi2 + xi3) * (xi2 + xi3);
  }

  // eval edge R
  if(fabs(xi2 - 1.0) > 1.0e-4) { // xi2 != 1.0
    in = Point3d(xi3/(xi3+xi1), 0.0, 0.0);
    m_edge_R->eval_by_xi(in, tmp);
    m_edge_R->diff_by_xi(in, 0, tmp_diff);

    dxyz_dxi1 = dxyz_dxi1 - tmp_diff * (xi1 + xi3) * (xi1 + xi3);
    dxyz_dxi2 = dxyz_dxi2 - tmp * 3.0 * (xi1 + xi3) * (xi1 + xi3) - tmp_diff * xi1 * (xi1 + xi3);
  }

  // eval at vertices

  dxyz_dxi1 = dxyz_dxi1 - m_vertex_100 * 3.0 * xi1 * xi1 + m_vertex_001 * 3.0 * xi3 * xi3;
  dxyz_dxi2 = dxyz_dxi2 - m_vertex_010 * 3.0 * xi2 * xi2 + m_vertex_001 * 3.0 * xi3 * xi3;

}
