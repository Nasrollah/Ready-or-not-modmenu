////////////////////////////////////////////////////////////////////////////////
//
//  File:  SegGeom.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/CrvSegGeom.h>
#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/GeomFactorsPUMI.h>

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdSegExp.h>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc

// --- scorec --- //
#include <CrvEdge2.h>
#include <CrvEdge3.h>
#include <CrvEdge4.h>
#include <PUMI_NektarppAdapt.h>
#include <pumi_mesh.h>

using namespace pumi;

namespace Nektar
{
    namespace SpatialDomains
    {
        CrvSegGeom::CrvSegGeom()
        {
            m_shapeType = LibUtilities::eSegment;
        }

        CrvSegGeom::CrvSegGeom(int id, const int coordim):
            SegGeom(id, coordim)
        {
        }

        CrvSegGeom::CrvSegGeom(
                int id,
                const int coordim,
                const PointGeomSharedPtr vertex[]):
            SegGeom(id, coordim, vertex)
        {
        }

        CrvSegGeom::CrvSegGeom(
                int id,
                const int coordim,
                const PointGeomSharedPtr vertex[],
                const CurveSharedPtr& curve):
            SegGeom(id, coordim, vertex, curve)
        {
        }

        // Kai -- 11/3/2014 
        CrvSegGeom::CrvSegGeom(
                int id,
                const int coordim,
                const PointGeomSharedPtr vertex[],
                CrvEdge * crv_edge):
            m_crvedge(crv_edge),
//            SegGeom(id, coordim, vertex)
            SegGeom()
        {

            m_shapeType = LibUtilities::eSegment;
            m_globalID =  m_eid = id;
            m_state = eNotFilled;
            m_coordim = coordim;

            if (coordim > 0)
            {
                int npts = m_crvedge->get_geom_order() + 1; // including end vertices, e.g., 3 pts for a quadratic edge 
                LibUtilities::PointsKey pkey(npts+1,LibUtilities::eGaussLobattoLegendre);
                const LibUtilities::BasisKey B(LibUtilities::eModified_A, npts, pkey);
                Array<OneD,NekDouble> tmp(npts);

                m_xmap = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
                SetUpCoeffs(m_xmap->GetNcoeffs());

//                double v0[3];
//                double v1[3];
//                double mid[3];

//                m_crvedge->get_edge_node_coords(0, mid);
//                m_crvedge->get_vertex_coords(v0, v1);

                for(int i = 0; i < m_coordim; ++i)
                {
                    // set up tmp to store edge nodes including end vertices           
                    for(int j = 0; j < npts; ++j) {
                      Point3d in_pt, out;
                      double lcoord = 1.0 - 1.0 / (npts - 1) * j;
                      in_pt = Point3d(lcoord, 0.0, 0.0);
                      m_crvedge->eval_by_xi(in_pt, out);
                      tmp[j] = out[i];
                    }
                    // Interpolate to GLL points
                    DNekMatSharedPtr I0;

                    LibUtilities::PointsKey fkey(npts,LibUtilities::ePolyEvenlySpaced);
                    I0 = LibUtilities::PointsManager()[fkey]->GetI(pkey);

                    NekVector<NekDouble> in (npts, tmp, eWrapper);
                    NekVector<NekDouble> out(npts+1);
                    out = (*I0)*in;

                    m_xmap->FwdTrans(out.GetPtr(), m_coeffs[i]);
                }
            }

            m_verts[0] = vertex[0];
            m_verts[1] = vertex[1];  

        }


        CrvSegGeom::CrvSegGeom(
                const int id,
                const PointGeomSharedPtr& vert1,
                const PointGeomSharedPtr& vert2):
            SegGeom(id, vert1, vert2)
        {
        }

        CrvSegGeom::CrvSegGeom(const CrvSegGeom &in) :
          SegGeom(in)
        {
        }

        CrvSegGeom::~CrvSegGeom()
        {
        }

        CrvEdge * CrvSegGeom::get_crv_edge()
        {
          return m_crvedge;
        }


        void CrvSegGeom::v_GenGeomFactors()
        {
//            SegGeom::v_GenGeomFactors();
//            return;

            int seg_id = GetGlobalID();
            if (m_geomFactorsState != ePtsFilled)
            {
/*                SpatialDomains::GeomType gType = eRegular;

                std::cout << "SegGeom::v_GenGeomFactors() -- " << std::endl;
                std::cout << "         gType set to eRegular" << std::endl;

                SegGeom::v_FillGeom();

                if(m_xmap->GetBasisNumModes(0)!=2)
                {
                    gType = eDeformed;
                    std::cout << "SegGeom::v_GenGeomFactors() -- " << std::endl;
                    std::cout << "         gType set to eDeformed" << std::endl;
                }
*/
              SpatialDomains::GeomType gType = eDeformed;
/*              PointGeomSharedPtr verts[2] = {
                             GetVertex(0),
                             GetVertex(1)
                                            };

              pMeshEnt pumi_verts[2];
              double coords[2][3];
              for(int i = 0; i < 2; ++i) {
                (verts[i])->GetCoords(coords[i][0], coords[i][1], coords[i][2]);
                pumi_verts[i] = PUMI_NektarppAdapt::Instance()->get_pumi_mesh_vtx(coords[i][0], coords[i][1], coords[i][2]);
              }
              std::vector<pVertex> pumi_vert_vec;
              for (int i = 0; i < 2; ++i)
                pumi_vert_vec.push_back(pumi_verts[i]);

              CrvEdge2 * bezseg = new CrvEdge2(pumi_vert_vec);
*/
              m_geomFactors = MemoryManager<GeomFactorsPUMI>::AllocateSharedPtr(
                  gType, m_coordim, m_xmap, m_coeffs, seg_id, m_crvedge);


//                m_geomFactors = MemoryManager<GeomFactors>::AllocateSharedPtr(
//                    gType, m_coordim, m_xmap, m_coeffs);
              m_geomFactorsState = ePtsFilled;
            }
        }


    }; //end of namespace
}; //end of namespace

