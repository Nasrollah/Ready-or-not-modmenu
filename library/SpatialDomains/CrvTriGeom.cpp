////////////////////////////////////////////////////////////////////////////////
//
//  File: CrvTriGeom.cpp
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

#include <SpatialDomains/CrvTriGeom.h>
#include <StdRegions/StdNodalTriExp.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/Interp.h>

#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/GeomFactorsPUMI.h>

// -- scorec -- //
#include "CrvTri1.h"
#include "CrvTri2.h"
#include "CrvTri3.h"
#include "CrvTri4C0.h"
#include "CrvTriBlending.h"
#include <PUMI_NektarppAdapt.h>
#include <pumi_mesh.h>
#include <vector>

using namespace pumi;

namespace Nektar
{
    namespace SpatialDomains
    {
        /**
         *
         */
        CrvTriGeom::CrvTriGeom()
        {
            m_shapeType = LibUtilities::eTriangle;
        }

        /**
         *
         */
        CrvTriGeom::CrvTriGeom(int id, const int coordim):
                                   TriGeom(id, coordim)
        {
        }


        /**
         *
         */
        CrvTriGeom::CrvTriGeom(const int id,
                const PointGeomSharedPtr verts[],
                const SegGeomSharedPtr edges[],
                const StdRegions::Orientation eorient[]):
                TriGeom(id, verts, edges, eorient)
        {
        }


        /**
         *
         */
        CrvTriGeom::CrvTriGeom(const int id, const SegGeomSharedPtr edges[],
                const StdRegions::Orientation eorient[], int order):
                TriGeom(id, edges, eorient), m_order(order)
        {
        }

        CrvTriGeom::CrvTriGeom(const int id,
                               const SegGeomSharedPtr edges[],
                               const StdRegions::Orientation eorient[],
                               CrvFace * in_crvface):
                TriGeom(id, edges, eorient),
                m_crvface(in_crvface)
        {
        }

        /**
         *
         */
        CrvTriGeom::CrvTriGeom(const int id,
                const SegGeomSharedPtr edges[],
                const StdRegions::Orientation eorient[],
                const CurveSharedPtr &curve) :
                TriGeom(id, edges, eorient, curve)
        {
        }


        /**
         *
         */
        CrvTriGeom::CrvTriGeom(const CrvTriGeom &in)
           : TriGeom(in)
        {
        }


        /**
         *
         */
        CrvTriGeom::~CrvTriGeom()
        {
        }

        CrvFace * CrvTriGeom::get_crv_tri()
        {
            return m_crvface;
        }
        /**
         * Set up GeoFac for this geometry using Coord quadrature distribution
         */
        void CrvTriGeom::v_GenGeomFactors()
        {
//            TriGeom::v_GenGeomFactors();
//            return;

            int tri_id = GetGlobalID();
//            std::cout << "The global ID of this CrvTri is " << tri_id << std::endl;
            if (m_geomFactorsState != ePtsFilled)
            {
                GeomType Gtype = eRegular;

                TriGeom::v_FillGeom();

                // check to see if expansions are linear
/*                for(int i = 0; i < m_coordim; ++i)
                {
                    if(m_xmap->GetBasisNumModes(0) != 2 ||
                       m_xmap->GetBasisNumModes(1) != 2)
                    {
                        Gtype = eDeformed;
                    }
                }
*/
                // skip the checking,
                // set to deformed directly
                // Kai -- 10/29/2014
                Gtype = eDeformed;
              // try to create a LagrangeTri object
              PointGeomSharedPtr verts[3] = {
                             GetVertex(0),
                             GetVertex(1),
                             GetVertex(2)
                                            };
              pMeshEnt pumi_verts[3];
              double coords[3][3];
              for(int i = 0; i < 3; ++i) {
                (verts[i])->GetCoords(coords[i][0], coords[i][1], coords[i][2]);
//              coords[i][0] = verts[i][0];  // x
//              coords[i][1] = verts[i][1];  // y
//              coords[i][2] = verts[i][2];  // z
                pumi_verts[i] = PUMI_NektarppAdapt::Instance()->get_pumi_mesh_vtx(coords[i][0], coords[i][1], coords[i][2]);
              }

/*
              std::cout << "Three vert coords: " << std::endl;
              std::cout << coords[0][0] << ", " << coords[0][1] << ", " << coords[0][2] << std::endl;
              std::cout << coords[1][0] << ", " << coords[1][1] << ", " << coords[1][2] << std::endl;
              std::cout << coords[2][0] << ", " << coords[2][1] << ", " << coords[2][2] << std::endl;
              std::cout << std::endl;
*/

              std::vector<pVertex> pumi_vert_vec;
              for (int i = 0; i < 3; ++i)
                pumi_vert_vec.push_back(pumi_verts[i]);

//              pMeshMdl pumi_mesh = PUMI_NektarppAdapt::Instance()->PumiMesh();
//              pPart pumi_part;
//              PUMI_Mesh_GetPart(pumi_mesh, 0, pumi_part);
              CrvFace * beztri = NULL;
              cout << "ORDER " << m_order << endl;
              switch (m_order){
                case 1:
                {
                  beztri = new CrvTri1(pumi_vert_vec);                  
                  break;
                }
                case 2:
                {
                  beztri = new CrvTri2(pumi_vert_vec);                  
                  break;
                }
                case 3:
                {
                  beztri = new CrvTri3(pumi_vert_vec);
                  break;
                }
                case 4:
                {
                  beztri = new CrvTri4C0(pumi_vert_vec);
                  break;
                }
                case 5:
                {
                  beztri = new CrvTriBlending(pumi_vert_vec,m_order);
                  break;
                }
              }
//              CrvFace * beztri = new CrvTri2(pumi_vert_vec);
//              CrvFace * beztri = new CrvTri3(pumi_vert_vec);
//              CrvFace * beztri = new CrvTri4C0(pumi_vert_vec);
              // CrvFace * beztri = new CrvTriBlending(pumi_vert_vec,m_order);
              std::cout << "interp error of CrvFace: " << beztri->interp_error() << std::endl;

              m_crvface = beztri;

/*
              std::cout << "PUMI compute: " << std::endl;

              double dxdxi = -0.5 * coords[0][0] + 0.5 * coords[1][0];
              double dydxi = -0.5 * coords[0][1] + 0.5 * coords[1][1];
              double dxdeta= -0.5 * coords[0][0] + 0.5 * coords[2][0];
              double dydeta= -0.5 * coords[0][1] + 0.5 * coords[2][1];
              std::cout << dxdxi << ", " << dydxi << std::endl;
              std::cout << dxdeta<< ", " << dydeta<< std::endl;
              std::cout << std::endl;
*/

                // creates a GeomFactorsPUMI object instead of original GeomFactors,  Kai -- 10-22-2014
                m_geomFactors = MemoryManager<GeomFactorsPUMI>::AllocateSharedPtr(
                    Gtype, m_coordim, m_xmap, m_coeffs, tri_id, beztri);

                m_geomFactorsState = ePtsFilled;
            }

        }

    }; //end of namespace
}; //end of namespace
