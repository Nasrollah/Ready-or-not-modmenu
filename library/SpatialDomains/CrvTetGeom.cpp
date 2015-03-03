////////////////////////////////////////////////////////////////////////////////
//
//  File: CrvTetGeom.cpp
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
//  Description: Curved Tetrahedral geometry information.
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/CrvTetGeom.h>

#include <SpatialDomains/Geometry1D.h>
//#include <StdRegions/StdCrvTetExp.h>
#include <SpatialDomains/SegGeom.h>
//#include <iostream>

#include <SpatialDomains/GeomFactorsPUMI.h>
// -- scorec -- //
#include "CrvTet2.h"
#include "CrvTetBlending.h"
#include <PUMI_NektarppAdapt.h>
#include <pumi_mesh.h>
#include <vector>

using namespace pumi;


namespace Nektar
{
    namespace SpatialDomains
    {
        
        CrvTetGeom::CrvTetGeom()
        {
        }
        
        CrvTetGeom::CrvTetGeom(const TriGeomSharedPtr faces[]) : 
            TetGeom(faces)
        {
        }
        
        CrvTetGeom::~CrvTetGeom()
        {  
        }

        CrvTet * CrvTetGeom::get_crv_tet()
        {
            return m_crvtet;
        }

        void CrvTetGeom::v_GenGeomFactors()
        {
            // create pumi geom factors using crv tet 2 
            int tet_id = GetGlobalID();
//            std::cout << "The global ID of this CrvTet is " << tet_id << std::endl;
            if (m_geomFactorsState != ePtsFilled)
            {
                GeomType Gtype = eDeformed;

                TetGeom::v_FillGeom();
                PointGeomSharedPtr verts[4] = {
                             GetVertex(0),
                             GetVertex(1),
                             GetVertex(2),
                             GetVertex(3)
                                              };
                pMeshEnt pumi_verts[4];
                double coords[4][3];
                for(int i = 0; i < 4; ++i) {
                  (verts[i])->GetCoords(coords[i][0], coords[i][1], coords[i][2]);
                  pumi_verts[i] = PUMI_NektarppAdapt::Instance()->get_pumi_mesh_vtx(coords[i][0], coords[i][1], coords[i][2]);
                }

                std::vector<pVertex> pumi_vert_vec;
                for (int i = 0; i < 4; ++i)
                  pumi_vert_vec.push_back(pumi_verts[i]);
                m_crvtet = new CrvTetBlending(pumi_vert_vec);
            

                // creates a GeomFactorsPUMI object instead of original GeomFactors,  Kai -- 11-19-2014
                m_geomFactors = MemoryManager<GeomFactorsPUMI>::AllocateSharedPtr(
                    Gtype, m_coordim, m_xmap, m_coeffs, tet_id, m_crvtet);

                m_geomFactorsState = ePtsFilled;

            } // end of if
        }
    }; //end of namespace
}; //end of namespace
