////////////////////////////////////////////////////////////////////////////////
//
//  File: QuadGeom.cpp
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

#include <SpatialDomains/QuadGeom.h>
#include <LibUtilities/Foundations/Interp.h>

#include <StdRegions/StdQuadExp.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/Curve.hpp>
#include <SpatialDomains/GeomFactors.h>


namespace Nektar
{
    namespace SpatialDomains
    {
        /**
         *
         */
        QuadGeom::QuadGeom()
        {
            m_shapeType = LibUtilities::eQuadrilateral;
        }


        /**
         *
         */
        QuadGeom::QuadGeom(int id, const int coordim):
            Geometry2D(coordim)
        {

            const LibUtilities::BasisKey B(LibUtilities::eModified_A, 2,
            LibUtilities::PointsKey(3,LibUtilities::eGaussLobattoLegendre));

            m_globalID = m_fid = id;

            m_xmap = MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(B,B);
            SetUpCoeffs(m_xmap->GetNcoeffs());
        }

        /**
         *
         */
        QuadGeom::QuadGeom(const int id,
                           const PointGeomSharedPtr verts[],
                           const SegGeomSharedPtr edges[],
                           const StdRegions::Orientation eorient[]):
            Geometry2D(verts[0]->GetCoordim()),
            m_fid(id)
        {
            m_globalID = id;
            m_shapeType = LibUtilities::eQuadrilateral;

            /// Copy the vert shared pointers.
            m_verts.insert(m_verts.begin(), verts, verts+QuadGeom::kNverts);

            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+QuadGeom::kNedges);


            for (int j=0; j<kNedges; ++j)
            {
                m_eorient[j] = eorient[j];
            }

            m_coordim = verts[0]->GetCoordim();
            ASSERTL0(m_coordim > 1,
                     "Cannot call function with dim == 1");

            int order0  = max(edges[0]->GetBasis(0)->GetNumModes(),
                              edges[2]->GetBasis(0)->GetNumModes());
            int points0 = max(edges[0]->GetBasis(0)->GetNumPoints(),
                              edges[2]->GetBasis(0)->GetNumPoints());
            int order1  = max(edges[1]->GetBasis(0)->GetNumModes(),
                              edges[3]->GetBasis(0)->GetNumModes());
            int points1 = max(edges[1]->GetBasis(0)->GetNumPoints(),
                              edges[3]->GetBasis(0)->GetNumPoints());

            const LibUtilities::BasisKey B0(LibUtilities::eModified_A, order0,
                  LibUtilities::PointsKey(points0,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey B1(LibUtilities::eModified_A, order1,
                  LibUtilities::PointsKey(points1,LibUtilities::eGaussLobattoLegendre));

            m_xmap = MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(B0,B1);
            SetUpCoeffs(m_xmap->GetNcoeffs());
        }


        /**
         *
         */
        QuadGeom::QuadGeom(const int id,
                           const SegGeomSharedPtr edges[],
                           const StdRegions::Orientation eorient[],
                           const CurveSharedPtr &curve) :
            Geometry2D(edges[0]->GetVertex(0)->GetCoordim()),
            m_fid(id)
        {
            int j;

            m_globalID = m_fid;

            m_shapeType = LibUtilities::eQuadrilateral;

            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+QuadGeom::kNedges);

            for(j=0; j <kNedges; ++j)
            {
                if(eorient[j] == StdRegions::eForwards)
                {
                    m_verts.push_back(edges[j]->GetVertex(0));
                }
                else
                {
                    m_verts.push_back(edges[j]->GetVertex(1));
                }
            }

            for (j=0; j<kNedges; ++j)
            {
                m_eorient[j] = eorient[j];
            }

            m_coordim = edges[0]->GetVertex(0)->GetCoordim();
            ASSERTL0(m_coordim > 1,
                "Cannot call function with dim == 1");

            int order0  = max(edges[0]->GetBasis(0)->GetNumModes(),
                              edges[2]->GetBasis(0)->GetNumModes());
            int points0 = max(edges[0]->GetBasis(0)->GetNumPoints(),
                              edges[2]->GetBasis(0)->GetNumPoints());
            int order1  = max(edges[1]->GetBasis(0)->GetNumModes(),
                              edges[3]->GetBasis(0)->GetNumModes());
            int points1 = max(edges[1]->GetBasis(0)->GetNumPoints(),
                              edges[3]->GetBasis(0)->GetNumPoints());

            const LibUtilities::BasisKey B0(LibUtilities::eModified_A, order0,
                  LibUtilities::PointsKey(points0,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey B1(LibUtilities::eModified_A, order1,
                  LibUtilities::PointsKey(points1,LibUtilities::eGaussLobattoLegendre));

            m_xmap = MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(B0,B1);
            SetUpCoeffs(m_xmap->GetNcoeffs());

            for(int i = 0; i < m_coordim; ++i)
            {
                int npts = curve->m_points.size();
                int nEdgePts = (int)sqrt(static_cast<NekDouble>(npts));
                Array<OneD,NekDouble> tmp(npts);
                LibUtilities::PointsKey curveKey(nEdgePts, curve->m_ptype);

                // Sanity checks:
                // - Curved faces should have square number of points;
                // - Each edge should have sqrt(npts) points.
                ASSERTL0(nEdgePts*nEdgePts == npts,
                         "NUMPOINTS should be a square number");
                
                for (j = 0; j < kNedges; ++j)
                {
                    ASSERTL0(edges[j]->GetXmap()->GetNcoeffs() == nEdgePts,
                             "Number of edge points does not correspond "
                             "to number of face points.");
                }
                
                for (j = 0; j < npts; ++j)
                {
                    tmp[j] = (curve->m_points[j]->GetPtr())[i];
                }
                
                // Interpolate curve points to GLL points
                Array<OneD, NekDouble> tmp2(points0*points1);
                LibUtilities::Interp2D(curveKey,curveKey,tmp,
                                       B0.GetPointsKey(),B1.GetPointsKey(),
                                       tmp2);
                
                // Forwards transform to get coefficient space.
                m_xmap->FwdTrans(tmp2, m_coeffs[i]);
            }
        }


        /**
         *
         */
        QuadGeom::QuadGeom(const int id,
                           const SegGeomSharedPtr edges[],
                           const StdRegions::Orientation eorient[]):
            Geometry2D(edges[0]->GetVertex(0)->GetCoordim()),
            m_fid(id)
        {
            int j;

            m_globalID = m_fid;
            m_shapeType = LibUtilities::eQuadrilateral;

            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+QuadGeom::kNedges);

            for(j=0; j <kNedges; ++j)
            {
                if(eorient[j] == StdRegions::eForwards)
                {
                    m_verts.push_back(edges[j]->GetVertex(0));
                }
                else
                {
                    m_verts.push_back(edges[j]->GetVertex(1));
                }
            }

            for (j=0; j<kNedges; ++j)
            {
                m_eorient[j] = eorient[j];
            }

            m_coordim = edges[0]->GetVertex(0)->GetCoordim();
            ASSERTL0(m_coordim > 1,
                "Cannot call function with dim == 1");

            int order0  = max(edges[0]->GetBasis(0)->GetNumModes(),
                              edges[2]->GetBasis(0)->GetNumModes());
            int points0 = max(edges[0]->GetBasis(0)->GetNumPoints(),
                              edges[2]->GetBasis(0)->GetNumPoints());
            int order1  = max(edges[1]->GetBasis(0)->GetNumModes(),
                              edges[3]->GetBasis(0)->GetNumModes());
            int points1 = max(edges[1]->GetBasis(0)->GetNumPoints(),
                              edges[3]->GetBasis(0)->GetNumPoints());

            const LibUtilities::BasisKey B0(LibUtilities::eModified_A, order0,
                  LibUtilities::PointsKey(points0,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey B1(LibUtilities::eModified_A, order1,
                  LibUtilities::PointsKey(points1,LibUtilities::eGaussLobattoLegendre));

            m_xmap = MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(B0,B1);
            SetUpCoeffs(m_xmap->GetNcoeffs());
        }


        /**
         *
         */
        QuadGeom::QuadGeom(const QuadGeom &in)
        {
            // From Geometry
            m_shapeType = in.m_shapeType;

            // From QuadFaceComponent
            m_fid = in.m_fid;
            m_globalID = m_fid;
            m_ownVerts = in.m_ownVerts;
            std::list<CompToElmt>::const_iterator def;
            for(def = in.m_elmtMap.begin(); def != in.m_elmtMap.end(); def++)
            {
                m_elmtMap.push_back(*def);
            }

			// From QuadGeom
			m_verts = in.m_verts;
			m_edges = in.m_edges;
            for (int i = 0; i < kNedges; i++)
            {
                m_eorient[i] = in.m_eorient[i];
            }
            m_ownData = in.m_ownData;
        }


        /**
         *
         */
        QuadGeom::~QuadGeom()
        {
        }


        /**
         *
         */
        NekDouble QuadGeom::GetCoord(const int i,
                                  const Array<OneD, const NekDouble> &Lcoord)
        {
            ASSERTL1(m_state == ePtsFilled,
                "Geometry is not in physical space");

            Array<OneD, NekDouble> tmp(m_xmap->GetTotPoints());
            m_xmap->BwdTrans(m_coeffs[i], tmp);

            return m_xmap->PhysEvaluate(Lcoord, tmp);
        }

        StdRegions::Orientation QuadGeom::GetFaceOrientation(
            const QuadGeom &face1,
            const QuadGeom &face2)
        {
            return GetFaceOrientation(face1.m_verts, face2.m_verts);
        }

        /** 
         * Calculate the orientation of face2 to face1 (note this is
         * not face1 to face2!). 
         */ 
        StdRegions::Orientation QuadGeom::GetFaceOrientation(
            const PointGeomVector &face1,
            const PointGeomVector &face2)
        {
            int i, j, vmap[4] = {-1,-1,-1,-1};
            NekDouble x, y, z, x1, y1, z1, cx = 0.0, cy = 0.0, cz = 0.0;
           
            // For periodic faces, we calculate the vector between the centre
            // points of the two faces. (For connected faces this will be
            // zero). We can then use this to determine alignment later in the
            // algorithm.
            for (i = 0; i < 4; ++i)
            {
                cx += (*face2[i])(0) - (*face1[i])(0);
                cy += (*face2[i])(1) - (*face1[i])(1);
                cz += (*face2[i])(2) - (*face1[i])(2);
            }
            cx /= 4;
            cy /= 4;
            cz /= 4;

            // Now construct a mapping which takes us from the vertices of one
            // face to the other. That is, vertex j of face2 corresponds to
            // vertex vmap[j] of face1.
            for (i = 0; i < 4; ++i)
            {
                x = (*face1[i])(0);
                y = (*face1[i])(1);
                z = (*face1[i])(2);
                for (j = 0; j < 4; ++j)
                {
                    x1 = (*face2[j])(0)-cx;
                    y1 = (*face2[j])(1)-cy;
                    z1 = (*face2[j])(2)-cz;
                    if (sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y)+(z1-z)*(z1-z)) < 1e-8)
                    {
                        vmap[j] = i;
                        break;
                    }
                }
            }

            // Use the mapping to determine the eight alignment options between
            // faces.
            if (vmap[1] == (vmap[0]+1) % 4)
            {
                switch (vmap[0])
                {
                    case 0:
                        return StdRegions::eDir1FwdDir1_Dir2FwdDir2;
                        break;
                    case 1:
                        return StdRegions::eDir1BwdDir2_Dir2FwdDir1;
                        break;
                    case 2:
                        return StdRegions::eDir1BwdDir1_Dir2BwdDir2;
                        break;
                    case 3:
                        return StdRegions::eDir1FwdDir2_Dir2BwdDir1;
                        break;
                }
            }
            else
            {
                switch (vmap[0])
                {
                    case 0:
                        return StdRegions::eDir1FwdDir2_Dir2FwdDir1;
                        break;
                    case 1:
                        return StdRegions::eDir1BwdDir1_Dir2FwdDir2;
                        break;
                    case 2:
                        return StdRegions::eDir1BwdDir2_Dir2BwdDir1;
                        break;
                    case 3:
                        return StdRegions::eDir1FwdDir1_Dir2BwdDir2;
                        break;
                }
            }

            ASSERTL0(false, "unable to determine face orientation");
            return StdRegions::eDir1FwdDir1_Dir2FwdDir2;
       }


        /**
         *
         */
        void QuadGeom::v_AddElmtConnected(int gvo_id, int locid)
        {
            CompToElmt ee(gvo_id,locid);
            m_elmtMap.push_back(ee);
        }


        /**
         *
         */
        int QuadGeom::v_NumElmtConnected() const
        {
            return int(m_elmtMap.size());
        }


        /**
         *
         */
        bool QuadGeom::v_IsElmtConnected(int gvo_id, int locid) const
        {
            std::list<CompToElmt>::const_iterator def;
            CompToElmt ee(gvo_id,locid);

            def = find(m_elmtMap.begin(),m_elmtMap.end(),ee);

            // Found the element connectivity object in the list
            if(def != m_elmtMap.end())
            {
                return(true);
            }

            return(false);
        }


        /**
         *
         */
        int QuadGeom::v_GetFid() const
        {
            return m_fid;
        }


        /**
         *
         */
        int QuadGeom::v_GetCoordim() const
        {
            return m_coordim;
        }


        /**
         *
         */
        const LibUtilities::BasisSharedPtr QuadGeom::v_GetBasis(const int i)
        {
            return m_xmap->GetBasis(i);
        }


        /**
         *
         */
        const LibUtilities::BasisSharedPtr QuadGeom::v_GetEdgeBasis(const int i)
        {
            ASSERTL1(i <= 3,"edge is out of range");

            if (i == 0 || i == 2)
            {
                return m_xmap->GetBasis(0);
            }
            else
            {
                return m_xmap->GetBasis(1);
            }
        }

        /**
         *
         */
        NekDouble QuadGeom::v_GetCoord(const int i, const Array<OneD, const NekDouble> &Lcoord)
	{
            return GetCoord(i,Lcoord);
	}


        /**
         * Set up GeoFac for this geometry using Coord quadrature distribution
         */
        void QuadGeom::v_GenGeomFactors()
        {
            if (m_geomFactorsState != ePtsFilled)
            {
                int i;
                GeomType Gtype = eRegular;

                QuadGeom::v_FillGeom();

                // We will first check whether we have a regular or deformed
                // geometry. We will define regular as those cases where the
                // Jacobian and the metric terms of the derivative are constants
                // (i.e. not coordinate dependent)

                // Check to see if expansions are linear
                // If not linear => deformed geometry
                for(i = 0; i < m_coordim; ++i)
                {
                    if((m_xmap->GetBasisNumModes(0) != 2)||
                       (m_xmap->GetBasisNumModes(1) != 2))
                    {
                        Gtype = eDeformed;
                    }
                }

                // For linear expansions, the mapping from standard to local
                // element is given by the relation:
                // x_i = 0.25 * [ ( x_i^A + x_i^B + x_i^C + x_i^D)       +
                //                (-x_i^A + x_i^B + x_i^C - x_i^D)*xi_1  +
                //                (-x_i^A - x_i^B + x_i^C + x_i^D)*xi_2  +
                //                ( x_i^A - x_i^B + x_i^C - x_i^D)*xi_1*xi_2 ]
                //
                // The jacobian of the transformation and the metric terms
                // dxi_i/dx_j, involve only terms of the form dx_i/dxi_j (both
                // for coordim == 2 or 3). Inspecting the formula above, it can
                // be appreciated that the derivatives dx_i/dxi_j will be
                // constant, if the coefficient of the non-linear term is zero.
                //
                // That is why for regular geometry, we require
                //
                //     x_i^A - x_i^B + x_i^C - x_i^D = 0
                //
                // or equivalently
                //
                //     x_i^A - x_i^B = x_i^D - x_i^C
                //
                // This corresponds to quadrilaterals which are paralellograms.
                if(Gtype == eRegular)
                {
                    for(i = 0; i < m_coordim; i++)
                    {
                        if( fabs( (*m_verts[0])(i) - (*m_verts[1])(i) +
                                  (*m_verts[2])(i) - (*m_verts[3])(i) )
                                > NekConstants::kNekZeroTol )
                        {
                            Gtype = eDeformed;
                            break;
                        }
                    }
                }

                m_geomFactors = MemoryManager<GeomFactors>::AllocateSharedPtr(
                    Gtype, m_coordim, m_xmap, m_coeffs);
                m_geomFactorsState = ePtsFilled;
            }
        }


        /**
         *
         */
        void QuadGeom::v_SetOwnData()
        {
            m_ownData = true;
        }


        /**
         * Note verts and edges are listed according to anticlockwise
         * convention but points in _coeffs have to be in array format from
         * left to right.
         */
        void QuadGeom::v_FillGeom()
        {
            // check to see if geometry structure is already filled
            if(m_state != ePtsFilled)
            {
                int i,j,k;
                int nEdgeCoeffs;

                Array<OneD, unsigned int> mapArray;
                Array<OneD, int>          signArray;

                for(i = 0; i < kNedges; i++)
                {
                    m_edges[i]->FillGeom();
                    m_xmap->GetEdgeToElementMap(i,m_eorient[i],
                                                mapArray,signArray);

                    nEdgeCoeffs = m_edges[i]->GetXmap()->GetNcoeffs();

                    for(j = 0; j < m_coordim; j++)
                    {
                        for(k = 0; k < nEdgeCoeffs; k++)
                        {
                            m_coeffs[j][mapArray[k]]
                                = signArray[k]*(m_edges[i]->GetCoeffs(j))[k];
                        }
                    }
                }

                m_state = ePtsFilled;
            }
        }
        
        /**
         *
         */
        NekDouble QuadGeom::v_GetLocCoords(const Array<OneD, const NekDouble> &coords, 
                                      Array<OneD,NekDouble> &Lcoords)
        {
            NekDouble resid = 0.0;
            if(GetMetricInfo()->GetGtype() == eRegular)
            { 
                NekDouble coords2 = (m_coordim == 3)? coords[2]: 0.0; 
                PointGeom dv1, dv2, norm, orth1, orth2;
                PointGeom xin(m_coordim,0,coords[0],coords[1],coords2);

                // Calculate edge vectors from 0-1 and 0-3 edges. 
                dv1.Sub(*m_verts[1],*m_verts[0]);
                dv2.Sub(*m_verts[3],*m_verts[0]);

                // Obtain normal to plane in which dv1 and dv2 lie
                norm.Mult(dv1,dv2);
                
                // Obtain vector which are normal to dv1 and dv2. 
                orth1.Mult(norm,dv1);
                orth2.Mult(norm,dv2);
                
                // Start with vector of desired points minus vertex_0
                xin -= *m_verts[0];

                // Calculate length using L/|dv1| = (x-v0).n1/(dv1.n1) for coordiante 1
                // Then rescale to [-1,1]. 
                Lcoords[0] = xin.dot(orth2)/dv1.dot(orth2);
                Lcoords[0] = 2*Lcoords[0]-1;
                Lcoords[1] = xin.dot(orth1)/dv2.dot(orth1);
                Lcoords[1] = 2*Lcoords[1]-1;
            }
            else
            {
                QuadGeom::v_FillGeom();                       
                
                // Determine nearest point of coords  to values in m_xmap
                int npts = m_xmap->GetTotPoints();
                Array<OneD, NekDouble> ptsx(npts), ptsy(npts);
                Array<OneD, NekDouble> tmpx(npts), tmpy(npts);

                m_xmap->BwdTrans(m_coeffs[0], ptsx);
                m_xmap->BwdTrans(m_coeffs[1], ptsy);
                
                const Array<OneD, const NekDouble> za = m_xmap->GetPoints(0);
                const Array<OneD, const NekDouble> zb = m_xmap->GetPoints(1);

                //guess the first local coords based on nearest point
                Vmath::Sadd(npts, -coords[0], ptsx,1,tmpx,1);
                Vmath::Sadd(npts, -coords[1], ptsy,1,tmpy,1);
                Vmath::Vmul (npts, tmpx,1,tmpx,1,tmpx,1);
                Vmath::Vvtvp(npts, tmpy,1,tmpy,1,tmpx,1,tmpx,1);
                          
                int min_i = Vmath::Imin(npts,tmpx,1);
                
                Lcoords[0] = za[min_i%za.num_elements()];
                Lcoords[1] = zb[min_i/za.num_elements()];

                // Perform newton iteration to find local coordinates 
                NewtonIterationForLocCoord(coords, ptsx, ptsy, Lcoords,resid);
            }
            return resid;
        }
            
            
        /**
         *
         */
        int QuadGeom::v_GetEid(int i) const
        {
            ASSERTL2((i >=0) && (i <= 3),"Edge id must be between 0 and 3");
            return m_edges[i]->GetEid();
        }


        /**
         *
         */
        int QuadGeom::v_GetVid(int i) const
        {
            ASSERTL2((i >=0) && (i <= 3),"Verted id must be between 0 and 3");
            return m_verts[i]->GetVid();
        }


        /**
         *
         */
        PointGeomSharedPtr QuadGeom::v_GetVertex(int i) const
        {
            ASSERTL2((i >=0) && (i <= 3),"Vertex id must be between 0 and 3");
            return m_verts[i];
        }


        /**
         *
         */
        const Geometry1DSharedPtr QuadGeom::v_GetEdge(int i) const
        {
            ASSERTL2((i >=0) && (i <= 3),"Edge id must be between 0 and 3");
            return m_edges[i];
        }


        /**
         *
         */
        StdRegions::Orientation QuadGeom::v_GetEorient(const int i) const
        {
            ASSERTL2((i >=0) && (i <= 3),"Edge id must be between 0 and 3");
            return m_eorient[i];
        }


        /**
         *
         */
        StdRegions::Orientation QuadGeom::v_GetCartesianEorient(const int i) const
        {
            ASSERTL2((i >=0) && (i <= 3),"Edge id must be between 0 and 3");
            if(i < 2)
            {
                return m_eorient[i];
            }
            else
            {
                if(m_eorient[i] == StdRegions::eForwards)
                {
                    return StdRegions::eBackwards;
                }
                else
                {
                    return StdRegions::eForwards;
                }
            }
        }


        /** 
         *
         */
        int QuadGeom::v_WhichEdge(SegGeomSharedPtr edge)
        {
            int returnval = -1;

            SegGeomVector::iterator edgeIter;
            int i;

            for (i=0,edgeIter = m_edges.begin(); edgeIter != m_edges.end(); ++edgeIter,++i)
            {
                if (*edgeIter == edge)
                {
                    returnval = i;
                    break;
                }
            }

            return returnval;
        }


        /**
         *
         */
        int QuadGeom::v_GetNumVerts() const
        {
            return kNverts;
        }


        /**
         *
         */
        int QuadGeom::v_GetNumEdges() const
        {
            return kNedges;
        }
 

        bool QuadGeom::v_ContainsPoint(
            const Array<OneD, const NekDouble> &gloCoord, NekDouble tol)
        {
            Array<OneD,NekDouble> locCoord(GetCoordim(),0.0);
            return v_ContainsPoint(gloCoord,locCoord,tol);

        }
        /**
         *
         */
        bool QuadGeom::v_ContainsPoint(const Array<OneD, const NekDouble> &gloCoord,
                                       Array<OneD, NekDouble> &stdCoord,
                                       NekDouble tol)
        {
            NekDouble resid;
            return v_ContainsPoint(gloCoord,stdCoord,tol,resid);
        }

        bool QuadGeom::v_ContainsPoint(const Array<OneD, const NekDouble> &gloCoord,
                                       Array<OneD, NekDouble> &stdCoord,
                                       NekDouble tol,
                                       NekDouble &resid)
        {
            ASSERTL1(gloCoord.num_elements() >= 2,
                     "Two dimensional geometry expects at least two coordinates.");
            
            resid = GetLocCoords(gloCoord, stdCoord);
            if (stdCoord[0] >= -(1+tol) && stdCoord[1] >= -(1+tol)
                && stdCoord[0] <= (1+tol) && stdCoord[1] <= (1+tol))
            {
                return true;
            }
            return false;
        }
    }; //end of namespace
}; //end of namespace
