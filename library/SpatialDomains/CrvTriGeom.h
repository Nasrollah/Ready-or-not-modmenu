////////////////////////////////////////////////////////////////////////////////
//
//  File:  CrvTriGeom.h
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
#ifndef NEKTAR_SPATIALDOMAINS_CRVTRIGEOM_H
#define NEKTAR_SPATIALDOMAINS_CRVTRIGEOM_H

#include <StdRegions/StdRegions.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <SpatialDomains/Geometry2D.h>
#include <SpatialDomains/PointGeom.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/TriGeom.h>

#include <SpatialDomains/GeomFactors.h>
#include <StdRegions/StdTriExp.h>

class CrvTri2;
class CrvTriBlending;
class CrvFace;

namespace Nektar
{
    namespace SpatialDomains
    {
        class  CrvTriGeom;
        class  SegGeom;
        struct Curve;

        typedef boost::shared_ptr<Curve> CurveSharedPtr;
        typedef boost::shared_ptr<SegGeom> SegGeomSharedPtr;
        typedef boost::shared_ptr<CrvTriGeom> CrvTriGeomSharedPtr;
        typedef std::vector< SegGeomSharedPtr > SegGeomVector;
        typedef std::vector< CrvTriGeomSharedPtr > CrvTriGeomVector;
        typedef std::vector< CrvTriGeomSharedPtr >::iterator CrvTriGeomVectorIter;
        typedef std::map<int, CrvTriGeomSharedPtr> CrvTriGeomMap;
        typedef std::map<int, CrvTriGeomSharedPtr>::iterator CrvTriGeomMapIter;

        class CrvTriGeom: public TriGeom
        {
            public:
                SPATIAL_DOMAINS_EXPORT CrvTriGeom();

                SPATIAL_DOMAINS_EXPORT CrvTriGeom(int id, const int coordim);

                SPATIAL_DOMAINS_EXPORT CrvTriGeom(
                        const int id,
                        const PointGeomSharedPtr verts[],
                        const SegGeomSharedPtr edges[],
                        const StdRegions::Orientation eorient[]);

                SPATIAL_DOMAINS_EXPORT CrvTriGeom(
                        const int id,
                        const SegGeomSharedPtr edges[],
                        const StdRegions::Orientation eorient[],
                        int order);

                /// added new ctor, Kai -- 11/6/2014
                SPATIAL_DOMAINS_EXPORT CrvTriGeom(
                        const int id,
                        const SegGeomSharedPtr edges[],
                        const StdRegions::Orientation eorient[],
                        CrvFace * in_crvface);

                SPATIAL_DOMAINS_EXPORT CrvTriGeom(
                        const int id,
                        const SegGeomSharedPtr edges[],
                        const StdRegions::Orientation eorient[],
                        const CurveSharedPtr &curve);

                SPATIAL_DOMAINS_EXPORT CrvTriGeom(const CrvTriGeom &in);

                SPATIAL_DOMAINS_EXPORT ~CrvTriGeom();
                // re-implement gen geom factors
                SPATIAL_DOMAINS_EXPORT virtual void v_GenGeomFactors();

                SPATIAL_DOMAINS_EXPORT CrvFace * get_crv_tri();

            private:

                CrvFace * m_crvface;
                int m_order;

        };

    }; //end of namespace SpatialDomains
}; //end of namespace Nektar

#endif
