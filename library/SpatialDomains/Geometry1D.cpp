////////////////////////////////////////////////////////////////////////////////
//
//  File: Geometry1D.cpp
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
//  Description:  1D geometry information
//
//
////////////////////////////////////////////////////////////////////////////////
#include <SpatialDomains/Geometry1D.h>

namespace Nektar
{
    namespace SpatialDomains
    {

        Geometry1D::Geometry1D()
        {
        }

        Geometry1D::Geometry1D(const int coordim):
            Geometry(coordim)
        {
        }

        Geometry1D::~Geometry1D()
        {
        }

        PointGeomSharedPtr Geometry1D::GetVertex(const int i) const
        {
            return v_GetVertex(i);
        }

        LibUtilities::ShapeType Geometry1D::DetShapeType() const
        {
            return v_DetShapeType();
        }

        int Geometry1D::GetEid() const
        {
            return v_GetEid();
        }


        int Geometry1D::v_GetShapeDim() const
        {
            return 1;
        }


        int Geometry1D::v_GetEid() const 
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for expansion type geometries");
            return 0;
        }

        PointGeomSharedPtr Geometry1D::v_GetVertex(const int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for expansion type geometries");
            PointGeomSharedPtr returnval;
            return returnval;
        }

        int Geometry1D::v_GetVid(int i) const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for expansion type geometries");
            return 0;
        }

        LibUtilities::ShapeType Geometry1D::v_DetShapeType() const
        {
            NEKERROR(ErrorUtil::efatal,
                     "This function is only valid for expansion type geometries");
            return LibUtilities::eNoShapeType;
        }
    }; //end of namespace
}; //end of namespace
