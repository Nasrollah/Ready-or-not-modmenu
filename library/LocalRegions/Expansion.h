///////////////////////////////////////////////////////////////////////////////
//
// File Expansion.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Header file for Expansion routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPANSION_H
#define EXPANSION_H

#include <StdRegions/StdExpansion.h>
#include <SpatialDomains/Geometry.h>
#include <SpatialDomains/GeomFactors.h>
#include <LocalRegions/LocalRegionsDeclspec.h>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <vector>
#include <map>

namespace Nektar
{
    namespace LocalRegions 
    {

        class Expansion;
        class MatrixKey;

        enum MetricType
        {
            MetricLaplacian00,
            MetricLaplacian01,
            MetricLaplacian02,
            MetricLaplacian11,
            MetricLaplacian12,
            MetricLaplacian22,
            MetricQuadrature
        };

        // type defines for use of PrismExp in a boost vector
        typedef boost::shared_ptr<Expansion> ExpansionSharedPtr;
        typedef boost::weak_ptr<Expansion> ExpansionWeakPtr;
        typedef std::vector< ExpansionSharedPtr > ExpansionVector;
        typedef std::vector< ExpansionSharedPtr >::iterator ExpansionVectorIter;
        typedef std::map<MetricType, Array<OneD, NekDouble> > MetricMap;

        class Expansion : virtual public StdRegions::StdExpansion
        {
            public:
                LOCAL_REGIONS_EXPORT Expansion(SpatialDomains::GeometrySharedPtr pGeom); // default constructor.
                LOCAL_REGIONS_EXPORT Expansion(const Expansion &pSrc); // copy constructor.
                LOCAL_REGIONS_EXPORT virtual ~Expansion();

                LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr GetLocMatrix(const LocalRegions::MatrixKey &mkey);

                LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr GetLocMatrix(const StdRegions::MatrixType mtype,
                            const StdRegions::ConstFactorMap &factors = StdRegions::NullConstFactorMap,
                            const StdRegions::VarCoeffMap &varcoeffs = StdRegions::NullVarCoeffMap);

                LOCAL_REGIONS_EXPORT SpatialDomains::GeometrySharedPtr GetGeom() const;

                LOCAL_REGIONS_EXPORT virtual const
                    SpatialDomains::GeomFactorsSharedPtr& v_GetMetricInfo() const;

                LOCAL_REGIONS_EXPORT DNekMatSharedPtr BuildTransformationMatrix(
                    const DNekScalMatSharedPtr &r_bnd, 
                    const StdRegions::MatrixType matrixType);

                LOCAL_REGIONS_EXPORT DNekMatSharedPtr BuildVertexMatrix(
                    const DNekScalMatSharedPtr &r_bnd);
			
                LOCAL_REGIONS_EXPORT void AddEdgeNormBoundaryInt(
                    const int                           edge,
                    const boost::shared_ptr<Expansion> &EdgeExp,
                    const Array<OneD, const NekDouble> &Fx,
                    const Array<OneD, const NekDouble> &Fy,
                          Array<OneD,       NekDouble> &outarray);
                LOCAL_REGIONS_EXPORT void AddEdgeNormBoundaryInt(
                    const int                           edge,
                    const boost::shared_ptr<Expansion> &EdgeExp,
                    const Array<OneD, const NekDouble> &Fn,
                          Array<OneD,       NekDouble> &outarray);
                LOCAL_REGIONS_EXPORT void AddFaceNormBoundaryInt(
                    const int                           face,
                    const boost::shared_ptr<Expansion> &FaceExp,
                    const Array<OneD, const NekDouble> &Fn,
                          Array<OneD,       NekDouble> &outarray);
                LOCAL_REGIONS_EXPORT void DGDeriv(
                    const int                                   dir,
                    const Array<OneD, const NekDouble>&         inarray,
                          Array<OneD, ExpansionSharedPtr>      &EdgeExp,
                          Array<OneD, Array<OneD, NekDouble> > &coeffs,
                          Array<OneD,             NekDouble>   &outarray);

            protected:
                SpatialDomains::GeometrySharedPtr  m_geom;
                SpatialDomains::GeomFactorsSharedPtr m_metricinfo;
                MetricMap m_metrics;

                void ComputeLaplacianMetric();
                void ComputeQuadratureMetric();

                virtual void v_MultiplyByQuadratureMetric(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray);

                virtual void v_ComputeLaplacianMetric() {};

                virtual void v_GetCoords(Array<OneD,NekDouble> &coords_1,
                                         Array<OneD,NekDouble> &coords_2,
                                         Array<OneD,NekDouble> &coords_3);

                virtual DNekScalMatSharedPtr v_GetLocMatrix(const LocalRegions::MatrixKey &mkey);

                virtual DNekMatSharedPtr v_BuildTransformationMatrix(
                    const DNekScalMatSharedPtr &r_bnd, 
                    const StdRegions::MatrixType matrixType);

                virtual DNekMatSharedPtr v_BuildVertexMatrix(
                    const DNekScalMatSharedPtr &r_bnd); 

                virtual void v_AddEdgeNormBoundaryInt(
                    const int                           edge,
                    const boost::shared_ptr<Expansion> &EdgeExp,
                    const Array<OneD, const NekDouble> &Fx,
                    const Array<OneD, const NekDouble> &Fy,
                          Array<OneD,       NekDouble> &outarray);
                virtual void v_AddEdgeNormBoundaryInt(
                    const int                           edge,
                    const boost::shared_ptr<Expansion> &EdgeExp,
                    const Array<OneD, const NekDouble> &Fn,
                          Array<OneD,       NekDouble> &outarray);
                virtual void v_AddFaceNormBoundaryInt(
                    const int                           face,
                    const boost::shared_ptr<Expansion> &FaceExp,
                    const Array<OneD, const NekDouble> &Fn,
                          Array<OneD,       NekDouble> &outarray);
                virtual void v_DGDeriv(
                    const int                                   dir,
                    const Array<OneD, const NekDouble>&         inarray,
                          Array<OneD, ExpansionSharedPtr>      &EdgeExp,
                          Array<OneD, Array<OneD, NekDouble> > &coeffs,
                          Array<OneD,             NekDouble>   &outarray);

            private:

        };
    } //end of namespace
} //end of namespace

#endif
