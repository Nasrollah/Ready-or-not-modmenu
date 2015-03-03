////////////////////////////////////////////////////////////////////////////////
//
//  File:  GeomFactors.h
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
//  Description: Geometric Factors base class
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_GEOMFACTORSPUMI_H
#define NEKTAR_SPATIALDOMAINS_GEOMFACTORSPUMI_H

#include <boost/unordered_set.hpp>

#include <LibUtilities/Foundations/Basis.h>
#include <SpatialDomains/SpatialDomains.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <SpatialDomains/GeomFactors.h>
#include <StdRegions/StdExpansion.h>
#include <StdRegions/StdRegions.hpp>

//#include <CrvTri2.h>
class CrvEnt; 
// ============================================================= //

// -- Kai, 10-22-2014
namespace Nektar
{
namespace SpatialDomains
{
    class GeomFactorsPUMI : public GeomFactors {
      public:
        /// Constructor for GeomFactorsPUMI class.
        GeomFactorsPUMI(const GeomType                              gtype,
                    const int                                   coordim,
                    const StdRegions::StdExpansionSharedPtr    &xmap,
                    const Array<OneD, Array<OneD, NekDouble> > &coords,
                    const int globalid,
                    CrvEnt * crvent);

        /// Copy constructor.
        GeomFactorsPUMI(const GeomFactorsPUMI &S);

        /// Tests if two GeomFactors classes are equal.
        SPATIAL_DOMAINS_EXPORT friend bool operator==(
            const GeomFactorsPUMI &lhs,
            const GeomFactorsPUMI &rhs);
        /// Destructor.
        SPATIAL_DOMAINS_EXPORT virtual ~GeomFactorsPUMI();
        /// re-impl ComputeDeriv
        SPATIAL_DOMAINS_EXPORT virtual DerivStorage ComputeDeriv(
                    const LibUtilities::PointsKeyVector &keyTgt) const;

        /// global id
        int m_globalId;
        CrvEnt * m_crvEnt;
    
    };
} // end of SpatialDomains
} // end of Nektar
#endif
