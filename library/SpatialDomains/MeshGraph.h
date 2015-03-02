////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshGraph.h
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
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_MESHGRAPH_H
#define NEKTAR_SPATIALDOMAINS_MESHGRAPH_H

#include <boost/unordered_map.hpp>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/FieldIO.h>

#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/TetGeom.h>
#include <SpatialDomains/PyrGeom.h>
#include <SpatialDomains/PrismGeom.h>
#include <SpatialDomains/HexGeom.h>

#include <SpatialDomains/Curve.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

class TiXmlDocument;

namespace Nektar
{
    namespace SpatialDomains
    {
        enum ExpansionType
        {
            eNoExpansionType,
            eModified,
            eModifiedQuadPlus1,
            eModifiedQuadPlus2,
            eOrthogonal,
            eGLL_Lagrange,
            eGLL_Lagrange_SEM,
            eGauss_Lagrange,
            eGauss_Lagrange_SEM,
            eFourier,
            eFourierSingleMode,
            eFourierHalfModeRe,
            eFourierHalfModeIm,
            eChebyshev,
            eFourierChebyshev,
            eChebyshevFourier,
            eFourierModified,
            eExpansionTypeSize
        };

        // Keep this consistent with the enums in ExpansionType.
        // This is used in the BC file to specify the expansion type.
        const std::string kExpansionTypeStr[] =
        {
            "NOTYPE",
            "MODIFIED",
            "MODIFIEDQUADPLUS1",
            "MODIFIEDQUADPLUS2",
            "ORTHOGONAL",
            "GLL_LAGRANGE",
            "GLL_LAGRANGE_SEM",
            "GAUSS_LAGRANGE",
            "GAUSS_LAGRANGE_SEM",
            "FOURIER",
            "FOURIERSINGLEMODE",
            "FOURIERHALFMODERE",
            "FOURIERHALFMODEIM",
            "CHEBYSHEV",
            "FOURIER-CHEBYSHEV",
            "CHEBYSHEV-FOURIER",
            "FOURIER-MODIFIED"
        };

        class InterfaceComponent;
        typedef boost::shared_ptr< InterfaceComponent > SharedInterfaceCompPtr;
        typedef std::vector< PointGeomSharedPtr >       PointGeomVector;
        typedef std::map<int, PointGeomSharedPtr>       PointGeomMap;
        typedef std::list< SharedInterfaceCompPtr >     InterfaceCompList;

        typedef boost::shared_ptr< GeometryVector >     Composite;
        typedef std::map<int, Composite>                CompositeMap;
        typedef std::map<int, Composite>::iterator      CompositeMapIter;
        typedef std::map<int, Composite>::const_iterator      CompositeMapConstIter;

        struct ElementEdge
        {
            GeometrySharedPtr m_Element;
            int m_EdgeIndx;
        };

        struct ElementFace
        {
            GeometrySharedPtr m_Element;
            int m_FaceIndx;
        };


        typedef boost::shared_ptr<ElementEdge> ElementEdgeSharedPtr;
        typedef std::vector<ElementEdgeSharedPtr> ElementEdgeVector;
        typedef boost::shared_ptr<ElementEdgeVector> ElementEdgeVectorSharedPtr;

        typedef boost::shared_ptr<ElementFace> ElementFaceSharedPtr;
        typedef std::vector<ElementFaceSharedPtr> ElementFaceVector;
        typedef boost::shared_ptr<ElementFaceVector> ElementFaceVectorSharedPtr;

        // set restriction on domain range for post-processing.
        struct DomainRange
        {
            bool                    m_doXrange;
            NekDouble               m_xmin;
            NekDouble               m_xmax;
            bool                    m_doYrange;
            NekDouble               m_ymin;
            NekDouble               m_ymax;
            bool                    m_doZrange;
            NekDouble               m_zmin;
            NekDouble               m_zmax;

            bool                    m_checkShape;
            LibUtilities::ShapeType m_shapeType;
        };

        typedef boost::shared_ptr<DomainRange> DomainRangeShPtr;
        static DomainRangeShPtr NullDomainRangeShPtr;

        struct Expansion
        {
            Expansion(GeometrySharedPtr geomShPtr,
                      const LibUtilities::BasisKeyVector basiskeyvec):
                m_geomShPtr(geomShPtr),
                m_basisKeyVector(basiskeyvec)
            {
            }

            GeometrySharedPtr             m_geomShPtr;
            LibUtilities::BasisKeyVector  m_basisKeyVector;
        };

        typedef boost::shared_ptr<Expansion> ExpansionShPtr;
        typedef std::map<int, ExpansionShPtr> ExpansionMap;
        typedef std::map<int, ExpansionShPtr>::iterator ExpansionMapIter;
        typedef std::map<int, ExpansionShPtr>::const_iterator ExpansionMapConstIter;

        typedef boost::shared_ptr<ExpansionMap> ExpansionMapShPtr;
        typedef std::map<std::string, ExpansionMapShPtr> ExpansionMapShPtrMap;
        typedef std::map<std::string, ExpansionMapShPtr>::iterator  ExpansionMapShPtrMapIter;


        typedef std::map<std::string, std::string> GeomInfoMap;

        /// Base class for a spectral/hp element mesh.
        class MeshGraph
        {
            public:
                SPATIAL_DOMAINS_EXPORT MeshGraph();

                SPATIAL_DOMAINS_EXPORT MeshGraph(
                        unsigned int meshDimension,
                        unsigned int spaceDimension);

                SPATIAL_DOMAINS_EXPORT MeshGraph(
                        const LibUtilities::SessionReaderSharedPtr &pSession,
                        const DomainRangeShPtr &rng = NullDomainRangeShPtr);


                SPATIAL_DOMAINS_EXPORT virtual ~MeshGraph();


                /* ---- Mesh Reading routines ---- */
                SPATIAL_DOMAINS_EXPORT static boost::shared_ptr<MeshGraph> Read(
                        const LibUtilities::SessionReaderSharedPtr &pSession,
                        DomainRangeShPtr &rng = NullDomainRangeShPtr);

                /// \todo Remove updated routine
                SPATIAL_DOMAINS_EXPORT static boost::shared_ptr<MeshGraph> Read(
                        const std::string& infilename,
                        bool pReadExpansions = true);


                /// Read will read the meshgraph vertices given a filename.
                SPATIAL_DOMAINS_EXPORT virtual void ReadGeometry(
                        const std::string& infilename);

                /// Read will read the meshgraph vertices given a TiXmlDocument.
                SPATIAL_DOMAINS_EXPORT virtual void ReadGeometry(
                        TiXmlDocument &doc);

                /// Read geometric information from a file.
                SPATIAL_DOMAINS_EXPORT void ReadGeometryInfo(
                        const std::string &infilename);

                /// Read geometric information from an XML document.
                SPATIAL_DOMAINS_EXPORT void ReadGeometryInfo(
                        TiXmlDocument &doc);

                /// Read the expansions given the XML file path.
                SPATIAL_DOMAINS_EXPORT void ReadExpansions(
                        const std::string &infilename);

                /// Read the expansions given the XML document reference.
                SPATIAL_DOMAINS_EXPORT void ReadExpansions(
                        TiXmlDocument &doc);

                SPATIAL_DOMAINS_EXPORT void ReadDomain(
                        TiXmlDocument &doc);

                SPATIAL_DOMAINS_EXPORT void ReadCurves(
                        TiXmlDocument &doc);

                SPATIAL_DOMAINS_EXPORT void ReadCurves(
                        std::string &infilename);

                /* ---- Helper functions ---- */
                /// Dimension of the mesh (can be a 1D curve in 3D space).
                inline int GetMeshDimension() const;

                /// Dimension of the space (can be a 1D curve in 3D space).
                inline int GetSpaceDimension() const;


                /* Range definitions for postprorcessing */
                SPATIAL_DOMAINS_EXPORT void SetDomainRange
                    (NekDouble xmin, NekDouble xmax,
                     NekDouble ymin = NekConstants::kNekUnsetDouble,
                     NekDouble ymax = NekConstants::kNekUnsetDouble,
                     NekDouble zmin = NekConstants::kNekUnsetDouble,
                     NekDouble zmax = NekConstants::kNekUnsetDouble);

                /// Check if goemetry is in range definition if activated
                bool CheckRange(Geometry2D &geom);

                /// Check if goemetry is in range definition if activated
                bool CheckRange(Geometry3D &geom);

                /* ---- Composites and Domain ---- */
                inline Composite GetComposite(int whichComposite) const;

                SPATIAL_DOMAINS_EXPORT GeometrySharedPtr GetCompositeItem(
                        int whichComposite,
                        int whichItem);

                SPATIAL_DOMAINS_EXPORT void GetCompositeList(
                        const std::string &compositeStr,
                        CompositeMap &compositeVector) const;

                inline const CompositeMap &GetComposites() const;

                inline const std::vector<CompositeMap> &GetDomain(void) const;

                inline const CompositeMap &GetDomain(int domain) const;


                /* ---- Expansions ---- */
                inline const ExpansionMap &GetExpansions();

                SPATIAL_DOMAINS_EXPORT const ExpansionMap &GetExpansions(
                        const std::string variable);

                SPATIAL_DOMAINS_EXPORT ExpansionShPtr GetExpansion(
                                                                   GeometrySharedPtr geom, const std::string variable = "DefaultVar");

                /// Sets expansions given field definitions
                SPATIAL_DOMAINS_EXPORT void SetExpansions(
                        std::vector<LibUtilities::FieldDefinitionsSharedPtr>
                                                                &fielddef);

                /// Sets expansions given field definition, quadrature points.
                SPATIAL_DOMAINS_EXPORT void SetExpansions(
                        std::vector<LibUtilities::FieldDefinitionsSharedPtr>
                                                                &fielddef,
                        std::vector< std::vector<LibUtilities::PointsType> >
                                                                &pointstype );

                /// Sets expansions to have equispaced points
                SPATIAL_DOMAINS_EXPORT void SetExpansionsToEvenlySpacedPoints(
                                                        int npoints = 0);

                /// This function sets the expansion #exp in map with entry #variable
                inline void SetExpansions(
                        const std::string variable,
                        ExpansionMapShPtr &exp);

                /// Sets the basis key for all expansions of the given shape.
                SPATIAL_DOMAINS_EXPORT void SetBasisKey(
                        LibUtilities::ShapeType shape,
                        LibUtilities::BasisKeyVector &keys,
                        std::string var = "DefaultVar");

                inline bool  SameExpansions(
                        const std::string var1,
                        const std::string var2);

                inline bool CheckForGeomInfo(std::string parameter);

                inline const std::string GetGeomInfo(std::string parameter);

                SPATIAL_DOMAINS_EXPORT static LibUtilities::BasisKeyVector
                                        DefineBasisKeyFromExpansionType(
                        GeometrySharedPtr in,
                        ExpansionType type,
                        const int order);

                SPATIAL_DOMAINS_EXPORT LibUtilities::BasisKeyVector
                                        DefineBasisKeyFromExpansionTypeHomo(
                        GeometrySharedPtr in,
                        ExpansionType type_x,
                        ExpansionType type_y,
                        ExpansionType type_z,
                        const int nummodes_x,
                        const int nummodes_y,
                        const int nummodes_z);


                /* ---- Manipulation of mesh ---- */
                inline int GetNvertices() const;

                inline PointGeomSharedPtr GetVertex(int id);
                /// Adds a vertex to the with the next available ID.
                SPATIAL_DOMAINS_EXPORT PointGeomSharedPtr AddVertex(
                        NekDouble x,
                        NekDouble y,
                        NekDouble z);

                /// \brief Adds an edge between two points.  If curveDefinition is
                /// null, then the edge is straight, otherwise it is curved according
                /// to the curveDefinition.
                SPATIAL_DOMAINS_EXPORT SegGeomSharedPtr AddEdge(PointGeomSharedPtr v0, PointGeomSharedPtr v1,
                    CurveSharedPtr curveDefinition = CurveSharedPtr());
                SPATIAL_DOMAINS_EXPORT SegGeomSharedPtr GetEdge(unsigned int id) { return m_segGeoms[id]; }

                SPATIAL_DOMAINS_EXPORT TriGeomSharedPtr  AddTriangle(SegGeomSharedPtr edges[], StdRegions::Orientation orient[]);
                SPATIAL_DOMAINS_EXPORT QuadGeomSharedPtr AddQuadrilateral(SegGeomSharedPtr edges[], StdRegions::Orientation orient[]);
                SPATIAL_DOMAINS_EXPORT TetGeomSharedPtr AddTetrahedron(TriGeomSharedPtr tfaces[TetGeom::kNtfaces]);
                SPATIAL_DOMAINS_EXPORT PyrGeomSharedPtr AddPyramid(TriGeomSharedPtr tfaces[PyrGeom::kNtfaces],
                    QuadGeomSharedPtr qfaces[PyrGeom::kNqfaces]);
                SPATIAL_DOMAINS_EXPORT PrismGeomSharedPtr AddPrism(TriGeomSharedPtr tfaces[PrismGeom::kNtfaces],
                    QuadGeomSharedPtr qfaces[PrismGeom::kNqfaces]);
                SPATIAL_DOMAINS_EXPORT HexGeomSharedPtr AddHexahedron(QuadGeomSharedPtr qfaces[HexGeom::kNqfaces]);

                SPATIAL_DOMAINS_EXPORT const CurveVector& GetCurvedEdges() const { return m_curvedEdges; }

                SPATIAL_DOMAINS_EXPORT const CurveVector& GetCurvedFaces() const { return m_curvedFaces; }
                // void AddExpansion(ExpansionShPtr expansion) { m_expansions[expansion->m_geomShPtr->GetGlobalID()] = expansion; }
                SPATIAL_DOMAINS_EXPORT const SegGeomMap& GetAllSegGeoms() const { return m_segGeoms; }
                SPATIAL_DOMAINS_EXPORT const TriGeomMap& GetAllTriGeoms() const { return m_triGeoms; }
                SPATIAL_DOMAINS_EXPORT const QuadGeomMap& GetAllQuadGeoms() const { return m_quadGeoms; }
                SPATIAL_DOMAINS_EXPORT const TetGeomMap& GetAllTetGeoms() const { return m_tetGeoms; }
                SPATIAL_DOMAINS_EXPORT const PyrGeomMap& GetAllPyrGeoms() const { return m_pyrGeoms; }
                SPATIAL_DOMAINS_EXPORT const PrismGeomMap& GetAllPrismGeoms() const { return m_prismGeoms; }
                SPATIAL_DOMAINS_EXPORT const HexGeomMap& GetAllHexGeoms() const { return m_hexGeoms; }

                /// Convenience method for ElVis.
                template<typename ElementType>
                const std::map<int, boost::shared_ptr<ElementType> >& GetAllElementsOfType() const;

            protected:
                LibUtilities::SessionReaderSharedPtr    m_session;
                PointGeomMap                            m_vertSet;
                InterfaceCompList                       m_iComps;

                CurveVector                             m_curvedEdges;
                CurveVector                             m_curvedFaces;

                SegGeomMap                              m_segGeoms;

                TriGeomMap                              m_triGeoms;
                QuadGeomMap                             m_quadGeoms;
                TetGeomMap                              m_tetGeoms;
                PyrGeomMap                              m_pyrGeoms;
                PrismGeomMap                            m_prismGeoms;
                HexGeomMap                              m_hexGeoms;

                int                                     m_meshDimension;
                int                                     m_spaceDimension;
                int                                     m_partition;
                bool                                    m_meshPartitioned;

                CompositeMap                            m_meshComposites;
                std::vector<CompositeMap>               m_domain;
                DomainRangeShPtr                        m_domainRange;

                ExpansionMapShPtrMap                    m_expansionMapShPtrMap;

                GeomInfoMap                             m_geomInfo;


                ExpansionMapShPtr    SetUpExpansionMap(void);
        };
        typedef boost::shared_ptr<MeshGraph> MeshGraphSharedPtr;


        /**
         *
         */
        inline int MeshGraph::GetMeshDimension(void) const
        {
            return m_meshDimension;
        }


        /**
         *
         */
        inline int MeshGraph::GetSpaceDimension(void) const
        {
            return m_spaceDimension;
        }


        /**
         *
         */
        inline Composite MeshGraph::GetComposite(int whichComposite) const
        {
            Composite returnval;
            ASSERTL0(m_meshComposites.find(whichComposite) != m_meshComposites.end(),
                    "Composite not found.");
            return m_meshComposites.find(whichComposite)->second;
        }


        /**
         *
         */
        inline const CompositeMap &MeshGraph::GetComposites() const
        {
            return m_meshComposites;
        }


        /**
         *
         */
        inline const std::vector<CompositeMap> &MeshGraph::GetDomain(void) const
        {
            return m_domain;
        }

        /**
         *
         */
        inline const CompositeMap &MeshGraph::GetDomain(const int domain) const
        {
            ASSERTL1(domain < m_domain.size(),"Request for domain which does not exist");
            return m_domain[domain];
        }


        /**
         *
         */
        inline const ExpansionMap &MeshGraph::GetExpansions()
        {
            std::string defstr = "DefaultVar";
            return GetExpansions(defstr);
        }


        /**
         *
         */
        void  MeshGraph::SetExpansions(const std::string variable, ExpansionMapShPtr &exp)
        {
            if(m_expansionMapShPtrMap.count(variable) != 0)
            {
                ASSERTL0(false,(std::string("Expansion field is already set for variable ") + variable).c_str());
            }
            else
            {
                m_expansionMapShPtrMap[variable] = exp;
            }
        }


        /**
         *
         */
        inline bool MeshGraph::SameExpansions(const std::string var1, const std::string var2)
        {
            ExpansionMapShPtr expVec1 = m_expansionMapShPtrMap.find(var1)->second;
            ExpansionMapShPtr expVec2 = m_expansionMapShPtrMap.find(var2)->second;

            if(expVec1.get() == expVec2.get())
            {
                return true;
            }

            return false;
        }


        /**
         *
         */
        inline bool MeshGraph::CheckForGeomInfo(std::string parameter)
        {
            return m_geomInfo.find(parameter) != m_geomInfo.end();
        }


        /**
         *
         */
        inline const std::string MeshGraph::GetGeomInfo(std::string parameter)
        {
            ASSERTL1(m_geomInfo.find(parameter) != m_geomInfo.end(),
                    "Parameter " + parameter + " does not exist.");
            return m_geomInfo[parameter];
        }


        /**
         *
         */
        inline int MeshGraph::GetNvertices() const
        {
            return int(m_vertSet.size());
        }


        /**
         *
         */
        inline PointGeomSharedPtr MeshGraph::GetVertex(int id)
        {
            PointGeomSharedPtr returnval;
            PointGeomMap::iterator x = m_vertSet.find(id);
            ASSERTL0(x != m_vertSet.end(),
                     "Vertex " + boost::lexical_cast<string>(id)
                     + " not found.");
            return x->second;
        }


        /**
         *
         */
        template<>
        inline const std::map<int, boost::shared_ptr<HexGeom> >& MeshGraph::GetAllElementsOfType() const
        {
            return GetAllHexGeoms();
        }


        /**
         *
         */
        template<>
        inline const std::map<int, boost::shared_ptr<PrismGeom> >& MeshGraph::GetAllElementsOfType() const
        {
            return GetAllPrismGeoms();
        }


        /**
         *
         */
        template<>
        inline const std::map<int, boost::shared_ptr<TetGeom> >& MeshGraph::GetAllElementsOfType() const
        {
            return GetAllTetGeoms();
        }


        /**
         *
         */
        template<>
        inline const std::map<int, boost::shared_ptr<PyrGeom> >& MeshGraph::GetAllElementsOfType() const
        {
            return GetAllPyrGeoms();
        }
    };
};

#endif

