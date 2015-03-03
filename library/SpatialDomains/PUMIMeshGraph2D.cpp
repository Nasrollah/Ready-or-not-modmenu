//////////////////////////////////////////////////////////////////////////////////
//  File:  MeshGraph2D.cpp
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

#include <SpatialDomains/PUMIMeshGraph2D.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/CrvTriGeom.h> // Kai 10-22-2014
#include <SpatialDomains/CrvSegGeom.h> // Kai 10-30-2014
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <tinyxml.h>

#include <PUMI_NektarppAdapt.h>
#include <pumi_mesh.h>
#include <CrvEdge2.h>
#include <CrvEdge3.h>
#include <CrvEdgeCAD.h>

using namespace pumi;

namespace Nektar
{
    namespace SpatialDomains
    {
        PUMIMeshGraph2D::PUMIMeshGraph2D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                 const DomainRangeShPtr &rng)
            : MeshGraph2D()
        {
            // load in pumi mesh
            PUMI_NektarppAdapt::Instance()->init(NULL, NULL);
            TiXmlElement* nektar_tag = pSession->GetElement("NEKTAR");
            TiXmlElement* pumi_model_tag = nektar_tag->FirstChildElement("GEOMODEL");
            if(pumi_model_tag)
            {
              const char *pumi_model_file = pumi_model_tag->Attribute("FILENAME");
              PUMI_NektarppAdapt::Instance()->load_model(pumi_model_file);
            }


            TiXmlElement* pumi_mesh_tag = nektar_tag->FirstChildElement("PUMIMESH");
            if(pumi_mesh_tag)
            {
              const char *pumi_mesh_file =  pumi_mesh_tag->Attribute("FILENAME");
              PUMI_NektarppAdapt::Instance()->load_mesh(pumi_mesh_file, NULL);
              PUMI_NektarppAdapt::Instance()->verify_mesh();
              // label entities
//              PUMI_NektarppAdapt::Instance()->read_or_create_label_for_mesh_entities();
              pMeshMdl pumi_mesh = PUMI_NektarppAdapt::Instance()->PumiMesh();

              pPart pumi_part;
              PUMI_Mesh_GetPart(pumi_mesh, 0, pumi_part);
              m_meshPartitioned = false;
              m_spaceDimension = 2;
              m_meshDimension = 2;
              // **** assemble Nektar mesh vertices
              pVertex pumi_vert;
              VIter viter = M_vertexIter(pumi_part);
              int indx = 0;
              while ( (pumi_vert = VIter_next(viter)) ) {
                double coords[3];
                V_coord(pumi_vert, coords);
                // get pumi vert id
                indx = EN_id(pumi_vert);
                // debug output
//                std::cout << "pumi vert id: " << indx << std::endl;
                PointGeomSharedPtr vert(MemoryManager<PointGeom>::AllocateSharedPtr(m_spaceDimension, indx, coords[0], coords[1], coords[2]));
                m_vertSet[indx] = vert;
                // indx++;
              }
              VIter_delete(viter);

              // **** assemble Nektar mesh edges
              pEdge pumi_edge;
              pVertex edge_verts[2];
              int edge_vert_ids[2];
              // reset index
              indx = 0;
              EIter eiter = M_edgeIter(pumi_part);              
              while ( (pumi_edge = EIter_next(eiter)) ) {
                EN_setID(pumi_edge, indx);
                // get edge end vertex entities
                edge_verts[0] = E_vertex(pumi_edge, 0);
                edge_verts[1] = E_vertex(pumi_edge, 1);
                // get vertex ides
                edge_vert_ids[0] = EN_id(edge_verts[0]);
                edge_vert_ids[1] = EN_id(edge_verts[1]);
                // get the Nektar vert entities by id
                PointGeomSharedPtr vertices[2] = {m_vertSet[edge_vert_ids[0]], m_vertSet[edge_vert_ids[1]]};
                pMeshEnt pumi_verts[2];
                double coords[2][3];
                for(int i = 0; i < 2; ++i) {
                  (vertices[i])->GetCoords(coords[i][0], coords[i][1], coords[i][2]);
                  pumi_verts[i] = PUMI_NektarppAdapt::Instance()->get_pumi_mesh_vtx(coords[i][0], coords[i][1], coords[i][2]);
                }
                std::vector<pVertex> pumi_vert_vec;
                for (int i = 0; i < 2; ++i)
                  pumi_vert_vec.push_back(pumi_verts[i]);
                CrvEdge * bezseg = new CrvEdge2(pumi_vert_vec);
//                CrvEdge * bezseg = new CrvEdgeCAD(pumi_vert_vec);

                // create Nektar edge object
                CrvSegGeomSharedPtr edge = MemoryManager<CrvSegGeom>::AllocateSharedPtr(indx, m_spaceDimension, vertices, bezseg);
//                SegGeomSharedPtr edge = MemoryManager<SegGeom>::AllocateSharedPtr(indx, m_spaceDimension, vertices);
                edge->SetGlobalID(indx);
                m_segGeoms[indx] = edge;
                indx++;
              }
              EIter_delete(eiter);

              // **** assemble Nektar mesh faces -- for 2d meshes, elements
              pFace pumi_face;
              pEdge face_edges[3];
              int face_edge_ids[3];
              // reset index
              indx = 0;
              FIter fiter = M_faceIter(pumi_part);
              while ( (pumi_face = FIter_next(fiter)) ) {
                EN_setID(pumi_face, indx);
                // get face bounding edge entities
                for (int i = 0; i < 3; ++i) {
                  face_edges[i] = F_edge(pumi_face, i);
                  face_edge_ids[i] = EN_id(face_edges[i]);
                }
                
                SegGeomSharedPtr edges[TriGeom::kNedges] =
                {
                    m_segGeoms[face_edge_ids[0]],
                    m_segGeoms[face_edge_ids[1]],
                    m_segGeoms[face_edge_ids[2]]
                };
                StdRegions::Orientation edgeorient[TriGeom::kNedges] =
                {
                    SegGeom::GetEdgeOrientation(*edges[0], *edges[1]),
                    SegGeom::GetEdgeOrientation(*edges[1], *edges[2]),
                    SegGeom::GetEdgeOrientation(*edges[2], *edges[0])
                };

                // now re-route tri geom creation to use CrvTriGeom
                CrvTriGeomSharedPtr crvtrigeom = MemoryManager<CrvTriGeom>::AllocateSharedPtr(indx, edges, edgeorient);
//                TriGeomSharedPtr crvtrigeom = MemoryManager<TriGeom>::AllocateSharedPtr(indx, edges, edgeorient);

                crvtrigeom->SetGlobalID(indx);

                m_triGeoms[indx] = crvtrigeom;

                indx++;
              }
              FIter_delete(fiter);

//              PUMI_NektarppAdapt::Instance()->AllocTriDerivMap(indx);
              // **** assemble Nektar Composites for BC and Domain refs
              Composite curVector(MemoryManager<GeometryVector>::AllocateSharedPtr());
              m_meshComposites[0] = curVector;
              // all tri elements => composite 0
              int num_faces = M_numFaces(pumi_part);
              for(int i = 0; i < num_faces; ++i)
                (m_meshComposites[0])->push_back(m_triGeoms[i]);
              // left is inflow
              Composite inflowVector(MemoryManager<GeometryVector>::AllocateSharedPtr());
              m_meshComposites[1] = inflowVector;
              // right is outflow
              Composite outflowVector(MemoryManager<GeometryVector>::AllocateSharedPtr());
              m_meshComposites[2] = outflowVector;
              // top and bottom is wall
              Composite wallVector(MemoryManager<GeometryVector>::AllocateSharedPtr());
              m_meshComposites[3] = wallVector;
              
              int wallid[2] = {14,22};
              int inflowid[3] ={6, 10, 25};
              int outflowid =18;
              
              EIter eiter2 = M_edgeIter(pumi_part);
              while (pEdge medge = EIter_next(eiter2)) {
                if(E_whatInType(medge) == Tedge) {
                  pGEdge gedge = (pGEdge)E_whatIn(medge);
                  int gedgeid = GEN_tag(gedge);
                  if(gedgeid == wallid[0] || gedgeid == wallid[1])
                    (m_meshComposites[3])->push_back(m_segGeoms[EN_id(medge)]);
                  if(gedgeid == inflowid[0] || gedgeid == inflowid[1] || gedgeid == inflowid[2])
                    (m_meshComposites[1])->push_back(m_segGeoms[EN_id(medge)]);
                  if(gedgeid == outflowid)
                    (m_meshComposites[2])->push_back(m_segGeoms[EN_id(medge)]);
                }
              }
              EIter_delete(eiter2);

              // Read the domain composites.
              // Parse the composites into a list.
              // Assumes theres only 1
              // TODO: GENERALIZE THIS
              CompositeMap fullDomain;
              GetCompositeList(0, fullDomain);
              m_domain.push_back(fullDomain);

              // use the original ReadExpansions function unchanged
              ReadExpansions(pSession->GetDocument());
            }
        }

        PUMIMeshGraph2D::~PUMIMeshGraph2D()
        {
        }
    }; //end of namespace
}; //end of namespace
