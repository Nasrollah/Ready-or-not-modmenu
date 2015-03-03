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

#include <SpatialDomains/PUMIMeshGraph3D.h>
#include <SpatialDomains/CrvTriGeom.h> // Kai 10-22-2014
#include <SpatialDomains/CrvSegGeom.h> // Kai 10-30-2014
#include <SpatialDomains/CrvTetGeom.h> // Kai 11-11-2014
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <tinyxml.h>

#include <PUMI_NektarppAdapt.h>
#include <pumi_mesh.h>
#include <CrvEdge2.h>
#include <CrvEdge3.h>
#include <CrvEdge4.h>
#include <CrvEdge5.h>
#include <CrvEdgeCAD.h>

using namespace pumi;

namespace Nektar
{
    namespace SpatialDomains
    {
        PUMIMeshGraph3D::PUMIMeshGraph3D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                 const DomainRangeShPtr &rng)
            : MeshGraph3D()
        {

            TiXmlElement* mesh = NULL;
            TiXmlElement* geom = NULL;
            TiXmlElement* master = NULL;    // Master tag within which all data is contained.

  
            // load in pumi mesh
            PUMI_NektarppAdapt::Instance()->init(NULL, NULL);
            master = pSession->GetElement("NEKTAR");

            geom = master->FirstChildElement("GEOMODEL");
            ASSERTL0(geom, "Unable to find GEOMODEL tag in file.");
            
            // get important faces

            const char *pumi_model_file = geom->Attribute("FILENAME");
            PUMI_NektarppAdapt::Instance()->load_model(pumi_model_file);

            mesh = master->FirstChildElement("PUMI");
            ASSERTL0(mesh, "Unable to find PUMI tag in file.");
            TiXmlAttribute *attr = mesh->FirstAttribute();

            while (attr)
            {
                std::string attrName(attr->Name());
                if (attrName == "ORDER")
                {
                    int err = attr->QueryIntValue(&m_CurveOrder);
                    ASSERTL1(err==TIXML_SUCCESS, "Unable to read PUMI curve order.");
                    break;
                }
                else if (attrName == "FILENAME")
                {
                    const char *pumi_mesh_file = mesh->Attribute("FILENAME");
                    PUMI_NektarppAdapt::Instance()->load_mesh(pumi_mesh_file, NULL);
                    PUMI_NektarppAdapt::Instance()->verify_mesh();
                }
                else
                {
                    std::string errstr("Unknown attribute: ");
                    errstr += attrName;
                    ASSERTL0(false, errstr.c_str());
                }

                // Get the next attribute.
                attr = attr->Next();
            }

            pMeshMdl pumi_mesh = PUMI_NektarppAdapt::Instance()->PumiMesh();

            pPart pumi_part;
            PUMI_Mesh_GetPart(pumi_mesh, 0, pumi_part);
            m_meshPartitioned = false;
            m_spaceDimension = 3;
            m_meshDimension = 3;

            // **** assemble Nektar mesh vertices
            pVertex pumi_vert;
            VIter viter = M_vertexIter(pumi_part);
            int indx = 0;
            while ( (pumi_vert = VIter_next(viter)) ){
              double coords[3];
              V_coord(pumi_vert, coords);
              // get pumi vert id
              indx = EN_id(pumi_vert);
              // debug output
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
            while ( (pumi_edge = EIter_next(eiter)) ){
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

              // create pumi curved edge object
              CrvEdge * bezseg = NULL;
              switch (m_CurveOrder){
                case 2:
                {
                  bezseg = new CrvEdge2(pumi_vert_vec);
                  break;
                }
                case 3:
                {
                  bezseg = new CrvEdge3(pumi_vert_vec);
                  break;
                }
                case 4:
                {
                  bezseg = new CrvEdge4(pumi_vert_vec);
                  break;
                }
                case 5:
                {
                  bezseg = new CrvEdge5(pumi_vert_vec);
                  break;
                }
                default:
                  ASSERTL1(false, "Curve Order not defined.");
                  break;
              }

              std::cout << "interp error of CrvEdge: " << bezseg->interp_error() << std::endl;

              // modified version of nektar edge with bezier geom
              CrvSegGeomSharedPtr edge = MemoryManager<CrvSegGeom>::AllocateSharedPtr(indx, m_spaceDimension, vertices, bezseg);

              edge->SetGlobalID(indx);
              m_segGeoms[indx] = edge;
              indx++;
            }
            EIter_delete(eiter);

            // **** go through tets for the first time **** //
            std::vector< std::vector<pVertex> > tet_verts;
            RIter riter2 = M_regionIter(pumi_part);
            while(pRegion rgn = RIter_next(riter2)) {
              std::vector<pVertex> ordered_verts;
              std::vector<int> ids(4);
              map<int, pMeshEnt> id_vert_map;
              std::vector<pMeshEnt> tmp_adj_verts;
              PUMI_MeshEnt_GetAdj(rgn, PUMI_VERTEX, 0, tmp_adj_verts);
              assert( 4 == tmp_adj_verts.size() );
              // init ids vector and id_vert map
              for(int ii = 0; ii < 4; ii++)
              {
               
                pVertex vert_ent = tmp_adj_verts[ii];
                int vid = EN_id(vert_ent);
                id_vert_map[vid] = vert_ent;
                ids[ii] = vid;
              }

              // sort ids vector in ascending order
              sort(ids.begin(), ids.end());

              for(int jj = 0; jj < 4; ++jj) {
                // put ordered vertex in vector
                ordered_verts.push_back(id_vert_map[ids[jj]]);
              }
              assert(ordered_verts.size() == 4);

              // make sure ordered_verts has ids sorted in ascending order
              int id0, id1;
              for(int i = 0; i < 3; ++i) {
                id0 = EN_id(ordered_verts[i]);
                id1 = EN_id(ordered_verts[i+1]);
                assert(id0 < id1);
              }

              // get vertex coords
              double vert_coords[4][3];
              for(int i = 0; i < 4; ++i)
                V_coord(ordered_verts[i], vert_coords[i]);

              // check volume positivity
              double volume = XYZ_volume(vert_coords);
              // swap first two verts if volume is found to be negative
              if(volume < 0.0) {
                swap(ordered_verts[0], ordered_verts[1]);
              }

              tet_verts.push_back(ordered_verts);

            }
            RIter_delete(riter2);

            // **** go through tet_verts to create faces **** //
            // reset index
            indx = 0;
            std::vector<pFace> processed_faces;
            int face_vert_map[4][3] = {
              {0, 1, 2},
              {0, 1, 3},
              {1, 2, 3},
              {0, 2, 3}
            };
            for(int i = 0 ; i < tet_verts.size() ; ++i) {
              std::vector<pVertex> ordered_verts = tet_verts[i];

              for(int iface = 0; iface < 4; ++iface) {
                  int vid0 = face_vert_map[iface][0];
                  int vid1 = face_vert_map[iface][1];
                  int vid2 = face_vert_map[iface][2];
                  pFace face = F_exist(Tvertex,
                                       ordered_verts[vid0],
                                       ordered_verts[vid1],
                                       ordered_verts[vid2],
                                       0);
                  assert(face);
                  if(std::find(processed_faces.begin(), processed_faces.end(), face) == processed_faces.end()) {
                      processed_faces.push_back(face);
                      EN_setID(face, indx);
                      // get edge 0, 1
                      pEdge edge0 = E_exist(ordered_verts[vid0], ordered_verts[vid1]);
                      assert(edge0);
                      int eid0 = EN_id(edge0);

                      // get edge 1, 2
                      pEdge edge1 = E_exist(ordered_verts[vid1], ordered_verts[vid2]);
                      assert(edge1);
                      int eid1 = EN_id(edge1);

                      // get edge 2, 0
                      pEdge edge2 = E_exist(ordered_verts[vid2], ordered_verts[vid0]);
                      assert(edge2);
                      int eid2 = EN_id(edge2);

                      SegGeomSharedPtr edges[TriGeom::kNedges] =
                      {
                          m_segGeoms[eid0],
                          m_segGeoms[eid1],
                          m_segGeoms[eid2]
                      };
                      StdRegions::Orientation edgeorient[TriGeom::kNedges] =
                      {
                          SegGeom::GetEdgeOrientation(*edges[0], *edges[1]),
                          SegGeom::GetEdgeOrientation(*edges[1], *edges[2]),
                          SegGeom::GetEdgeOrientation(*edges[2], *edges[0])
                      };

                      // now re-route tri geom creation to use CrvTriGeom
//                        std::cout << "pumi face id: " << indx << std::endl;

                      // original version of nektar tri face
//                      TriGeomSharedPtr crvtrigeom = MemoryManager<TriGeom>::AllocateSharedPtr(indx, edges, edgeorient);

                      // modified version of nektar crv tri face
                      CrvTriGeomSharedPtr crvtrigeom = MemoryManager<CrvTriGeom>::AllocateSharedPtr(indx, edges, edgeorient);

                      crvtrigeom->SetGlobalID(indx);
                      m_triGeoms[indx] = crvtrigeom;
                      indx++;

                  }// end if
              }// end for       

            } // end for

            /**** go through tet_verts again to assemble tets ****/
            // reset indx
            indx = 0;
            for(int j = 0 ; j < tet_verts.size() ; ++j) {
              std::vector<pVertex> ordered_verts = tet_verts[j];
              TriGeomSharedPtr tfaces[4];
              // face node; 0, 1, 2
              pFace face0 = F_exist(Tvertex, ordered_verts[0], ordered_verts[1], ordered_verts[2], 0);
              assert(face0);
              int fid = EN_id(face0);
              tfaces[0] = m_triGeoms[fid];
              // face node; 0, 1, 3
              face0 = F_exist(Tvertex, ordered_verts[0], ordered_verts[1], ordered_verts[3], 0);
              fid = EN_id(face0);
              tfaces[1] = m_triGeoms[fid];
              // face node; 1, 2, 3
              face0 = F_exist(Tvertex, ordered_verts[1], ordered_verts[2], ordered_verts[3], 0);
              assert(face0);
              fid = EN_id(face0);
              tfaces[2] = m_triGeoms[fid];
              // face node; 0, 2, 3
              face0 = F_exist(Tvertex, ordered_verts[0], ordered_verts[2], ordered_verts[3], 0);
              assert(face0);
              fid = EN_id(face0);
              tfaces[3] = m_triGeoms[fid];

              // create tet
              CrvTetGeomSharedPtr tetgeom(MemoryManager<CrvTetGeom>::AllocateSharedPtr(tfaces));
              tetgeom->SetGlobalID(indx);

              m_tetGeoms[indx] = tetgeom;
              PopulateFaceToElMap(tetgeom, 4);
              indx++;
            }

            // TODO FIX ME
            // **** assemble Nektar Composites for BC and Domain refs
            Composite curVector(MemoryManager<GeometryVector>::AllocateSharedPtr());
            m_meshComposites[0] = curVector;
            // all tet elements => composite 0
            int num_rgns = M_numRegions(pumi_part);
            for(int i = 0; i < num_rgns; ++i)
              (m_meshComposites[0])->push_back(m_tetGeoms[i]);

            // left is inflow
            Composite inflowVector(MemoryManager<GeometryVector>::AllocateSharedPtr());
            m_meshComposites[1] = inflowVector;
            // right is outflow
            Composite outflowVector(MemoryManager<GeometryVector>::AllocateSharedPtr());
            m_meshComposites[2] = outflowVector;
            // front side wall
            Composite frontVector(MemoryManager<GeometryVector>::AllocateSharedPtr());
            m_meshComposites[3] = frontVector;
            // back side wall
            Composite backVector(MemoryManager<GeometryVector>::AllocateSharedPtr());
            m_meshComposites[4] = backVector;
            // top wall
            Composite topVector(MemoryManager<GeometryVector>::AllocateSharedPtr());
            m_meshComposites[5] = topVector;
            // bottom wall
            Composite bottomVector(MemoryManager<GeometryVector>::AllocateSharedPtr());
            m_meshComposites[6] = bottomVector;
/*
            int inflowid = 82;
            int outflowid =78;
            int topid = 24;
            int bottomid = 42;
            int frontid = 80;
            int backid = 76;
*/
            int inflowids[3] = {114, 112, 104};
            int outflowid = 108;
            int topid = 34;
            int bottomid = 56;
            int frontid = 110;
            int backid = 106;
            FIter fiter2 = M_faceIter(pumi_part);
            while (pFace mface = FIter_next(fiter2)) {
              if(F_whatInType(mface) == Tface) {
                pGFace gface = (pGFace)F_whatIn(mface);
                int gfaceid = GEN_tag(gface);
                if(gfaceid == inflowids[0] ||
                   gfaceid == inflowids[1] ||
                   gfaceid == inflowids[2]) {
                  (m_meshComposites[1])->push_back(m_triGeoms[EN_id(mface)]);
                }
                else if(gfaceid == outflowid) {
                  (m_meshComposites[2])->push_back(m_triGeoms[EN_id(mface)]);
                }
                else if(gfaceid == frontid) {
                  (m_meshComposites[3])->push_back(m_triGeoms[EN_id(mface)]);
                }
                else if(gfaceid == backid) {
                  (m_meshComposites[4])->push_back(m_triGeoms[EN_id(mface)]);
                }
                else if(gfaceid == topid) {
                  (m_meshComposites[5])->push_back(m_triGeoms[EN_id(mface)]);
                }
                else if(gfaceid == bottomid) {
                  (m_meshComposites[6])->push_back(m_triGeoms[EN_id(mface)]);
                }
                else {
                  std::cout << "Undefined boundary condition face id" << std::endl;
                  exit(1);
                }
              }
            }
            FIter_delete(fiter2);
            // Read the domain composites.
            // Parse the composites into a list.
            // Assumes theres only 1
            // TODO: GENERALIZE THIS
            CompositeMap one;
            one[0] = m_meshComposites[0];
            m_domain.push_back(one);

            // use the original ReadExpansions function unchanged
            ReadExpansions(pSession->GetDocument());

            
        }

        PUMIMeshGraph3D::~PUMIMeshGraph3D()
        {
        }
    }; //end of namespace
}; //end of namespace
