#include <cstdlib>
#include <math.h>

#include <StdRegions/StdExpansion3D.h>
#include <LocalRegions/HexExp.h>
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/PyrExp.h>
#include <LocalRegions/TetExp.h>

#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include <SpatialDomains/MeshComponents.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::LocalRegions;
using namespace Nektar::StdRegions;
using namespace Nektar::SpatialDomains;

/// Defines a solution which excites all modes in the Tet expansion.
NekDouble Tet_sol( NekDouble x, NekDouble y, NekDouble z,
                   int order1, int order2, int order3);

/// Defines a solution which excites all modes in the Prism expansion.
NekDouble Prism_sol( NekDouble x, NekDouble y, NekDouble z,
                     int order1, int order2, int order3);

/// Defines a solution which excites all modes in the Hex expansion.
NekDouble Hex_sol( NekDouble x, NekDouble y, NekDouble z,
                   int order1, int order2, int order3,
                   LibUtilities::BasisType btype1,
                   LibUtilities::BasisType btype2,
                   LibUtilities::BasisType btype3);

/// Generates a Hex geometry based on the command-line parameters.
SpatialDomains::HexGeomSharedPtr CreateHexGeom(int argc, char *argv[]);

/// Generates a Hex geometry based on the command-line parameters.
SpatialDomains::PrismGeomSharedPtr CreatePrismGeom(int argc, char *argv[]);

/// Generates a Pyramid geometry based on the command-line parameters.
SpatialDomains::PyrGeomSharedPtr CreatePyrGeom(int argc, char *argv[]);

/// Generates a Tet geometry based on the command-line parameters.
SpatialDomains::TetGeomSharedPtr CreateTetGeom(int argc, char *argv[]);

/// modification to deal with exact solution. Return 1 if integer < 0
static double  pow_loc(const double val, const int i)
{
    return (i < 0)? 1.0: pow(val,i);
}

/// This routine projects a polynomial or trigonmetric functions which
/// has energy in all modes of the expansions and reports and error
int main(int argc, char *argv[]){
    int           i;

    int           order1,order2,order3, nq1,nq2,nq3;
    PointsType    Qtype1,Qtype2,Qtype3;
    BasisType     btype1,btype2,btype3;

    ShapeType regionshape;
    StdExpansion* E = NULL;
    Array<OneD, NekDouble> sol;

    if(argc < 23)
    {
        fprintf(stderr,"Usage: StdProject2D RegionShape Type1 Type2 Type3 "
                       "order1 order2 order3 nq1 nq2 nq3 x1 y1 z1 x2 y2 z2 "
                       "x3 y3 z3 [x4 y4 z4...]\n");
        fprintf(stderr,"Where RegionShape is an integer value which "
                       "dictates the region shape:\n");
        fprintf(stderr,"\t Tetrahedron   = 5\n");
        fprintf(stderr,"\t Pyramid       = 6\n");
        fprintf(stderr,"\t Prism         = 7\n");
        fprintf(stderr,"\t Hexahedron    = 8\n");


        fprintf(stderr,"Where type is an integer value which "
                       "dictates the basis as:\n");

        fprintf(stderr,"\t Ortho_A    = 1\n");
        fprintf(stderr,"\t Ortho_B    = 2\n");
        fprintf(stderr,"\t Ortho_C    = 3\n");
        fprintf(stderr,"\t Modified_A = 4\n");
        fprintf(stderr,"\t Modified_B = 5\n");
        fprintf(stderr,"\t Modified_C = 6\n");
        fprintf(stderr,"\t Fourier    = 7\n");
        fprintf(stderr,"\t Lagrange   = 8\n");
        fprintf(stderr,"\t Legendre   = 9\n");
        fprintf(stderr,"\t Chebyshev  = 10\n");

        exit(1);
    }

    regionshape = (LibUtilities::ShapeType) atoi(argv[1]);

    // Check to see if 3D region
    if (regionshape != LibUtilities::eTetrahedron &&
        regionshape != LibUtilities::ePyramid     &&
        regionshape != LibUtilities::ePrism       &&
        regionshape != LibUtilities::eHexahedron)
    {
        NEKERROR(ErrorUtil::efatal,"This shape is not a 3D region");
    }

    int btype1_val  = atoi(argv[2]);
    int btype2_val  = atoi(argv[3]);
    int btype3_val  = atoi(argv[4]);
    btype1          = static_cast<BasisType>(btype1_val);
    btype2          = static_cast<BasisType>(btype2_val);
    btype3          = static_cast<BasisType>(btype3_val);

    // Check to see that correct Expansions are used
    switch(regionshape)
    {
        case LibUtilities::eTetrahedron:
            if((btype1 == eOrtho_B) || (btype1 == eOrtho_C)
               || (btype1 == eModified_B) || (btype1 == eModified_C))
            {
                NEKERROR(ErrorUtil::efatal,
                         "Basis 1 cannot be of type Ortho_B, Ortho_C, Modified_B "
                         "or Modified_C");
            }
            if((btype2 == eOrtho_A) || (btype2 == eOrtho_C)
               || (btype2 == eModified_A) || (btype2 == eModified_C))
            {
                NEKERROR(ErrorUtil::efatal,
                         "Basis 2 cannot be of type Ortho_A, Ortho_C, Modified_A "
                         "or Modified_C");
            }
            if((btype3 == eOrtho_A) || (btype3 == eOrtho_B)
               || (btype3 == eModified_A) || (btype3 == eModified_B))
            {
                NEKERROR(ErrorUtil::efatal,
                         "Basis 3 cannot be of type Ortho_A, Ortho_B, Modified_A "
                         "or Modified_B");
            }
            break;
        case LibUtilities::ePyramid:
            if((btype1 == eOrtho_B) || (btype1 == eOrtho_C)
               || (btype1 == eModified_B) || (btype1 == eModified_C))
            {
                NEKERROR(ErrorUtil::efatal,
                         "Basis 1 cannot be of type Ortho_B, Ortho_C, Modified_B "
                         "or Modified_C");
            }
            if((btype2 == eOrtho_B) || (btype2 == eOrtho_C)
               || (btype2 == eModified_B) || (btype2 == eModified_C))
            {
                NEKERROR(ErrorUtil::efatal,
                         "Basis 2 cannot be of type Ortho_B, Ortho_C, Modified_B "
                         "or Modified_C");
            }
            if((btype3 == eOrtho_A) || (btype3 == eOrtho_B)
               || (btype3 == eModified_A) || (btype3 == eModified_B))
            {
                NEKERROR(ErrorUtil::efatal,
                         "Basis 3 cannot be of type Ortho_A, Ortho_B, Modified_A "
                         "or Modified_B");
            }
            break;
        case LibUtilities::ePrism:
            if((btype1 == eOrtho_B) || (btype1 == eOrtho_C)
               || (btype1 == eModified_B) || (btype1 == eModified_C))
            {
                NEKERROR(ErrorUtil::efatal,
                         "Basis 1 cannot be of type Ortho_B, Ortho_C, Modified_B "
                         "or Modified_C");
            }
            if((btype2 == eOrtho_B) || (btype2 == eOrtho_C)
               || (btype2 == eModified_B) || (btype2 == eModified_C))
            {
                NEKERROR(ErrorUtil::efatal,
                         "Basis 2 cannot be of type Ortho_B, Ortho_C, Modified_B "
                         "or Modified_C");
            }
            if((btype3 == eOrtho_A) || (btype3 == eOrtho_C)
               || (btype3 == eModified_A) || (btype3 == eModified_C))
            {
                NEKERROR(ErrorUtil::efatal,
                         "Basis 3 cannot be of type Ortho_A, Ortho_C, Modified_A "
                         "or Modified_C");
            }
            break;
        case LibUtilities::eHexahedron:
            if((btype1 == eOrtho_B) || (btype1 == eOrtho_C)
               || (btype1 == eModified_B) || (btype1 == eModified_C))
            {
                NEKERROR(ErrorUtil::efatal,
                         "Basis 1 is for 2 or 3D expansions");
            }

            if((btype2 == eOrtho_B) || (btype2 == eOrtho_C)
               || (btype2 == eModified_B) || (btype2 == eModified_C))
            {
                NEKERROR(ErrorUtil::efatal,
                         "Basis 2 is for 2 or 3D expansions");
            }

            if((btype3 == eOrtho_B) || (btype3 == eOrtho_C)
               || (btype3 == eModified_B) || (btype3 == eModified_C))
            {
                NEKERROR(ErrorUtil::efatal,
                         "Basis 3 is for 2 or 3D expansions");
            }
            break;
        default:
            ASSERTL0(false, "Not a 3D expansion.");
            break;
    }

    order1 =   atoi(argv[5]);
    order2 =   atoi(argv[6]);
    order3 =   atoi(argv[7]);
    nq1    =   atoi(argv[8]);
    nq2    =   atoi(argv[9]);
    nq3    =   atoi(argv[10]);

    sol = Array<OneD, NekDouble>(nq1*nq2*nq3);

    if(btype1 != LibUtilities::eFourier)
    {
        Qtype1 = LibUtilities::eGaussLobattoLegendre;
    }
    else
    {
        Qtype1 = LibUtilities::eFourierEvenlySpaced;
    }

    if(btype2 != LibUtilities::eFourier)
    {
        if (regionshape == LibUtilities::eTetrahedron) {
            Qtype2 = LibUtilities::eGaussRadauMAlpha1Beta0;
        }
        else
        {
            Qtype2 = LibUtilities::eGaussLobattoLegendre;
        }
    }
    else
    {
        Qtype2 = LibUtilities::eFourierEvenlySpaced;
    }

    if(btype3 != LibUtilities::eFourier)
    {
        if (regionshape == LibUtilities::eTetrahedron ||
            regionshape == LibUtilities::ePyramid)
        {
            Qtype3 = LibUtilities::eGaussRadauMAlpha2Beta0;
        }
        else if (regionshape == LibUtilities::ePrism)
        {
            Qtype3 = LibUtilities::eGaussRadauMAlpha1Beta0;
        }
        else
        {
            Qtype3 = LibUtilities::eGaussLobattoLegendre;
        }
    }
    else
    {
        Qtype3 = LibUtilities::eFourierEvenlySpaced;
    }

    //-----------------------------------------------
    // Define a 3D expansion based on basis definition
    const LibUtilities::PointsKey Pkey1(nq1,Qtype1);
    const LibUtilities::PointsKey Pkey2(nq2,Qtype2);
    const LibUtilities::PointsKey Pkey3(nq3,Qtype3);
    const LibUtilities::BasisKey  Bkey1(btype1,order1,Pkey1);
    const LibUtilities::BasisKey  Bkey2(btype2,order2,Pkey2);
    const LibUtilities::BasisKey  Bkey3(btype3,order3,Pkey3);
    Array<OneD,NekDouble> x = Array<OneD,NekDouble>(nq1*nq2*nq3);
    Array<OneD,NekDouble> y = Array<OneD,NekDouble>(nq1*nq2*nq3);
    Array<OneD,NekDouble> z = Array<OneD,NekDouble>(nq1*nq2*nq3);


    switch(regionshape)
    {
    case LibUtilities::eTetrahedron:
        {
            SpatialDomains::TetGeomSharedPtr geom = CreateTetGeom(argc, argv);
            E = new LocalRegions::TetExp(Bkey1, Bkey2, Bkey3, geom);
            E->GetCoords(x,y,z);

            //----------------------------------------------
            // Define solution to be projected
            for(i = 0; i < nq1*nq2*nq3; ++i)
            {
                sol[i]  = Tet_sol(x[i],y[i],z[i],order1,order2,order3);
            }
            //----------------------------------------------
        }
        break;
    case LibUtilities::ePyramid:
        {
            SpatialDomains::PyrGeomSharedPtr geom = CreatePyrGeom(argc, argv);
            E = new LocalRegions::PyrExp(Bkey1, Bkey2, Bkey3, geom);
            E->GetCoords(x,y,z);

            //----------------------------------------------
            // Define solution to be projected (same as tetrahedron space)
            for(i = 0; i < nq1*nq2*nq3; ++i)
            {
                sol[i]  = Tet_sol(x[i],y[i],z[i],order1,order2,order3);
            }
            //----------------------------------------------
        }
        break;
    case LibUtilities::ePrism:
        {
            SpatialDomains::PrismGeomSharedPtr geom = CreatePrismGeom(argc, argv);
            E = new LocalRegions::PrismExp(Bkey1, Bkey2, Bkey3, geom);
            E->GetCoords(x,y,z);
            
            //----------------------------------------------
            // Define solution to be projected
            for(i = 0; i < nq1*nq2*nq3; ++i)
            {
                sol[i]  = Prism_sol(x[i],y[i],z[i],order1,order2,order3);
            }
            //----------------------------------------------
        }
        break;
    case LibUtilities::eHexahedron:
        {
            SpatialDomains::HexGeomSharedPtr geom = CreateHexGeom(argc, argv);
            E = new LocalRegions::HexExp(Bkey1, Bkey2, Bkey3, geom);
            E->GetCoords(x,y,z);

            //----------------------------------------------
            // Define solution to be projected
            for(i = 0; i < nq1*nq2*nq3; ++i)
            {
                sol[i]  = Hex_sol(x[i],y[i],z[i],order1,order2,order3,
                                  btype1,btype2,btype3);
            }
            //---------------------------------------------
        }
        break;
        default:
            ASSERTL0(false, "Not a 3D expansion.");
            break;
    }

    //---------------------------------------------
    // Project onto Expansion
    Array<OneD, NekDouble> coeffs(E->GetNcoeffs());
    Array<OneD, NekDouble> phys  (nq1*nq2*nq3);
    E->FwdTrans(sol, coeffs);
    //---------------------------------------------

    //-------------------------------------------
    // Backward Transform Solution to get projected values
    E->BwdTrans(coeffs, phys);
    //-------------------------------------------

    //--------------------------------------------
    // Calculate L_inf error
    cout << "L infinity error: " << E->Linf(phys, sol) << endl;
    cout << "L 2 error:        " << E->L2  (phys, sol) << endl;
    //--------------------------------------------

    return 0;
}


NekDouble Tet_sol(NekDouble x, NekDouble y, NekDouble z,
                  int order1, int order2, int order3){
    int k, l, m;
    NekDouble sol = 0;

    for(k = 0; k < order1; ++k)
    {
        for(l = 0; l < order2-k; ++l)
        {
            for(m = 0; m < order3-k-l; ++m)
            {
                sol += pow(x,k)*pow(y,l)*pow(z,m);
            }
        }
    }

    return sol;
}


NekDouble Prism_sol(NekDouble x, NekDouble y, NekDouble z,
                    int order1, int order2, int order3)
{
    int k, l, m;
    NekDouble sol = 0;
    
    for(k = 0; k < order1; ++k)
    {
        for(l = 0; l < order2; ++l)
        {
            for(m = 0; m < order3-k; ++m)
            {
                sol += pow(x,k)*pow(y,l)*pow(z,m);
            }
        }
    }

    return sol;
}


NekDouble Hex_sol(NekDouble x, NekDouble y, NekDouble z,
                  int order1, int order2, int order3,
                  LibUtilities::BasisType btype1,
                  LibUtilities::BasisType btype2,
                  LibUtilities::BasisType btype3)
{
    int i,j,k;
    NekDouble sol = 0.0;

    int  Nx = (btype1 == LibUtilities::eFourier ? order1/2 : order1);
    int  Ny = (btype2 == LibUtilities::eFourier ? order2/2 : order2);
    int  Nz = (btype3 == LibUtilities::eFourier ? order3/2 : order3);
    bool Fx = (btype1 == LibUtilities::eFourier);
    bool Fy = (btype2 == LibUtilities::eFourier);
    bool Fz = (btype3 == LibUtilities::eFourier);
    NekDouble a;

    for (i = 0; i < Nx; ++i)
    {
        for (j = 0; j < Ny; ++j)
        {
            for (k = 0; k < Nz; ++k)
            {
                a  = (Fx ? sin(M_PI*i*x) + cos(M_PI*i*x) : pow_loc(x,i));
                a *= (Fy ? sin(M_PI*j*y) + cos(M_PI*j*y) : pow_loc(y,j));
                a *= (Fz ? sin(M_PI*k*z) + cos(M_PI*k*z) : pow_loc(z,k));
                sol += a;
            }
        }
    }
    return sol;
}


SpatialDomains::HexGeomSharedPtr CreateHexGeom(int argc, char *argv[])
{
    if (argc != 35)
    {
        cout << "Insufficient points for a hex!" << endl;
    }

    // /////////////////////////////////////////////////////////////////////
    // Set up Hexahedron vertex coordinates
    // PointGeom (const int coordim, const int vid, double x,
    //   double y, double z)
    const int nVerts = 8;
    const double point[][3] = {
        {atof(argv[11]),atof(argv[12]),atof(argv[13])},
        {atof(argv[14]),atof(argv[15]),atof(argv[16])},
        {atof(argv[17]),atof(argv[18]),atof(argv[19])},
        {atof(argv[20]),atof(argv[21]),atof(argv[22])},
        {atof(argv[23]),atof(argv[24]),atof(argv[25])},
        {atof(argv[26]),atof(argv[27]),atof(argv[28])},
        {atof(argv[29]),atof(argv[30]),atof(argv[31])},
        {atof(argv[32]),atof(argv[33]),atof(argv[34])}
    };

    // Populate the list of verts
    PointGeomSharedPtr verts[8];
    const int three = 3;

    for( int i = 0; i < nVerts; ++i ) {
        verts[i] = MemoryManager<PointGeom>
            ::AllocateSharedPtr( three,  i,   point[i][0],
                                 point[i][1], point[i][2] );
    }

    // /////////////////////////////////////////////////////////////////////
    // Set up Hexahedron Edges
    // SegGeom (int id, const int coordim), EdgeComponent(id, coordim)
    const int nEdges = 12;
    const int vertexConnectivity[][2] = {
        {0,1}, {1,2}, {2,3}, {0,3}, {0,4}, {1,5},
        {2,6}, {3,7}, {4,5}, {5,6}, {6,7}, {4,7}
    };

    // Populate the list of edges
    SegGeomSharedPtr edges[nEdges];
    for( int i = 0; i < nEdges; ++i ) {
        PointGeomSharedPtr vertsArray[2];
        for( int j = 0; j < 2; ++j ) {
            vertsArray[j] = verts[vertexConnectivity[i][j]];
        }
        edges[i] = MemoryManager<SegGeom>
            ::AllocateSharedPtr( i, three, vertsArray);
    }

    // ////////////////////////////////////////////////////////////////////
    // Set up Hexahedron faces
    const int nFaces = 6;
    const int edgeConnectivity[][4] = {
        {0,1,2,3}, {0,5,8,4}, {1,6,9,5},
        {2,7,10,6}, {3,7,11,4}, {8,9,10,11}
    };
    const bool   isEdgeFlipped[][4] = {
        {0,0,0,1}, {0,0,1,1}, {0,0,1,1},
        {0,0, 1,1}, {0,0, 1,1}, {0,0, 0, 1}
    };

    // Populate the list of faces
    QuadGeomSharedPtr faces[nFaces];
    for( int i = 0; i < nFaces; ++i ) {
        SegGeomSharedPtr edgeArray[4];
        Orientation eorientArray[4];
        for( int j = 0; j < 4; ++j ) {
            edgeArray[j]    = edges[edgeConnectivity[i][j]];
            eorientArray[j] = isEdgeFlipped[i][j] ? eBackwards : eForwards;
        }
        faces[i] = MemoryManager<QuadGeom>
            ::AllocateSharedPtr( i, edgeArray, eorientArray);
    }

    SpatialDomains::HexGeomSharedPtr geom =
        MemoryManager<SpatialDomains::HexGeom>::AllocateSharedPtr(faces);
    geom->SetOwnData();

    return geom;
}

SpatialDomains::PrismGeomSharedPtr CreatePrismGeom(int argc, char *argv[])
{
    if (argc != 29)
    {
        cout << "Insufficient points for a prism!" << endl;
    }
    
    // /////////////////////////////////////////////////////////////////////
    // Set up Prism vertex coordinates
    // PointGeom (const int coordim, const int vid, double x,
    //   double y, double z)
    const int nVerts = 6;
    const double point[][3] = {
        {atof(argv[11]),atof(argv[12]),atof(argv[13])},
        {atof(argv[14]),atof(argv[15]),atof(argv[16])},
        {atof(argv[17]),atof(argv[18]),atof(argv[19])},
        {atof(argv[20]),atof(argv[21]),atof(argv[22])},
        {atof(argv[23]),atof(argv[24]),atof(argv[25])},
        {atof(argv[26]),atof(argv[27]),atof(argv[28])}
    };

    // Populate the list of verts
    PointGeomSharedPtr verts[nVerts];
    const int three = 3;
    
    for( int i = 0; i < nVerts; ++i ) {
        verts[i] = MemoryManager<PointGeom>
            ::AllocateSharedPtr( three,  i,   point[i][0],
                                 point[i][1], point[i][2] );
    }

    // /////////////////////////////////////////////////////////////////////
    // Set up Prism Edges
    // SegGeom (int id, const int coordim), EdgeComponent(id, coordim)
    const int nEdges = 9;
    const int vertexConnectivity[][2] = {
        {0,1}, {1,2}, {3,2}, {0,3}, {0,4}, 
        {1,4}, {2,5}, {3,5}, {4,5}
    };
    
    // Populate the list of edges
    SegGeomSharedPtr edges[nEdges];
    for( int i = 0; i < nEdges; ++i ) {
        PointGeomSharedPtr vertsArray[2];
        for( int j = 0; j < 2; ++j ) {
            vertsArray[j] = verts[vertexConnectivity[i][j]];
        }
        edges[i] = MemoryManager<SegGeom>
            ::AllocateSharedPtr( i, three, vertsArray);
    }
    
    // ////////////////////////////////////////////////////////////////////
    // Set up Prism faces
    const int nFaces  = 5;
    const int edgeConnectivity[][4] = {
        {0,1,2,3}, 
        {0,5,4,-1}, // Triangular face
        {1,6,8,5},
        {2,6,7,-1}, // Triangular face 
        {3,7,8,4}
    };
    const bool   isEdgeFlipped[][4] = {
        {0,0,1,1}, 
        {0,0,1,0}, 
        {0,0,1,1},
        {0,0,1,0}, 
        {0,0,1,1}
    };
    
    // Populate the list of faces
    Geometry2DSharedPtr faces[5];

    for (int i = 0; i < nFaces; ++i) {
        if (i == 1 || i == 3) 
        {
            SegGeomSharedPtr edgeArray[3];
            Orientation eorientArray[3];
            
            for (int j = 0; j < 3; ++j) {
                edgeArray[j]    = edges[edgeConnectivity[i][j]];
                eorientArray[j] = isEdgeFlipped[i][j] ? eBackwards : eForwards;
            }
            
            faces[i] = MemoryManager<TriGeom>
                ::AllocateSharedPtr( i, edgeArray, eorientArray);
        }
        else
        {
            SegGeomSharedPtr edgeArray[4];
            Orientation eorientArray[4];
            
            for (int j = 0; j < 4; ++j) {
                edgeArray[j]    = edges[edgeConnectivity[i][j]];
                eorientArray[j] = isEdgeFlipped[i][j] ? eBackwards : eForwards;
            }
            faces[i] = MemoryManager<QuadGeom>
                ::AllocateSharedPtr( i, edgeArray, eorientArray);
        }
    }

    SpatialDomains::PrismGeomSharedPtr geom =
        MemoryManager<SpatialDomains::PrismGeom>::AllocateSharedPtr(faces);
    geom->SetOwnData();

    return geom;
}

SpatialDomains::PyrGeomSharedPtr CreatePyrGeom(int argc, char *argv[])
{
    if (argc != 26)
    {
        cout << "Insufficient points for a pyramid!" << endl;
    }
    
    // /////////////////////////////////////////////////////////////////////
    // Set up Pyramid vertex coordinates
    // PointGeom (const int coordim, const int vid, double x,
    //   double y, double z)
    const int nVerts = 6;
    const double point[][3] = {
        {atof(argv[11]),atof(argv[12]),atof(argv[13])},
        {atof(argv[14]),atof(argv[15]),atof(argv[16])},
        {atof(argv[17]),atof(argv[18]),atof(argv[19])},
        {atof(argv[20]),atof(argv[21]),atof(argv[22])},
        {atof(argv[23]),atof(argv[24]),atof(argv[25])}
    };

    // Populate the list of verts
    PointGeomSharedPtr verts[nVerts];
    const int three = 3;
    
    for( int i = 0; i < nVerts; ++i ) {
        verts[i] = MemoryManager<PointGeom>
            ::AllocateSharedPtr( three,  i,   point[i][0],
                                 point[i][1], point[i][2] );
    }

    // /////////////////////////////////////////////////////////////////////
    // Set up Pyramid Edges
    // SegGeom (int id, const int coordim), EdgeComponent(id, coordim)
    const int nEdges = 8;
    const int vertexConnectivity[][2] = {
        {0,1}, {1,2}, {3,2}, {0,3},
        {0,4}, {1,4}, {2,4}, {3,4}
    };

    // Populate the list of edges
    SegGeomSharedPtr edges[nEdges];
    for( int i = 0; i < nEdges; ++i ) {
        PointGeomSharedPtr vertsArray[2];
        for( int j = 0; j < 2; ++j ) {
            vertsArray[j] = verts[vertexConnectivity[i][j]];
        }
        edges[i] = MemoryManager<SegGeom>
            ::AllocateSharedPtr( i, three, vertsArray);
    }
    
    // ////////////////////////////////////////////////////////////////////
    // Set up Pyramid faces
    const int nFaces  = 5;
    const int edgeConnectivity[][4] = {
        {0,1,2,3}, 
        {0,5,4,-1}, // Triangular face
        {1,6,5,-1}, // Triangular face
        {2,6,7,-1}, // Triangular face 
        {3,7,4,-1}  // Triangular face
    };
    const bool   isEdgeFlipped[][4] = {
        {0,0,1,1}, 
        {0,0,1,0}, 
        {0,0,1,0},
        {0,0,1,0}, 
        {0,0,1,0}
    };
    
    // Populate the list of faces
    Geometry2DSharedPtr faces[5];

    for (int i = 0; i < nFaces; ++i) {
        if (i > 0)
        {
            SegGeomSharedPtr edgeArray[3];
            Orientation eorientArray[3];
            
            for (int j = 0; j < 3; ++j) {
                edgeArray[j]    = edges[edgeConnectivity[i][j]];
                eorientArray[j] = isEdgeFlipped[i][j] ? eBackwards : eForwards;
            }
            
            faces[i] = MemoryManager<TriGeom>
                ::AllocateSharedPtr( i, edgeArray, eorientArray);
        }
        else
        {
            SegGeomSharedPtr edgeArray[4];
            Orientation eorientArray[4];
            
            for (int j = 0; j < 4; ++j) {
                edgeArray[j]    = edges[edgeConnectivity[i][j]];
                eorientArray[j] = isEdgeFlipped[i][j] ? eBackwards : eForwards;
            }
            faces[i] = MemoryManager<QuadGeom>
                ::AllocateSharedPtr( i, edgeArray, eorientArray);
        }
    }

    SpatialDomains::PyrGeomSharedPtr geom =
        MemoryManager<SpatialDomains::PyrGeom>::AllocateSharedPtr(faces);
    geom->SetOwnData();

    return geom;
}

SpatialDomains::TetGeomSharedPtr CreateTetGeom(int argc, char *argv[])
{
    if (argc != 23)
    {
        cout << "Insufficient points for Tet!" << endl;
    }

    const int nVerts = 4;
    const double point[][3] = {
        {atof(argv[11]),atof(argv[12]),atof(argv[13])},
        {atof(argv[14]),atof(argv[15]),atof(argv[16])},
        {atof(argv[17]),atof(argv[18]),atof(argv[19])},
        {atof(argv[20]),atof(argv[21]),atof(argv[22])}
    };

    // ///////////////////////////////////////////////////////////////////////
    // Populate the list of verts
    // PointGeom (const int coordim, const int vid, double x, double y,
    //    double z)
    PointGeomSharedPtr verts[4];
    const int three = 3;
    for(int i=0; i < nVerts; ++i){
        verts[i] =  MemoryManager<PointGeom>
            ::AllocateSharedPtr(three,i,point[i][0],point[i][1],point[i][2]);
    }

    // /////////////////////////////////////////////////////////////////////
    // Set up Tetrahedron Edges
    // SegGeom (int id, const int coordim), EdgeComponent(id, coordim)
    const int nEdges = 6;
    const int vertexConnectivity[][2] = {
        {0,1},{1,2},{0,2},{0,3},{1,3},{2,3}
    };


    // Populate the list of edges
    SegGeomSharedPtr edges[nEdges];
    for(int i=0; i < nEdges; ++i){
        PointGeomSharedPtr vertsArray[2];
        for(int j=0; j<2; ++j){
            vertsArray[j] = verts[vertexConnectivity[i][j]];
        }
        edges[i] = MemoryManager<SegGeom>
            ::AllocateSharedPtr(i, three, vertsArray);
    }

    // //////////////////////////////////////////////////////////////
    // Set up Tetrahedron faces
    const int nFaces = 4;
    const int edgeConnectivity[][3] = {
        {0,1,2}, {0,4,3}, {1,5,4}, {2,5,3}
    };
    const bool   isEdgeFlipped[][3] = {
        {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}
    };

    // Populate the list of faces
    TriGeomSharedPtr faces[nFaces];
    for(int i=0; i < nFaces; ++i){
        SegGeomSharedPtr edgeArray[3];
        Orientation eorientArray[3];
        for(int j=0; j < 3; ++j){
            edgeArray[j] = edges[edgeConnectivity[i][j]];
            eorientArray[j] = isEdgeFlipped[i][j] ? eBackwards : eForwards;
        }
        faces[i] = MemoryManager<TriGeom>
            ::AllocateSharedPtr(i, edgeArray, eorientArray);
    }

    SpatialDomains::TetGeomSharedPtr geom =
        MemoryManager<SpatialDomains::TetGeom>::AllocateSharedPtr(faces);
    geom->SetOwnData();

    return geom;
}

