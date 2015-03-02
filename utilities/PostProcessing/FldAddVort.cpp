#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    int i,j;

    if(argc != 4)
    {
        fprintf(stderr,"Usage: FldAddVort  meshfile infld outfld\n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);


    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-3]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(vSession);//meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // Import field file.
    string fieldfile(argv[argc-2]);
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    LibUtilities::Import(fieldfile,fielddef,fielddata);
    bool useFFT = false;
    bool dealiasing = false;
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int expdim  = graphShPt->GetMeshDimension();
    int nfields = fielddef[0]->m_fields.size();
    int addfields = (nfields == 4)? 3:1;
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields + addfields);

    switch(expdim)
    {
    case 1:
        {
            ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 2,"Quasi-3D approach is only set up for 1 or 2 homogeneous directions");

            if(fielddef[0]->m_numHomogeneousDir == 1)
            {
                MultiRegions::ExpList2DHomogeneous1DSharedPtr Exp2DH1;

                // Define Homogeneous expansion
                //int nplanes = fielddef[0]->m_numModes[1];
                int nplanes;
                vSession->LoadParameter("HomModesZ",nplanes,fielddef[0]->m_numModes[1]);

                // choose points to be at evenly spaced points at
                // nplanes points
                const LibUtilities::PointsKey Pkey(nplanes,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[1],nplanes,Pkey);
                NekDouble ly = fielddef[0]->m_homogeneousLengths[0];

                Exp2DH1 = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,ly,useFFT,dealiasing,graphShPt);
                Exp[0] = Exp2DH1;

                for(i = 1; i < nfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::AllocateSharedPtr(*Exp2DH1);
                }
            }
            else if(fielddef[0]->m_numHomogeneousDir == 2)
            {
                MultiRegions::ExpList3DHomogeneous2DSharedPtr Exp3DH2;

                // Define Homogeneous expansion
                //int nylines = fielddef[0]->m_numModes[1];
                //int nzlines = fielddef[0]->m_numModes[2];

				int nylines;
				int nzlines;
				vSession->LoadParameter("HomModesY",nylines,fielddef[0]->m_numModes[1]);
				vSession->LoadParameter("HomModesZ",nzlines,fielddef[0]->m_numModes[2]);

                // choose points to be at evenly spaced points at
                // nplanes points
                const LibUtilities::PointsKey PkeyY(nylines,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  BkeyY(fielddef[0]->m_basis[1],nylines,PkeyY);

                const LibUtilities::PointsKey PkeyZ(nzlines,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  BkeyZ(fielddef[0]->m_basis[2],nzlines,PkeyZ);

                NekDouble ly = fielddef[0]->m_homogeneousLengths[0];
                NekDouble lz = fielddef[0]->m_homogeneousLengths[1];

                Exp3DH2 = MemoryManager<MultiRegions::ExpList3DHomogeneous2D>::AllocateSharedPtr(vSession,BkeyY,BkeyZ,ly,lz,useFFT,dealiasing,graphShPt);
                Exp[0] = Exp3DH2;

                for(i = 1; i < nfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous2D>::AllocateSharedPtr(*Exp3DH2);
                }
            }
            else
            {
                MultiRegions::ExpList1DSharedPtr Exp1D;
                Exp1D = MemoryManager<MultiRegions::ExpList1D>
                                                        ::AllocateSharedPtr(vSession,graphShPt);
                Exp[0] = Exp1D;
                for(i = 1; i < nfields + addfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList1D>
                                                        ::AllocateSharedPtr(*Exp1D);
                }
            }
        }
        break;
    case 2:
        {
            ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 1,"NumHomogeneousDir is only set up for 1");

            if(fielddef[0]->m_numHomogeneousDir == 1)
            {
                MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;

                // Define Homogeneous expansion
                //int nplanes = fielddef[0]->m_numModes[2];

				int nplanes;
				vSession->LoadParameter("HomModesZ",nplanes,fielddef[0]->m_numModes[2]);

                // choose points to be at evenly spaced points at
                // nplanes points
                const LibUtilities::PointsKey Pkey(nplanes,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[2],nplanes,Pkey);
                NekDouble lz = fielddef[0]->m_homogeneousLengths[0];

                Exp3DH1 = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,lz,useFFT,dealiasing,graphShPt);
                Exp[0] = Exp3DH1;

                for(i = 1; i < nfields + addfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(*Exp3DH1);

                }
            }
            else
            {
                MultiRegions::ExpList2DSharedPtr Exp2D;
                Exp2D = MemoryManager<MultiRegions::ExpList2D>
                                                        ::AllocateSharedPtr(vSession,graphShPt);
                Exp[0] =  Exp2D;

                for(i = 1; i < nfields + addfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList2D>
                        ::AllocateSharedPtr(*Exp2D);
                }
            }
        }
        break;
    case 3:
        {
            MultiRegions::ExpList3DSharedPtr Exp3D;
            Exp3D = MemoryManager<MultiRegions::ExpList3D>
                                                    ::AllocateSharedPtr(vSession,graphShPt);
            Exp[0] =  Exp3D;

            for(i = 1; i < nfields + addfields; ++i)
            {
                Exp[i] = MemoryManager<MultiRegions::ExpList3D>
                    ::AllocateSharedPtr(*Exp3D);
            }
        }
        break;
    default:
        ASSERTL0(false,"Expansion dimension not recognised");
        break;
    }
    //----------------------------------------------

    //----------------------------------------------
    // Copy data from field file
    for(j = 0; j < nfields; ++j)
    {
        for(int i = 0; i < fielddata.size(); ++i)
        {
            Exp[j]->ExtractDataToCoeffs(fielddef [i],
                                        fielddata[i],
                                        fielddef [i]->m_fields[j],
                                        Exp[j]->UpdateCoeffs());
        }
        Exp[j]->BwdTrans(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
    }
    //----------------------------------------------

    //----------------------------------------------
    // Compute gradients of fields
    ASSERTL0(nfields >= 3, "Need two fields (u,v) to add reentricity");
    int nq = Exp[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble> > grad(nfields*nfields);
    Array<OneD, Array<OneD, NekDouble> > outfield(addfields);

    for(i = 0; i < nfields*nfields; ++i)
    {
        grad[i] = Array<OneD, NekDouble>(nq);
    }

    for(i = 0; i < addfields; ++i)
    {
        outfield[i] = Array<OneD, NekDouble>(nq);
    }

    // Calculate Gradient & Vorticity
    if(nfields == 3)
    {
        for(i = 0; i < nfields; ++i)
        {
            Exp[i]->PhysDeriv(Exp[i]->GetPhys(), grad[i*nfields],grad[i*nfields+1]);
        }
        // W_z = Vx - Uy
        Vmath::Vsub(nq,grad[1*nfields+0],1,grad[0*nfields+1],1,outfield[0],1);
    }
    else
    {
        for(i = 0; i < nfields; ++i)
        {
            Exp[i]->PhysDeriv(Exp[i]->GetPhys(), grad[i*nfields],grad[i*nfields+1],grad[i*nfields+2]);
        }

        // W_x = Wy - Vz
        Vmath::Vsub(nq,grad[2*nfields+1],1,grad[1*nfields+2],1,outfield[0],1);
        // W_y = Uz - Wx
        Vmath::Vsub(nq,grad[0*nfields+2],1,grad[2*nfields+0],1,outfield[1],1);
        // W_z = Vx - Uy
        Vmath::Vsub(nq,grad[1*nfields+0],1,grad[0*nfields+1],1,outfield[2],1);
    }

    for (i = 0; i < addfields; ++i)
    {
        Exp[nfields + i]->FwdTrans(outfield[i], Exp[nfields+i]->UpdateCoeffs());
    }

    //-----------------------------------------------
    // Write solution to file with additional computed fields
    string   out(argv[argc-1]);
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = Exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    vector<string > outname;

    if(addfields == 1)
    {
        outname.push_back("W_z");
    }
    else
    {
        outname.push_back("W_x");
        outname.push_back("W_y");
        outname.push_back("W_z");
    }

    for(j = 0; j < nfields + addfields; ++j)
    {
        for(i = 0; i < FieldDef.size(); ++i)
        {
            if (j >= nfields)
            {
                FieldDef[i]->m_fields.push_back(outname[j-nfields]);
            }
            else
            {
                FieldDef[i]->m_fields.push_back(fielddef[i]->m_fields[j]);
            }
            Exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
        }
    }
    LibUtilities::Write(out, FieldDef, FieldData);
    //-----------------------------------------------

    return 0;
}

