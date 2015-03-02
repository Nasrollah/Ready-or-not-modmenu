///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysStaticCond.cpp
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
// Description: Implementation to linear solver using single-
//              or multi-level static condensation
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysStaticCond.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LinearAlgebra/StorageSmvBsr.hpp>
#include <LibUtilities/LinearAlgebra/SparseDiagBlkMatrix.hpp>
#include <LibUtilities/LinearAlgebra/SparseUtils.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysStaticCond
         *
         * Solves a linear system using single- or multi-level static
         * condensation.
         */

        /**
         * For a matrix system of the form @f[
         * \left[ \begin{array}{cc}
         * \boldsymbol{A} & \boldsymbol{B}\\
         * \boldsymbol{C} & \boldsymbol{D}
         * \end{array} \right]
         * \left[ \begin{array}{c} \boldsymbol{x_1}\\ \boldsymbol{x_2}
         * \end{array}\right]
         * = \left[ \begin{array}{c} \boldsymbol{y_1}\\ \boldsymbol{y_2}
         * \end{array}\right],
         * @f]
         * where @f$\boldsymbol{D}@f$ and
         * @f$(\boldsymbol{A-BD^{-1}C})@f$ are invertible, store and assemble
         * a static condensation system, according to a given local to global
         * mapping. #m_linSys is constructed by AssembleSchurComplement().
         * @param   mKey        Associated matrix key.
         * @param   pLocMatSys  LocalMatrixSystem
         * @param   locToGloMap Local to global mapping.
         */
        GlobalLinSysStaticCond::GlobalLinSysStaticCond(
            const GlobalLinSysKey                &pKey,
            const boost::weak_ptr<ExpList>       &pExpList,
            const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
                : GlobalLinSys(pKey, pExpList, pLocToGloMap),
                  m_locToGloMap (pLocToGloMap)
        {
        }

        void GlobalLinSysStaticCond::v_InitObject()
        {
            // Allocate memory for top-level structure
            SetupTopLevel(m_locToGloMap);

            // Construct this level
            Initialise(m_locToGloMap);
        }
        
        /**
         *
         */
        GlobalLinSysStaticCond::~GlobalLinSysStaticCond()
        {

        }
        
        
        /**
         *
         */
        void GlobalLinSysStaticCond::v_Solve(
            const Array<OneD, const NekDouble> &in,
                  Array<OneD,       NekDouble> &out,
            const AssemblyMapSharedPtr         &pLocToGloMap,
            const Array<OneD, const NekDouble> &dirForcing)
        {
            bool dirForcCalculated = (bool) dirForcing.num_elements();
            bool atLastLevel       = pLocToGloMap->AtLastLevel();
            int  scLevel           = pLocToGloMap->GetStaticCondLevel();
            
            int nGlobDofs          = pLocToGloMap->GetNumGlobalCoeffs();
            int nGlobBndDofs       = pLocToGloMap->GetNumGlobalBndCoeffs();
            int nDirBndDofs        = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobHomBndDofs    = nGlobBndDofs - nDirBndDofs;
            int nLocBndDofs        = pLocToGloMap->GetNumLocalBndCoeffs();
            int nIntDofs           = pLocToGloMap->GetNumGlobalCoeffs()
                - nGlobBndDofs;

            Array<OneD, NekDouble> F = m_wsp + 2*nLocBndDofs;
            Array<OneD, NekDouble> tmp;
            if(nDirBndDofs && dirForcCalculated)
            {
                Vmath::Vsub(nGlobDofs,in.get(),1,dirForcing.get(),1,F.get(),1);
            }
            else
            {
                Vmath::Vcopy(nGlobDofs,in.get(),1,F.get(),1);
            }
            
            NekVector<NekDouble> F_HomBnd(nGlobHomBndDofs,tmp=F+nDirBndDofs,
                                          eWrapper);
            NekVector<NekDouble> F_GlobBnd(nGlobBndDofs,F,eWrapper);
            NekVector<NekDouble> F_LocBnd(nLocBndDofs,0.0);
            NekVector<NekDouble> F_Int(nIntDofs,tmp=F+nGlobBndDofs,eWrapper);
            
            NekVector<NekDouble> V_GlobBnd(nGlobBndDofs,out,eWrapper);
            NekVector<NekDouble> V_GlobHomBnd(nGlobHomBndDofs,
                                              tmp=out+nDirBndDofs,
                                              eWrapper);
            NekVector<NekDouble> V_Int(nIntDofs,tmp=out+nGlobBndDofs,eWrapper);
            NekVector<NekDouble> V_LocBnd(nLocBndDofs,m_wsp,eWrapper);
            
            NekVector<NekDouble> V_GlobHomBndTmp(nGlobHomBndDofs,0.0);

            // set up normalisation factor for right hand side on first SC level
            DNekScalBlkMatSharedPtr sc = v_PreSolve(scLevel, F_GlobBnd);

            if(nGlobHomBndDofs)
            {
                // construct boundary forcing
                if( nIntDofs  && ((!dirForcCalculated) && (atLastLevel)) )
                {
                    DNekScalBlkMat &BinvD      = *m_BinvD;
                    DNekScalBlkMat &SchurCompl = *sc;

                    // include dirichlet boundary forcing
                    pLocToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                    V_LocBnd = BinvD*F_Int + SchurCompl*V_LocBnd;
                }
                else if((!dirForcCalculated) && (atLastLevel))
                {
                    // include dirichlet boundary forcing
                    DNekScalBlkMat &SchurCompl = *sc;
                    pLocToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                    V_LocBnd = SchurCompl*V_LocBnd;
                }
                else
                {
                    DNekScalBlkMat &BinvD      = *m_BinvD;
                    V_LocBnd = BinvD*F_Int;
                }
                
                pLocToGloMap->AssembleBnd(V_LocBnd,V_GlobHomBndTmp,
                                          nDirBndDofs);
                F_HomBnd = F_HomBnd - V_GlobHomBndTmp;

                // Transform from original basis to low energy
                v_BasisTransform(F, nDirBndDofs);

                // For parallel multi-level static condensation some
                // processors may have different levels to others. This
                // routine receives contributions to partition vertices from
                // those lower levels, whilst not sending anything to the
                // other partitions, and includes them in the modified right
                // hand side vector.
                int lcLevel = pLocToGloMap->GetLowestStaticCondLevel();
                if(atLastLevel && scLevel < lcLevel)
                {
                    // If this level is not the lowest level across all
                    // processes, we must do dummy communication for the
                    // remaining levels
                    Array<OneD, NekDouble> tmp(nGlobBndDofs);
                    for (int i = scLevel; i < lcLevel; ++i)
                    {
                        Vmath::Fill(nGlobBndDofs, 0.0, tmp, 1);
                        pLocToGloMap->UniversalAssembleBnd(tmp);
                        Vmath::Vcopy(nGlobHomBndDofs,
                                     tmp.get()+nDirBndDofs,          1,
                                     V_GlobHomBndTmp.GetPtr().get(), 1);
                        F_HomBnd = F_HomBnd - V_GlobHomBndTmp;
                    }
                }

                // solve boundary system
                if(atLastLevel)
                {
                    Array<OneD, NekDouble> pert(nGlobBndDofs,0.0);
                    NekVector<NekDouble>   Pert(nGlobBndDofs,pert,eWrapper);

                    // Solve for difference from initial solution given inout;
                    SolveLinearSystem(
                        nGlobBndDofs, F, pert, pLocToGloMap, nDirBndDofs);

                    // Transform back to original basis
                    v_BasisInvTransform(pert);

                    // Add back initial conditions onto difference
                    Vmath::Vadd(nGlobHomBndDofs,&out[nDirBndDofs],1,
                                &pert[nDirBndDofs],1,&out[nDirBndDofs],1);
                }
                else
                {
                    m_recursiveSchurCompl->Solve(F,
                                V_GlobBnd.GetPtr(),
                                pLocToGloMap->GetNextLevelLocalToGlobalMap());
                }
            }

            // solve interior system
            if(nIntDofs)
            {
                DNekScalBlkMat &invD  = *m_invD;

                if(nGlobHomBndDofs || nDirBndDofs)
                {
                    DNekScalBlkMat &C     = *m_C;

                    if(dirForcCalculated && nDirBndDofs)
                    {
                        pLocToGloMap->GlobalToLocalBnd(V_GlobHomBnd,V_LocBnd,
                                                      nDirBndDofs);
                    }
                    else
                    {
                        pLocToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                    }
                    F_Int = F_Int - C*V_LocBnd;
                }

                V_Int = invD*F_Int;
            }
        }


        /**
         * If at the last level of recursion (or the only level in the case of
         * single-level static condensation), assemble the Schur complement.
         * For other levels, in the case of multi-level static condensation,
         * the next level of the condensed system is computed.
         * @param   pLocToGloMap    Local to global mapping.
         */
        void GlobalLinSysStaticCond::v_Initialise(
                const boost::shared_ptr<AssemblyMap>& pLocToGloMap)
        {
            int nLocalBnd = m_locToGloMap->GetNumLocalBndCoeffs();
            int nGlobal = m_locToGloMap->GetNumGlobalCoeffs();
            m_wsp = Array<OneD, NekDouble>(2*nLocalBnd + nGlobal, 0.0);

            if (pLocToGloMap->AtLastLevel())
            {
                v_AssembleSchurComplement(pLocToGloMap);
            }
            else
            {
                ConstructNextLevelCondensedSystem(
                        pLocToGloMap->GetNextLevelLocalToGlobalMap());
            }
        }

        int GlobalLinSysStaticCond::v_GetNumBlocks()
        {
            return m_schurCompl->GetNumberOfBlockRows();
        }

        /**
         * For the first level in multi-level static condensation, or the only
         * level in the case of single-level static condensation, allocate the
         * condensed matrices and populate them with the local matrices
         * retrieved from the expansion list.
         * @param
         */
        void GlobalLinSysStaticCond::SetupTopLevel(
                const boost::shared_ptr<AssemblyMap>& pLocToGloMap)
        {
            int n;
            int n_exp = m_expList.lock()->GetNumElmts();

            const Array<OneD,const unsigned int>& nbdry_size
                    = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
            const Array<OneD,const unsigned int>& nint_size
                    = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();

            // Setup Block Matrix systems
            MatrixStorage blkmatStorage = eDIAGONAL;
            m_schurCompl = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nbdry_size, blkmatStorage);
            m_BinvD      = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nint_size , blkmatStorage);
            m_C          = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nbdry_size, blkmatStorage);
            m_invD       = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nint_size , blkmatStorage);

            for(n = 0; n < n_exp; ++n)
            {
                if (m_linSysKey.GetMatrixType() ==
                        StdRegions::eHybridDGHelmBndLam)
                {
                    DNekScalMatSharedPtr loc_mat
                        = GlobalLinSys::v_GetBlock(
                            m_expList.lock()->GetOffset_Elmt_Id(n));
                    m_schurCompl->SetBlock(n,n,loc_mat);
                }
                else
                {
                    DNekScalBlkMatSharedPtr loc_schur
                        = GlobalLinSys::v_GetStaticCondBlock(
                            m_expList.lock()->GetOffset_Elmt_Id(n));
                    DNekScalMatSharedPtr t;
                    m_schurCompl->SetBlock(n, n, t = loc_schur->GetBlock(0,0));
                    m_BinvD     ->SetBlock(n, n, t = loc_schur->GetBlock(0,1));
                    m_C         ->SetBlock(n, n, t = loc_schur->GetBlock(1,0));
                    m_invD      ->SetBlock(n, n, t = loc_schur->GetBlock(1,1));
                }
            }
        }

        /**
         *
         */
        void GlobalLinSysStaticCond::ConstructNextLevelCondensedSystem(
                        const AssemblyMapSharedPtr& pLocToGloMap)
        {
            int i,j,n,cnt;
            DNekScalBlkMatSharedPtr blkMatrices[4];

            // Create temporary matrices within an inner-local scope to ensure
            // any references to the intermediate storage is lost before
            // the recursive step, rather than at the end of the routine.
            // This allows the schur complement matrix from this level to be
            // disposed of in the next level after use without being retained
            // due to lingering shared pointers.
            {

                const Array<OneD,const unsigned int>& nBndDofsPerPatch
                                = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
                const Array<OneD,const unsigned int>& nIntDofsPerPatch
                                = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();

                // STEP 1:
                // Based upon the schur complement of the the current level we
                // will substructure this matrix in the form
                //      --     --
                //      | A   B |
                //      | C   D |
                //      --     --
                // All matrices A,B,C and D are (diagonal) blockmatrices.
                // However, as a start we will use an array of DNekMatrices as
                // it is too hard to change the individual entries of a
                // DNekScalBlkMatSharedPtr.

                // In addition, we will also try to ensure that the memory of
                // the blockmatrices will be contiguous. This will probably
                // enhance the efficiency
                // - Calculate the total number of entries in the blockmatrices
                int nPatches  = pLocToGloMap->GetNumPatches();
                int nEntriesA = 0; int nEntriesB = 0;
                int nEntriesC = 0; int nEntriesD = 0;

                for(i = 0; i < nPatches; i++)
                {
                    nEntriesA += nBndDofsPerPatch[i]*nBndDofsPerPatch[i];
                    nEntriesB += nBndDofsPerPatch[i]*nIntDofsPerPatch[i];
                    nEntriesC += nIntDofsPerPatch[i]*nBndDofsPerPatch[i];
                    nEntriesD += nIntDofsPerPatch[i]*nIntDofsPerPatch[i];
                }

                // Now create the DNekMatrices and link them to the memory
                // allocated above
                Array<OneD, DNekMatSharedPtr> substructuredMat[4]
                    = {Array<OneD, DNekMatSharedPtr>(nPatches),  //Matrix A
                       Array<OneD, DNekMatSharedPtr>(nPatches),  //Matrix B
                       Array<OneD, DNekMatSharedPtr>(nPatches),  //Matrix C
                       Array<OneD, DNekMatSharedPtr>(nPatches)}; //Matrix D

                // Initialise storage for the matrices. We do this separately
                // for each matrix so the matrices may be independently
                // deallocated when no longer required.
                Array<OneD, NekDouble> storageA(nEntriesA,0.0);
                Array<OneD, NekDouble> storageB(nEntriesB,0.0);
                Array<OneD, NekDouble> storageC(nEntriesC,0.0);
                Array<OneD, NekDouble> storageD(nEntriesD,0.0);

                Array<OneD, NekDouble> tmparray;
                PointerWrapper wType = eWrapper;
                int cntA = 0;
                int cntB = 0;
                int cntC = 0;
                int cntD = 0;

                for(i = 0; i < nPatches; i++)
                {
                    // Matrix A
                    tmparray = storageA+cntA;
                    substructuredMat[0][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nBndDofsPerPatch[i],
                                                        nBndDofsPerPatch[i],
                                                        tmparray, wType);
                    // Matrix B
                    tmparray = storageB+cntB;
                    substructuredMat[1][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nBndDofsPerPatch[i],
                                                        nIntDofsPerPatch[i],
                                                        tmparray, wType);
                    // Matrix C
                    tmparray = storageC+cntC;
                    substructuredMat[2][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nIntDofsPerPatch[i],
                                                        nBndDofsPerPatch[i],
                                                        tmparray, wType);
                    // Matrix D
                    tmparray = storageD+cntD;
                    substructuredMat[3][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nIntDofsPerPatch[i],
                                                        nIntDofsPerPatch[i],
                                                        tmparray, wType);

                    cntA += nBndDofsPerPatch[i] * nBndDofsPerPatch[i];
                    cntB += nBndDofsPerPatch[i] * nIntDofsPerPatch[i];
                    cntC += nIntDofsPerPatch[i] * nBndDofsPerPatch[i];
                    cntD += nIntDofsPerPatch[i] * nIntDofsPerPatch[i];
                }

                // Then, project SchurComplement onto
                // the substructured matrices of the next level
                DNekScalBlkMatSharedPtr SchurCompl  = m_schurCompl;
                DNekScalMatSharedPtr schurComplSubMat;
                int       schurComplSubMatnRows;
                Array<OneD, const int>       patchId, dofId;
                Array<OneD, const unsigned int>      isBndDof;
                Array<OneD, const NekDouble> sign;
                NekDouble scale;

                for(n = cnt = 0; n < SchurCompl->GetNumberOfBlockRows(); ++n)
                {
                    schurComplSubMat      = SchurCompl->GetBlock(n,n);
                    schurComplSubMatnRows = schurComplSubMat->GetRows();
                    
                    scale = SchurCompl->GetBlock(n,n)->Scale();
                    Array<OneD, NekDouble> schurSubMat
                        = SchurCompl->GetBlock(n,n)->GetOwnedMatrix()->GetPtr();
                    
                    patchId  = pLocToGloMap->GetPatchMapFromPrevLevel()
                               ->GetPatchId() + cnt;
                    dofId    = pLocToGloMap->GetPatchMapFromPrevLevel()
                               ->GetDofId()   + cnt;
                    isBndDof = pLocToGloMap->GetPatchMapFromPrevLevel()
                               ->IsBndDof() + cnt;
                    sign     = pLocToGloMap->GetPatchMapFromPrevLevel()
                               ->GetSign() + cnt;

                    // Set up  Matrix;
                    for(i = 0; i < schurComplSubMatnRows; ++i)
                    {
                        int pId = patchId[i];
                        Array<OneD, NekDouble> subMat0
                            = substructuredMat[0][pId]->GetPtr();
                        Array<OneD, NekDouble> subMat1
                            = substructuredMat[1][patchId[i]]->GetPtr();
                        Array<OneD, NekDouble> subMat2
                            = substructuredMat[2][patchId[i]]->GetPtr();
                        Array<OneD, NekDouble> subMat3
                            = substructuredMat[3][patchId[i]]->GetPtr();
                        int subMat0rows = substructuredMat[0][pId]->GetRows();
                        int subMat1rows = substructuredMat[1][pId]->GetRows();
                        int subMat2rows = substructuredMat[2][pId]->GetRows();
                        int subMat3rows = substructuredMat[3][pId]->GetRows();
                        
                        if(isBndDof[i])
                        {
                            for(j = 0; j < schurComplSubMatnRows; ++j)
                            {
                                ASSERTL0(patchId[i]==patchId[j],
                                         "These values should be equal");
                                
                                if(isBndDof[j])
                                {
                                    subMat0[dofId[i]+dofId[j]*subMat0rows] +=
                                        sign[i]*sign[j]*(
                                            scale*schurSubMat[
                                                i+j*schurComplSubMatnRows]);
                                }
                                else
                                {
                                    subMat1[dofId[i]+dofId[j]*subMat1rows] +=
                                        sign[i]*sign[j]*(
                                            scale*schurSubMat[
                                                i+j*schurComplSubMatnRows]);
                                }
                            }
                        }
                        else
                        {
                            for(j = 0; j < schurComplSubMatnRows; ++j)
                            {
                                ASSERTL0(patchId[i]==patchId[j],
                                         "These values should be equal");
                                
                                if(isBndDof[j])
                                {
                                    subMat2[dofId[i]+dofId[j]*subMat2rows] +=
                                        sign[i]*sign[j]*(
                                            scale*schurSubMat[
                                                i+j*schurComplSubMatnRows]);
                                }
                                else
                                {
                                    subMat3[dofId[i]+dofId[j]*subMat3rows] +=
                                        sign[i]*sign[j]*(
                                            scale*schurSubMat[
                                                i+j*schurComplSubMatnRows]);
                                }
                            }
                        }
                    }
                    cnt += schurComplSubMatnRows;
                }

                // STEP 2: condense the system
                // This can be done elementally (i.e. patch per patch)
                for(i = 0; i < nPatches; i++)
                {
                    if(nIntDofsPerPatch[i])
                    {
                        Array<OneD, NekDouble> subMat0
                            = substructuredMat[0][i]->GetPtr();
                        Array<OneD, NekDouble> subMat1
                            = substructuredMat[1][i]->GetPtr();
                        Array<OneD, NekDouble> subMat2
                            = substructuredMat[2][i]->GetPtr();
                        int subMat0rows = substructuredMat[0][i]->GetRows();
                        int subMat1rows = substructuredMat[1][i]->GetRows();
                        int subMat2rows = substructuredMat[2][i]->GetRows();
                        int subMat2cols = substructuredMat[2][i]->GetColumns();

                        // 1. D -> InvD
                        substructuredMat[3][i]->Invert();
                        // 2. B -> BInvD
                        (*substructuredMat[1][i]) = (*substructuredMat[1][i])*
                            (*substructuredMat[3][i]);
                        // 3. A -> A - BInvD*C (= schurcomplement)
                        // (*substructuredMat[0][i]) =(*substructuredMat[0][i])-
                        // (*substructuredMat[1][i])*(*substructuredMat[2][i]);
                        // Note: faster to use blas directly
                        Blas::Dgemm('N','N', subMat1rows, subMat2cols,
                                    subMat2rows, -1.0, &subMat1[0], subMat1rows,
                                    &subMat2[0], subMat2rows, 1.0, &subMat0[0],
                                    subMat0rows);
                    }
                }

                // STEP 3: fill the block matrices. However, do note that we
                // first have to convert them to a DNekScalMat in order to be
                // compatible with the first level of static condensation
                const Array<OneD,const unsigned int>& nbdry_size
                    = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
                const Array<OneD,const unsigned int>& nint_size
                    = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();
                MatrixStorage blkmatStorage = eDIAGONAL;

                blkMatrices[0] = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nbdry_size, blkmatStorage);
                blkMatrices[1] = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nint_size , blkmatStorage);
                blkMatrices[2] = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nbdry_size, blkmatStorage);
                blkMatrices[3] = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nint_size , blkmatStorage);

                DNekScalMatSharedPtr tmpscalmat;
                for(i = 0; i < nPatches; i++)
                {
                    for(j = 0; j < 4; j++)
                    {
                        tmpscalmat= MemoryManager<DNekScalMat>
                            ::AllocateSharedPtr(1.0,substructuredMat[j][i]);
                        blkMatrices[j]->SetBlock(i,i,tmpscalmat);
                    }
                }
            }

            // We've finished with the Schur complement matrix passed to this
            // level, so return the memory to the system. The Schur complement
            // matrix need only be retained at the last level. Save the other
            // matrices at this level though.
            m_schurCompl.reset();

            m_recursiveSchurCompl = v_Recurse(
                m_linSysKey, m_expList, blkMatrices[0], blkMatrices[1],
                blkMatrices[2], blkMatrices[3], pLocToGloMap);
        }
    }
}
