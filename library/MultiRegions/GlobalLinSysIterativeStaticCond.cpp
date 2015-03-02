///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysIterativeStaticCond.cpp
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

#include <MultiRegions/GlobalLinSysIterativeStaticCond.h>
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
         * @class GlobalLinSysIterativeStaticCond
         *
         * Solves a linear system iteratively using single- or multi-level
         * static condensation.
         */

        /**
         * Registers the class with the Factory.
         */
        string GlobalLinSysIterativeStaticCond::className
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "IterativeStaticCond",
                    GlobalLinSysIterativeStaticCond::create,
                    "Iterative static condensation.");

        string GlobalLinSysIterativeStaticCond::className2
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "IterativeMultiLevelStaticCond",
                    GlobalLinSysIterativeStaticCond::create,
                    "Iterative multi-level static condensation.");


        std::string GlobalLinSysIterativeStaticCond::storagedef = 
            LibUtilities::SessionReader::RegisterDefaultSolverInfo(
                "LocalMatrixStorageStrategy",
                "Sparse");
        std::string GlobalLinSysIterativeStaticCond::storagelookupIds[3] = {
            LibUtilities::SessionReader::RegisterEnumValue(
                "LocalMatrixStorageStrategy",
                "Contiguous",
                MultiRegions::eContiguous),
            LibUtilities::SessionReader::RegisterEnumValue(
                "LocalMatrixStorageStrategy",
                "Non-contiguous",
                MultiRegions::eNonContiguous),
            LibUtilities::SessionReader::RegisterEnumValue(
                "LocalMatrixStorageStrategy",
                "Sparse",
                MultiRegions::eSparse),
        };

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
        GlobalLinSysIterativeStaticCond::GlobalLinSysIterativeStaticCond(
            const GlobalLinSysKey                &pKey,
            const boost::weak_ptr<ExpList>       &pExpList,
            const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
            : GlobalLinSys          (pKey, pExpList, pLocToGloMap),
              GlobalLinSysIterative (pKey, pExpList, pLocToGloMap),
              GlobalLinSysStaticCond(pKey, pExpList, pLocToGloMap)
        {
            ASSERTL1((pKey.GetGlobalSysSolnType()==eIterativeStaticCond)||
                     (pKey.GetGlobalSysSolnType()==eIterativeMultiLevelStaticCond),
                     "This constructor is only valid when using static "
                     "condensation");
            ASSERTL1(pKey.GetGlobalSysSolnType()
                        == pLocToGloMap->GetGlobalSysSolnType(),
                     "The local to global map is not set up for the requested "
                     "solution type");
        }

        
        /**
         *
         */
        GlobalLinSysIterativeStaticCond::GlobalLinSysIterativeStaticCond(
            const GlobalLinSysKey                &pKey,
            const boost::weak_ptr<ExpList>       &pExpList,
            const DNekScalBlkMatSharedPtr         pSchurCompl,
            const DNekScalBlkMatSharedPtr         pBinvD,
            const DNekScalBlkMatSharedPtr         pC,
            const DNekScalBlkMatSharedPtr         pInvD,
            const boost::shared_ptr<AssemblyMap> &pLocToGloMap,
            const PreconditionerSharedPtr         pPrecon)
            : GlobalLinSys          (pKey, pExpList, pLocToGloMap),
              GlobalLinSysIterative (pKey, pExpList, pLocToGloMap),
              GlobalLinSysStaticCond(pKey, pExpList, pLocToGloMap)
        {
            m_schurCompl  = pSchurCompl;
            m_BinvD       = pBinvD;
            m_C           = pC;
            m_invD        = pInvD;
            m_precon      = pPrecon;
        }


        void GlobalLinSysIterativeStaticCond::v_InitObject()
        {
            MultiRegions::PreconditionerType pType
                = m_locToGloMap->GetPreconType();
            std::string PreconType
                = MultiRegions::PreconditionerTypeMap[pType];
            m_precon = GetPreconFactory().CreateInstance(
                PreconType,GetSharedThisPtr(),m_locToGloMap);

            // Allocate memory for top-level structure
            SetupTopLevel(m_locToGloMap);

            // Setup Block Matrix systems
            int n, n_exp = m_expList.lock()->GetNumElmts();

            MatrixStorage blkmatStorage = eDIAGONAL;
            const Array<OneD,const unsigned int>& nbdry_size
                    = m_locToGloMap->GetNumLocalBndCoeffsPerPatch();

            m_S1Blk      = MemoryManager<DNekScalBlkMat>
                ::AllocateSharedPtr(nbdry_size, nbdry_size , blkmatStorage);

            // Preserve original matrix in m_S1Blk
            for (n = 0; n < n_exp; ++n)
            {
                DNekScalMatSharedPtr mat = m_schurCompl->GetBlock(n, n);
                m_S1Blk->SetBlock(n, n, mat);
            }

            // Build preconditioner
            m_precon->BuildPreconditioner();

            // Do transform of Schur complement matrix
            for (n = 0; n < n_exp; ++n)
            {
                if (m_linSysKey.GetMatrixType() !=
                        StdRegions::eHybridDGHelmBndLam)
                {
                    DNekScalMatSharedPtr mat = m_S1Blk->GetBlock(n, n);
                    DNekScalMatSharedPtr t = m_precon->TransformedSchurCompl(
                        m_expList.lock()->GetOffset_Elmt_Id(n), mat);
                    m_schurCompl->SetBlock(n, n, t);
                }
            }

            // Construct this level
            Initialise(m_locToGloMap);
        }
        
        /**
         *
         */
        GlobalLinSysIterativeStaticCond::~GlobalLinSysIterativeStaticCond()
        {
            
        }

        DNekScalBlkMatSharedPtr GlobalLinSysIterativeStaticCond::
            v_GetStaticCondBlock(unsigned int n)
        {
            DNekScalBlkMatSharedPtr schurComplBlock;
            int  scLevel           = m_locToGloMap->GetStaticCondLevel();
            DNekScalBlkMatSharedPtr sc = scLevel == 0 ? m_S1Blk : m_schurCompl;
            DNekScalMatSharedPtr    localMat = sc->GetBlock(n,n);
            unsigned int nbdry    = localMat->GetRows();
            unsigned int nblks    = 1;
            unsigned int esize[1] = {nbdry};

            schurComplBlock = MemoryManager<DNekScalBlkMat>
                ::AllocateSharedPtr(nblks, nblks, esize, esize);
            schurComplBlock->SetBlock(0, 0, localMat);

            return schurComplBlock;
        }

        /**
         * Assemble the schur complement matrix from the block matrices stored
         * in #m_blkMatrices and the given local to global mapping information.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysIterativeStaticCond::v_AssembleSchurComplement(
            const AssemblyMapSharedPtr pLocToGloMap)
        {
            int i,j,n,cnt,gid1,gid2;
            NekDouble sign1,sign2;

            bool doGlobalOp = m_expList.lock()->GetGlobalOptParam()->
                DoGlobalMatOp(m_linSysKey.GetMatrixType());

            // Set up unique map
            v_UniqueMap();

            // Build precon again if we in multi-level static condensation (a
            // bit of a hack)
            if (m_linSysKey.GetGlobalSysSolnType()==eIterativeMultiLevelStaticCond)
            {
                MultiRegions::PreconditionerType pType
                    = m_locToGloMap->GetPreconType();
                std::string PreconType
                    = MultiRegions::PreconditionerTypeMap[pType];
                m_precon = GetPreconFactory().CreateInstance(
                    PreconType,GetSharedThisPtr(),m_locToGloMap);
                m_precon->BuildPreconditioner();
            }

            if (!doGlobalOp)
            {
                PrepareLocalSchurComplement();
                return;
            }

            int nBndDofs  = pLocToGloMap->GetNumGlobalBndCoeffs();
            int NumDirBCs = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            unsigned int rows = nBndDofs - NumDirBCs;
            unsigned int cols = nBndDofs - NumDirBCs;

            // COO sparse storage to assist in assembly
            COOMatType gmat_coo;

            // Get the matrix storage structure
            // (whether to store only one triangular part, if symmetric)
            MatrixStorage matStorage = eFULL;

            // assemble globally
            DNekScalMatSharedPtr loc_mat;
            int loc_lda;
            for(n = cnt = 0; n < m_schurCompl->GetNumberOfBlockRows(); ++n)
            {
                loc_mat = m_schurCompl->GetBlock(n,n);
                loc_lda = loc_mat->GetRows();

                // Set up  Matrix;
                for(i = 0; i < loc_lda; ++i)
                {
                    gid1  = pLocToGloMap->GetLocalToGlobalBndMap (cnt + i)
                                                                    - NumDirBCs;
                    sign1 = pLocToGloMap->GetLocalToGlobalBndSign(cnt + i);

                    if(gid1 >= 0)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2  = pLocToGloMap->GetLocalToGlobalBndMap(cnt+j)
                                                                 - NumDirBCs;
                            sign2 = pLocToGloMap->GetLocalToGlobalBndSign(cnt+j);

                            if (gid2 >= 0)
                            {
                                gmat_coo[std::make_pair(gid1,gid2)] +=
                                    sign1*sign2*(*loc_mat)(i,j);
                            }
                        }
                    }
                }
                cnt += loc_lda;
            }

            DNekSmvBsrDiagBlkMat::SparseStorageSharedPtrVector
                sparseStorage (1);

            BCOMatType partMat;
            convertCooToBco(rows, cols, 1, gmat_coo, partMat);

            sparseStorage[0] =
                 MemoryManager<DNekSmvBsrDiagBlkMat::StorageType>::
                    AllocateSharedPtr(rows, cols, 1, partMat, matStorage );

            // Create block diagonal matrix
            m_sparseSchurCompl = MemoryManager<DNekSmvBsrDiagBlkMat>::
                                            AllocateSharedPtr(sparseStorage);
        }


        /**
         * Populates sparse block-diagonal schur complement matrix from
         * the block matrices stored in #m_blkMatrices.
         */
        void GlobalLinSysIterativeStaticCond::PrepareLocalSchurComplement()
        {
            LocalMatrixStorageStrategy storageStrategy =
                m_expList.lock()->GetSession()->
                    GetSolverInfoAsEnum<LocalMatrixStorageStrategy>(
                                       "LocalMatrixStorageStrategy");

            switch(storageStrategy)
            {
                case MultiRegions::eContiguous:
                case MultiRegions::eNonContiguous:
                {
                    size_t storageSize = 0;
                    int nBlk           = m_schurCompl->GetNumberOfBlockRows();

                    m_scale = Array<OneD, NekDouble> (nBlk, 1.0);
                    m_rows  = Array<OneD, unsigned int> (nBlk, 0U);

                    // Determine storage requirements for dense blocks.
                    for (int i = 0; i < nBlk; ++i)
                    {
                        m_rows[i]    = m_schurCompl->GetBlock(i,i)->GetRows();
                        m_scale[i]   = m_schurCompl->GetBlock(i,i)->Scale();
                        storageSize += m_rows[i] * m_rows[i];
                    }

                    // Assemble dense storage blocks.
                    DNekScalMatSharedPtr loc_mat;
                    m_denseBlocks.resize(nBlk);
                    double *ptr = 0;

                    if (MultiRegions::eContiguous == storageStrategy)
                    {
                        m_storage.resize    (storageSize);
                        ptr = &m_storage[0];
                    }

                    for (unsigned int n = 0; n < nBlk; ++n)
                    {
                        loc_mat = m_schurCompl->GetBlock(n,n);

                        if (MultiRegions::eContiguous == storageStrategy)
                        {
                            int loc_lda      = loc_mat->GetRows();
                            int blockSize    = loc_lda * loc_lda;
                            m_denseBlocks[n] = ptr;
                            for(int i = 0; i < loc_lda; ++i)
                            {
                                for(int j = 0; j < loc_lda; ++j)
                                {
                                    ptr[j*loc_lda+i] = (*loc_mat)(i,j);
                                }
                            }
                            ptr += blockSize;
                            GlobalLinSys::v_DropStaticCondBlock(
                                m_expList.lock()->GetOffset_Elmt_Id(n));
                        }
                        else
                        {
                            m_denseBlocks[n] = loc_mat->GetRawPtr();
                        }
                    }
                    break;
                }
                case MultiRegions::eSparse:
                {
                    DNekScalMatSharedPtr loc_mat;
                    int loc_lda;
                    int blockSize = 0;
                    
                    // First run through to split the set of local matrices into
                    // partitions of fixed block size, and count number of local
                    // matrices that belong to each partition.
                    std::vector<std::pair<int,int> > partitions;
                    for(int n = 0; n < m_schurCompl->GetNumberOfBlockRows(); ++n)
                    {
                        loc_mat = m_schurCompl->GetBlock(n,n);
                        loc_lda = loc_mat->GetRows();

                        ASSERTL1(loc_lda >= 0,
                                 boost::lexical_cast<std::string>(n) + "-th "
                                 "matrix block in Schur complement has "
                                 "rank 0!");

                        if (blockSize == loc_lda)
                        {
                            partitions[partitions.size()-1].first++;
                        }
                        else
                        {
                            blockSize = loc_lda;
                            partitions.push_back(make_pair(1,loc_lda));
                        }
                    }

                    MatrixStorage matStorage = eFULL;

                    // Create a vector of sparse storage holders
                    DNekSmvBsrDiagBlkMat::SparseStorageSharedPtrVector
                            sparseStorage (partitions.size());

                    for (int part = 0, n = 0; part < partitions.size(); ++part)
                    {
                        BCOMatType partMat;

                        for(int k = 0; k < partitions[part].first; ++k, ++n)
                        {
                            loc_mat = m_schurCompl->GetBlock(n,n);
                            loc_lda = loc_mat->GetRows();

                            ASSERTL1(loc_lda == partitions[part].second,
                                     boost::lexical_cast<std::string>(n) + "-th"
                                     " matrix block in Schur complement has "
                                     "unexpected rank");

                            partMat[make_pair(k,k)] = BCOEntryType(
                                loc_lda*loc_lda, loc_mat->GetRawPtr());

                            GlobalLinSys::v_DropStaticCondBlock(
                                m_expList.lock()->GetOffset_Elmt_Id(n));
                        }

                        sparseStorage[part] =
                        MemoryManager<DNekSmvBsrDiagBlkMat::StorageType>::
                            AllocateSharedPtr(
                                partitions[part].first, partitions[part].first,
                                partitions[part].second, partMat, matStorage );
                    }

                    // Create block diagonal matrix
                    m_sparseSchurCompl = MemoryManager<DNekSmvBsrDiagBlkMat>::
                                            AllocateSharedPtr(sparseStorage);

                    break;
                }
                default:
                    ErrorUtil::NekError("Solver info property \
                        LocalMatrixStorageStrategy takes values \
                        Contiguous, Non-contiguous and Sparse");
            }
        }

        /**
         *
         */
        void GlobalLinSysIterativeStaticCond::v_DoMatrixMultiply(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            int nLocal = m_locToGloMap->GetNumLocalBndCoeffs();
            int nDir = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            bool doGlobalOp = m_expList.lock()->GetGlobalOptParam()->
                    DoGlobalMatOp(m_linSysKey.GetMatrixType());

            if(doGlobalOp)
            {
                // Do matrix multiply globally
                Array<OneD, NekDouble> in  = pInput  + nDir;
                Array<OneD, NekDouble> out = pOutput + nDir;

                m_sparseSchurCompl->Multiply(in,out);
                m_locToGloMap->UniversalAssembleBnd(pOutput, nDir);
            }
            else if (m_sparseSchurCompl)
            {
                // Do matrix multiply locally using block-diagonal sparse matrix
                Array<OneD, NekDouble> tmp = m_wsp + nLocal;

                m_locToGloMap->GlobalToLocalBnd(pInput, m_wsp);
                m_sparseSchurCompl->Multiply(m_wsp,tmp);
                m_locToGloMap->AssembleBnd(tmp, pOutput);
            }
            else
            {
                // Do matrix multiply locally, using direct BLAS calls
                m_locToGloMap->GlobalToLocalBnd(pInput, m_wsp);
                int i, cnt;
                Array<OneD, NekDouble> tmpout = m_wsp + nLocal;
                for (i = cnt = 0; i < m_denseBlocks.size(); cnt += m_rows[i], ++i)
                {
                    const int rows = m_rows[i];
                    Blas::Dgemv('N', rows, rows,
                                m_scale[i], m_denseBlocks[i], rows, 
                                m_wsp.get()+cnt, 1, 
                                0.0, tmpout.get()+cnt, 1);
                }
                m_locToGloMap->AssembleBnd(tmpout, pOutput);
            }
        }

        void GlobalLinSysIterativeStaticCond::v_UniqueMap()
        {
            m_map = m_locToGloMap->GetGlobalToUniversalBndMapUnique();
        }

        DNekScalBlkMatSharedPtr GlobalLinSysIterativeStaticCond::v_PreSolve(
            int                     scLevel,
            NekVector<NekDouble>   &F_GlobBnd)
        {
            if (scLevel == 0)
            {
                Set_Rhs_Magnitude(F_GlobBnd);
                return m_S1Blk;
            }
            else
            {
                return m_schurCompl;
            }
        }

        void GlobalLinSysIterativeStaticCond::v_BasisTransform(
            Array<OneD, NekDouble>& pInOut,
            int                     offset)
        {
            m_precon->DoTransformToLowEnergy(pInOut, offset);
        }

        void GlobalLinSysIterativeStaticCond::v_BasisInvTransform(
            Array<OneD, NekDouble>& pInOut)
        {
            m_precon->DoTransformFromLowEnergy(pInOut);
        }

        GlobalLinSysStaticCondSharedPtr GlobalLinSysIterativeStaticCond::v_Recurse(
            const GlobalLinSysKey                &mkey,
            const boost::weak_ptr<ExpList>       &pExpList,
            const DNekScalBlkMatSharedPtr         pSchurCompl,
            const DNekScalBlkMatSharedPtr         pBinvD,
            const DNekScalBlkMatSharedPtr         pC,
            const DNekScalBlkMatSharedPtr         pInvD,
            const boost::shared_ptr<AssemblyMap> &l2gMap)
        {
            GlobalLinSysIterativeStaticCondSharedPtr sys = MemoryManager<
                GlobalLinSysIterativeStaticCond>::AllocateSharedPtr(
                    mkey, pExpList, pSchurCompl, pBinvD, pC, pInvD, l2gMap,
                    m_precon);
            sys->Initialise(l2gMap);
            return sys;
        }
    }
}
