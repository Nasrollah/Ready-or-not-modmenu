///////////////////////////////////////////////////////////////////////////////
//
// File FilterEnergy.h
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
// Description: Outputs solution fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_COMPRESSIBLEFLOWSOLVER_FILTERS_FILTERENERGY_H
#define NEKTAR_COMPRESSIBLEFLOWSOLVER_FILTERS_FILTERENERGY_H

#include <SolverUtils/Filters/FilterEnergyBase.h>

namespace Nektar
{
    namespace SolverUtils
    {
        class FilterEnergy : public FilterEnergyBase
        {
        public:
            friend class MemoryManager<FilterEnergy>;

            /// Creates an instance of this class
            static FilterSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const std::map<std::string, std::string> &pParams) {
                FilterSharedPtr p = MemoryManager<FilterEnergy>
                    ::AllocateSharedPtr(pSession, pParams);
                return p;
            }

            ///Name of the class
            static std::string className;

            SOLVER_UTILS_EXPORT FilterEnergy(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const std::map<std::string, std::string> &pParams);
            SOLVER_UTILS_EXPORT ~FilterEnergy();

        protected:
            virtual void v_GetVelocity(
                const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
                const int i,
                Array<OneD, NekDouble> &velocity);
            virtual Array<OneD, NekDouble> v_GetDensity(
                const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields);
        };
    }
}

#endif /* NEKTAR_SOLVERUTILS_FILTERS_FILTERCHECKPOINT_H */
