///////////////////////////////////////////////////////////////////////////////
//
// File FilterThresholdMax.cpp
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
// Description: Outputs time when solution first exceeds a threshold value.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Filters/FilterThresholdMax.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string FilterThresholdMax::className = GetFilterFactory().RegisterCreatorFunction("ThresholdMax", FilterThresholdMax::create);

        FilterThresholdMax::FilterThresholdMax(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::map<std::string, std::string> &pParams) :
            Filter(pSession)
        {
            ASSERTL0(pParams.find("ThresholdValue") != pParams.end(),
                     "Missing parameter 'ThresholdValue'.");
            m_thresholdValue = atof(pParams.find("ThresholdValue")->second.c_str());
            ASSERTL0(pParams.find("InitialValue") != pParams.end(),
                     "Missing parameter 'InitialValue'.");
            m_initialValue = atof(pParams.find("InitialValue")->second.c_str());

            if (pParams.find("StartTime") != pParams.end())
            {
                m_startTime = atof(pParams.find("StartTime")->second.c_str());
            }

            m_outputFile = pSession->GetSessionName() + "_max.fld";
            if (pParams.find("OutputFile") != pParams.end())
            {
                m_outputFile = pParams.find("OutputFile")->second;
            }

            m_thresholdVar = 0;
            if (pParams.find("ThresholdVar") != pParams.end())
            {
                std::string var = pParams.find("ThresholdVar")->second.c_str();
                std::vector<string> varlist = pSession->GetVariables();
                std::vector<string>::const_iterator x;
                ASSERTL0((x=std::find(varlist.begin(), varlist.end(), var)) != varlist.end(),
                         "Specified variable " + var +
                         " in ThresholdMax filter is not available.");
                m_thresholdVar = x - varlist.begin();
            }

            m_fld = MemoryManager<LibUtilities::FieldIO>::AllocateSharedPtr(pSession->GetComm());

        }

        FilterThresholdMax::~FilterThresholdMax()
        {

        }

        void FilterThresholdMax::v_Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            m_threshold = Array<OneD, NekDouble> (pFields[m_thresholdVar]->GetNpoints(), m_initialValue);
        }

        void FilterThresholdMax::v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            if (time < m_startTime)
            {
                return;
            }

            int i;
            NekDouble timestep = pFields[m_thresholdVar]->GetSession()->GetParameter("TimeStep");

            for (i = 0; i < pFields[m_thresholdVar]->GetNpoints(); ++i)
            {
                if (m_threshold[i] < timestep && pFields[m_thresholdVar]->GetPhys()[i] > m_thresholdValue)
                {
                    m_threshold[i] = time;
                }
            }
        }

        void FilterThresholdMax::v_Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                = pFields[m_thresholdVar]->GetFieldDefinitions();
            std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

            Array<OneD, NekDouble> vCoeffs(pFields[0]->GetNcoeffs());
            pFields[m_thresholdVar]->FwdTrans_IterPerExp(m_threshold, vCoeffs);

            // copy Data into FieldData and set variable
            for(int i = 0; i < FieldDef.size(); ++i)
            {
                // Could do a search here to find correct variable
                FieldDef[i]->m_fields.push_back("m");
                pFields[m_thresholdVar]->AppendFieldData(FieldDef[i], FieldData[i], vCoeffs);
            }

            m_fld->Write(m_outputFile,FieldDef,FieldData);

        }

        bool FilterThresholdMax::v_IsTimeDependent()
        {
            return true;
        }
    }
}
