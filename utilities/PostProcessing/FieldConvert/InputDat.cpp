////////////////////////////////////////////////////////////////////////////////
//
//  File: InputDat.cpp
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
//  Description: Read tecplot dat for isontours in 2D triangular FE block format
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include <tinyxml.h>

#include "InputDat.h"

namespace Nektar
{
    namespace Utilities
    {

        ModuleKey InputDat::m_className[1] = {
            GetModuleFactory().RegisterCreatorFunction(
                    ModuleKey(eInputModule, "dat"),
                    InputDat::create,
                    "Reads Tecplot dat file for FE block triangular format."),
        };

        /**
         * @brief Set up InputDat object.
         *
         */
        InputDat::InputDat(FieldSharedPtr f) : InputModule(f)
        {
            m_allowedFiles.insert("dat");
            f->m_fieldPts = MemoryManager<FieldPts>::AllocateSharedPtr();
        }

        InputDat::~InputDat()
        {
        }

        /**
         *
         */
        void InputDat::Process(po::variables_map &vm)
        {

            if(m_f->m_verbose)
            {
                cout << "Processing input dat file" << endl;
            }

            string      line, word, tag;
            std::ifstream datFile;
            stringstream s;

            // Open the file stream.
            string fname = m_f->m_inputfiles["dat"][0];


            datFile.open(fname.c_str());
            if (!datFile.good())
            {
                cerr << "Error opening file: " << fname << endl;
                abort();
            }

            // read variables
            while (!datFile.eof())
            {
                getline(datFile, line);

                if(line.find("VARIABLES") != string::npos)
                {
                    vector<string> variables;
                    std::size_t  pos = line.find('=');
                    pos++;

                    // note this expects a comma separated list but
                    // does not work for white space separated lists!
                    bool valid = ParseUtils::GenerateOrderedStringVector(
                                        line.substr(pos).c_str(), variables);
                    ASSERTL0(valid,"Unable to process list of field variable in "
                             " VARIABLES list:  "+ line.substr(pos));


                    // currently assum there are x y and z coordinates
                    for(int i = 3; i < variables.size(); ++i)
                    {
                        m_f->m_fieldPts->m_fields.push_back(variables[i]);
                    }

                    break;
                }
            }

            // set up basic parameters
            m_f->m_fieldPts->m_ptsDim  = 3;
            m_f->m_fieldPts->m_nFields = m_f->m_fieldPts->m_fields.size();
            m_f->m_fieldPts->m_pts     = Array<OneD, Array<OneD, NekDouble> > (
                                                m_f->m_fieldPts->m_nFields +
                                                m_f->m_fieldPts->m_ptsDim);
            m_f->m_fieldPts->m_ptype  = ePtsTriBlock;

            // read zones
            while (!datFile.eof())
            {
                getline(datFile, line);

                if((line.find("ZONE") != string::npos)||
                   (line.find("Zone") != string::npos)||
                   (line.find("zone") != string::npos))
                {
                    ReadTecplotFEBlockZone(datFile,line);
                }
            }

            datFile.close();
        }

        void InputDat::ReadTecplotFEBlockZone(
                std::ifstream   &datFile,
                string          &line)
        {

            ASSERTL0(line.find("FEBlock") != string::npos,
                     "Routine only set up for FEBLock format");
            ASSERTL0(line.find("ET") != string::npos,
                     "Routine only set up TRIANLES");

            // read the number of nodes

            stringstream s;
            string tag;
            int start,end;

            s.clear();
            s.str(line);
            tag = s.str();

            // read the number of vertices
            start = tag.find("N=");
            end   = tag.find_first_of(',',start);
            int nvert = atoi(tag.substr(start+2,end).c_str());

            // read the number of elements
            start = tag.find("E=");
            end   = tag.find_first_of(',',start);
            int nelmt = atoi(tag.substr(start+2,end).c_str());


            // set-up or extend m_pts array;
            int norigpts  = m_f->m_fieldPts->m_pts[0].num_elements();
            int totfields = m_f->m_fieldPts->m_pts.num_elements();
            Array<OneD, Array<OneD, NekDouble> > origpts(totfields);
            for(int i = 0; i < totfields; ++i)
            {
                origpts[i] = m_f->m_fieldPts->m_pts[i];
                m_f->m_fieldPts->m_pts[i] =
                                    Array<OneD, NekDouble>(norigpts + nvert);
            }

            NekDouble value;
            for(int n = 0; n < totfields; ++n)
            {

                for(int i = 0; i < norigpts; ++i)
                {
                    m_f->m_fieldPts->m_pts[n][i] = origpts[n][i];
                }
                for(int i = 0; i < nvert; ++i)
                {
                    datFile >> value;
                    m_f->m_fieldPts->m_pts[n][norigpts+i] = value;
                }
            }

            // read connectivity and add to list
            int intvalue;
            Array<OneD, int> conn(3*nelmt);
            for(int i = 0; i < 3*nelmt; ++i)
            {
                datFile >> intvalue;
                intvalue -=1; // decrement intvalue by 1 for c array convention
                conn[i] = norigpts + intvalue;
            }
            m_f->m_fieldPts->m_ptsConn.push_back(conn);

            getline(datFile, line);
        }
    }
}

