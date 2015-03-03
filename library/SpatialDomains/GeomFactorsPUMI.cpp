////////////////////////////////////////////////////////////////////////////////
//
//  File: GeomFactors.cpp
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
//  Description: Geometric factors base class.
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/GeomFactorsPUMI.h>
#include <LibUtilities/Foundations/Interp.h>

#include <PUMI_NektarppAdapt.h>
#include <pumi_mesh.h>
#include <CrvTri2.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace pumi;

namespace Nektar
{
    namespace SpatialDomains
    {
        // constructor
        GeomFactorsPUMI::GeomFactorsPUMI(
                const GeomType                              gtype,
                const int                                   coordim,
                const StdRegions::StdExpansionSharedPtr    &xmap,
                const Array<OneD, Array<OneD, NekDouble> > &coords,
                const int globalid,
                CrvEnt * crvent) :
            m_globalId(globalid),
            m_crvEnt(crvent),
            GeomFactors(gtype, coordim, xmap, coords)
        {
            CheckIfValid(); 
        }

        // copy constructor
        GeomFactorsPUMI::GeomFactorsPUMI(const GeomFactorsPUMI &S) :
            GeomFactors(S)
        {
          m_globalId = S.m_globalId;
          m_crvEnt = S.m_crvEnt;
        }

        // dtor
        GeomFactorsPUMI::~GeomFactorsPUMI()
        {
        }

        // assignment
        bool operator==(const GeomFactorsPUMI &lhs, const GeomFactorsPUMI &rhs)
        {
            if(!(lhs.m_type == rhs.m_type))
            {
                return false;
            }

            if(!(lhs.m_expDim == rhs.m_expDim))
            {
                return false;
            }

            if(!(lhs.m_coordDim == rhs.m_coordDim))
            {
                return false;
            }

            const Array<OneD, const NekDouble> jac_lhs =
                            lhs.ComputeJac(lhs.m_xmap->GetPointsKeys());
            const Array<OneD, const NekDouble> jac_rhs =
                            rhs.ComputeJac(rhs.m_xmap->GetPointsKeys());
            if(!(jac_lhs == jac_rhs))
            {
                return false;
            }

            return true;
        }

        DerivStorage GeomFactorsPUMI::ComputeDeriv(
                const LibUtilities::PointsKeyVector &keyTgt) const
        {
            ASSERTL1(keyTgt.size() == m_expDim,
                     "Dimension of target point distribution does not match "
                     "expansion dimension.");
            
//            std::cout << "GeomFactorsPUMI::ComputeDeriv" << std::endl;
            int i = 0, j = 0;
            int nqtot_map      = 1;
            int nqtot_tbasis   = 1;
            DerivStorage deriv = DerivStorage(m_expDim);
            DerivStorage d_map = DerivStorage(m_expDim);
            LibUtilities::PointsKeyVector map_points(m_expDim);

            // Allocate storage and compute number of points
            for (i = 0; i < m_expDim; ++i)
            {
                map_points[i]  = m_xmap->GetBasis(i)->GetPointsKey();
                nqtot_map     *= map_points[i].GetNumPoints();
                nqtot_tbasis  *= keyTgt[i].GetNumPoints();
                deriv[i] = Array<OneD, Array<OneD,NekDouble> >(m_coordDim);
                d_map[i] = Array<OneD, Array<OneD,NekDouble> >(m_coordDim);
            }

            // Calculate local derivatives
            for(i = 0; i < m_coordDim; ++i)
            {
                Array<OneD, NekDouble> tmp(nqtot_map);
                // Transform from coefficient space to physical space
                m_xmap->BwdTrans(m_coords[i], tmp);

                // Allocate storage and take the derivative (calculated at the
                // points as specified in 'Coords')
                for (j = 0; j < m_expDim; ++j)
                {
                    d_map[j][i] = Array<OneD,NekDouble>(nqtot_map);
                    deriv[j][i] = Array<OneD,NekDouble>(nqtot_tbasis);
                    m_xmap->StdPhysDeriv(j, tmp, d_map[j][i]);
                }
            }

            for (i = 0; i < m_coordDim; ++i)
            {
                // Interpolate the derivatives:
                // - from the points as defined in the mapping ('Coords')
                // - to the points we at which we want to know the metrics
                //   ('tbasis')
                bool same = true;
                for (j = 0; j < m_expDim; ++j)
                {
                    same = same && (map_points[j] == keyTgt[j]);
                }
                if( same )
                {
                    for (j = 0; j < m_expDim; ++j)
                    {
                        deriv[j][i] = d_map[j][i];
                    }
                }
                else
                {
                    for (j = 0; j < m_expDim; ++j)
                    {
                        Interp(map_points, d_map[j][i], keyTgt, deriv[j][i]);
                    }
                }
            }

            int tmpNumPts = keyTgt[0].GetNumPoints();
            int tmpNumPts1 = 0;
            int tmpNumPts2 = 0;
            Array<OneD, NekDouble> keyPts(tmpNumPts);
            Array<OneD, NekDouble> keyPts1(tmpNumPts1);
            Array<OneD, NekDouble> keyPts2(tmpNumPts2);

            LibUtilities::PointsManager()[keyTgt[0]]->GetPoints(keyPts);

            if(m_crvEnt->get_exp_dim() > 1){
              tmpNumPts1 = keyTgt[1].GetNumPoints();
              keyPts1 = Array<OneD, NekDouble>(tmpNumPts1);
              LibUtilities::PointsManager()[keyTgt[1]]->GetPoints(keyPts1);
            }

            if(m_crvEnt->get_exp_dim() > 2){
              tmpNumPts2 = keyTgt[2].GetNumPoints();
              keyPts2 = Array<OneD, NekDouble>(tmpNumPts2);
              LibUtilities::PointsManager()[keyTgt[2]]->GetPoints(keyPts2);
            }

/*
            std::cout << "GenGeomFactor: deriv: " << std::endl;
            std::cout << "Original Values: " << std::endl;
            std::cout << "m_expDim: " << m_expDim << std::endl;
            std::cout << "m_coordDim: " << m_coordDim << std::endl;
            std::cout << "nqtot_tbasis: " << nqtot_tbasis << std::endl;
            for(int i = 0; i < m_expDim; ++i)
              for(int j = 0; j < m_coordDim; ++j)
                for(int k = 0; k < nqtot_tbasis; ++k)
                  std::cout << deriv[i][j][k] << std::endl;
*/
//            std::cout << "New Values Computed by PUMI: " << std::endl;

//            std::cout << "printing out control points of beztri : " << std::endl;
//              m_beztriPtr->print();
 
        if(m_crvEnt->get_exp_dim() == 1) {  // input crv ent is a mesh edge
          for(int j = 0; j < tmpNumPts; ++j) // eta1 dir
          {
            Point3d in_eta_3d = Point3d(keyPts[j], 0.0, 0.0);
            Point3d dxyz_deta1;
            int dir = 0;
            m_crvEnt->diff_by_eta(in_eta_3d, dir, dxyz_deta1);
//            std::cout << "BezSeg: dxyz_deta1: (" << dxyz_deta1[0] <<", " << dxyz_deta1[1] << ", " << dxyz_deta1[2] << ")" << std::endl;
            deriv[0][0][j] = dxyz_deta1[0];
            deriv[0][1][j] = dxyz_deta1[1];
            deriv[0][2][j] = dxyz_deta1[2];
          }

        }

        if(m_crvEnt->get_exp_dim() == 3) {  // input is a region
          for(int i = 0; i < tmpNumPts2; ++i) // eta3 dir
          {  
            for(int j = 0; j < tmpNumPts1; ++j) // eta2 dir
            {
              for(int k = 0; k < tmpNumPts; ++k) // eta1 dir             
              {
                Point3d dxyz_deta1(0.0, 0.0, 0.0);
                Point3d dxyz_deta2(0.0, 0.0, 0.0);
                Point3d dxyz_deta3(0.0, 0.0, 0.0);
                // tensor/collapsed coords
                Point3d ceta = Point3d(keyPts[k], keyPts1[j], keyPts2[i]);
                // convert to standard cartesian coords
                Point3d eta = Point3d(0.0, 0.0, 0.0);
                eta[0] = (1.0 + ceta[0]) * (1.0 - ceta[2]) * (1.0 - ceta[1]) / 4.0 - 1.0;
                eta[1] = (1.0 + ceta[1]) * (1.0 - ceta[2]) / 2.0 - 1.0;
                eta[2] = ceta[2];

                int dir[3] = {0, 1, 2};
                m_crvEnt->diff_by_eta(eta, dir[0], dxyz_deta1);
                m_crvEnt->diff_by_eta(eta, dir[1], dxyz_deta2);
                m_crvEnt->diff_by_eta(eta, dir[2], dxyz_deta3);
//                std::cout << "BezTet: dxyz_deta1: (" << dxyz_deta1[0] <<", " << dxyz_deta1[1] << ", " << dxyz_deta1[2] << ")" << std::endl;
//                std::cout << "        dxyz_deta2: (" << dxyz_deta2[0] <<", " << dxyz_deta2[1] << ", " << dxyz_deta2[2] << ")" << std::endl;
//                std::cout << "        dxyz_deta3: (" << dxyz_deta3[0] <<", " << dxyz_deta3[1] << ", " << dxyz_deta3[2] << ")" << std::endl;
                int qindex = i*tmpNumPts1*tmpNumPts + j*tmpNumPts + k;
                deriv[0][0][qindex] = dxyz_deta1[0];
                deriv[0][1][qindex] = dxyz_deta1[1];
                deriv[0][2][qindex] = dxyz_deta1[2];
                deriv[1][0][qindex] = dxyz_deta2[0];
                deriv[1][1][qindex] = dxyz_deta2[1];
                deriv[1][2][qindex] = dxyz_deta2[2];
                deriv[2][0][qindex] = dxyz_deta3[0];
                deriv[2][1][qindex] = dxyz_deta3[1];
                deriv[2][2][qindex] = dxyz_deta3[2];
              }
            }
          }
        }

        if(m_crvEnt->get_exp_dim() == 2) {  // input crv ent is a mesh face
          for(int i = 0; i < tmpNumPts1; ++i) // eta2 dir
          {
            for(int j = 0; j < tmpNumPts; ++j) // eta1 dir             
            {
//              Point2d in_eta = Point2d(keyPts[j], keyPts1[i]);
              Point3d dxyz_deta1, dxyz_deta2;
              //m_beztriPtr->eval_deriv1_by_eta(in_eta, dxyz_deta1, dxyz_deta2);
              // tensor/collapsed coords
              Point2d eta = Point2d(keyPts[j], keyPts1[i]);
              // convert to standard cartesian coords
              Point2d cart = Point2d((0.5 * (1.0 + eta[0]) * (1.0 - eta[1]) - 1.0), 
                                      eta[1]);
              // input
              Point3d in_eta_3d = Point3d(cart[0], cart[1], 0.0);
              int dir[2] = {0, 1};
              m_crvEnt->diff_by_eta(in_eta_3d, dir[0], dxyz_deta1);
              m_crvEnt->diff_by_eta(in_eta_3d, dir[1], dxyz_deta2);
//              std::cout << "BezTri: dxyz_deta1: (" << dxyz_deta1[0] <<", " << dxyz_deta1[1] << ", " << dxyz_deta1[2] << ")" << std::endl;
//              std::cout << "        dxyz_deta2: (" << dxyz_deta2[0] <<", " << dxyz_deta2[1] << ", " << dxyz_deta2[2] << ")" << std::endl;
              deriv[0][0][i*tmpNumPts + j] = dxyz_deta1[0];
              deriv[0][1][i*tmpNumPts + j] = dxyz_deta1[1];
              if (m_coordDim == 3)
                  deriv[0][2][i*tmpNumPts + j] = dxyz_deta1[2];
              deriv[1][0][i*tmpNumPts + j] = dxyz_deta2[0];
              deriv[1][1][i*tmpNumPts + j] = dxyz_deta2[1];
              if (m_coordDim == 3)
                  deriv[1][2][i*tmpNumPts + j] = dxyz_deta2[2];
//            std::cout << deriv[0][0][k] << std::endl;
//            std::cout << deriv[0][1][k] << std::endl;
//            std::cout << deriv[1][0][k] << std::endl;
//            std::cout << deriv[1][1][k] << std::endl;
//            std::cout << std::endl;
            }
          }
        }
/*
           if(m_expDim == 2 && m_coordDim == 2){
              for(int k = 0; k < nqtot_tbasis; ++k){
              std::cout << "deriv[][]["<<k<<"]" << std::endl;
              std::cout << deriv[0][0][k] << std::endl;
              std::cout << deriv[0][1][k] << std::endl;
              std::cout << deriv[1][0][k] << std::endl;
              std::cout << deriv[1][1][k] << std::endl;
              std::cout << std::endl;

              }
            }
*/
            return deriv;
        }

    }; // end of SpatialDomains
}; // end of Nektar
