// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2015 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_EDGE_TOPOLOGY_H
#define HEDRA_EDGE_TOPOLOGY_H

#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <vector>

namespace hedra
{
    // Initialize Edges and their topological relations
    
    //input:
    //  D  eigen int vector     #F by 1 - face degrees
    //  F  eigen int matrix     #F by max(D) - vertex indices in face
  
    // Output:
    // EV   #E by 2, Stores the edge description as pair of indices to vertices
    // FE : #F by max(D), Stores the Face-Edge relation
    // EF : #E by 2: Stores the Edge-Face relation

    template <typename DerivedF>
    IGL_INLINE void hedra_edge_topology(const Eigen::PlainObjectBase<DerivedF>& D,
                                          const Eigen::PlainObjectBase<DerivedF>& F,
                                          Eigen::MatrixXi& EV,
                                          Eigen::MatrixXi& gFE,
                                          Eigen::MatrixXi& gEF)
    {
        // Only needs to be edge-manifold
        std::vector<std::vector<int> > ETT;
        for(int f=0;f<gD.rows();++f)
            for (int i=0;i<gD(f);++i)
            {
                // v1 v2 f vi
                int v1 = gF(f,i);
                int v2 = gF(f,(i+1)%gD(f));
                if (v1 > v2) std::swap(v1,v2);
                std::vector<int> r(4);
                r[0] = v1; r[1] = v2;
                r[2] = f;  r[3] = i;
                ETT.push_back(r);
            }
        std::sort(ETT.begin(),ETT.end());
        
        // count the number of edges (assume manifoldness)
        int En = 1; // the last is always counted
        for(unsigned i=0;i<ETT.size()-1;++i)
            if (!((ETT[i][0] == ETT[i+1][0]) && (ETT[i][1] == ETT[i+1][1])))
                ++En;
        
        gEV = Eigen::MatrixXi::Constant((int)(En),2,-1);
        gFE = Eigen::MatrixXi::Constant((int)(gF.rows()),(int)(gF.cols()),-1);
        gEF = Eigen::MatrixXi::Constant((int)(En),2,-1);
        En = 0;
        
        for(unsigned i=0;i<ETT.size();++i)
        {
            if (i == ETT.size()-1 ||
                !((ETT[i][0] == ETT[i+1][0]) && (ETT[i][1] == ETT[i+1][1]))
                )
            {
                // Border edge
                std::vector<int>& r1 = ETT[i];
                gEV(En,0)     = r1[0];
                gEV(En,1)     = r1[1];
                gEF(En,0)    = r1[2];
                gFE(r1[2],r1[3]) = En;
            }
            else
            {
                std::vector<int>& r1 = ETT[i];
                std::vector<int>& r2 = ETT[i+1];
                gEV(En,0)     = r1[0];
                gEV(En,1)     = r1[1];
                gEF(En,0)    = r1[2];
                gEF(En,1)    = r2[2];
                gFE(r1[2],r1[3]) = En;
                gFE(r2[2],r2[3]) = En;
                ++i; // skip the next one
            }
            ++En;
        }
        
        // Sort the relation EF, accordingly to EV
        // the first one is the face on the left of the edge
        
        for(unsigned i=0; i<gEF.rows(); ++i)
        {
            int fid = gEF(i,0);
            bool flip = true;
            // search for edge EV.row(i)
            for (unsigned j=0; j<gD(fid); ++j)
            {
                if ((gF(fid,j) == gEV(i,0)) && (gF(fid,(j+1)%gD(fid)) == gEV(i,1)))
                    flip = false;
            }
            
            if (flip)
            {
                int tmp = gEF(i,0);
                gEF(i,0) = gEF(i,1);
                gEF(i,1) = tmp;
            }
        }

    }
}


#endif
