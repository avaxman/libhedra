// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_DUAL_MESH_H
#define HEDRA_DUAL_MESH_H
#include <igl/igl_inline.h>
#include <igl/setdiff.h>
#include <hedra/polygonal_face_centers.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
    // returns a dual mesh such the the face baycenters become new vertices, and also boundary edges
    // Inputs:
    //  V  eigen double matrix     #V by 3 - vertex coordinates
    //  D  eigen int vector        #F by 1 - face degrees
    //  F  eigen int matrix        #F by max(D) - vertex indices in face
    //  EV eigen int matrix     #E by 2 - map from edges to end vertices
    //  EF eigen int matrix     #E by 2 - map from edges to adjacent faces
    //  innerEdges eigen int vector list of inner edges indices into EV/EF
    //
    // Outputs:
    //  dualV  eigen double matrix  new vertices
    //  dualD  eigen int vector    new valences
    //  dualF eigen int matrix     new faces
    IGL_INLINE bool dual_mesh(const Eigen::MatrixXd& V,
                              const Eigen::VectorXi& D,
                              const Eigen::MatrixXi& F,
                              const Eigen::MatrixXi& EV,
                              const Eigen::MatrixXi& EF,
                              const Eigen::MatrixXi& FE;
                              const Eigen::VectorXi innerEdges,
                              Eigen::MatrixXd& dualV,
                              Eigen::VectorXi& dualD,
                              Eigen::VectorXi& dualF)
    {
        using namespace Eigen;
        MatrixXd faceCenters;
        polygonal_face_centers(V,D,F,faceCenters);
        
        VectorXi allIndices, boundEdges, stub;
        allIndices=VectorXi::LinSpaced(EV.rows(), 0, Ev.rows()-1);
        setdiff(allIndices,innerEdges, boundEdges, stub);
        
        MatrixXd midBoundEdges(boundEdges.rows(),3);
        for (int i=0;i<boundEdges.rows();i++)
            midBoundEdges.row(i)<<(V.row(EV(boundEdges(i),0))+V.row(EV(boundEdges(i),1)))/2.0;
        
        
        //adding the boundEdges vertices as "faces"
        MatrixXi extEF=EF;
        for (int i=0;i<boundEdges.rows();i++)
            extEF(boundEdges(i),1)=F.rows()+i;
        
        dualV.resize(F.rows()+boundEdges.rows(),3);
        dualV.block(0,0,faceCenters.rows(),3)=faceCenters;
        dualV.block(faceCenters.rows(),0,midBoundEdges.rows(),3)=midBoundEdges;
        
        vector<vector<int > > VE;
        for (int i=0;i<EV.rows();i++){
            VE(EV(i,0)).push_back(i));
            VE(EV(i,1)).push_back(i));
        }
        
        dualD.resize(V.rows());
        for (int i=0;i<dualD.rows();i++)
            dualD(i)=VE[i].size();
        
        dualF.resize(dualD.rows(),dualD.maxCoeff()+1);
        
        //creating faces
        for (int i=0;i<VE.rows();i++){
            int currEdge=VE[i][0];
            dualF(i,j)=EF(currEdge,0);
            for (int j=1;j<VE[i].size();i++){
                int nextVertex=(EF(currEdge,0)==F(i,j-1) ? EF(currEdge,1) : EF(currEdge,0));
                dualF(i,j)=EF(currEdge,nextVertex);
                //looking for next edge
                for (int k=0;k<VE[0].size();k++)
                    if (((EF(VE[i][k],0)==nextVertex)||(EF(VE[i][k],1)==nextVertex))&&(VE[0][k]!=currEdge))
                        currEdge=k;
            }
        }
        return true;
    }
}


#endif


