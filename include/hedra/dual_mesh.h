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
#include <hedra/dcel.h>
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
                              const Eigen::MatrixXi& EFi,
                              const Eigen::MatrixXi& FE,
                              const Eigen::VectorXi& innerEdges,
                              Eigen::MatrixXd& dualV,
                              Eigen::VectorXi& dualD,
                              Eigen::MatrixXi& dualF) 
    {
        using namespace Eigen;
        polygonal_face_centers(V,D,F,dualV);
    
        VectorXi allIndices, boundEdges, stub;
        allIndices=VectorXi::LinSpaced(EV.rows(), 0, EV.rows()-1);
        igl::setdiff(allIndices,innerEdges, boundEdges, stub);
        
        VectorXi boundVertexMask=VectorXi::Constant(V.rows(),0);
        for (int i=0;i<boundEdges.rows();i++){
            boundVertexMask(EV(boundEdges(i),0))=1;
            boundVertexMask(EV(boundEdges(i),1))=1;
        }
        /*MatrixXd midBoundEdges(boundEdges.rows(),3);
        for (int i=0;i<boundEdges.rows();i++)
            midBoundEdges.row(i)<<(V.row(EV(boundEdges(i),0))+V.row(EV(boundEdges(i),1)))/2.0;*/
        
        
        Eigen::VectorXi VH;
        Eigen::MatrixXi EH,FH;
        Eigen::VectorXi HV,HE,HF,nextH,prevH,twinH;
        hedra::dcel(D, F, EV, EF, EFi, innerEdges,VH,EH,FH,HV,HE,HF,nextH,prevH,twinH);
        
        cout<<"EF: "<<EF<<endl;
        cout<<"HF: "<<HF<<endl;
        
        //traversing one rings to collect faces
        int currDualFace=0;
        vector<vector<int> > dualFList(V.rows()-boundVertexMask.sum());
        for (int i=0;i<VH.rows();i++){
            if (boundVertexMask(i))
                continue;  //at the moment, not making dual faces for boundary vertices
            int beginH=VH(i);
            int currH=beginH;
            do{
                dualFList[currDualFace].push_back(HF(currH));
                currH=twinH(prevH(currH));
            }while((currH!=beginH)&&(currH!=-1));
            currDualFace++;
        }
        
        dualD.resize(dualFList.size());
        for (int i=0;i<dualD.rows();i++)
            dualD(i)=dualFList[i].size();
        
        dualF=MatrixXi::Constant(dualD.rows(),dualD.maxCoeff()+1,-1);
        
        //putting vertices in random order
        for (int i=0;i<dualFList.size();i++)
            for (int j=0;j<dualFList[i].size();j++)
                dualF(i,j)=dualFList[i][j];
        
        cout<<"dualD: "<<dualD<<endl;
        cout<<"dualF: "<<dualF<<endl;
        return true;
    }
}


#endif


