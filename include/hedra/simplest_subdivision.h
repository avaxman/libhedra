// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_SIMPLEST_SUBDIVISION_H
#define HEDRA_SIMPLEST_SUBDIVISION_H
#include <igl/igl_inline.h>
#include <hedra/polygonal_face_centers.h>
#include <hedra/dual_mesh.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
    // returns a mesh after vertex insertion, which creates a mesh by connected midedges inside each face and around each vertex
    // Inputs:
    //  V  eigen double matrix     #V by 3 - vertex coordinates
    //  D  eigen int vector        #F by 1 - face degrees
    //  F  eigen int matrix        #F by max(D) - vertex indices in face
    //  FE eign int matrix         #F by max(D) - edges by order in face
    
    // Outputs:
    //  newV  eigen double matrix  new vertices
    //  newD  eigen int vector    new valences
    //  newF eigen int matrix     new faces
    IGL_INLINE bool simplest_subdivision(const Eigen::MatrixXd& V,
                                         const Eigen::VectorXi& D,
                                         const Eigen::MatrixXi& F,
                                         const Eigen::MatrixXi& EV,
                                         const Eigen::MatrixXi& EF,
                                         const Eigen::MatrixXi& EFi,
                                         const Eigen::MatrixXi& FE,
                                         const Eigen::VectorXi& innerEdges,
                                         Eigen::MatrixXd& newV,
                                         Eigen::VectorXi& newD,
                                         Eigen::MatrixXi& newF)
    {
        using namespace Eigen;
        
        newV.resize(EV.rows(),3);
        for (int i=0;i<EV.rows();i++)
            newV.row(i)<<(V.row(EV(i,0))+V.row(EV(i,1)))/2.0;
        
        VectorXi allIndices, boundEdges, stub;
        allIndices=VectorXi::LinSpaced(EV.rows(), 0, EV.rows()-1);
        igl::setdiff(allIndices,innerEdges, boundEdges, stub);
        
        VectorXi boundVertexMask=VectorXi::Constant(V.rows(),0);
        for (int i=0;i<boundEdges.rows();i++){
            boundVertexMask(EV(boundEdges(i),0))=1;
            boundVertexMask(EV(boundEdges(i),1))=1;
        }
        
        Eigen::VectorXi VH;
        Eigen::MatrixXi EH,FH;
        Eigen::VectorXi HV,HE,HF,nextH,prevH,twinH;
        hedra::dcel(D, F, EV, EF, EFi, innerEdges,VH,EH,FH,HV,HE,HF,nextH,prevH,twinH);
        
        vector< vector<int> > newFList(V.rows()+FE.rows()-boundVertexMask.sum());
        for (int i=0;i<D.rows();i++)
            for (int j=0;j<D(i);j++)
                newFList[i].push_back(FE(i,j));
        
        
        //traversing one rings to collect faces
        int currNewFace=FE.rows();
        //cout<<"VH: "<<VH<<endl;
        for (int i=0;i<VH.rows();i++){
            if (boundVertexMask(i))
                continue;  //at the moment, not making faces for boundary vertices
            int beginH=VH(i);
            int currH=beginH;
            do{
                newFList[currNewFace].push_back(HE(currH));
                //cout<<"HE(currH): "<<HE(currH)<<endl;
                currH=twinH(prevH(currH));
            }while((currH!=beginH)&&(currH!=-1));
            for (int j=0;j<newFList[currNewFace].size();j++)
                cout<<newFList[currNewFace][j]<<endl;
            currNewFace++;
        }
        
        newD.resize(newFList.size());
        for (int i=0;i<newFList.size();i++)
            newD(i)=newFList[i].size();
        
        //cout<<"newD: "<<newD<<endl;
        
        newF.resize(newFList.size(), newD.maxCoeff());
        for (int i=0;i<newD.rows();i++)
            for (int j=0;j<newD(i);j++)
                newF(i,j)=newFList[i][j];
        
        //cout<<"newF: "<<newF<<endl;
        

        return true;
    }
}


#endif


