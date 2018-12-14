// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_DUAL_TRUNCATION_H
#define HEDRA_DUAL_TRUNCATION_H
#include <igl/igl_inline.h>
#include <hedra/polygonal_face_centers.h>
#include <hedra/linear_vi_subdivision.h>
#include <hedra/moebius_vi_subdivision.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
  // returns a mesh after dual truncation, which is basically vertex insertion in the barycenter of each face, connected with all other vertices
  // Inputs:
  //  V  eigen double matrix     #V by 3 - vertex coordinates
  //  D  eigen int vector        #F by 1 - face degrees
  //  F  eigen int matrix        #F by max(D) - vertex indices in face
  
  // Outputs:
  //  newV  eigen double matrix  new vertices
  //  newD  eigen int vector    new valences
  //  newF eigen int matrix     new faces
  
  
  IGL_INLINE bool dual_truncation(const Eigen::MatrixXd& V,
                                  const Eigen::VectorXi& D,
                                  const Eigen::MatrixXi& F,
                                  OneRingSubdivisionData& sd,
                                  Eigen::MatrixXd& fineV,
                                  Eigen::VectorXi& fineD,
                                  Eigen::MatrixXi& fineF)
  {
    
    
    using namespace Eigen;
    
    sd.setup(V,D,F);
    
    Eigen::MatrixXd fineFacePoints(F.rows(),3);
    
    MatrixXd candidateFacePoints(F.rows(), D.maxCoeff()*3);
    
    //canonical embedding candidate points
    for (int i=0;i<V.rows();i++){
      MatrixXd origStarVertices(sd.vertexValences(i),3);
      for (int j=0;j<sd.vertexValences(i);j++)
        origStarVertices.row(j)=V.row(sd.starVertices(i,j));
      
      MatrixXd canonStarVertices;
      canonStarVertices=sd.original2Canonical(i,origStarVertices);
      
      RowVector3d canonCenter;
      canonCenter=sd.original2Canonical(i,V.row(i));
      
      
      //face points
      int numRingFaces=sd.vertexValences(i)-sd.isBoundaryVertex(i);  //one less face for boundary vertices
      MatrixXd canonFacePoints(numRingFaces,3);
      for (int j=0;j<numRingFaces;j++){
        MatrixXd origFaceVertices(D(sd.ringFaces(i,j)),3);
        int currVertexinFace=-1;
        for (int k=0;k<D(sd.ringFaces(i,j));k++){
          origFaceVertices.row(k)=V.row(F(sd.ringFaces(i,j),k));
          if (F(sd.ringFaces(i,j),k)==i)
            currVertexinFace=k;
        }
        
        //std::cout<<"origFaceVertices: "<<origFaceVertices<<std::endl;
        MatrixXd canonFaceVertices;
        canonFaceVertices=sd.original2Canonical(i,origFaceVertices);
        //std::cout<<"canonFaceVertices: "<<canonFaceVertices<<std::endl;
        canonFacePoints.row(j)=canonFaceVertices.colwise().mean();
        //std::cout<<"canonFacePoints.row(j): "<<canonFacePoints.row(j)<<std::endl;
        
        //std::cout<<"sd.canonical2Original(i,canonFacePoints.row(j)): "<<sd.canonical2Original(i,canonFacePoints.row(j))<<std::endl;
        
        //Lifting to candidate points
        candidateFacePoints.block(sd.ringFaces(i,j), 3*currVertexinFace,1,3)=sd.canonical2Original(i,canonFacePoints.row(j));  //not entirely optimized
      }
      
    }
    
    //Blending face points from candidates
    fineFacePoints = sd.facePointBlend(candidateFacePoints);
    
    fineV.conservativeResize(V.rows()+fineFacePoints.rows(),3);
    fineV<<V, fineFacePoints;
    
    fineD=VectorXi::Constant(D.sum(),3);
    fineF.conservativeResize(D.sum(),3);
    
    int currNewFace=0;
    for (int i=0;i<D.rows();i++)
      for (int j=0;j<D(i);j++)
        fineF.row(currNewFace++)<<F(i,j), F(i,(j+1)%D(i)), V.rows()+i;
    
    return true;
  }
  
  
  
  IGL_INLINE bool dual_truncation(const Eigen::MatrixXd& V,
                                  const Eigen::VectorXi& D,
                                  const Eigen::MatrixXi& F,
                                  const int& st,
                                  Eigen::MatrixXd& fineV,
                                  Eigen::VectorXi& fineD,
                                  Eigen::MatrixXi& fineF)
  {
    using namespace Eigen;
    using namespace std;
    
    switch (st){
      case hedra::LINEAR_SUBDIVISION: {
        hedra::LinearVISubdivisionData lsd;
        dual_truncation(V, D,F,lsd, fineV, fineD, fineF);
        break;
      }
      case hedra::CANONICAL_MOEBIUS_SUBDIVISION: {
        hedra::MoebiusVISubdivisionData msd;
        dual_truncation(V, D,F,msd, fineV, fineD, fineF);
        break;
      }
      default: return false;
    }
    
    return true;
  }
}


#endif


