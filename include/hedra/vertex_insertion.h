// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_VERTEX_INSERTION_H
#define HEDRA_VERTEX_INSERTION_H
#include <igl/igl_inline.h>
#include <hedra/polygonal_face_centers.h>
#include <hedra/subdivision_basics.h>
#include <hedra/linear_vi_subdivision.h>
#include <hedra/moebius_vi_subdivision.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
  // returns a mesh after vertex insertion, which is basically vertex insertion in the barycenter of each face, connected with all midedges
  // Inputs:
  //  V  eigen double matrix     #V by 3 - vertex coordinates
  //  D  eigen int vector        #F by 1 - face degrees
  //  F  eigen int matrix        #F by max(D) - vertex indices in face
  //  FE eign int matrix         #F by max(D) - edges by order in face
  
  // Outputs:
  //  newV  eigen double matrix  new vertices
  //  newD  eigen int vector    new valences
  //  newF eigen int matrix     new faces
  IGL_INLINE bool vertex_insertion(const Eigen::MatrixXd& V,
                                   const Eigen::VectorXi& D,
                                   const Eigen::MatrixXi& F,
                                   OneRingSubdivisionData& sd,
                                   Eigen::MatrixXd& fineV,
                                   Eigen::VectorXi& fineD,
                                   Eigen::MatrixXi& fineF)
  {
    
    
    using namespace Eigen;
    
    sd.setup(V,D,F);
    
    Eigen::MatrixXd fineVertexPoints(V.rows(),3);
    Eigen::MatrixXd fineFacePoints(F.rows(),3);
    Eigen::MatrixXd fineEdgePoints(sd.EV.rows(),3);
    
    
    MatrixXd candidateFacePoints(F.rows(), D.maxCoeff()*3);
    MatrixXd candidateEdgePoints(sd.EV.rows(), 6);
    
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
      
      //std::cout<<"candidateFacePoints: "<<candidateFacePoints<<std::endl;
      
      //candidate edge points
      //MatrixXd canonEdgePoints(sd.vertexValences(i),3);
      /*for (int j=sd.isBoundaryVertex(i);j<sd.vertexValences(i)-sd.isBoundaryVertex(i);j++)
        canonEdgePoints.row(j)=(canonCenter+canonStarVertices.row(j)+canonFacePoints.row(j)+canonFacePoints.row((j+sd.vertexValences(i)-1)%sd.vertexValences(i)))/4.0;*/
      
      MatrixXd canonEdgePoints = sd.canonicalEdgePoints(i, canonCenter, canonStarVertices, canonFacePoints);
      
      /*if (isBoundaryVertex(i)){
       canonEdgePoints.row(0)=(canonCenter+canonStarVertices.row(0))/2.0;
       canonEdgePoints.row(vertexValences(i)-1)=(canonCenter+canonStarVertices.row(vertexValences(i)-1))/2.0;
       }*/
      
      //Lifting to candidate points
      MatrixXd localCandidateEdgePoints;
      localCandidateEdgePoints=sd.canonical2Original(i,canonEdgePoints);
      for (int j=sd.isBoundaryVertex(i);j<sd.vertexValences(i)-sd.isBoundaryVertex(i);j++){
        int onEdge=sd.HE(sd.starHalfedges(i,j));
        // std::cout<<"sd.starHalfedges(i,j): "<<sd.starHalfedges(i,j)<<std::endl;
        // std::cout<<"onEdge: "<<onEdge<<std::endl;
        int inEdge = (sd.EV(onEdge,0)==i ? 0 : 1);
        candidateEdgePoints.block(onEdge,3*inEdge,1,3)=localCandidateEdgePoints.row(j);  //WRONG!!!
      }
      
      //vertex points
      RowVector3d canonFineCenter;
      if (!sd.isBoundaryVertex(i)){
        canonFineCenter = sd.innerVertexCanonicalBlend(canonCenter,canonEdgePoints,canonFacePoints);
       
      } //boundary will be assigned later
      
      fineVertexPoints.row(i)=sd.canonical2Original(i,canonFineCenter);
    }
    
    //Blending face points from candidates
    fineFacePoints = sd.facePointBlend(candidateFacePoints);


    //Blending edge points from candidates, and boundary edge points from boundary curves.
    
    Eigen::MatrixXd a(sd.EH.rows(),3), b(sd.EH.rows(),3), c(sd.EH.rows(),3), d(sd.EH.rows(),3);
    
    
    for (int i=0;i<sd.EH.rows();i++){
      int currH;
      if (sd.EH(i,0)==-1)
        currH=sd.EH(i,1);
      else currH=sd.EH(i,0);
      int ix2=sd.HV(currH);
      int ix3=sd.HV(sd.nextH(currH));
      
      Eigen::RowVector3d x2; x2<<V.row(ix2);
      Eigen::RowVector3d x3; x3<<V.row(ix3);
      
      if (sd.twinH(currH)!=-1){  //inner edge
        a.row(i)=x2;
        d.row(i)=x3;
        if (ix2==sd.EV(i,0)){
          b.row(i)=candidateEdgePoints.block(i,0,1,3);
          c.row(i)=candidateEdgePoints.block(i,3,1,3);
        } else {
          c.row(i)=candidateEdgePoints.block(i,0,1,3);
          b.row(i)=candidateEdgePoints.block(i,3,1,3);
        }
        continue;
      }
      
      bool x1Ear=(sd.twinH(sd.prevH(currH))==-1);
      bool x4Ear=(sd.twinH(sd.nextH(currH))==-1);
      Eigen::RowVector3d x1, x4;
      
      if (!x1Ear){
        int otherBoundH=currH;
        do{
          otherBoundH=sd.twinH(sd.prevH(otherBoundH));
        }while(sd.twinH(sd.prevH(otherBoundH))!=-1);
        int ix1=sd.HV(sd.prevH(otherBoundH));
        x1=V.row(ix1);
      }
      
      if (!x4Ear){
        int otherBoundH=currH;
        do{
          otherBoundH=sd.twinH(sd.nextH(otherBoundH));
        }while(sd.twinH(sd.nextH(otherBoundH))!=-1);
        int ix4=sd.HV(sd.nextH(sd.nextH(otherBoundH)));
        x4=V.row(ix4);
      }
      
      if (x1Ear)
        x1=sd.threePointsExtrapolation(x4, x3, x2);
      if (x4Ear)
        x4=sd.threePointsExtrapolation(x1, x2, x3);
      
      a.row(i)=x1;
      b.row(i)=x2;
      c.row(i)=x3;
      d.row(i)=x4;
    }
    
    //std::cout<<"a,b,c,d: "<<a<<b<<c<<d<<std::endl;
    fineEdgePoints=sd.fourPointsInterpolation(a,b,c,d);
    //std::cout<<"candidateEdgePoints: "<<candidateEdgePoints<<std::endl;
    //std::cout<<"fineEdgePoints: "<<fineEdgePoints<<std::endl;
    
    
    //Blending vertex boundary points
    for (int i=0;i<sd.EH.rows();i++){
      int currH;
      if (sd.EH(i,0)==-1)
        currH=sd.EH(i,1);
      else currH=sd.EH(i,0);
      
      if (sd.twinH(currH)!=-1)
        continue;  //not a boundary
      
      int ix2=sd.HV(currH);
      int ix3=sd.HV(sd.nextH(currH));
      RowVector3d x2; x2<<V.row(ix2);
      RowVector3d x3; x3<<V.row(ix3);
      bool x1Ear=(sd.twinH(sd.prevH(currH))==-1);
      RowVector3d x1;
      if (!x1Ear){
        int otherBoundH=currH;
        do{
          otherBoundH=sd.twinH(sd.prevH(otherBoundH));
        }while(sd.twinH(sd.prevH(otherBoundH))!=-1);
        int ix1=sd.HV(sd.prevH(otherBoundH));
        x1=V.row(ix1);
        
        fineVertexPoints.row(ix2)=sd.boundaryVertexPoint(x1,fineEdgePoints.row(sd.HE(sd.prevH(otherBoundH))),x2, fineEdgePoints.row(i), x3);
      } else {
        fineVertexPoints.row(ix2) = V.row(ix2);
      }
    }
    
    std::cout<<"(fineVertexPoints-V).lpNorm<Infinity>()" <<(fineVertexPoints-V).lpNorm<Infinity>()<<std::endl;
    
    fineV.conservativeResize(fineVertexPoints.rows()+fineFacePoints.rows()+fineEdgePoints.rows(),3);
    fineV<<fineVertexPoints, fineEdgePoints, fineFacePoints;
    int numNewFaces=D.sum();
    fineD=VectorXi::Constant(numNewFaces, 4);
    fineF.conservativeResize(numNewFaces,4);
    int currNewFace=0;
    for (int i=0;i<D.rows();i++){
      for (int j=0;j<D(i);j++){
        fineF.row(currNewFace++)<<F(i,j),
        V.rows()+sd.FE(i,j),
        V.rows()+fineEdgePoints.rows()+i,
        V.rows()+sd.FE(i,(j+D(i)-1)%D(i));
      }
    }
    
    return true;
  }
  
  
  //user version
  IGL_INLINE bool vertex_insertion(const Eigen::MatrixXd& V,
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
        vertex_insertion(V, D,F,lsd, fineV, fineD, fineF);
        break;
      }
      case hedra::CANONICAL_MOEBIUS_SUBDIVISION: {
        hedra::MoebiusVISubdivisionData msd;
        vertex_insertion(V, D,F,msd, fineV, fineD, fineF);
        break;
      }
      default: return false;
    }
    
    return true;
  }
  
  
}


#endif


