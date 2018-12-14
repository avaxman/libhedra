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
#include <hedra/linear_vi_subdivision.h>
#include <hedra/moebius_vi_subdivision.h>
#include <hedra/subdivision_basics.h>
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
                            OneRingSubdivisionData& sd,
                            Eigen::MatrixXd& dualV,
                            Eigen::VectorXi& dualD,
                            Eigen::MatrixXi& dualF)
  {
    using namespace Eigen;
    using namespace std;
    
    sd.setup(V,D,F);
    
    Eigen::MatrixXd dualFacePoints(F.rows(),3);
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
      
      MatrixXd canonEdgePoints = sd.canonicalEdgePoints(i, canonCenter, canonStarVertices, canonFacePoints);
      
      /*if (isBoundaryVertex(i)){
       canonEdgePoints.row(0)=(canonCenter+canonStarVertices.row(0))/2.0;
       canonEdgePoints.row(vertexValences(i)-1)=(canonCenter+canonStarVertices.row(vertexValences(i)-1))/2.0;
       }*/
      
      //Lifting to candidate points
      MatrixXd localCandidateEdgePoints;
      localCandidateEdgePoints=sd.canonical2Original(i,canonEdgePoints);
      for (int j=0;j<sd.vertexValences(i);j++){
        int onEdge=sd.HE(sd.starHalfedges(i,j));
        // std::cout<<"sd.starHalfedges(i,j): "<<sd.starHalfedges(i,j)<<std::endl;
        // std::cout<<"onEdge: "<<onEdge<<std::endl;
        int inEdge = (sd.EV(onEdge,0)==i ? 0 : 1);
        candidateEdgePoints.block(onEdge,3*inEdge,1,3)=localCandidateEdgePoints.row(j);  //WRONG!!!
      }
    }
    
    //Blending face points from candidates
    dualFacePoints = sd.facePointBlend(candidateFacePoints);
    
    //cellecting only boundary points
    vector<RowVector3d> boundaryEdgePointList;
    VectorXi boundaryVertexIndices(sd.EV.rows());
    int boundaryCounter=0;
    for (int i=0;i<sd.EV.rows();i++)
      if ((sd.isBoundaryVertex[sd.EV(i,0)])&&(sd.isBoundaryVertex[sd.EV(i,1)])){
        MatrixXd stub =sd.fourPointsInterpolation(V.row(sd.EV(i,0)),candidateEdgePoints.block(i,0,1,3), candidateEdgePoints.block(i,3,1,3), V.row(sd.EV(i,1)));  //WASTEFUL
        boundaryEdgePointList.push_back(stub.row(0));
        boundaryVertexIndices(i)=boundaryCounter++;
      }
    
   
    //finding ears
    VectorXi vertices2Ears(V.rows()); vertices2Ears.setConstant(-1.0);
    vector<int> earsList;
    for (int i=0;i<V.rows();i++){
      if (!sd.isBoundaryVertex[i])
        continue;
      
      if (sd.vertexValences(i)<3){
        earsList.push_back(i);
        vertices2Ears(i)=earsList.size()-1;
      }
    }
    
    
    VectorXi allIndices, boundEdges, stub;
    allIndices=VectorXi::LinSpaced(sd.EV.rows(), 0, sd.EV.rows()-1);
    igl::setdiff(allIndices,sd.innerEdges, boundEdges, stub);
    
    VectorXi boundVertexMask=VectorXi::Constant(V.rows(),0);
    for (int i=0;i<boundEdges.rows();i++){
      boundVertexMask(sd.EV(boundEdges(i),0))=1;
      boundVertexMask(sd.EV(boundEdges(i),1))=1;
    }
    
    
    //constructing the topology
    dualV.conservativeResize(dualFacePoints.rows()+boundaryEdgePointList.size()+earsList.size(),3);
    dualV.block(0,0,dualFacePoints.rows(),3)=dualFacePoints;
    for (int i=0;i<boundaryEdgePointList.size();i++)
      dualV.row(dualFacePoints.rows()+i)=boundaryEdgePointList[i];
    
    for (int i=0;i<earsList.size();i++)
      dualV.row(dualFacePoints.rows()+boundaryEdgePointList.size()+i)=V.row(earsList[i]);
    
    vector<int> emptyIntVector;
    //traversing one rings to collect faces
    int currDualFace=0;
    vector<vector<int> > dualFList;
    for (int i=0;i<sd.VH.rows();i++){
      
      dualFList.push_back(emptyIntVector);
      int beginH=sd.VH(i);
      int currH=beginH;
      
      //reseting to first boundary edge (or something arbitrary otherwise)
      while ((sd.twinH(currH)!=-1)){
        currH=sd.nextH(sd.twinH(currH));
        if (currH==beginH) break;
      }
      
      beginH=currH;
      
      do{
        if (sd.twinH(currH)==-1)
          dualFList[currDualFace].push_back(F.rows()+boundaryVertexIndices(sd.HE(currH)));
        
        dualFList[currDualFace].push_back(sd.HF(currH));
        if (sd.twinH(sd.prevH(currH))==-1)  //before last edge
          dualFList[currDualFace].push_back(F.rows()+boundaryVertexIndices(sd.HE(sd.prevH(currH))));
        if (sd.twinH(sd.prevH(beginH))==-1)  //an ear
          dualFList[currDualFace].push_back(F.rows()+boundaryEdgePointList.size()+vertices2Ears(i));
        currH=sd.twinH(sd.prevH(currH));
      }while((currH!=beginH)&&(currH!=-1));
      currDualFace++;
    }
    
    dualD.conservativeResize(dualFList.size());
    for (int i=0;i<dualD.rows();i++)
      dualD(i)=dualFList[i].size();
    
    dualF=MatrixXi::Constant(dualD.rows(),dualD.maxCoeff()+1,-1);
    
    //putting vertices in random order
    for (int i=0;i<dualFList.size();i++)
      for (int j=0;j<dualFList[i].size();j++)
        dualF(i,j)=dualFList[i][j];
    
    return true;
  }
  
  
  
  //user version
  IGL_INLINE bool dual_mesh(const Eigen::MatrixXd& V,
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
        dual_mesh(V, D,F,lsd, fineV, fineD, fineF);
        break;
      }
      case hedra::CANONICAL_MOEBIUS_SUBDIVISION: {
        hedra::MoebiusVISubdivisionData msd;
        dual_mesh(V, D,F,msd, fineV, fineD, fineF);
        break;
      }
      default: return false;
    }
    
    return true;
  }
}


#endif


