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
#include "polygonal_face_centers.h"
#include "linear_vi_subdivision.h"
#include "moebius_vi_subdivision.h"
#include "subdivision_basics.h"
#include "dcel.h"
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
                            Eigen::MatrixXi& dualF,
                            const bool clipBoundary)
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
        
        MatrixXd canonFaceVertices;
        canonFaceVertices=sd.original2Canonical(i,origFaceVertices);
        canonFacePoints.row(j)=canonFaceVertices.colwise().mean();
        
        
        //Lifting to candidate points
        candidateFacePoints.block(sd.ringFaces(i,j), 3*currVertexinFace,1,3)=sd.canonical2Original(i,canonFacePoints.row(j));  //not entirely optimized
      }
      
      MatrixXd canonEdgePoints = sd.canonicalEdgePoints(i, canonCenter, canonStarVertices, canonFacePoints);
      
      //Lifting to candidate points
      MatrixXd localCandidateEdgePoints;
      localCandidateEdgePoints=sd.canonical2Original(i,canonEdgePoints);
      for (int j=0;j<sd.vertexValences(i);j++){
        int onEdge=sd.HE(sd.starHalfedges(i,j));
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
    
    //addendum: killing dual faces (all faces that are a boundary)
    if (clipBoundary){
      MatrixXi EV, FE, EF, EFi;
      MatrixXd FEs;
      VectorXi innerEdges;
      hedra::polygonal_edge_topology(dualD, dualF, EV, FE, EF, EFi, FEs, innerEdges);
      VectorXi isBoundaryFace=VectorXi::Zero(dualF.rows());
      for (int i=0;i<EF.rows();i++){
        if (EF(i,0)==-1)
          isBoundaryFace(EF(i,1))=1;
        if (EF(i,1)==-1)
          isBoundaryFace(EF(i,0))=1;
      }
      cout<<"Clipping "<<isBoundaryFace.sum()<<" boundary faces"<<endl;
      MatrixXi newDualF(dualF.rows()-isBoundaryFace.sum(),dualF.cols());
      VectorXi newDualD(dualF.rows()-isBoundaryFace.sum());
      int fCount=0;
      for (int i=0;i<dualF.rows();i++){
        if (!isBoundaryFace(i)){
          newDualD(fCount)=dualD(i);
          newDualF.row(fCount++)=dualF.row(i);
        }
      }
      
      dualF=newDualF;
      dualD=newDualD;
      //removing isolated vertices
      VectorXi isUnrefVertex=VectorXi::Constant(dualV.rows(),1);
      for (int i=0;i<dualD.rows();i++)
        for (int j=0;j<dualD(i);j++)
          isUnrefVertex(dualF(i,j))=0;
      
      int vCount=0;
      VectorXi transVertices=VectorXi::Constant(dualV.rows(),-1);
      for (int i=0;i<dualV.rows();i++)
        if (!isUnrefVertex(i))
          transVertices(i)=vCount++;
      
      for (int i=0;i<dualD.rows();i++)
        for (int j=0;j<dualD(i);j++)
          dualF(i,j)=transVertices(dualF(i,j));
      
      MatrixXd newDualV(vCount, 3);
      for (int i=0;i<dualV.rows();i++)
        if (!isUnrefVertex(i))
          newDualV.row(transVertices(i))=dualV.row(i);
      
      dualV = newDualV;
      
    }
    
    
    return true;
  }
  
  
  
  //user version
  IGL_INLINE bool dual_mesh(const Eigen::MatrixXd& V,
                            const Eigen::VectorXi& D,
                            const Eigen::MatrixXi& F,
                            const int& st,
                            Eigen::MatrixXd& fineV,
                            Eigen::VectorXi& fineD,
                            Eigen::MatrixXi& fineF,
                            const bool clipBoundary)
  
  {
    using namespace Eigen;
    using namespace std;
    
    switch (st){
      case hedra::LINEAR_SUBDIVISION: {
        hedra::LinearVISubdivisionData lsd;
        dual_mesh(V, D,F,lsd, fineV, fineD, fineF,clipBoundary);
        break;
      }
      case hedra::CANONICAL_MOEBIUS_SUBDIVISION: {
        hedra::MoebiusVISubdivisionData msd;
        dual_mesh(V, D,F,msd, fineV, fineD, fineF,clipBoundary);
        break;
      }
      default: return false;
    }
    
    return true;
  }
}


#endif


