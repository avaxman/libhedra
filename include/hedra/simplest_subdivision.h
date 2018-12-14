// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_SIMPLEST_SUBDIVISION_H
#define HEDRA_SIMPLEST_SUBDIVISION_H
#include <igl/igl_inline.h>
#include <hedra/polygonal_face_centers.h>
#include <hedra/subdivision_basics.h>
#include <hedra/moebius_simplest_subdivision.h>
#include <hedra/linear_simplest_subdivision.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
  
  IGL_INLINE bool simplest_subdivision(const Eigen::MatrixXd& V,
                                       const Eigen::VectorXi& D,
                                       const Eigen::MatrixXi& F,
                                       OneRingSubdivisionData& sd,
                                       Eigen::MatrixXd& fineV,
                                       Eigen::VectorXi& fineD,
                                       Eigen::MatrixXi& fineF)
  {
    
    using namespace Eigen;
    using namespace std;
    sd.setup(V,D,F);
    
    Eigen::MatrixXd fineEdgePoints(sd.EV.rows(),3);
    Eigen::MatrixXd candidateEdgePoints(sd.EV.rows(), 6);
    
    //canonical embedding candidate points
    for (int i=0;i<V.rows();i++){
      MatrixXd origStarVertices(sd.vertexValences(i),3);
      for (int j=0;j<sd.vertexValences(i);j++)
        origStarVertices.row(j)=V.row(sd.starVertices(i,j));
      
      MatrixXd canonStarVertices;
      canonStarVertices=sd.original2Canonical(i,origStarVertices);
      
      RowVector3d canonCenter;
      canonCenter=sd.original2Canonical(i,V.row(i));
      
      //candidate edge points
      MatrixXd canonEdgePoints(sd.vertexValences(i),3);
      for (int j=0;j<sd.vertexValences(i);j++)
        canonEdgePoints.row(j)=(canonCenter+canonStarVertices.row(j))/2.0;
      
      //Lifting to candidate points
      MatrixXd localCandidateEdgePoints;
      localCandidateEdgePoints=sd.canonical2Original(i,canonEdgePoints);
      for (int j=0;j<sd.vertexValences(i);j++){
        int onEdge=sd.HE(sd.starHalfedges(i,j));
        int inEdge = (sd.EV(onEdge,0)==i ? 0 : 1);
        candidateEdgePoints.block(onEdge,3*inEdge,1,3)=localCandidateEdgePoints.row(j);
      }
    }
    
    
    //blending edge points
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
      
      a.row(i)=x2;
      d.row(i)=x3;
      if (ix2==sd.EV(i,0)){
        b.row(i)=candidateEdgePoints.block(i,0,1,3);
        c.row(i)=candidateEdgePoints.block(i,3,1,3);
      } else {
        c.row(i)=candidateEdgePoints.block(i,0,1,3);
        b.row(i)=candidateEdgePoints.block(i,3,1,3);
      }
    }
    
    fineEdgePoints=sd.fourPointsInterpolation(a,b,c,d);
    
    //constructing the topology
    fineV=fineEdgePoints;
  
    vector<int> emptyIntVector;
    vector< vector<int> > newFList;
    
    //new faces from old faces
    for (int i=0;i<D.rows();i++){
      newFList.push_back(emptyIntVector);
      for (int j=0;j<D(i);j++)
        newFList[i].push_back(sd.FE(i,j));
    }
    
    for (int i=0;i<sd.VH.rows();i++){
      if (sd.isBoundaryVertex[i])
        continue;  //at the moment, not making faces for boundary vertices
      int beginH=sd.VH(i);
      int currH=beginH;
      newFList.push_back(emptyIntVector);
      do{
        newFList[newFList.size()-1].push_back(sd.HE(currH));
        currH=sd.twinH(sd.prevH(currH));
      }while((currH!=beginH)&&(currH!=-1));
    }
    
    fineD.conservativeResize(newFList.size());
    for (int i=0;i<fineD.rows();i++)
      fineD(i)=newFList[i].size();
    
    fineF=MatrixXi::Constant(fineD.rows(),fineD.maxCoeff()+1,-1);
    
    //putting vertices in random order
    for (int i=0;i<newFList.size();i++)
      for (int j=0;j<newFList[i].size();j++)
        fineF(i,j)=newFList[i][j];
    
    return true;
  }
  
  
  IGL_INLINE bool simplest_subdivision(const Eigen::MatrixXd& V,
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
        hedra::LinearSimplestSubdivisionData lsd;
        simplest_subdivision(V, D,F,lsd, fineV, fineD, fineF);
        break;
      }
      case hedra::CANONICAL_MOEBIUS_SUBDIVISION: {
        hedra::MoebiusSimplestSubdivisionData msd;
        simplest_subdivision(V, D,F,msd, fineV, fineD, fineF);
        break;
      }
      default: return false;
    }
    
    return true;
  }
}


#endif


