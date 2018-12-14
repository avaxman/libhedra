// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_SUBDIVISION_BASICS_H
#define HEDRA_SUBDIVISION_BASICS_H
#include <igl/igl_inline.h>
#include <hedra/vertex_valences.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
  
  const int LINEAR_SUBDIVISION=0;
  const int CANONICAL_MOEBIUS_SUBDIVISION=1;
  
  IGL_INLINE bool vertex_stars(const Eigen::MatrixXi& EV,
                               const Eigen::VectorXi& VH,
                               const Eigen::MatrixXi& EH,
                               const Eigen::MatrixXi& FH,
                               const Eigen::VectorXi& HV,
                               const Eigen::VectorXi& HE,
                               const Eigen::VectorXi& HF,
                               const Eigen::VectorXi& nextH,
                               const Eigen::VectorXi& prevH,
                               const Eigen::VectorXi& twinH,
                               Eigen::VectorXi& vertexValences,
                               Eigen::MatrixXi& starVertices,
                               Eigen::MatrixXi& starHalfedges,
                               Eigen::MatrixXi& ringFaces,
                               Eigen::VectorXi& isBoundaryVertex)
  {
    
    using namespace Eigen;
    hedra::vertex_valences(EV,vertexValences);
    
    int maxValence =vertexValences.maxCoeff();
    starVertices.conservativeResize(vertexValences.rows(),maxValence);
    starHalfedges.conservativeResize(vertexValences.rows(),maxValence);
    ringFaces.conservativeResize(vertexValences.rows(),maxValence);
    isBoundaryVertex.conservativeResize(vertexValences.rows());
    
    for (int i=0;i<vertexValences.rows();i++){
      int beginH=VH(i);
      int currH=beginH;
      isBoundaryVertex(i)=1;
      
      int currCounter=0;
      while ((twinH(currH)!=-1)){
        currH=nextH(twinH(currH));
        if (currH==beginH) {isBoundaryVertex(i)=0; break;}
      }
      
      beginH=currH;
      
      do{
        starVertices(i,currCounter)=HV(nextH(currH));
        ringFaces(i,currCounter)=HF(currH);
        starHalfedges(i,currCounter++)=currH;
        if(twinH(prevH(currH))==-1){  //last edge on the boundary should be accounted for
          starVertices(i,currCounter)=HV(prevH(currH));
          starHalfedges(i,currCounter++)=prevH(currH);
        }
        currH=twinH(prevH(currH));
      }while ((beginH!=currH)&&(currH!=-1));
    }
    
    return true;
  }
  
  class OneRingSubdivisionData{
  public:
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXi D;
    Eigen::MatrixXi EV, FE, EFi, EF;
    Eigen::MatrixXd FEs;
    Eigen::VectorXi VH, innerEdges;
    Eigen::MatrixXi EH;
    Eigen::MatrixXi FH;
    Eigen::VectorXi HV;
    Eigen::VectorXi HE;
    Eigen::VectorXi HF;
    Eigen::VectorXi nextH;
    Eigen::VectorXi prevH;
    Eigen::VectorXi twinH;
    Eigen::VectorXi vertexValences;
    Eigen::MatrixXi starVertices;
    Eigen::MatrixXi starHalfedges;
    Eigen::MatrixXi ringFaces;
    Eigen::VectorXi isBoundaryVertex;
    
    virtual void setup(const Eigen::MatrixXd&, const Eigen::VectorXi&, const Eigen::MatrixXi&)=0;
    virtual Eigen::MatrixXd original2Canonical(const int, const Eigen::MatrixXd&)=0;
    virtual Eigen::MatrixXd canonical2Original(const int, const Eigen::MatrixXd&)=0;
    virtual Eigen::MatrixXd threePointsExtrapolation(const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&)=0;
    virtual Eigen::MatrixXd fourPointsInterpolation(const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&)=0;
    virtual Eigen::MatrixXd facePointBlend(const Eigen::MatrixXd&)=0;
    //virtual Eigen::MatrixXd EdgePointBlend(const Eigen::MatrixXd&){}
    virtual Eigen::RowVector3d innerVertexCanonicalBlend(const Eigen::RowVector3d&, const Eigen::MatrixXd&, const Eigen::MatrixXd&)=0;
    virtual Eigen::MatrixXd boundaryVertexPoint(const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&)=0;
    
    virtual Eigen::MatrixXd canonicalEdgePoints(const int v0, const Eigen::RowVector3d& canonCenter,const Eigen::MatrixXd& canonStarVertices,
                                               const Eigen::MatrixXd& canonFacePoints)=0;
    
    //virtual OneRingSubdivisionData()=0;
    //virtual ~OneRingSubdivisionData()=0;
    
  };
  
  
}


#endif


