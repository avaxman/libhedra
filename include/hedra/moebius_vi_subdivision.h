// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_MOEBIUS_VI_SUBDIVISION_H
#define HEDRA_MOEBIUS_VI_SUBDIVISION_H
#include <igl/igl_inline.h>
#include <hedra/quaternionic_operations.h>
#include <hedra/vertex_valences.h>
#include <hedra/moebius_refinement.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
  
  class MoebiusVISubdivisionData:public OneRingSubdivisionData{
  public:
    Eigen::MatrixXd centers;
    Eigen::VectorXd radii;
    
    void setup(const Eigen::MatrixXd& _V, const Eigen::VectorXi& _D, const Eigen::MatrixXi& _F){
      V=_V; D=_D; F=_F;
      hedra::polygonal_edge_topology(D, F, EV, FE, EF, EFi, FEs, innerEdges);
      hedra::dcel(D, F, EV,EF,EFi,innerEdges,VH,EH,FH,HV,HE,HF,nextH, prevH,twinH);
      hedra::vertex_stars(EV,VH,EH,FH,HV,HE,HF,nextH,prevH,twinH,vertexValences,starVertices,starHalfedges,ringFaces,isBoundaryVertex);
      
      //creating Canonical Moebius forms
      //first creating the q_{ki}^{-1} form
      centers.conservativeResize(V.rows(),3);
      radii.conservativeResize(V.rows());
      for (int i=0;i<V.rows();i++){
        Eigen::MatrixXd qce; qce.resize(vertexValences(i),4);
        
        
        for (int j=0;j<vertexValences(i);j++){
          Eigen::RowVector4d qki; qki<<0.0,V.row(i)-V.row(starVertices(i,j));
          qce.row(j)<<QInv(qki);
        }
        
        //cout<<"qce after initial inversion: "<<qce<<endl;
        
        if (vertexValences(i)<=2){  //an ear; the canonical embedding is not well defined anyhow
          centers.row(i).setZero();
          
          radii(i)=2.0/((V.row(starVertices(i,0))-V.row(i)).norm()+(V.row(starVertices(i,1))-V.row(i)).norm());
          std::cout<<"trying to compute canonical embedding for an ear!!"<<std::endl;
          return;
        }
        
        
        //computing best fit circle
        Eigen::RowVector3d currCenter;
        
        hedra::average_circle(qce.block(0,1,qce.rows(),3), isBoundaryVertex(i), currCenter, radii(i));
        centers.row(i)=currCenter;
      }
    }
    //there is no special canonical forms for linear subdivision
    Eigen::MatrixXd original2Canonical(const int v0, const Eigen::MatrixXd& origPoints){
      Eigen::MatrixXd qop1(origPoints.rows(),4);
      Eigen::MatrixXd qop2(origPoints.rows(),4);
      
      //inverting B->T
      for (int i=0;i<origPoints.rows();i++){
        Eigen::RowVector4d qki; qki<<0.0,V.row(v0)-origPoints.row(i);
        qop1.row(i)<<QInv(qki);
      }
      
      //cout<<"qop1: "<<qop1<<endl;
      
      Eigen::RowVector4d qc; qc<<0.0,centers.row(v0);
      
      //inverting T->T'
      for (int i=0;i<origPoints.rows();i++){
        if ((origPoints.row(i)-V.row(v0)).squaredNorm()<10e-10)
          qop2.row(i)<<0.0,centers.row(v0);
        else qop2.row(i)<<radii(v0)*radii(v0)*QInv(qc-qop1.row(i))+qc;
      }
      
      //cout<<"qop2: "<<qop1<<endl;
      return (qop2.block(0,1,qop2.rows(),3));
    }
    Eigen::MatrixXd canonical2Original(const int v0, const Eigen::MatrixXd& origPoints){
      Eigen::MatrixXd qop1(origPoints.rows(),4);
      Eigen::MatrixXd qop(origPoints.rows(),4);
      
      Eigen::RowVector4d qc; qc<<0.0,centers.row(v0);
      Eigen::RowVector4d qv0; qv0<<0.0,V.row(v0);
      for (int i=0;i<origPoints.rows();i++){
        Eigen::RowVector4d qep; qep<<0.0, origPoints.row(i);
        if ((qep-qc).squaredNorm()>10e-10){
          qep=qc-radii(v0)*radii(v0)*QInv(qep-qc);
          
          qop.row(i)=qv0-QInv(qep);
        } else qop.row(i)=qv0;
      }
      
      return (qop.block(0,1,qop.rows(),3));
      
    }
    Eigen::MatrixXd threePointsExtrapolation(const Eigen::MatrixXd& va, const Eigen::MatrixXd& vb , const Eigen::MatrixXd& vc){
      using namespace Eigen;
      MatrixXd qa(va.rows(),4); qa<<VectorXd::Zero(va.rows()),va;
      MatrixXd qb(vb.rows(),4); qb<<VectorXd::Zero(vb.rows()),vb;
      MatrixXd qc(vc.rows(),4); qc<<VectorXd::Zero(vc.rows()),vc;
      MatrixXd abcd=MatrixXd::Zero(va.rows(),4);
      abcd.col(0).setConstant(-1.0/3.0);
      
      MatrixXd Unit=MatrixXd::Zero(qa.rows(),4); Unit.col(0).setOnes();
      
      RowVector4d d=QMultN(QInvN(Unit+QMultN(QMultN(qc-qb,QInvN(qb-qa)),abcd)),qc+QMultN(QMultN(qc-qb,QInvN(qb-qa)),QMultN(abcd,qa)));
      
      return d.block(0,1,d.rows(),3);
      
    }
    Eigen::MatrixXd boundaryEdgePoint(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, const Eigen::MatrixXd& c, const Eigen::MatrixXd& d){
      return moebius_four_points_blend(a,b,c,d);
    }
    
    Eigen::MatrixXd fourPointsInterpolation(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, const Eigen::MatrixXd& c, const Eigen::MatrixXd& d){
      return moebius_four_points_blend(a,b,c,d);
    }
    Eigen::MatrixXd boundaryVertexPoint(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, const Eigen::MatrixXd& p, const Eigen::MatrixXd& c, const Eigen::MatrixXd& d){
      return p;
    }
    Eigen::MatrixXd facePointBlend(const Eigen::MatrixXd& candidateFacePoints)
    {
      using namespace Eigen;
      MatrixXd fineFacePoints=MatrixXd::Zero(F.rows(),3);
      for (int i=0;i<D.rows();i++){
        
        //finding weirdest cross-ratio in terms of angle
        double nonCircularity=0.0;
        int startIndex=0;
        for (int j=0;j<D(i);j++){
          RowVector4d qa; qa<<0.0,V.row(F(i,j));
          RowVector4d qb; qb<<0.0,candidateFacePoints.block(i, 3*j,1,3);
          RowVector4d qc; qc<<0.0,candidateFacePoints.block(i, 3*((j+D(i)/2)%D(i)),1,3) ;
          RowVector4d qd; qd<<0.0,V.row(((j+D(i)/2)%D(i)));
          /*std::cout<<"qa: "<<qa<<std::endl;
           std::cout<<"qb: "<<qb<<std::endl;
           std::cout<<"qc: "<<qc<<std::endl;
           std::cout<<"qd: "<<qd<<std::endl;*/
          RowVector4d cabd=QMult(QMult(qa-qc, QInv(qb-qa)), QMult((qd-qb), QInv(qc-qd)));
          double currNonCirularity=abs(cabd(1)/cabd.norm()-1.0);
          if (currNonCirularity >= nonCircularity){
            nonCircularity =currNonCirularity;
            startIndex=i;
          }
        }
        
        MatrixXd oppositePoints(D(i),3);
        for (int j=0;j<D(i);j++)
          oppositePoints.row(j)=moebius_four_points_blend(V.row(F(i,j)), candidateFacePoints.block(i, 3*j,1,3), candidateFacePoints.block(i, 3*((j+D(i)/2)%D(i)),1,3), V.row(F(i,(j+D(i)/2)%D(i))));
        
        //std::cout<<"oppositePoints: "<<oppositePoints<<std::endl;
        
        //sequentially computing six points
        MatrixXd seqFacePoint=oppositePoints.row(1);
        for (int j=startIndex;j<D(i)+startIndex;j++)
          seqFacePoint=moebius_six_points_blend(V.row(F(i,j%D(i))), seqFacePoint,  V.row(F(i,(j+D(i)/2)%D(i))), V.row(F(i,(j+1)%D(i))), oppositePoints.row((j+1)%D(i)),  V.row(F(i,(j+1+D(i)/2)%D(i))));
        
        fineFacePoints.row(i)= seqFacePoint;
      }
      return fineFacePoints;
    }
    
    Eigen::MatrixXd EdgePointBlend(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, const Eigen::MatrixXd& c, const Eigen::MatrixXd& d){
      return moebius_four_points_blend(a,b,c,d);
    }
    
    Eigen::RowVector3d innerVertexCanonicalBlend(const Eigen::RowVector3d& canonCenter,
                                                 const Eigen::MatrixXd& canonEdgePoints,
                                                 const Eigen::MatrixXd& canonFacePoints)
    {
      
      return canonCenter;
    }
    
    Eigen::MatrixXd canonicalEdgePoints(const int v0,
                                        const Eigen::RowVector3d& canonCenter,
                                        const Eigen::MatrixXd& canonStarVertices,
                                        const Eigen::MatrixXd& canonFacePoints)
    {
      Eigen::MatrixXd canonEdgePoints(vertexValences(v0),3);
      for (int j=0;j<vertexValences(v0);j++)
        canonEdgePoints.row(j)=(canonCenter+canonStarVertices.row(j))/2.0;
      
      return canonEdgePoints;
      
    }
    
    MoebiusVISubdivisionData(){}
    ~MoebiusVISubdivisionData(){}
    
  };
  
}


#endif


