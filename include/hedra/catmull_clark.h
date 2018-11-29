// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_CATMULL_CLARK_H
#define HEDRA_CATMULL_CLARK_H
#include <igl/igl_inline.h>
#include <hedra/polygonal_face_centers.h>
#include <hedra/dual_mesh.h>
#include <hedra/dcel.h>
#include <hedra/quaternionic_operations.h>
//#include <hedra/linear_subdivision_basics.h>
//#include <hedra/moebius_subdivision_basics.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
  
  typedef enum {LINEAR_SUBDIVISION, CANONICAL_MOEBIUS_SUBDIVISION} SubdivisionType;
  void vertex_valences(const Eigen::MatrixXi& EV,
                       Eigen::VectorXi& vertexValences)
  {
    vertexValences.conservativeResize(EV.maxCoeff()+1);
    vertexValences.setZero();
    for (int i=0;i<EV.rows();i++){
      vertexValences(EV(i,0))++;
      vertexValences(EV(i,1))++;
    }
  }
  
  
  void average_circle(const Eigen::MatrixXd& p,
                      const bool& isBoundary,
                      Eigen::RowVector3d& c,
                      double& r)
  {
    
    if (isBoundary)
      c=(p.row(0)+p.row(p.rows()-1))/2.0;
    else
      c=p.colwise().mean();
    r=(p.rowwise()-c).rowwise().norm().mean();
    return;
    
  }
  
  //Computes the point p on the edge jm with combinatorial edges ij, jk, jm, mn, ml
  IGL_INLINE Eigen::MatrixXd moebius_six_points_blend(const Eigen::MatrixXd& vi,
                                                      const Eigen::MatrixXd& vj,
                                                      const Eigen::MatrixXd& vk,
                                                      const Eigen::MatrixXd& vl,
                                                      const Eigen::MatrixXd& vm,
                                                      const Eigen::MatrixXd& vn)
  {
    using namespace Eigen;
    MatrixXd qi(vi.rows(),4); qi<<VectorXd::Zero(vi.rows()),vi;
    MatrixXd qj(vi.rows(),4); qj<<VectorXd::Zero(vj.rows()),vj;
    MatrixXd qk(vi.rows(),4); qk<<VectorXd::Zero(vk.rows()),vk;
    MatrixXd ql(vi.rows(),4); ql<<VectorXd::Zero(vl.rows()),vl;
    MatrixXd qm(vi.rows(),4); qm<<VectorXd::Zero(vm.rows()),vm;
    MatrixXd qn(vi.rows(),4); qn<<VectorXd::Zero(vn.rows()),vn;
    
    VectorXd qjmNorm=(qj-qm).rowwise().norm();
    
    double wj=0.5;
    MatrixXd Unit=MatrixXd::Zero(qi.rows(),4); Unit.col(0).setOnes();
    
    MatrixXd mijl=QMultN(QMultN(qm-qi, QInvN(qi-qj)), QMultN((qj-ql), QInvN(ql-qm)));
    MatrixXd mijn=QMultN(QMultN(qm-qi, QInvN(qi-qj)), QMultN((qj-qn), QInvN(qn-qm)));
    MatrixXd mkjn=QMultN(QMultN(qm-qk, QInvN(qk-qj)), QMultN((qj-qn), QInvN(qn-qm)));
    
    MatrixXd smijl=QExpN(QLogN(mijl)*0.5);
    MatrixXd smkjn=QExpN(QLogN(mkjn)*0.5);;
    MatrixXd r2=QMultN(QMultN(QInvN(smijl),mijn),QInv(smkjn));
    MatrixXd mijp=QMultN(smijl,QExpN(QLogN(r2)*0.5))*((1.0-wj)/wj);
    
    MatrixXd p=QMult(QInv(QMult(QMult(qi-qj, QInv(qm-qi)), -mijp)+Unit),QMult(QMult(qi-qj, QInv(qm-qi)),QMult(-mijp,qm))+qj);
    for (int i=0;i<vi.rows();i++)
      if (qjmNorm(i)<10e-6)
        p.row(i)<<0.0,vj.row(i);
    
    return p.block(0,1,p.rows(),3);
    
  }
  
  //computes the midpoint p to form the series a-b-p-c-d
  IGL_INLINE Eigen::MatrixXd moebius_four_points_blend(const Eigen::MatrixXd& va,
                                                       const Eigen::MatrixXd& vb,
                                                       const Eigen::MatrixXd& vc,
                                                       const Eigen::MatrixXd& vd)
  {
    using namespace Eigen;
    double wb=0.5;
    
    VectorXd normab=(va-vb).rowwise().squaredNorm();
    VectorXd normcd=(vc-vd).rowwise().squaredNorm();
    MatrixXd qa(va.rows(),4); qa<<VectorXd::Zero(va.rows()),va;
    MatrixXd qb(vb.rows(),4); qb<<VectorXd::Zero(va.rows()),vb;
    MatrixXd qc(vc.rows(),4); qc<<VectorXd::Zero(va.rows()),vc;
    MatrixXd qd(vd.rows(),4); qd<<VectorXd::Zero(va.rows()),vd;
    
    MatrixXd cabd=QMultN(QMultN(qa-qc, QInvN(qb-qa)), QMultN((qd-qb), QInvN(qc-qd)));
    
    RowVector4d minusOneQuat; minusOneQuat<<-1.0,0.0,0.0,0.0;
    
    MatrixXd cabp=-QExpN(QLogN(cabd)*0.5)*((1.0-wb)/wb);
    for (int i=0;i<va.rows();i++)
      if ((cabd.row(i)-minusOneQuat).squaredNorm()<10e-4){  //-1,0,0,0 square root not well defined
        RowVector3d faceNormal=QMult(qa.row(i)-qc.row(i), QInv(qb.row(i)-qa.row(i))).tail(3).normalized();
        cabp.row(i)<<0,faceNormal;
      }
    
    MatrixXd Unit=MatrixXd::Zero(qa.rows(),4); Unit.col(0).setOnes();
    
    MatrixXd p=QMultN(QInvN(QMultN(QMultN(qb-qa, QInvN(qa-qc)), -cabp)+Unit),QMultN(QMultN(qb-qa, QInvN(qa-qc)),QMultN(-cabp,qc))+qb);
    
    for (int i=0;i<va.rows();i++){
      if (normab(i)<10e-14)
        p.row(i)<<0.0,va.row(i);
      if (normcd(i)<10e-14)
        p.row(i)<<0.0,vd.row(i);
      
      
    }
    return p.block(0,1,p.rows(),3);
  }
  
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
    
    /*std::cout<<"vertexValences: "<<vertexValences<<std::endl;
    std::cout<<"starVertices: "<<starVertices<<std::endl;
    std::cout<<"ringFaces: "<<ringFaces<<std::endl;
    std::cout<<"starHalfedges: "<<starHalfedges<<std::endl;*/
    
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
    
    //virtual OneRingSubdivisionData()=0;
    //virtual ~OneRingSubdivisionData()=0;
    
  };
  
  
  class LinearCCSubdivisionData:public OneRingSubdivisionData{
  public:
    Eigen::MatrixXd centers;
    Eigen::VectorXd radii;
    
    void setup(const Eigen::MatrixXd& _V, const Eigen::VectorXi& _D, const Eigen::MatrixXi& _F){
      V=_V; D=_D; F=_F;
      hedra::polygonal_edge_topology(D, F, EV, FE, EF, EFi, FEs, innerEdges);
      hedra::dcel(D, F, EV,EF,EFi,innerEdges,VH,EH,FH,HV,HE,HF,nextH, prevH,twinH);
      hedra::vertex_stars(EV,VH,EH,FH,HV,HE,HF,nextH,prevH,twinH,vertexValences,starVertices,starHalfedges,ringFaces,isBoundaryVertex);
    }
    //there is no special canonical forms for linear subdivision
    Eigen::MatrixXd original2Canonical(const int, const Eigen::MatrixXd& origPoints){return origPoints;}
    Eigen::MatrixXd canonical2Original(const int, const Eigen::MatrixXd& origPoints){return origPoints;}
    Eigen::MatrixXd threePointsExtrapolation(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b , const Eigen::MatrixXd& c){
      return c*2.0-b;
    }
    Eigen::MatrixXd boundaryEdgePoint(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, const Eigen::MatrixXd& c, const Eigen::MatrixXd& d){
      return (b+c)/2.0;
    }
    Eigen::MatrixXd boundaryVertexPoint(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, const Eigen::MatrixXd& p, const Eigen::MatrixXd& c, const Eigen::MatrixXd& d){
      return ((a+d).array()/8.0+p.array()*3.0/4.0);
    }
    
    Eigen::MatrixXd fourPointsInterpolation(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, const Eigen::MatrixXd& c, const Eigen::MatrixXd& d){
      return (b+c)/2.0;
    }
    Eigen::MatrixXd facePointBlend(const Eigen::MatrixXd& candidateFacePoints)
    {
      
      Eigen::MatrixXd fineFacePoints=Eigen::MatrixXd::Zero(F.rows(),3);
      for (int i=0;i<D.rows();i++){
        for (int j=0;j<D(i);j++)
          fineFacePoints.row(i)+=candidateFacePoints.block(i,3*j,1,3);
        
        fineFacePoints.row(i)/=(double)D(i);
      }
      return fineFacePoints;
    }
    
    /*Eigen::MatrixXd EdgePointBlend(const Eigen::MatrixXd& candidateEdgePoints)
     {
     //these points are anyhow equal and this is also the boundary formula from the vertex stars function, so no special treatment needed.
     return(candidateEdgePoints.block(0,0,candidateEdgePoints.rows(),3)+candidateEdgePoints.block(0,3,candidateEdgePoints.rows(),3))/2.0;
     }*/
    
    Eigen::RowVector3d innerVertexCanonicalBlend(const Eigen::RowVector3d& canonCenter,
                                                 const Eigen::MatrixXd& canonEdgePoints,
                                                 const Eigen::MatrixXd& canonFacePoints)
    {
      
      Eigen::RowVector3d F = canonFacePoints.colwise().mean();
      Eigen::RowVector3d E = canonEdgePoints.colwise().mean();
      return((F+E*4.0-F*2.0+(double)(canonEdgePoints.rows()-3)*canonCenter)/(double)canonEdgePoints.rows());
    }
    
    LinearCCSubdivisionData(){}
    ~LinearCCSubdivisionData(){}
    
  };
  
  class MoebiusCCSubdivisionData:public OneRingSubdivisionData{
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
      return moebius_four_points_blend(a,b,c,d);
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
      
      Eigen::RowVector3d F = canonFacePoints.colwise().mean();
      Eigen::RowVector3d E = canonEdgePoints.colwise().mean();
      return((F+E*4.0-F*2.0+(double)(canonEdgePoints.rows()-3)*canonCenter)/(double)canonEdgePoints.rows());
    }
    
    MoebiusCCSubdivisionData(){}
    ~MoebiusCCSubdivisionData(){}

  };
  
  
  
  
  
  
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
        
        MatrixXd canonFaceVertices;
        canonFaceVertices=sd.original2Canonical(i,origFaceVertices);
        canonFacePoints.row(j)=canonFaceVertices.colwise().mean();
        //std::cout<<"canonFacePoints.row(j): "<<canonFacePoints.row(j)<<std::endl;
        
        //Lifting to candidate points
        candidateFacePoints.block(sd.ringFaces(i,j), 3*currVertexinFace,1,3)=sd.canonical2Original(i,canonFacePoints.row(j));  //not entirely optimized
      }
      
      //std::cout<<"candidateFacePoints: "<<candidateFacePoints<<std::endl;
      
      //candidate edge points
      MatrixXd canonEdgePoints(sd.vertexValences(i),3);
      for (int j=sd.isBoundaryVertex(i);j<sd.vertexValences(i)-sd.isBoundaryVertex(i);j++)
        canonEdgePoints.row(j)=(canonCenter+canonStarVertices.row(j)+canonFacePoints.row(j)+canonFacePoints.row((j+sd.vertexValences(i)-1)%sd.vertexValences(i)))/4.0;
      
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
        /*RowVector3d F = canonFacePoints.colwise().mean();
         RowVector3d E = canonEdgePoints.colwise().mean();
         canonFineCenter=(F+E*4.0-F*2.0+(double)(vertexValences(i)-3)*canonCenter)/(double)vertexValences(i);*/
      } else {  //this will get assigned later
        //This will get overri
        /*if (vertexValences(i)>2)
         canonFineCenter=canonCenter*3.0/4.0+(canonStarVertices.row(0)+canonStarVertices.row(canonStarVertices.rows()-1))/8.0;
         else  //an ear
         canonFineCenter=canonCenter;*/
      }
      
      fineVertexPoints.row(i)=sd.canonical2Original(i,canonFineCenter);
    }
    
    //std::cout<<"candidateFacePoints: "<<candidateFacePoints<<std::endl;
    //std::cout<<"candidateEdgePoints: "<<candidateEdgePoints<<std::endl;
    
    //Blending face points from candidates
    fineFacePoints = sd.facePointBlend(candidateFacePoints);
    
    //Blending edge points from candidates, and boundary edge points from boundary curves.
    
    Eigen::MatrixXd a(sd.EH.rows(),3), b(sd.EH.rows(),3), c(sd.EH.rows(),3), d(sd.EH.rows(),3);
    
    //default values for inner edges
    /*b=candidateEdgePoints.block(0,0,candidateEdgePoints.rows(),3);
     c=candidateEdgePoints.block(0,3,candidateEdgePoints.rows(),3);*/
    
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
    
    fineEdgePoints=sd.fourPointsInterpolation(a,b,c,d);
    
    
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
  
  IGL_INLINE bool catmull_clark(const Eigen::MatrixXd& V,
                                const Eigen::VectorXi& D,
                                const Eigen::MatrixXi& F,
                                const hedra::SubdivisionType st,
                                Eigen::MatrixXd& fineV,
                                Eigen::VectorXi& fineD,
                                Eigen::MatrixXi& fineF)
  
  {
    using namespace Eigen;
    using namespace std;
    
    
    
    switch (st){
      case hedra::LINEAR_SUBDIVISION: {
        hedra::LinearCCSubdivisionData lsd;
        vertex_insertion(V, D,F,lsd, fineV, fineD, fineF);
        break;
      }
      case hedra::CANONICAL_MOEBIUS_SUBDIVISION: {
        hedra::MoebiusCCSubdivisionData msd;
        vertex_insertion(V, D,F,msd, fineV, fineD, fineF);
        break;
      }
      default: return false;
    }
    
    return true;
  }
  
  
}


#endif


