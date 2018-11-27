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
//#include <hedra/linear_subdivision_basics.h>
//#include <hedra/moebius_subdivision_basics.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
  
  void computeBestFitCircle(const MatrixXd& p, const bool isBoundary, RowVector3d& c, double& r)
  {
    
    if (isBoundary)
      c=(p.row(0)+p.row(p.rows()-1))/2.0;
    else
      c=p.colwise().mean();
    r=(p.rowwise()-c).rowwise().norm().mean();
    return;
  
  }
  
  struct CanonicalFunction{
    
    cf.original2canonical = hedra::trivial_original2canonical;
    cf.canonical2original = hedra::trivial_canonical2original;
    
  };

  
  
  //Computes the point p on the edge jm with combinatorial edges ij, jk, jm, mn, ml
  IGL_INLINE Eigen::MatrixXd moebiusSixPointsBlend(const Eigen::MatrixXd& vi,
                                                   const Eigen::MatrixXd& vj,
                                                   const Eigen::MatrixXd& vk,
                                                   const Eigen::MatrixXd& vl,
                                                   const Eigen::MatrixXd& vm,
                                                   const Eigen::MatrixXd& vn)
  {
    using namespace Eigen;
    MatrixXd qi(vi.rows(),4); qi<<0.0,vi;
    MatrixXd qj(vi.rows(),4); qj<<0.0,vj;
    MatrixXd qk(vi.rows(),4); qk<<0.0,vk;
    MatrixXd ql(vi.rows(),4); ql<<0.0,vl;
    MatrixXd qm(vi.rows(),4); qm<<0.0,vm;
    MatrixXd qn(vi.rows(),4); qn<<0.0,vn;
    
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
  IGL_INLINE Eigen::MatrixXd moebiusFourPointsBlend(const Eigen::MatrixXd& va,
                                                    const Eigen::MatrixXd& vb,
                                                    const Eigen::MatrixXd& vc,
                                                    const Eigen::MatrixXd& vd)
  {
    double wb=0.5;
    
    VectorXd normab=(va-vb).rowwise().squaredNorm();
    VectorXd normcd=(vc-vd).rowwise().squaredNorm();
    
    
    using namespace Eigen;
    MatrixXd qa(vi.rows(),4); qa<<0.0,va;
    MatrixXd qb(vi.rows(),4); qb<<0.0,vb;
    MatrixXd qc(vi.rows(),4); qc<<0.0,vc;
    MatrixXd qd(vi.rows(),4); qd<<0.0,vd;
    
    MatrixXd cabd=QMult(QMultN(qa-qc, QInvN(qb-qa)), QMultN((qd-qb), QInvN(qc-qd)));
    
    RowVector4d minusOneQuat; minusOneQuat<<-1.0,0.0,0.0,0.0;
    
    MatrixXd cabp=-QExpN(QLogN(cabd)*0.5)*((1.0-wb)/wb);
    for (int i=0;i<vi.rows();i++)
      if ((cabd.row(i)-minusOneQuat).squaredNorm()<10e-4){  //-1,0,0,0 square root not well defined
        RowVector3d faceNormal=QMult(qa.row(i)-qc.row(i), QInv(qb.row(i)-qa.row(i))).tail(3).normalized();
        cabp.row(i)<<0,faceNormal;
      }
    
    MatrixXd Unit=MatrixXd::Zero(qi.rows(),4); Unit.col(0).setOnes();
    
    MatrixXd p=QMultN(QInvN(QMultN(QMultN(qb-qa, QInvN(qa-qc)), -cabp)+Unit),QMultN(QMultN(qb-qa, QInvN(qa-qc)),QMultN(-cabp,qc))+qb);
    
    for (int i=0;i<vi.rows();i++){
      if (normab(i)<10e-14)
        p.row(i)<<0.0,va.row(i);
      if (normcd(i)<10e-14)
        p.row(i)<<0.0,vd.row(i);
      
      return p.block(0,1,p.rows(),3);
    }
    
  
  IGL_INLINE bool vertex_stars(const Eigen::VectorXi& D,
                               const Eigen::MatrixXi& F,
                               const Eigen::MatrixXi& EV,
                               const Eigen::MatrixXi& EF,
                               const Eigen::MatrixXi& EFi,
                               const Eigen::MatrixXi& FE,
                               const Eigen::VectorXi& innerEdges,
                               const Eigen::VectorXi& VH,
                               const Eigen::MatrixXi& EH,
                               const Eigen::MatrixXi& FH,
                               const Eigen::VectorXi& HV,
                               const Eigen::VectorXi& HE,
                               const Eigen::VectorXi& HF,
                               const Eigen::VectorXi& nextH,
                               const Eigen::VectorXi& prevH,
                               const Eigen::VectorXi& twinH,
                               Eigen::MatrixXi& starVertices,
                               Eigen::MatrixXi& starHalfedges,
                               Eigen::MatrixXi& ringFaces){
    
    using namespace Eigen;
    
    VectorXi vertexValences;
    hedra::vertex_valences(D,F,vertexValences);
    
    int maxValence =vertexValences.maxCoeff();
    starVertices.conservativeResize(V.rows(),maxValence);
    starHalfedges.conservativeResize(V.rows(),maxValence);
    ringFaces.conservativeResize(V.rows(),maxValence);
    
    VectorXi isBoundaryVertex(V.rows());
    
    for (int i=0;i<F.maxCoeff().rows();i++){
      int beginH=VH(i);
      int currH=beginH;
      isBoundaryVertex(i)=1;
    
      int currCounter;
      while ((twinH(currH)!=-1)){
        currH=nextH(twinH(currH));
        if (currH==beginH) {isBoundaryVertex(i)=0; break;}
      }
      
      beginH=currH;
      
      do{
        starVertices(i,currCounter)=HV(nextH(currH));
        starFaces(i,currCounter)=HF(currH);
        starHalfedges(i,currCounter++)=HE(currH);
        if(twinH(prevH(currH))==-1){  //last edge on the boundary should be accounted for
          starVertices(i,currCounter)=HV(prevH(currH));
          starHalfedges(i,currCounter++)=HE(prevH(currH));
        }
        currH=twinH(prevH(currH));
      }while ((beginH!=currH)&&(currH!=-1));
    }
  }
  
  
  
  
  IGL_INLINE bool face_point_linear_average(const Eigen::MatrixXd& V,
                                            const Eigen::VectorXi& D,
                                            const Eigen::MatrixXi& F,
                                            const Eigen::MatrixXd& candidateFacePoints,
                                            Eigen::MatrixXd fineFacePoints)
  {
    
    fineFacePoints=MatrixXd::Zeros(F.rows(),3);
    for (int i=0;i<D.rows();i++){
      for (int j=0;j<D(i);j++)
        fineFacePoints.row(i)+=candidateFacePoints.block(i,3*j,1,3);
      
      fineFacePoints.row(i)/=(double)D(i);
    }
  }
  
  IGL_INLINE bool edge_point_linear_average(const Eigen::MatrixXd& V,
                                            const Eigen::VectorXi& D,
                                            const Eigen::MatrixXi& F,
                                            const Eigen::MatrixXi& EV,
                                            const Eigen::VectorXi& innerEdges,
                                            const Eigen::MatrixXd& candidateEdgePoints,
                                            const Eigen::MatrixXd& fineEdgePoints)
  {
    //these points are anyhow equal and this is also the boundary formula from the vertex stars function, so no special treatment needed.
    fineEdgePoints = (candidateEdgePoints.block(0,0,candidateEdgePoints.rows(),3)+candidateEdgePoints.block(0,3,candidateEdgePoints.rows(),3))/2.0;
  }
  
  IGL_INLINE bool face_point_moebius_average(const Eigen::MatrixXd& V,
                                             const Eigen::VectorXi& D,
                                             const Eigen::MatrixXi& F,
                                             const Eigen::MatrixXd& candidateFacePoints,
                                             Eigen::MatrixXd fineFacePoints)
  {
    
    using namespace Eigen;
    fineFacePoints=MatrixXd::Zeros(F.rows(),3);
    for (int i=0;i<D.rows();i++){
      
      //finding weirdest cross-ratio in terms of angle
      double nonCircularity=0.0;
      int startIndex=0;
      for (int j=0;j<D(i);j++){
        RowVector4d qa; qa<<0.0,V.row(F(i,j));
        RowVector4d qb; qb<<0.0,candidateFacePoints.block(F(i,j), 3*j,1,3);
        RowVector4d qc; qc<<0.0,candidateFacePoints.block(F(i,j), 3*((j+D(i)/2)%D(i)),1,3) ;
        RowVector4d qd; qd<<0.0,vertexPoints.row(((j+D(i)/2)%D(i)));
        RowVector4d cabd=QMult(QMult(qa-qc, QInv(qb-qa)), QMult((qd-qb), QInv(qc-qd)));
        double currNonCirularity=abs(cabd(1)/cabd.norm()-1.0);
        if (currNonCirularity >= nonCircularity){
          nonCircularity =currNonCirularity;
          startIndex=i;
        }
      }
      
      MatrixXd oppositePoints(D(i),3);
      for (int j=0;j<D(i);j++)
        oppositePoints.row(j)=moebiusFourPointsBlend(V.row(F(i,j)), candidateFacePoints.block(F(i,j), 3*j,1,3), candidateFacePoints.block(F(i,j), 3*((j+D(i)/2)%D(i)),1,3), vertexPoints.row(((j+D(i)/2)%D(i))));
      
      //sequentially computing six points
      MatrixXd seqFacePoint=oppositePoints.row(1);
      for (int j=startIndex;j<D(i)+startIndex;j++)
        seqFacePoint=moebiusSixPointsBlend(V.row(F(i,j%N)), seqFacePoint,  V.row(F(i,j+D(i)/2)%N), V.row(F(i,(j+1)%N)), oppositePoints.row((j+1)%N),  V.row(F(i,(j+1+D(i)/2)%N)));
                                                                                                        
      return seqFacePoint;
    }
  }
  
  IGL_INLINE bool edge_point_moebius_average(const Eigen::MatrixXd& V,
                                             const Eigen::VectorXi& D,
                                             const Eigen::MatrixXi& F,
                                             const Eigen::MatrixXi& EV,
                                             const Eigen::VectorXi& innerEdges,
                                             const Eigen::MatrixXd& candidateEdgePoints,
                                             const Eigen::MatrixXd& fineEdgePoints)
  {
    
    MatrixXd a(EV.rows(),3), b(EV.rows(),3) c(EV.rows(),3), d(EV.rows(),3);
    
    //default values for inner edges
    b=candidateEdgePoints.block(0,0,candidateEdgePoints.rows(),3);
    c=candidateEdgePoints.block(0,3,candidateEdgePoints.rows(),3);
    
    //for inner vertices: V->candidateedgepoints->V blending (already encoded from the vertex stars function)
    //for boundary vertices: like a 4-points curve (need to be modified as follows)
    
    for (int i=0;i<EH.rows();i++){
      int currH;
      if (EH(i,0)==-1)
        currH=EH(i,1);
      else currH=EH(i,0);
      int ix2=HV(currH);
      int ix3=HV(nextH(currH));
      
      RowVector3d x2; x2<<V.row(ix2);
      RowVector3d x3; x3<<V.row(ix3);
      
      if (twinH(currH)!=-1){  //inner edge
        a.row(i)=V.row(EV(i,0));
        d.row(i)=V.row(EV(i,1));
        continue;
      }
      
      bool x1Ear=(twinH(prevH(currH))==-1);
      bool x4Ear=(twinH(nextH(currH))==-1);
      RowVector3d x1, x4;
      
      if (!x1Ear){
        int otherBoundH=currH;
        do{
          otherBoundH=twinH(prevH(otherBoundH));
        }while(twinH[level](prevH(otherBoundH))!=-1);
        int ix1=HV(prevH(otherBoundH));
        x1=V.row(ix1);
      }
      
      if (!x4Ear){
        int otherBoundH=currH;
        do{
          otherBoundH=twinH(nextH[level](otherBoundH));
        }while(twinH(nextH(otherBoundH))!=-1);
        int ix4=HV(nextH(nextH(otherBoundH)));
        x4=V.row(ix4);
      }
      
      if (x1Ear)
        x1=moebiusExtrapolation(x4, x3, x2);
      if (x4Ear)
        x4=moebiusExtrapolation(x1, x2, x3);
      
      
      a.row(i)=x1;
      d.row(i)=x4;
    }
    
    fineEdgePoints=moebiusFourPointsBlend(a,b,c,d);
    
  }
  
  IGL_INLINE bool boundary_vertex_moebius_average(const Eigen::MatrixXd& V,
                                                  const Eigen::VectorXi& VH,
                                                  const Eigen::MatrixXi& EH,
                                                  const Eigen::MatrixXi& FH,
                                                  const Eigen::VectorXi& HV,
                                                  const Eigen::VectorXi& HE,
                                                  const Eigen::VectorXi& HF,
                                                  const Eigen::VectorXi& nextH,
                                                  const Eigen::VectorXi& prevH,
                                                  const Eigen::VectorXi& twinH,
                                                  const Eigen::MatrixXd edgePoints
                                                  Eigen::MatrixXi& fineV){
    
    for (int i=0;i<EH.rows();i++){
      int currH;
      if (EH(i,0)==-1)
        currH=EH(i,1);
      else currH=EH(i,0);
      
      if (twinH(currH)!=-1)
        continue;  //not a boundary
      
      int ix2=HV(currH);
      int ix3=HV(nextH[level](currH));
      RowVector3d x2; x2<<V.row(ix2);
      RowVector3d x3; x3<<V.row(ix3);
      bool x1Ear=(twinH(prevH(currH))==-1);
      RowVector3d x1;
      if (!x1Ear){
        int otherBoundH=currH;
        do{
          otherBoundH=twinH(prevH(otherBoundH));
        }while(twinH[level](prevH(otherBoundH))!=-1);
        int ix1=HV(prevH(otherBoundH));
        x1=V.row(ix1);
        
        fineV.row(ix2)=moebiusFourPointsBlend(x1,edgePoints.row(HE(prevH(otherBoundH))), edgePoints.row(i), x3);
      } else {
        fineV.row(ix2) = V.row(ix2);
      }
    }
    
  }
                                                                                                          
                                                                                                          
  
  
  IGL_INLINE bool canonical_star_subdivision_cc(const Eigen::MatrixXd& V,
                                 const Eigen::VectorXi& D,
                                 const Eigen::MatrixXi& F,
                                 const Eigen::MatrixXi& EV,
                                 const Eigen::MatrixXi& EF,
                                 const Eigen::MatrixXi& EFi,
                                 const Eigen::MatrixXi& FE,
                                 const Eigen::VectorXi& innerEdges,
                                 const Eigen::VectorXi& VH,
                                 const Eigen::MatrixXi& EH,
                                 const Eigen::MatrixXi& FH,
                                 const Eigen::VectorXi& HV,
                                 const Eigen::VectorXi& HE,
                                 const Eigen::VectorXi& HF,
                                 const Eigen::VectorXi& nextH,
                                 const Eigen::VectorXi& prevH,
                                 const Eigen::VectorXi& twinH,
                                 const hedra::CanonicalData cf,
                                 Eigen::MatrixXd& fineVertexPoints,
                                 Eigen::MatrixXd& candidateEdgePoints,
                                 Eigen::MatrixXd& candidateFacePoints){
    
    using namespace Eigen;
    
    hedra::vertex_stars(D,F,VH, EH, FH, HV, HE, HF,starVertices,starHalfedges,ringFaces, isBoundaryVertex);
    cf.setup_canonical_embeddings(V, vertexValences,starVertices,starHalfedges,ringFaces, isBoundaryVertex);
    
    fineVertexPoints.conservativeResize(V.rows(),3);
    candidateEdgePoints.conservativeResize(EV.rows(),2*3);  //two per edge
    candidateFacePoints.conservativeResize(F.rows(),D.maxCoeff*3);  //D(f) per face f
    
    for (int i=0;i<V.rows();i++){
      MatrixXd origStarVertices(vertexValences(i)));
      for (int j=0;j<vertexValences(i);j++)
        origStarVertices.row(j)=V.row(starVertices(i,j));
      
      MatrixXd canonStarVertices;
      cf.original2Canonical(origStarVertices, canonStarVertices);
      
      RowVector3d canonCenter;
      cf.original2Canonical(V.row(i), canonCenter);
      
      //face points
      int numRingFaces=vertexValences(i)-isBoundarVertex(i);  //one less face for boundary vertices
      MatrixXd canonFacePoints(numRingFaces,1);
      for (int j=0;j<numRingFaces;j++){
        MatrixXd origFaceVertices(D(ringFaces(i,j)),3);
        int currVertexinFace=-1;
        for (int k=0;k<D(ringFaces(i,j));k++){
          origFaceVertices.row(k)=V.row(F(ringFaces(i,j),k));
          if (F(ringFaces(i,j),k)==i)
            currVertexinFace=k;
        }
        
        MatrixXd canonFaceVertices;
        cf.original2Canonical(origFaceVertices, canonFaceVertices);
        canonFacePoints.row(j)=canonFaceVertices.colwise().mean();
      
        //Lifting to candidate points
        candidateFacePoints.block(ringFaces(i,j), 3*currVertexinFace,1,3)=cf.canonical2Original(canonFacePoints.row(j));  //not entirely optimized
      }
      
      //edge points
      MatrixXd canonEdgePoints(vertexValences(i),3);
      for (int j=isBoundaryVertex(i);j<vertexValences(i)-isBoundaryVertex(i);j++)
        canonEdgePoints.row(i)=(canonCenter+canonStarVertices.row(j)+canonFacePoints.row(j)+canonFacePoints.row((j+vertexValences(i)-1)%vertexValences(i)))/4.0;
      
      if (isBoundaryVertex(i)){
        canonEdgePoints.row(0)=(canonCenter+canonStarVertices.row(0))/2.0;
        canonEdgePoints.row(edgePoints.rows()-1)=(canonCenter+canonStarVertices.row(vertexValences(i)-1))/2.0;
      }
      
      //Lifting to candidate points
      MatrixXd localCandidateEdgePoints;
      cf.canonical2Original(canonEdgePoints,localCandidateEdgePoints);
      for (int j=0;j<vertexValences(i);j++)
        candidateEdgePoints.block(i,3*j,1,3)=localCandidateEdgePoints.row(j);  //WRONG!!!
      
      //vertex points
      RowVector3d canonFineCenter;
      if (!isBoundaryVertex){
        RowVector3d F = canonFacePoints.colwise().mean();
        RowVector3d E = canonEdgePoints.colwise().mean();
        canonFineCentral=(F+E*4.0-F*2.0+(double)(vertexValences(i)-3)*canonCenter)/(double)vertexValences(i);
      } else {
        if (vr.rows()>2)
         canonFineCenter=canonCenter*3.0/4.0+(canonStarVertices.row(0)+canonStarVertices.row(canonStarVertices.rows()-1))/8.0;
        else  //an ear
          canonFineCenter=canonCenter;
      }
      
      fineVertexPoints.row(i)=cf.canonical2Original(canonFineCenter);
      
    }
    
  }
  
  IGL_INLINE bool canonical_star_subdivision_cc_moebius(const Eigen::MatrixXd& V,
                                         const Eigen::VectorXi& D,
                                         const Eigen::MatrixXi& F,
                                         const Eigen::MatrixXi& EV,
                                         const Eigen::MatrixXi& EF,
                                         const Eigen::MatrixXi& EFi,
                                         const Eigen::MatrixXi& FE,
                                         const Eigen::VectorXi& innerEdges,
                                         const Eigen::VectorXi& VH,
                                         const Eigen::MatrixXi& EH,
                                         const Eigen::MatrixXi& FH,
                                         const Eigen::VectorXi& HV,
                                         const Eigen::VectorXi& HE,
                                         const Eigen::VectorXi& HF,
                                         const Eigen::VectorXi& nextH,
                                         const Eigen::VectorXi& prevH,
                                         const Eigen::VectorXi& twinH,
                                         Eigen::MatrixXd& fineVertexPoints,
                                         Eigen::MatrixXd& candidateEdgePoints,
                                         Eigen::MatrixXd& candidateFacePoints)
  {
    canonicalFunction cf;
    cf.original2canonical = hedra::moebius_original2canonical;
    cf.canonical2original = hedra::moebius_canonical2original;
    canonical_star_subdivision_cc(V,D,F,EV,EF,EFi,innerEdges,VH,EH,FH,HV,HE,HF,nextH,prevH,twinH, cf, fineVertexPoints, candidateEdgePoints, candidateFacePoints);
  }
  
  IGL_INLINE bool canonical_star_subdivision_cc_linear(const Eigen::MatrixXd& V,
                                        const Eigen::VectorXi& D,
                                        const Eigen::MatrixXi& F,
                                        const Eigen::MatrixXi& EV,
                                        const Eigen::MatrixXi& EF,
                                        const Eigen::MatrixXi& EFi,
                                        const Eigen::MatrixXi& FE,
                                        const Eigen::VectorXi& innerEdges,
                                        const Eigen::VectorXi& VH,
                                        const Eigen::MatrixXi& EH,
                                        const Eigen::MatrixXi& FH,
                                        const Eigen::VectorXi& HV,
                                        const Eigen::VectorXi& HE,
                                        const Eigen::VectorXi& HF,
                                        const Eigen::VectorXi& nextH,
                                        const Eigen::VectorXi& prevH,
                                        const Eigen::VectorXi& twinH,
                                        Eigen::MatrixXd& fineVertexPoints,
                                        Eigen::MatrixXd& candidateEdgePoints,
                                        Eigen::MatrixXd& candidateFacePoints)
  {
    canonicalFunction cf;
    cf.original2canonical = hedra::trivial_original2canonical;
    cf.canonical2original = hedra::trivial_canonical2original;
    canonical_star_subdivision_cc(V,D,F,EV,EF,EFi,innerEdges,VH,EH,FH,HV,HE,HF,nextH,prevH,twinH, cf, fineVertexPoints, candidateEdgePoints, candidateFacePoints);
  }
  
  IGL_INLINE bool catmull_clark(const Eigen::MatrixXd& V,
                                const Eigen::VectorXi& D,
                                const Eigen::MatrixXi& F,
                                const Eigen::MatrixXi& EV,
                                const Eigen::MatrixXi& EF,
                                const Eigen::MatrixXi& EFi,
                                const Eigen::MatrixXi& FE,
                                const Eigen::VectorXi& innerEdges,
                                const hedra::SubdivisionType st
                                Eigen::MatrixXd& newV,
                                Eigen::VectorXi& newD,
                                Eigen::MatrixXi& newF)
  
  {
    using namespace Eigen;
    using namespace std;
    
    switch (st){
      case hedra::LINEAR_SUBDIVISION: vertex_insertion(V, D,F, EV, EF, EFi, FE, innerEdges, &vertex_star_cc_linear, &edge_point_linear_average, &face_point_linear_average, fineV, fineD, fineF);
      case hedra::CANONICAL_MOEBIUS_SUBDIVISION: vertex_insertion(V, D,F, EV, EF, EFi, FE, innerEdges, &vertex_star_cc_moebius, &edge_point_blend_moebius_cc, &face_point_moebius_average, fineV, fineD, fineF);
      default: return false
    }
    
    return true;
  }
  
  
  IGL_INLINE bool vertex_insertion(const Eigen::MatrixXd& V,
                                   const Eigen::VectorXi& D,
                                   const Eigen::MatrixXi& F,
                                   const Eigen::MatrixXi& EV,
                                   const Eigen::MatrixXi& EF,
                                   const Eigen::MatrixXi& EFi,
                                   const Eigen::MatrixXi& FE,
                                   const Eigen::VectorXi& innerEdges,
                                   const std::function<RowVector3d(const MatrixXd&, const MatrixXi&, const MatrixXi&, const MatrixXi&,, const MatrixXi&,, const MatrixXi&,, const VectorXi&)>& vertexStarCandidatePointsFunc,
                                   const std::function<RowVector3d(const MatrixXd&, const MatrixXd&, const VectorXd&, MatrixXd&)>& facePointBlendFunc,
                                   const std::function<RowVector3d(const RowVector3d&, const MatrixXd&, const bool, const RowVector3d&, const double, const MatrixXd&, MatrixXd&)>& edgePointBlendFunc,
                                   Eigen::MatrixXd& fineV,
                                   Eigen::VectorXi& fineD,
                                   Eigen::MatrixXi& fineF)
  {
    
    
    
    Eigen::MatrixXd candidateFacePoints(F.rows(), D.maxCoeff*3);
    Eigen::MatrixXd candidateEdgePoints(EV.rows(), 6);
    
    Eigen::MatrixXd fineVertexPoints(V.rows(),3);
    hedra::dcel(D,F,EV, EF,EFi,innerEdges,VH,EH,FH,HV,HE,HF,nextH,prevH,twinH);
    
    vertexStarCandidatePointsFunc(V,D,F,EV,EF,EFi,innerEdges,VH,EH,FH,HV,HE,HF,nextH,prevH,twinH, fineVertexPoints, candidateEdgePoints, candidateFacePoints);
    facePointBlendFunc(V,D,F, candidateFacePoints, fineFacePoints);
    edgePointBlendFunc(V,D,F,EV, innerEdges, candidateEdgePoints, fineEdgePoints);
    boundaryVertexBlendFunc(V,VH,EH,FH,HV,HE,HF,nextH,prevH,twinH, fineEdgePoints,fineVertexPoints);
    
    cout<<"(fineVertexPoints-V).lpNorm<Infinity>()" <<(fineVertexPoints-V).lpNorm<Infinity>()<<endl;
    
    fineV.conservativeResize(fineVertexPoints.rows()+fineFacePoints.rows()+fineEdgePoints.rows(),3);
    fineV<<fineVertexPoints, fineEdgePoints, fineFacePoints;
    int numNewFaces=D[level].sum();
    fineD=VectorXi::Constant(numNewFaces, 4);
    fineF.conservativeResize(numNewFaces,4);
    int currNewFace=0;
    for (int i=0;i<D.rows();i++){
      for (int j=0;j<D(i);j++){
        F[level+1].row(currNewFace++)<<F[level](i,j),
        V.rows()+FE[level](i,j),
        V.rows()+edgePoints.rows()+i,
        V.rows()+FE[level](i,(j+D[level](i)-1)%D[level](i));
      }
    }
    
    return true;
  }
  
}


#endif


