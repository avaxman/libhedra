// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_MOEBIUS_REFINEMENT_H
#define HEDRA_MOEBIUS_REFINEMENT_H
#include <igl/igl_inline.h>
#include <hedra/quaternionic_operations.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
 
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
    
    /* std::cout<<"qa: "<<qa<<std::endl;
     std::cout<<"qb: "<<qb<<std::endl;
     std::cout<<"qc: "<<qc<<std::endl;
     std::cout<<"qd: "<<qd<<std::endl;*/
    
    MatrixXd cabd=QMultN(QMultN(qa-qc, QInvN(qb-qa)), QMultN((qd-qb), QInvN(qc-qd)));
    
    //std::cout<<"cabd: "<<cabd<<std::endl;
    
    RowVector4d minusOneQuat; minusOneQuat<<-1.0,0.0,0.0,0.0;
    
    MatrixXd cabp=-QExpN(QLogN(cabd)*0.5)*((1.0-wb)/wb);
    //std::cout<<"cabp: "<<cabp<<std::endl;
    for (int i=0;i<va.rows();i++)
      if ((cabd.row(i)-minusOneQuat).squaredNorm()<10e-4){  //-1,0,0,0 square root not well defined
        RowVector3d faceNormal=QMult(qa.row(i)-qc.row(i), QInv(qb.row(i)-qa.row(i))).tail(3).normalized();
        cabp.row(i)<<0,faceNormal;
      }
    
    MatrixXd Unit=MatrixXd::Zero(qa.rows(),4); Unit.col(0).setOnes();
    
    MatrixXd p=QMultN(QInvN(QMultN(QMultN(qb-qa, QInvN(qa-qc)), cabp)+Unit),QMultN(QMultN(qb-qa, QInvN(qa-qc)),QMultN(cabp,qc))+qb);
    
    //MatrixXd p1=QMultN(QInvN(QMultN(QMultN(qb-qa, QInvN(qa-qc)), cabp)+Unit),QMultN(QMultN(qb-qa, QInvN(qa-qc)),QMultN(cabp,qc))+qb);
    //MatrixXd p2=QMultN(QInvN(QMultN(QMultN(qb-qa, QInvN(qa-qc)), -cabp)+Unit),QMultN(QMultN(qb-qa, QInvN(qa-qc)),QMultN(-cabp,qc))+qb);
    
    //std::cout<<"p: "<<p<<std::endl;
    
    for (int i=0;i<va.rows();i++){
      if (normab(i)<10e-14)
        p.row(i)<<0.0,va.row(i);
      else if (normcd(i)<10e-14)
        p.row(i)<<0.0,vd.row(i);
      /*else {
       RowVector4d midEdge=qb.row(i)*wb+qc.row(i)*(1.0-wb);
       double dist1=(p1.row(i)-midEdge).squaredNorm();
       double dist2=(p2.row(i)-midEdge).squaredNorm();
       
       if (dist1>dist2)
       p.row(i)=p2.row(i);
       else
       p.row(i)=p1.row(i);
       }*/
      
    }
    return p.block(0,1,p.rows(),3);
  }
  
}


#endif


