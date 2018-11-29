// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef HEDRA_QUATERNIONIC_OPERATIONS
#define HEDRA_QUATERNIONIC_OPERATIONS

#include <iostream>
#include <Eigen/Dense>


inline Eigen::RowVector4d QConj(const Eigen::RowVector4d& q)
{
    Eigen::RowVector4d newq(q.rows(), q.cols());
    newq<<q(0), -q.tail(3);
    return newq;
}

inline Eigen::MatrixXd QConjN(const Eigen::MatrixXd& q)
{
  Eigen::MatrixXd newq(q.rows(), q.cols());
  newq<<q.col(0), -q.block(0,1,q.rows(),3);
  return newq;
}


inline Eigen::RowVector4d QMult(const Eigen::RowVector4d& q1, const Eigen::RowVector4d& q2)
{
    Eigen::RowVector4d newq;
    double r1=q1(0);
    double r2=q2(0);
    Eigen::RowVector3d v1=q1.tail(3);
    Eigen::RowVector3d v2=q2.tail(3);
    newq<<r1*r2-v1.dot(v2), r1*v2+r2*v1+v1.cross(v2);
    return newq;
}


/*{
  Eigen::MatrixXd newq;
  Eigen::VectorXd r1=q1.col(0);
  Eigen::VectorXd r2=q2.col(0);
  Eigen::MatrixXd v1=q1.block(0,1,q1.rows(),3);
  Eigen::MatrixXd v2=q2.block(0,1,q2.rows(),3);
  newq<<r1.array()*r2.array()-(v1.array()*v2.array()).rowwise().sum(), v2.rowwise()*r1.array()+v1.rowwise()*r2.array()+v1.rowwise().cross(v2);
  return newq;
}*/

inline Eigen::MatrixXd QMultN(const Eigen::MatrixXd& q1, const Eigen::MatrixXd& q2)
{
  Eigen::MatrixXd newq(q1.rows(),4);
  Eigen::VectorXd r1=q1.col(0);
  Eigen::VectorXd r2=q2.col(0);
  Eigen::MatrixXd v1=q1.block(0,1,q1.rows(), 3);
  Eigen::MatrixXd v2=q2.block(0,1,q2.rows(), 3);
  Eigen::MatrixXd r1mat; r1mat.resize(r1.rows(),3); r1mat<<r1, r1, r1;
  Eigen::MatrixXd r2mat; r2mat.resize(r2.rows(),3); r2mat<<r2, r2, r2;
  newq.col(0)=r1.cwiseProduct(r2)-(v1.cwiseProduct(v2)).rowwise().sum();
  Eigen::MatrixXd v1cv2; v1cv2.resize(v1.rows(),3);
  for (int i=0;i<v1.rows();i++){
    Eigen::Vector3d vv1=v1.row(i);
    Eigen::Vector3d vv2=v2.row(i);
    v1cv2.row(i)=vv1.cross(vv2);
    
  }
  newq.block(0,1,newq.rows(),newq.cols()-1)=v1.cwiseProduct(r2mat)+v2.cwiseProduct(r1mat)+v1cv2;
  return newq;
}

inline Eigen::RowVector4d QInv(const Eigen::RowVector4d& q)
{
    return(QConj(q)/q.squaredNorm());
}

inline Eigen::MatrixXd QInvN(const Eigen::MatrixXd& q)
{
  return(QConjN(q).cwiseQuotient(q.rowwise().squaredNorm().replicate(1,4)));
}

inline Eigen::MatrixXd QLogN(const Eigen::MatrixXd& q)
{
    Eigen::VectorXd nq=q.rowwise().norm();
    Eigen::VectorXd nv=q.block(0,1,q.rows(),q.cols()-1).rowwise().norm();
    Eigen::VectorXd acosqnq=acos((q.col(0).cwiseQuotient(nq)).array()).matrix().cwiseQuotient(nv);
    Eigen::MatrixXd acosmat(acosqnq.rows(),3); acosmat<<acosqnq, acosqnq, acosqnq;
    Eigen::MatrixXd logq(q.rows(),q.cols());
    logq<<log(nq.array()), q.block(0,1,q.rows(),q.cols()-1).cwiseProduct(acosmat);
    for (int i=0;i<logq.rows();i++)
        if (nv(i)<10e-6)
            logq.row(i)<<log(nq(i)),0.0,0.0,0.0;
    
    //cout<<"logq: "<<logq<<endl;
    return logq;
}

inline Eigen::MatrixXd QExpN(const Eigen::MatrixXd& q)
{
    Eigen::VectorXd nv=q.block(0,1,q.rows(),q.cols()-1).rowwise().norm();
    Eigen::VectorXd exp1=exp(q.col(0).array());
    Eigen::MatrixXd exp1mat(exp1.rows(),4); exp1mat<<exp1,exp1,exp1,exp1;
    Eigen::MatrixXd expq(q.rows(),q.cols());
    Eigen::VectorXd sinnv=(sin(nv.array()).matrix()).cwiseQuotient(nv);
    Eigen::MatrixXd sinnvmat(sinnv.rows(),3); sinnvmat<<sinnv,sinnv,sinnv;
    expq<<cos(nv.array()), q.block(0,1,q.rows(),q.cols()-1).cwiseProduct(sinnvmat);
    expq=expq.cwiseProduct(exp1mat);
    for (int i=0;i<expq.rows();i++)
        if (nv(i)<10e-6)
            expq.row(i)<<exp1(i),0.0,0.0,0.0;
    
    return expq;
}

inline void Quat2Coords(const Eigen::MatrixXd& QV, Eigen::MatrixXd& V)
{
    V.conservativeResize(QV.rows(), 3);
    V=QV.block(0,1,QV.rows(),QV.cols()-1);
}

inline void Coords2Quat(const Eigen::MatrixXd& V, Eigen::MatrixXd& QV)
{
    QV.conservativeResize(V.rows(),4); QV.setZero();
    QV.block(0,1,QV.rows(),QV.cols()-1)=V;
}

inline void quat2SparseMatrix(const Eigen::VectorXi& RowIndices,
                              const Eigen::VectorXi& ColIndices,
                              const Eigen::MatrixXd& Values,
                              Eigen::SparseMatrix<double>& Mat,
                              const int m,
                              const int n)
{
    /*[  rq, -vqx, -vqy, -vqz]
     [ vqx,   rq, -vqz,  vqy]
     [ vqy,  vqz,   rq, -vqx]
     [ vqz, -vqy,  vqx,   rq]*/
    using namespace Eigen;
    using namespace std;
    Mat.resize(4*m,4*n);
    vector<Triplet<double> > RealTris(16*RowIndices.size());
    for (int i=0;i<RowIndices.size();i++){
        
        //r=rq*rp-<vq,vp>
        RealTris[16*i]=Triplet<double>(4*RowIndices(i),      4*ColIndices(i),     Values(i,0));
        RealTris[16*i+1]=Triplet<double>(4*RowIndices(i)+1,    4*ColIndices(i)+1,   Values(i,0));
        RealTris[16*i+2]=Triplet<double>(4*RowIndices(i)+2,    4*ColIndices(i)+2,   Values(i,0));
        RealTris[16*i+3]=Triplet<double>(4*RowIndices(i)+3,    4*ColIndices(i)+3,   Values(i,0));
        
        //v=rq*vp+vp*vq+(vq x vp)
        RealTris[16*i+4]=Triplet<double>(4*RowIndices(i),      4*ColIndices(i)+1,   -Values(i,1));
        RealTris[16*i+5]=Triplet<double>(4*RowIndices(i)+1,    4*ColIndices(i),     Values(i,1));
        RealTris[16*i+6]=Triplet<double>(4*RowIndices(i)+2,    4*ColIndices(i)+3,   -Values(i,1));
        RealTris[16*i+7]=Triplet<double>(4*RowIndices(i)+3,    4*ColIndices(i)+2,   Values(i,1));
        
        RealTris[16*i+8]=Triplet<double>(4*RowIndices(i),      4*ColIndices(i)+2,   -Values(i,2));
        RealTris[16*i+9]=Triplet<double>(4*RowIndices(i)+1,    4*ColIndices(i)+3,   Values(i,2));
        RealTris[16*i+10]=Triplet<double>(4*RowIndices(i)+2,    4*ColIndices(i),     Values(i,2));
        RealTris[16*i+11]=Triplet<double>(4*RowIndices(i)+3,    4*ColIndices(i)+1,   -Values(i,2));
        
        RealTris[16*i+12]=Triplet<double>(4*RowIndices(i),      4*ColIndices(i)+3,   -Values(i,3));
        RealTris[16*i+13]=Triplet<double>(4*RowIndices(i)+1,    4*ColIndices(i)+2,   -Values(i,3));
        RealTris[16*i+14]=Triplet<double>(4*RowIndices(i)+2,    4*ColIndices(i)+1,   Values(i,3));
        RealTris[16*i+15]=Triplet<double>(4*RowIndices(i)+3,    4*ColIndices(i),     Values(i,3));
    }
    
    Mat.setFromTriplets(RealTris.begin(), RealTris.end());
}

inline Eigen::SparseMatrix<double> quatConjMat(int m){
    Eigen::SparseMatrix<double> Mat(4*m,4*m);
    std::vector<Eigen::Triplet<double> > Tris(4*m);
    for (int i=0;i<m;i++){
        Tris[4*i]=Eigen::Triplet<double>(4*i,4*i,1.0);
        Tris[4*i+1]=Eigen::Triplet<double>(4*i+1,4*i+1,-1.0);
        Tris[4*i+2]=Eigen::Triplet<double>(4*i+2,4*i+2,-1.0);
        Tris[4*i+3]=Eigen::Triplet<double>(4*i+3,4*i+3,-1.0);
    }
    
    Mat.setFromTriplets(Tris.begin(), Tris.end());
    return Mat;
    
}



/*inline Eigen::RowVector4d SLERP(const Eigen::RowVector4d& q1, const Eigen::RowVector4d& q2, const double t){
    return QMult(q1, QExp(QLog(QMult(QInv(q1), q2))*t));
}
*/


#endif
