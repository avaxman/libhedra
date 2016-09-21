// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef HEDRA_QUATERNION_OPS
#define HEDRA_QUATERNION_OPS

#include <iostream>
#include <Eigen/Dense>


inline RowVector4d QConj(const RowVector4d& q)
{
    RowVector4d newq(q.rows(), q.cols());
    newq<<q(0), -q.tail(3);
    return newq;
}


inline RowVector4d QMult1(const RowVector4d& q1, RowVector4d& q2)
{
    RowVector4d newq;
    double r1=q1(0);
    double r2=q2(0);
    RowVector3d v1=q1.tail(3);
    RowVector3d v2=q2.tail(3);
    newq<<r1*r2-v1.dot(v2), r1*v2+r2*v1+v1.cross(v2);
    return newq;
}

inline RowVector3d QInv(const RowVector3d& q)
{
    return(QConj(q)/q.squaredNorm());
}


inline MatrixXd QLog(const MatrixXd& q)
{
    VectorXd nq=q.rowwise().norm();
    VectorXd nv=q.block(0,1,q.rows(),q.cols()-1).rowwise().norm();
    VectorXd acosqnq=acos((q.col(0).cwiseQuotient(nq)).array()).matrix().cwiseQuotient(nv);
    MatrixXd acosmat(acosqnq.rows(),3); acosmat<<acosqnq, acosqnq, acosqnq;
    MatrixXd logq(q.rows(),q.cols());
    logq<<log(nq.array()), q.block(0,1,q.rows(),q.cols()-1).cwiseProduct(acosmat);
    for (int i=0;i<logq.rows();i++)
        if (nv(i)<10e-6)
            logq.row(i)<<log(nq(i)),0.0,0.0,0.0;
    
    //cout<<"logq: "<<logq<<endl;
    return logq;
}

inline MatrixXd QExp(const MatrixXd& q)
{
    VectorXd nv=q.block(0,1,q.rows(),q.cols()-1).rowwise().norm();
    VectorXd exp1=exp(q.col(0).array());
    MatrixXd exp1mat(exp1.rows(),4); exp1mat<<exp1,exp1,exp1,exp1;
    MatrixXd expq(q.rows(),q.cols());
    VectorXd sinnv=(sin(nv.array()).matrix()).cwiseQuotient(nv);
    MatrixXd sinnvmat(sinnv.rows(),3); sinnvmat<<sinnv,sinnv,sinnv;
    expq<<cos(nv.array()), q.block(0,1,q.rows(),q.cols()-1).cwiseProduct(sinnvmat);
    expq=expq.cwiseProduct(exp1mat);
    for (int i=0;i<expq.rows();i++)
        if (nv(i)<10e-6)
            expq.row(i)<<exp1(i),0.0,0.0,0.0;
    
    return expq;
}


inline RowVector4d SLERP(const RowVector4d& q1, const RowVector4d& q2, const double t){
    return QMult(q1, QExp(QLog(QMult(QInv(q1), q2))*t));
}



#endif
