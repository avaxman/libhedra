// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef HEDRA_QUATERNIONIC_DERIVATIVES_H
#define HEDRA_QUATERNIONIC_DERIVATIVES_H
#include <hedra/quaternionic_operations.h>
#include <Eigen/Core>
#include <string>
#include <vector>


namespace hedra {
    
    //deriving an expression a*X*b (or a*conj(X)*b) by X.
    
    inline void quatDerivativeIndices(Eigen::VectorXi& Rows,
                                      Eigen::VectorXi& Cols,
                                      const int CurrTriPos,
                                      const Eigen::Vector4i TriSkips,
                                      const int Row,
                                      const int Col)
    {
        for (int i=0;i<4;i++)
            for (int j=0;j<4;j++){
                Rows(CurrTriPos+TriSkips(i)+j)=Row+i;
                Cols(CurrTriPos+TriSkips(i)+j)=Col+j;
            }
    }
    
    
    inline void quatDerivativeValues(Eigen::VectorXd& Values,
                                     const int CurrTriPos,
                                     const Eigen::Vector4i TriSkips,
                                     const Eigen::RowVector4d& LeftCoeff,
                                     const Eigen::RowVector4d& RightCoeff,
                                     const bool isConj,
                                     const bool add)
    {
        //[  ra, -vax, -vay, -vaz]
        //[ vax,   ra, -vaz,  vay]
        //[ vay,  vaz,   ra, -vax]
        //[ vaz, -vay,  vax,   ra]
        
        Eigen::Matrix4d a; a<<LeftCoeff(0), -LeftCoeff(1), -LeftCoeff(2), -LeftCoeff(3),
        LeftCoeff(1),  LeftCoeff(0), -LeftCoeff(3), LeftCoeff(2),
        LeftCoeff(2),  LeftCoeff(3), LeftCoeff(0), -LeftCoeff(1),
        LeftCoeff(3), -LeftCoeff(2), LeftCoeff(1), LeftCoeff(0);
        
        //cout<<"a: "<<a<<endl;
        
        //[  rb, -vbx, -vby, -vbz]
        //[ vbx,   rb,  vbz, -vby]
        //[ vby, -vbz,   rb,  vbx]
        //[ vbz,  vby, -vbx,   rb]
        
        Eigen::Matrix4d b; b<<RightCoeff(0), -RightCoeff(1), -RightCoeff(2), -RightCoeff(3),
        RightCoeff(1),  RightCoeff(0),  RightCoeff(3), -RightCoeff(2),
        RightCoeff(2), -RightCoeff(3),  RightCoeff(0),  RightCoeff(1),
        RightCoeff(3),  RightCoeff(2), -RightCoeff(1),  RightCoeff(0);
        
        //cout<<"b: "<<b<<endl;
        
        if (!isConj){
            
            //not conjugate
            //[ ra*rb - vax*vbx - vay*vby - vaz*vbz, vay*vbz - rb*vax - ra*vbx - vaz*vby, vaz*vbx - rb*vay - vax*vbz - ra*vby, vax*vby - rb*vaz - ra*vbz - vay*vbx]
            //[ ra*vbx + rb*vax + vay*vbz - vaz*vby, ra*rb - vax*vbx + vay*vby + vaz*vbz, ra*vbz - rb*vaz - vax*vby - vay*vbx, rb*vay - ra*vby - vax*vbz - vaz*vbx]
            //[ ra*vby + rb*vay - vax*vbz + vaz*vbx, rb*vaz - ra*vbz - vax*vby - vay*vbx, ra*rb + vax*vbx - vay*vby + vaz*vbz, ra*vbx - rb*vax - vay*vbz - vaz*vby]
            //[ ra*vbz + rb*vaz + vax*vby - vay*vbx, ra*vby - rb*vay - vax*vbz - vaz*vbx, rb*vax - ra*vbx - vay*vbz - vaz*vby, ra*rb + vax*vbx + vay*vby - vaz*vbz]
            
            
            Eigen::Matrix4d ab=a*b;
            
            //cout<<"ab: "<<ab<<endl;
            
            for (int i=0;i<4;i++)
                for (int j=0;j<4;j++)
                    if (add)
                        Values(CurrTriPos+TriSkips(i)+j)+=ab(i,j);
                    else
                        Values(CurrTriPos+TriSkips(i)+j)=ab(i,j);
            
            
        } else {
            
            Eigen::Matrix4d ConjMat; ConjMat<<1.0,0.0,0.0,0.0,
            0.0,-1.0,0.0,0.0,
            0.0,0.0,-1.0,0.0,
            0.0,0.0,0.0,-1.0;
            Eigen::Matrix4d abc=a*b*ConjMat;
            
            //cout<<"abc:"<<abc<<endl;
            
            for (int i=0;i<4;i++)
                for (int j=0;j<4;j++)
                    if (add)
                        Values(CurrTriPos+TriSkips(i)+j)+=abc(i,j);
                    else
                        Values(CurrTriPos+TriSkips(i)+j)=abc(i,j);
            
        }
    }
}



#endif /* HEDRA_QUATERNIONIC_DERIVATIVES_H */
