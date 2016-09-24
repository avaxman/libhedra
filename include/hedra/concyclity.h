// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_CONCYCLITY_H
#define HEDRA_CONCYCLITY_H
#include <igl/igl_inline.h>
#include <hedra/quat_cross_ratio.h>
#include <Eigen/Core>
#include <vector>
#include <cmath> 


namespace hedra
{
    // Computes the concyclity of a face, which is pi-argument of the cross-ratio. That is cos(angle)=-Re(cr)/norm(cr)
    
    // Inputs:
    //  V           eigen double matrix     #V by 3 - mesh coordinates
    //  D           eigen int matrix        #F by 1 - face degrees
    //  F           eigen int matrix        #F by max(D)
    // Outputs:
    //  planarity   eigen double matix      #F by 1
    IGL_INLINE bool concyclity(const Eigen::MatrixXd& V,
                               const Eigen::VectorXi& D,
                               const Eigen::MatrixXi& F,
                               Eigen::VectorXd& concyclity)
    {
        using namespace Eigen;
        concyclity.resize(D.size());
        
        MatrixXi quadIndices(D.sum(),4);
        int currIndex=0;
        for (int i=0;i<D.size();i++)
            for (int j=0;j<D(i);j++)
                quadIndices.row(currIndex++)<<F(i,j),F(i,(j+1)%D(i)), F(i,(j+2)%D(i)), F(i,(j+3)%D(i));
        
        MatrixXd cr;
        quat_cross_ratio(V,quadIndices, cr);
        VectorXd realPart=cr.col(0);
        VectorXd absPart=cr.rowwise().norm();
        VectorXd crAngles=acos(-realPart.cwiseQuotient(absPart).array());
    
        currIndex=0;
        for (int i=0;i<D.rows();i++){
            if (D(i)>=4)
                concyclity(i)=sqrt((crAngles.segment(currIndex, D(i)).squaredNorm()/(double)D(i)));
            else
                concyclity(i)=0.0;  //TODO: avoid the computation altogether
            currIndex+=D(i);
        }
        //converting to degrees
        concyclity.array()*=180.0/M_PI;
        return true;
    }
}

    


#endif


