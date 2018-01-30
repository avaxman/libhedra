// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_PLANARITY_H
#define HEDRA_PLANARITY_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <vector>
#include <cmath> 


namespace hedra
{
    // Computes the planarity of every polygonal face. This is the RMSE of a vector where each element is the planarity of a consecutive four vertices (for quads it's just the planarity then.
    
    // Inputs:
    //  V           eigen double matrix     #V by 3 - mesh coordinates
    //  D           eigen int matrix        #F by 1 - face degrees
    //  F           eigen int matrix        #F by max(D)
    // Outputs:
    //  planarity   eigen double matix      #F by 1
    IGL_INLINE bool planarity(const Eigen::MatrixXd& V,
                              const Eigen::VectorXi& D,
                              const Eigen::MatrixXi& F,
                              Eigen::VectorXd& planarity)
    {
        using namespace Eigen;
        planarity.resize(D.size());
        
        for (int i=0;i<D.rows();i++){
            Eigen::VectorXd quadPlanarities(D(i));
            for (int j=0;j<D(i);j++){
                RowVector3d v1=V.row(F(i,j));
                RowVector3d v2=V.row(F(i,(j+1)%D(i)));
                RowVector3d v3=V.row(F(i,(j+2)%D(i)));
                RowVector3d v4=V.row(F(i,(j+3)%D(i)));
                RowVector3d diagCross=(v3-v1).cross(v4-v2);
                double denom = diagCross.norm()*(((v3-v1).norm()+(v4-v2).norm())/2);
                if (fabs(denom)<1e-8)
                    //degenerate quad is still planar
                    quadPlanarities(j) = 0;
                else
                    quadPlanarities(j) = (diagCross.dot(v2-v1)/denom);  //percentage
            }
            planarity(i)=100.0*sqrt(quadPlanarities.squaredNorm()/(double)D(i));
        }
        return true;
    }
}

    


#endif


