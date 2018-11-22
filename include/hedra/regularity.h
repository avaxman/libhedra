// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_REGULARITY_H
#define HEDRA_REGULARITY_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <vector>
#include <cmath> 


namespace hedra
{
    // Computes the regularity of every polygonal face. This is done by taking the coefficient of variation of both face angles and lengths independently, and returning their RMSE in percentage. The coefficient of variation is the (std dev)/mean, so the errors are normalized.
    //the faces are not assumed to be polyhedral (planar), but there is a mild assumption that nearby corner normals are well-defined to have angles less than PI between them.
    // Inputs:
    //  V           eigen double matrix     #V by 3 - mesh coordinates
    //  D           eigen int matrix        #F by 1 - face degrees
    //  F           eigen int matrix        #F by max(D)
    // Outputs:
    //  regularity   eigen double matix      #F by 1
    IGL_INLINE bool regularity(const Eigen::MatrixXd& V,
                              const Eigen::VectorXi& D,
                              const Eigen::MatrixXi& F,
                              Eigen::VectorXd& regularity)
    {
        using namespace Eigen;
        regularity.resize(D.size());
        
        for (int i=0;i<D.rows();i++){
            VectorXd lengths(D(i));
            VectorXd angles(D(i));
            //finding the minimal-coordinate vertex which is convex by definition and taking its normal as seed.
            
            Eigen::MatrixXd VFace(D(i),3);
            for (int j=0;j<D(i);j++)
                VFace.row(j)=V.row(F(i,j));
            
            int startVertex;
            for (int dim=0;dim<3;dim++){
                if (VFace.col(dim).maxCoeff(&startVertex)-VFace.col(dim).minCoeff(&startVertex)>10e-4)
                    break;
            }
            //TODO: assert if non has breaked, which means degenerate polygon
            RowVector3d v1=VFace.row((startVertex+D(i)-1)%D(i));
            RowVector3d v2=VFace.row(startVertex);
            RowVector3d v3=VFace.row((startVertex+1)%D(i));
            RowVector3d prevNormal=((v3-v2).cross(v1-v2)).normalized();
            //TODO: what if this is collinear?
            
            int currVertex=startVertex;
            int vecIndex=0;
            do{
                RowVector3d v1=VFace.row((currVertex+D(i)-1)%D(i));
                RowVector3d v2=VFace.row(currVertex);
                RowVector3d v3=VFace.row((currVertex+1)%D(i));
                RowVector3d currNormal=prevNormal;
                if (((v3-v2).cross(v1-v2)).norm()>10e-6)
                    currNormal=((v3-v2).cross(v1-v2)).normalized();
                
                //signing normal according to prevnormal (this is where the "norma smoothness" assumption comes in)
                currNormal*=(currNormal.dot(prevNormal)>0 ? 1.0 : -1.0);
                double sinangle=((v3-v2).cross(v1-v2)).dot(currNormal);
                double cosangle=(v3-v2).dot(v1-v2);
                lengths(vecIndex)=(v2-v1).norm();
                angles(vecIndex++)=atan2(sinangle, cosangle);
                currVertex=(currVertex+1)%D(i);
            }while(currVertex!=startVertex);
            
            double meanl=lengths.mean();
            double meana=angles.mean();
            double stddevl=((lengths.array()-lengths.mean())/(double)lengths.size()).matrix().norm();
            double stddeva=((angles.array()-angles.mean())/(double)angles.size()).matrix().norm();
            double covl=stddevl/meanl;
            double cova=stddeva/meana;
            regularity(i)=100.0*sqrt((covl*covl+cova*cova)/2);
        }
        return true;
    }
}

    


#endif


