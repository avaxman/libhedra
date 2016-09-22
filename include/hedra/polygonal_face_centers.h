// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_POLYGONAL_FACE_CENTERS_H
#define HEDRA_POLYGONAL_FACE_CENTERS_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
    // Computes barycenter of polygonal faces.
    // Inputs:
    //  V  eigen double matrix  #v by 3 vertex coordinates
    //  D  eigen int vector     #F by 1 - face degrees
    //  F  eigen int matrix     #F by max(D) - vertex indices in face
    // Outputs:
    //  faceCenters eigen double matrix #F by 3 face barycenter coordinates
    IGL_INLINE bool polygonal_face_centers(const Eigen::MatrixXd& V,
                                           const Eigen::VectorXi& D,
                                           const Eigen::MatrixXi& F,
                                           Eigen::MatrixXd& faceCenters)
    {
        using namespace Eigen;
        faceCenters=MatrixXd::Zero(F.rows(),3);
        for (int i=0;i<D.rows();i++){
            for (int j=0;j<D(i);j++)
                faceCenters.row(i)+=V.row(F(i,j));
                
            faceCenters.row(i)/=(double)D(i);
        }
        
        return true;
    }
}


#endif


