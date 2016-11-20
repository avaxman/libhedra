// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_POLYHEDRAL_FACE_NORMALS_H
#define HEDRA_POLYHEDRAL_FACE_NORMALS_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
    // Computes normals to (near)-polyhedral and convex faces. the assumption of ``nearness'' is to have a coherent normals for all corners which are similar enough and co-oriented, so that a simple nondegenerate average works.
    // Inputs:
    //  V           eigen double matrix     #V by 3 vertex coordinates
    //  D           eigen int vector        #F by 1 - face degrees
    //  F           eigen int matrix        #F by max(D) - vertex indices in face
    // Outputs:
    //  faceNormals eigen double matrix     #F by 3 face normals
    IGL_INLINE bool polyhedral_face_normals(const Eigen::MatrixXd& V,
                                            const Eigen::VectorXi& D,
                                            const Eigen::MatrixXi& F,
                                            Eigen::MatrixXd& faceNormals)
    {
        using namespace Eigen;
        faceNormals.resize(D.rows(),3);
        for (int i=0;i<D.rows();i++){
            RowVector3d faceNormal; faceNormal<<0.0,0.0,0.0;
            for (int j=0;j<D(i);j++){
                RowVector3d vn=V.row(F(i,(j+D(i)-1)%D(i)));
                RowVector3d v0=V.row(F(i,j));
                RowVector3d v1=V.row(F(i,(j+1)%D(i)));
                if (((v1-v0).cross(vn-v0)).norm()>10e-6)
                    faceNormal=faceNormal+((v1-v0).cross(vn-v0)).normalized();
            }
            
            faceNormals.row(i)=faceNormal.normalized();
        }
        
        return true;
    }
}


#endif


