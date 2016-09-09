// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_TRIANGULATE_MESH_H
#define HEDRA_TRIANGULATE_MESH_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
    // returns a triangulated version of a polygonal mesh without adding steiner points (only diagonals)
    //TODO: smarter diagonalization with no self intersections.
    // Inputs:
    //  D  eigen int vector     #F by 1 - face degrees
    //  F  eigen int matrix     #F by max(D) - vertex indices in face
    // Outputs:
    //  T  eigen int matrix     #T by 3 - triangles (T is actually sum(D)-rows(D))
    //  TF  eigen int vector    #T by 1 - the original polygonal face in F for each triangle
    IGL_INLINE bool triangulate_mesh(const Eigen::VectorXi& D,
                                     const Eigen::MatrixXi& F,
                                     Eigen::MatrixXi& T,
                                     Eigen::VectorXi& TF)
    {
        using namespace std;
        vector<Eigen::Vector3i> NewTriangles;
        vector<int> RawTF;
        
        for (int i=0;i<D.rows();i++){
        //triangulating the face greedily
            for (int CurrIndex=1;CurrIndex<D(i)-1;CurrIndex++){
                Eigen::Vector3i NewFace;
                NewFace<<F(i,0),F(i,CurrIndex),F(i,CurrIndex+1);
                RawTF.push_back(i);
                NewTriangles.push_back(NewFace);
            }
        }
        
        T.resize(NewTriangles.size(),3);
        TF.resize(RawTF.size());
        for (int i=0;i<NewTriangles.size();i++){
            T.row(i)=NewTriangles[i];
            TF(i)=RawTF[i];
        }
        
        return true;
    }
}


#endif


