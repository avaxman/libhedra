// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2015 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_POLYGONAL_READ_OFF_H
#define HEDRA_POLYGONAL_READ_OFF_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
    // reads mesh from an ascii OFF file of a polygonal mesh
    // Inputs:
    //   str  path to .off file
    // Outputs:
    //  V  eigen double matrix  #V by 3 - vertex coordinates
    //  D  eigen int vector     #F by 1 - face degrees
    //  F  eigen int matrix     #F by max(D) - vertex indices in face
    IGL_INLINE bool polygonal_read_OFF(const std::string str,
                                       Eigen::MatrixXd& V,
                                       Eigen::VectorXi& D,
                                       Eigen::MatrixXi& F)
    {
        
        using namespace std;
        ifstream FileHandle;
        FileHandle.open(str);
        if (!FileHandle.is_open())
            return false;
        int NumofVertices, NumofFaces, NumofEdges;
        char OFFString[6];
        vector<vector<int> > RawFaces;
        
        FileHandle>>OFFString>>NumofVertices>>NumofFaces>>NumofEdges;
        V.resize(NumofVertices,3);
        RawFaces.resize(NumofFaces);
        D.resize(NumofFaces,1);
        for (int i=0;i<NumofVertices;i++)
            FileHandle>>V(i,0)>>V(i,1)>>V(i,2);
        
        for (int i=0;i<NumofFaces;i++){
            FileHandle>>D(i);
            RawFaces[i].resize(D(i));
            for (int j=0;j<D(i);j++)
                FileHandle>>RawFaces[i][j];
        }
        
        F.resize(NumofFaces,D.maxCoeff());
        F.setConstant(-1);  //to "don't care" vertices
        for (int i=0;i<NumofFaces;i++)
            for (int j=0;j<RawFaces[i].size();j++)
                F(i,j)=RawFaces[i][j];
        
        FileHandle.close();
        return true;
    }
}


#endif


