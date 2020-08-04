// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_AVERAGE_ONTO_VERTICES_H
#define HEDRA_AVERAGE_ONTO_VERTICES_H

#include <Eigen/Dense>

namespace hedra
{   
    // average_onto_vertices
    // Move a field defined on faces to vertices by averaging
    //
    // Input:
    // V, D, F: mesh
    // S: field defined on faces
    // 
    // Output:
    // SF: field defined on vertices
    template<typename DerivedV,typename DerivedF,typename DerivedD,typename DerivedS>
    void average_onto_vertices(const Eigen::MatrixBase<DerivedV> &V,
                                      const Eigen::MatrixBase<DerivedF> &F,
                                      const Eigen::MatrixBase<DerivedD> &D,
                                      const Eigen::MatrixBase<DerivedS> &S,
                                      Eigen::MatrixBase<DerivedS> &SV)
    {
        SV = DerivedS::Zero(V.rows(),S.cols());
        Eigen::Matrix<typename DerivedF::Scalar,Eigen::Dynamic,1> COUNT(V.rows());
        COUNT.setZero();
        for (int i = 0; i <F.rows(); ++i) {
            for (int j = 0; j<D(i); ++j) {
                SV.row(F(i,j)) += S.row(i);
                COUNT[F(i,j)] ++;
            }
        }
        for (int i = 0; i <V.rows(); ++i)
            SV.row(i) /= COUNT[i];
    };
}

#endif
