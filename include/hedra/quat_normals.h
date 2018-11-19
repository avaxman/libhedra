// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_QUAT_NORMALS_H
#define HEDRA_QUAT_NORMALS_H
#include <igl/igl_inline.h>
#include <hedra/quaternionic_operations.h>
#include <Eigen/Core>
#include <vector>
#include <cmath> 


namespace hedra
{
  // Computes the quaternionic cross ratio of quadruplets, that are given as imaginary quaternions, representing 3D coordinates.
  
  // Inputs:
  //  V           eigen double matrix     #V by 3 - coordinates
  //  Q           eigen int matrix        #Q by 4 - quadruplets of indices into V
  // Outputs:
  //  cr          eigen double matix      #Q by 4 - quaternion (r,vx,vy,vz) cross ratio per quadruplet.
  IGL_INLINE bool quat_normals(const Eigen::MatrixXd& Vq,
                               const Eigen::MatrixXi& FaceTriads,
                               Eigen::MatrixXd& FN)
  {
    using namespace Eigen;
    FN.resize(FaceTriads.rows(),4);
    for (int i=0;i<FaceTriads.rows();i++){
      RowVector4d qi=Vq.row(FaceTriads(i,0));
      RowVector4d qj=Vq.row(FaceTriads(i,1));
      RowVector4d qk=Vq.row(FaceTriads(i,2));
      
      FN.row(i)=QMult(qj-qi, QInv(qk-qj));
    }
  }
  
}




#endif


