// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_AVERAGE_ONTO_FACE_H
#define HEDRA_AVERAGE_ONTO_FACE_H

#include <Eigen/Core>


namespace hedra
{
  // average_onto_faces
  // Move a field defined on vertices to faces by averaging
  //
  // Input:
  // V, D, F: mesh
  // S: field defined on vertices
  // 
  // Output:
  // SF: field defined on faces
  template <typename DerivedF, typename DerivedD, typename DerivedS, typename DerivedSF>
  IGL_INLINE void average_onto_faces(
    const Eigen::MatrixBase<DerivedF> &F,
    const Eigen::MatrixBase<DerivedD> &D,
    const Eigen::MatrixBase<DerivedS> &S,
    Eigen::PlainObjectBase<DerivedSF> &SF)
   {
      SF.setConstant(F.rows(),S.cols(),0);
      for (int i = 0; i <F.rows(); ++i)
      {
          for (int j = 0; j<D(i); ++j)
              SF.row(i) += S.row(F(i,j));
          SF.row(i) /= D(i);
      }
   }
}

#endif
