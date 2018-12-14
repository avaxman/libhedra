// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_VERTEX_VALENCES_H
#define HEDRA_VERTEX_VALENCES_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
  // computes vertex valences
  IGL_INLINE void vertex_valences(const Eigen::MatrixXi& EV,
                                  Eigen::VectorXi& vertexValences)
  {
    vertexValences.conservativeResize(EV.maxCoeff()+1);
    vertexValences.setZero();
    for (int i=0;i<EV.rows();i++){
      vertexValences(EV(i,0))++;
      vertexValences(EV(i,1))++;
    }
  }
}


#endif


