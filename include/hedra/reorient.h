// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2020 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_REORIENT_H
#define HEDRA_REORIENT_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
  // reorients mesh consistently
  IGL_INLINE void reorient(const Eigen::VectorXi& D,
                           Eigen::MatrixXi& F
                           int seedFace = 0)
  {
    
    VectorXi visited=VectorXi::Zero(D.size());
    std::queue<int> 
    while !
    vertexValences.conservativeResize(EV.maxCoeff()+1);
    vertexValences.setZero();
    for (int i=0;i<EV.rows();i++){
      vertexValences(EV(i,0))++;
      vertexValences(EV(i,1))++;
    }
  }
}


#endif


