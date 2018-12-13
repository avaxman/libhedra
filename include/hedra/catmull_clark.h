// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_CATMULL_CLARK_H
#define HEDRA_CATMULL_CLARK_H
#include <igl/igl_inline.h>
#include <hedra/polygonal_face_centers.h>
#include <hedra/dcel.h>
#include <hedra/vertex_valences.h>
#include <hedra/vertex_insertion.h>
#include <hedra/subdivision_basics.h>
#include <hedra/linear_cc_subdivision.h>
#include <hedra/moebius_cc_subdivision.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{

  
  IGL_INLINE bool catmull_clark(const Eigen::MatrixXd& V,
                                const Eigen::VectorXi& D,
                                const Eigen::MatrixXi& F,
                                const int& st,
                                Eigen::MatrixXd& fineV,
                                Eigen::VectorXi& fineD,
                                Eigen::MatrixXi& fineF)
  
  {
    using namespace Eigen;
    using namespace std;
    
    switch (st){
      case hedra::LINEAR_SUBDIVISION: {
        hedra::LinearCCSubdivisionData lsd;
        vertex_insertion(V, D,F,lsd, fineV, fineD, fineF);
        break;
      }
      case hedra::CANONICAL_MOEBIUS_SUBDIVISION: {
        hedra::MoebiusCCSubdivisionData msd;
        vertex_insertion(V, D,F,msd, fineV, fineD, fineF);
        break;
      }
      default: return false;
    }
    
    return true;
  }
  
  
}


#endif


