// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_COPYLEFT_CGAL_EXTRACT_MESH_H
#define HEDRA_COPYLEFT_CGAL_EXTRACT_MESH_H
#include <igl/igl_inline.h>

#include <Eigen/Dense>
#include <CGAL/Arrangement_2.h>
#include <vector>

namespace hedra
{
  namespace copyleft
  {
    namespace cgal
    {
      
      IGL_INLINE void extract_mesh(const Eigen::MatrixXd& V,
                                   const Eigen::MatrixXi& F,
                                   const Eigen::MatrixXd& TC,
                                   const Eigen::MatrixXi& FTC,
                                   Eigen::MatrixXd& newV,
                                   Eigen::VectorXi& newD,
                                   Eigen::MatrixXi& newF){
        
      }
    }
  }
}

#endif


