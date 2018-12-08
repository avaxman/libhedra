// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_POLYGONAL_LINES_H
#define HEDRA_POLYGONAL_LINES_H

#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <igl/igl_inline.h>
#include <igl/avg_edge_length.h>
#include <hedra/visualization_schemes.h>
#include <hedra/line_cylinders.h>

namespace hedra
{
  
  //creating meshes for polygonal lines, as curves of cylinders
  // Input:
  //  V:      #V X 3 vertex coordinates.
  //  F:      #F X 3 mesh triangles.
  //  EV:     #E x 2 edges to vertices indices
  //  seams:  which edges of EV are seams
  //  width:  the width of the seam lines
  //  res:    the resolution of the cylinder bases.
  // Output:
  //  VSeams:          The vertices of the seam cylinders.
  //  FSeams:          The faces of the seam cylinders.
  //  CSeams:         The colors of the seam cylinders (per face).
  void IGL_INLINE polygonal_edge_lines(const Eigen::MatrixXd &V,
                                       const Eigen::MatrixXi &F,
                                       const Eigen::MatrixXi &EV,
                                       double width,
                                       int res,
                                       Eigen::MatrixXd &VPolyLines,
                                       Eigen::MatrixXi &FPolyLines,
                                       Eigen::MatrixXd &CPolyLines)
  {
    
    Eigen::MatrixXd P1(EV.rows(),3), P2(EV.rows(),3);
    for (int i=0;i<EV.rows();i++){
      P1.row(i)=V.row(EV(i,0));
      P2.row(i)=V.row(EV(i,1));
    }
    
    hedra::line_cylinders(P1, P2, width, hedra::default_edge_color().replicate(P1.rows(),1), res, VPolyLines, FPolyLines, CPolyLines);
  }
  
  
  //A version with default width and resolution.
  void IGL_INLINE polygonal_edge_lines(const Eigen::MatrixXd &V,
                                       const Eigen::MatrixXi &F,
                                       const Eigen::MatrixXi &T,
                                       const Eigen::MatrixXi &EV,
                                       Eigen::MatrixXd &VPolyLines,
                                       Eigen::MatrixXi &FPolyLines,
                                       Eigen::MatrixXd &CPolyLines,
                                       const double lengthRatio=1.25)
  {
    double l = lengthRatio*igl::avg_edge_length(V, T);
    polygonal_edge_lines(V,F,EV, l/25.0,6, VPolyLines, FPolyLines, CPolyLines);
  }
  
}

#endif


