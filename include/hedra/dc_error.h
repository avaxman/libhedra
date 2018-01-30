// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_DC_ERROR_H
#define HEDRA_DC_ERROR_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <hedra/quat_cross_ratio.h>
#include <vector>
#include <cmath> 


namespace hedra
{
  // Computes the discrete conformal error between meshes on a set of predefined four vertices. The error is computed as abs(log(|crt|-log(|crs|)) (crt - target cross ratio of four points, and crs - source cross ratio.
  
  //Both meshes must agree on combinatorical structure.
  
  
  // Inputs:
  //  Vs           eigen double matrix     #V by 3 - source mesh coordinates
  //  Vt           eigen double matrix     #V by 3 - target mesh coordinates
  //  D           eigen int matrix        #F by 1 - face degrees
  //  F           eigen int matrix        #F by max(D) - face indices
  //  I4          eigen int matrix        #I by 4 - representing four indices into V of vertices on which the cross ratios are measured
  // Outputs:
  //  dcError     eigen double vector      #EV by 1 -
  IGL_INLINE bool dc_error(const Eigen::MatrixXd& Vs,
                           const Eigen::MatrixXd& Vt,
                           const Eigen::VectorXi& D,
                           const Eigen::MatrixXi& F,
                           const Eigen::MatrixXi& EV,
                           const Eigen::MatrixXi& EF,
                           const Eigen::MatrixXi& EFi,
                           const Eigen::VectorXi& innerEdges,
                           Eigen::VectorXd& dcError)
  {
    using namespace Eigen;
    
    MatrixXi I4(innerEdges.rows(),4);
    for (int i=0;i<innerEdges.rows();i++){
      int f=EF(innerEdges(i),0);
      int g=EF(innerEdges(i),1);
      int vi=EV(innerEdges(i),0);
      int vk=EV(innerEdges(i),1);
      int vj=F(g,(EFi(innerEdges(i),1)+2)%D(g));
      int vl=F(f,(EFi(innerEdges(i),0)+2)%D(f));
      I4.row(i)<<vi,vj,vk,vl;
    }
    
    dcError.conservativeResize(EV.rows());
    
    MatrixXd crs, crt;
    quat_cross_ratio(Vs,I4, crs);
    quat_cross_ratio(Vt,I4, crt);
  
    VectorXd lcrs=crs.rowwise().norm();
    VectorXd lcrt=crt.rowwise().norm();
    
    dcError.conservativeResize(EV.rows()); dcError.setZero();
    for (int i=0;i<innerEdges.rows();i++)
      dcError(innerEdges(i))=abs(log(lcrs(i))-log(lcrt(i)));
    
    return true;
  }
}






#endif


