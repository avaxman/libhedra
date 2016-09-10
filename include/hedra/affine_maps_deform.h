// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_AFFINE_MAPS_DEFORM_H
#define HEDRA_AFFINE_MAPS_DEFORM_H

#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <vector>

//This file implements the deformation algorithm described in:

//Amir Vaxman
//Modeling Polyhedral Meshes with Affine Maps
//Computer Graphics Forum (Proc. SGP) 31(5), 2012

namespace hedra
{
    // Deform a polyhedral mesh with a single affine map per face
    
    
    //precomputation of the necessary matrices
    
    //input:
    //  D  eigen int vector     #F by 1 - face degrees
    //  F  eigen int matrix     #F by max(D) - vertex indices in face
    //  EF eigen int matrix     #E by 2 - map from edges to adjacent faces
    //  EV eigen int matrix     #E by 2 - map from edges to end vertices
    //  h eigen int vector     #constraint vertex indices (handles)
    //  bendFactor double       #the relative similarty between affine maps on adjacent faces
  
    // Output:
    // E eigen double sparse matrix    Energy matrix
    // C eigen double sparse matrix    Constraint matrix

    IGL_INLINE void affine_maps_precompute(const Eigen::VectorXi& D,
                                           const Eigen::MatrixXi& F,
                                           const Eigen::MatrixXi& EF,
                                           const Eigen::MatrixXi& EV,
                                           const Eigen::VectorXi& h,
                                           const double bendFactor;
                                           Eigen::sparsematrix<double>& E,
                                           Eigen::sparsematrix<double>& C);
    {
        
    }
    
    
    //COmputing the deformation.
    //Prerequisite: affine_maps_precompute is called, and
    //qh input values are matching those in h.
    
    //input:
    // E eigen double sparse matrix     Energy matrix
    // C eigen double sparse matrix     Constraint matrix
    // qh eigen double matrix           h by 3 new handle positions
    
    //output:
    // q eigen double matrix            V by 3 new vertex positions (note: include handles)
    // A eigen double matrix            F by 9 affine maps (column major)
    
    IGL_INLINE void affine_maps_deform(const Eigen::sparsematrix<double>& E,
                                       const Eigen::sparsematrix<double>& C,
                                       const Eigen::MatrixXd qh,
                                       Eigen::MatrixXd q,
                                       Eigen::MatrixXd A);
    {
        
    }
}


#endif
