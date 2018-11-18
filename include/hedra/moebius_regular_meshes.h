// This file is part of libhedra, a library for polygonal mesh processing
//
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_MOEBIUS_REGULAR_MESHES_H
#define HEDRA_MOEBIUS_REGULAR_MESHES_H

#include <hedra/CeresMRSolver.h>

namespace hedra
{
  
  struct MoebiusRegularData{
    
    double FullCRFactor;
    double FullFNFactor;
    double UnitFactor;
    double HFactor;
    double VecCRFactor;
    double LengthCRFactor;
    double LengthFNFactor;
    double PhaseCRFactor;
    double PhaseFNFactor;
    
    Eigen::MatrixXi QuadVertexIndices;  //a list of all four consecutive vertices, corresponding to edge-mesh, that on which the cross ratio is defined
    Eigen::MatrixXi QuadFaceIndices;  //in case of non-triangular faces that are optimized to be cocircular
    Eigen::MatrixXi FaceTriads;  //a list of face-consecutive triads on which corner normals are defined
    
    Eigen::VectorXi D;  //face degrees
    Eigen::MatrixXi F;
    Eigen::MatrixXi EV,EF,FE, EFi;
    Eigen::MatrixXi ExtEV;  //extended with diagonals
    Eigen::MatrixXd FEs;
    Eigen::VectorXi innerEdges;  //indices of edges which are non-boundary
    Eigen::VectorXi vertexValences;
    Eigen::VectorXi boundaryVertices;
    Eigen::VectorXi boundaryMask;
    Eigen::MatrixXi cornerF;  //enumerating corners for FN
    
    Eigen::MatrixXi cornerPairs;  //neighboring corner pairs (f, if, g, ig, v)
    Eigen::MatrixXi oneRings;  //entries into corner pairs of one rings in progression (to compute willmore energy) #V x max(valence)
    
    Eigen::MatrixXi oneRingVertices;  //the one rings of every single vertex
    
    //the necessary prescription intrinsics for the mesh
    Eigen::VectorXd patternCRLengths;
    Eigen::VectorXd patternCRAngles;
    Eigen::VectorXd patternFNLengths;
    Eigen::VectorXd patternFNAngles;
    Eigen::VectorXd patternFaceCRLengths;
    Eigen::VectorXd patternFaceCRAngles;
    
    Eigen::VectorXd prescribedLengths;  //the lengths that the mesh should have. Initialized to the original lengths, but can be loaded.
    Eigen::VectorXd presFNLengths;
    Eigen::VectorXd presFNAngles;  //the resulting FN lengthd and angles
    
    Eigen::MatrixXd VOrig, QOrig;    //original vertex positions
    
    //variables for minimality
    Eigen::MatrixXd VCR;  //#V - variable ideal cross ratio vectors computed by the solver
    Eigen::MatrixXd FN;   //#F - variable ideal normals computed by the solver
    
    //of positional handles
    Eigen::VectorXi constIndices;
    Eigen::MatrixXd quatConstPoses;
    Eigen::VectorXi constMask;
    
    //edge-based cross ratios
    Eigen::MatrixXd origECR;
    Eigen::MatrixXd deformECR;
    
    //corner based normals
    Eigen::MatrixXd origCFN;
    Eigen::MatrixXd deformCFN;
    
    //(per-vertex) Moebius regularity
    Eigen::VectorXd origMR;
    Eigen::VectorXd deformMR;
    
    //(per-face) Euclidean regular energy
    Eigen::VectorXd origER;
    Eigen::VectorXd deformER;
    
    //Per vertex general Willmore energy
    Eigen::VectorXd origW;
    Eigen::VectorXd deformW;
    
    Eigen::VectorXd convErrors; //last process convergence errors
    
    //optimization operators
    CeresMRSolver CSolver;
    
  };
  
  IGL_INLINE bool setup_moebius_regular(const Eigen::MatrixXd& VOrig,
                                        const Eigen::VectorXi& D,
                                        const Eigen::MatrixXi& F,
                                        const Eigen::MatrixXi& EV,
                                        const Eigen::MatrixXi& FE,
                                        const Eigen::MatrixXi& EF,
                                        const Eigen::MatrixXi& EFi,
                                        const Eigen::MatrixXd& FEs,
                                        const Eigen::VectorXi& innerEdges,
                                        const Eigen::VectorXi& constIndices,
                                        MoebiusRegularData& MRData){
    
  }
  
  
  IGL_INLINE bool compute_moebius_regular(const MoebiusRegularData& MRData,
                                          const double MRCoeff,
                                          const double ERCoeff,
                                          const Eigen::MatrixXd& constPoses,
                                          Eigen::MatrixXd& VRegular)
  {
    
  }
  
}
#endif
