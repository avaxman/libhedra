// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_QUAT_MOEBIUS_DEFORM_H
#define HEDRA_QUAT_MOEBIUS_DEFORM_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <hedra/Moebius3DCornerVarsTraits.h>
#include <hedra/EigenSolverWrapper.h>
#include <hedra/LMSolver.h>
#include <hedra/quaternionic_operations.h>

using namespace Eigen;
using namespace std;


namespace hedra{
  
  
  struct QuatMoebiusData{
    Eigen::MatrixXi T;     //triangulated mesh (only for comparisons and visualization)
    Eigen::VectorXi TF;
    Eigen::MatrixXi edgeT;  //the edge mesh for MC/IAP error etc. visualization
    
    Eigen::MatrixXi F, edgeTE;
    Eigen::VectorXi D;
    Eigen::MatrixXi extEV, EV, EF, FE, EFi;
    Eigen::MatrixXd FEs;
    Eigen::VectorXi innerEdges;
    
    bool isExactDC;
    bool isExactIAP;
    double rigidityFactor;
    
    //Raw 3D positions
    MatrixXd origV;    //original vertex positions
    MatrixXd deformV;  //current deformed mesh
    
    //Complex Positions
    MatrixXd origVq;
    MatrixXd deformVq;
    
    //corner variables matrices
    MatrixXd deformX;
    VectorXi cornerOffset;  //where does every face begin in the corner list
    
    MatrixXi faceCornerPairs;   //faces (f,g, i,k) and vertices around every edge for compatibility and smoothness, like in the paper
    MatrixXi cornerPairs;  //neighboring corner pairs across edges for AMAP/MC/IAP comparisons
    
    //for inversion control and other usages
    //SparseMatrix<Complex> d0;
    
    //of positional handles
    VectorXi constIndices;
    MatrixXd constPoses;
    
    //optimization operators
    hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > > deformLinearSolver;
    hedra::optimization::Moebius3DCornerVarsTraits deformTraits;
    hedra::optimization::LMSolver<hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > >,hedra::optimization::Moebius3DCornerVarsTraits> deformSolver;
    
  };
  

  IGL_INLINE void quat_moebius_setup(const Eigen::MatrixXd& V,
                                        const Eigen::VectorXi& D,
                                        const Eigen::MatrixXi& F,
                                        const Eigen::MatrixXi& TF,
                                        const Eigen::MatrixXi& EV,
                                        const Eigen::MatrixXi& EF,
                                        const Eigen::MatrixXi& EFi,
                                        const Eigen::MatrixXi& FE,
                                        const Eigen::MatrixXd& FEs,
                                        const Eigen::VectorXi& innerEdges,
                                        struct QuatMoebiusData& mdata)
  {
    mdata.F=F;
    mdata.D=D;
    mdata.origV=V;
    mdata.deformV=V;
    mdata.TF=TF;
    mdata.EV=EV;
    mdata.EF=EF;
    mdata.EFi=EFi;
    mdata.FE=FE;
    mdata.FEs=FEs;
    mdata.innerEdges=innerEdges;
    
    //complex representations
    Coords2Quat(mdata.origV,mdata.origVq);
    mdata.deformVq=mdata.origVq;
    
    //creating full edge list
    vector<pair<int,int>>  diagonals;
    for (int i=0;i<F.rows();i++)
      for (int j=0;j<D(i);j++)
        for (int k=j+2;k<D(i);k++)
          diagonals.push_back(pair<int,int>(F(i,j),F(i,k)));
    
    mdata.extEV.resize(EV.rows()+diagonals.size(),2);
    mdata.extEV.block(0,0,EV.rows(),2)=EV;
    for (int i=0;i<diagonals.size();i++)
      mdata.extEV.row(EV.rows()+i)<<diagonals[i].first, diagonals[i].second;
    
    //reseting deformation result
    mdata.deformX.resize(D.sum(),4);
    mdata.deformX.setZero();
    mdata.deformX.col(0).setConstant(1.0);
    mdata.cornerOffset.resize(F.rows());
    mdata.cornerOffset(0)=0;
    
    for (int i=1;i<D.rows();i++)
      mdata.cornerOffset(i)=mdata.cornerOffset(i-1)+D(i-1);
    
    MatrixXi adjCorners(V.rows(),12);
    VectorXi valences(V.rows()); valences.setZero();
    for (int i=0;i<D.rows();i++){
      for (int j=0;j<D(i);j++){
        adjCorners(F(i,j),valences(F(i,j)))=mdata.cornerOffset(i)+j;
        valences(F(i,j))++;
      }
    }
    
    vector<pair<int,int> > cornerPairList;
    for (int i=0;i<V.rows();i++)
      for (int j=0;j<valences(i)-1;j++)
        cornerPairList.push_back(pair<int,int>(adjCorners(i,j),adjCorners(i,(j+1)%valences(i))));
    
    mdata.cornerPairs.resize(cornerPairList.size(),2);
    for (int i=0;i<cornerPairList.size();i++)
      mdata.cornerPairs.row(i)<<cornerPairList[i].first, cornerPairList[i].second;
    
  }
  
  IGL_INLINE void quat_moebius_precompute(const Eigen::VectorXi& h,
                                          const bool isExactDC,
                                          const bool isExactIAP,
                                          const double rigidityFactor,
                                          struct QuatMoebiusData& mdata){
    
    mdata.constIndices=h;
    mdata.isExactDC=isExactDC;
    mdata.isExactIAP=isExactIAP;
    mdata.rigidityFactor=rigidityFactor;
    mdata.deformTraits.init(mdata.origV, mdata.D, mdata.F, isExactDC, mdata.constIndices,rigidityFactor);
    mdata.deformSolver.init(&mdata.deformLinearSolver, &mdata.deformTraits);
    
  }
  

  IGL_INLINE void quat_moebius_deform(struct QuatMoebiusData& mdata,
                                         const Eigen::MatrixXd& qh,
                                         const int numIterations,
                                         Eigen::MatrixXd& q){
    
    //Coords2Quat(qh, mdata.quatConstPoses);
    
    //feeding initial solution as the previous one
    mdata.deformTraits.constPoses=mdata.constPoses;
    mdata.deformTraits.currPositions=mdata.deformV;
    mdata.deformTraits.currX=mdata.deformX;

    mdata.deformSolver.solve(true);
    mdata.deformV=mdata.deformTraits.finalPositions;
    mdata.deformX=mdata.deformTraits.finalX;
    
    //Complex2Coords(mdata.deformVc, mdata.deformV);
    q=mdata.deformV;
  }
}


#endif
