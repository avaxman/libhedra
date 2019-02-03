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
//#include <hedra/Moebius3DCornerVarsTraits.h>
#include <hedra/EigenSolverWrapper.h>
//#include <hedra/LMSolver.h>
#include <hedra/quaternionic_operations.h>
#include <hedra/CeresQuatDeformSolver.h>

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
    //VectorXi cornerOffset;  //where does every face begin in the corner list
    
    //MatrixXi faceCornerPairs;   //faces (f,g, i,k) and vertices around every edge for compatibility and smoothness, like in the paper
    //MatrixXi cornerPairs;  //neighboring corner pairs across edges for AMAP/MC/IAP comparisons
    
    //for inversion control and other usages
    //SparseMatrix<Complex> d0;
    
    //of positional handles
    VectorXi constIndices;
    MatrixXd constPoses;
    
    //optimization operators
    //hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > > deformLinearSolver;
    //hedra::optimization::Moebius3DCornerVarsTraits deformTraits;
    //hedra::optimization::LMSolver<hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > >,hedra::optimization::Moebius3DCornerVarsTraits> deformSolver;
    CeresQMDSolver solver;
    
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
                                        struct QuatMoebiusData& qmdata)
  {
    qmdata.F=F;
    qmdata.D=D;
    qmdata.origV=V;
    qmdata.deformV=V;
    qmdata.TF=TF;
    qmdata.EV=EV;
    qmdata.EF=EF;
    qmdata.EFi=EFi;
    qmdata.FE=FE;
    qmdata.FEs=FEs;
    qmdata.innerEdges=innerEdges;
    
    qmdata.origV=V;
    
    //creating full edge list
    vector<pair<int,int>>  diagonals;
    for (int i=0;i<F.rows();i++)
      for (int j=0;j<D(i);j++)
        for (int k=j+2;k<D(i);k++)
          diagonals.push_back(pair<int,int>(F(i,j),F(i,k)));
    
    qmdata.extEV.resize(EV.rows()+diagonals.size(),2);
    qmdata.extEV.block(0,0,EV.rows(),2)=EV;
    for (int i=0;i<diagonals.size();i++)
      qmdata.extEV.row(EV.rows()+i)<<diagonals[i].first, diagonals[i].second;
    
    //reseting deformation result
    /*qmdata.deformX.resize(D.sum(),4);
    qmdata.deformX.setZero();
    qmdata.deformX.col(0).setConstant(1.0);
    qmdata.cornerOffset.resize(F.rows());
    qmdata.cornerOffset(0)=0;
    
    for (int i=1;i<D.rows();i++)
      qmdata.cornerOffset(i)=mdata.cornerOffset(i-1)+D(i-1);
    
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
    
    qmdata.cornerPairs.resize(cornerPairList.size(),2);
    for (int i=0;i<cornerPairList.size();i++)
      qmdata.cornerPairs.row(i)<<cornerPairList[i].first, cornerPairList[i].second;*/
    
  }
  
  IGL_INLINE void quat_moebius_precompute(const Eigen::VectorXi& h,
                                          struct QuatMoebiusData& qmdata){
    
    /*qmdata.constIndices=h;
    qmdata.isExactDC=isExactDC;
    qmdata.isExactIAP=isExactIAP;
    qmdata.rigidityFactor=rigidityFactor;*/
    qmdata.solver.init(qmdata.origV, qmdata.D, qmdata.F,qmdata.extEV);
    qmdata.solver.set_constant_handles(h);
    //mdata.deformSolver.init(&mdata.deformLinearSolver, &mdata.deformTraits);
    
  }
  

  IGL_INLINE void quat_moebius_deform(struct QuatMoebiusData& qmdata,
                                      const double AMAPFactor,
                                      const double rigidityFactor,
                                      const bool isExactDC,
                                      const Eigen::MatrixXd& qh,
                                      const bool outputProgress,
                                      Eigen::MatrixXd& deformV){
    
    //Coords2Quat(qh, mdata.quatConstPoses);
    
    for (int i=0;i<qmdata.origV.rows();i++)
      for (int j=0;j<3;j++)
        qmdata.solver.currSolution[3*i+j]=qmdata.deformV(i,j);
    
    for (int i=0;i<qmdata.deformX.rows();i++)
      for (int j=0;j<4;j++)
        qmdata.solver.currSolution[3*qmdata.deformVq.rows()+4*i+j]=qmdata.deformX(i,j);
    
    for (int i=0;i<qmdata.constIndices.size();i++)
      for (int j=0;j<3;j++)
        qmdata.solver.currSolution[3*qmdata.constIndices(i)+j]=qh(i,j);
    
    
    //For now just a big DCFactor
    qmdata.solver.solve(AMAPFactor, rigidityFactor, (isExactDC ? 1000.0 : 0.0), outputProgress);
    
    for (int i=0;i<qmdata.origV.rows();i++)
      qmdata.deformV.row(i)<<qmdata.solver.currSolution[3*i],qmdata.solver.currSolution[3*i+1],qmdata.solver.currSolution[3*i+2];
    
    for (int i=0;i<qmdata.deformX.rows();i++){
      qmdata.deformX.row(i)<<qmdata.solver.currSolution[3*qmdata.origV.rows()+4*i],
      qmdata.solver.currSolution[3*qmdata.origV.rows()+4*i+1],
      qmdata.solver.currSolution[3*qmdata.origV.rows()+4*i+2],
      qmdata.solver.currSolution[3*qmdata.origV.rows()+4*i+3];
    }
  
  }
}


#endif
