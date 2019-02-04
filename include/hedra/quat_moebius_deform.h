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
    
    //vertex variables
    MatrixXd deformY;
    
    //of positional handles
    VectorXi constIndices;
    MatrixXd constPoses;
    
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
    qmdata.deformV=V;
    qmdata.deformY.resize(V.rows(),4);
    qmdata.deformY.setZero();
    qmdata.deformY.col(0).setConstant(1.0);
    
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
    
  }
  
  IGL_INLINE void quat_moebius_precompute(const Eigen::VectorXi& h,
                                          struct QuatMoebiusData& qmdata){
    
    qmdata.solver.init(qmdata.origV, qmdata.D, qmdata.F,qmdata.extEV);
    qmdata.constIndices=h;
    qmdata.solver.set_constant_handles(h);
    
    
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
    
    for (int i=0;i<qmdata.deformY.rows();i++)
      for (int j=0;j<4;j++)
        qmdata.solver.currSolution[3*qmdata.deformV.rows()+4*i+j]=qmdata.deformY(i,j);
    
    for (int i=0;i<qmdata.constIndices.size();i++)
      for (int j=0;j<3;j++)
        qmdata.solver.currSolution[3*qmdata.constIndices(i)+j]=qh(i,j);
    
    
    //For now just a big DCFactor
    qmdata.solver.solve(AMAPFactor, rigidityFactor, (isExactDC ? 1000.0 : 0.0), outputProgress);
    
    for (int i=0;i<qmdata.origV.rows();i++)
      qmdata.deformV.row(i)<<qmdata.solver.currSolution[3*i],qmdata.solver.currSolution[3*i+1],qmdata.solver.currSolution[3*i+2];
    
    deformV = qmdata.deformV;
    
    for (int i=0;i<qmdata.deformY.rows();i++){
      qmdata.deformY.row(i)<<qmdata.solver.currSolution[3*qmdata.origV.rows()+4*i],
      qmdata.solver.currSolution[3*qmdata.origV.rows()+4*i+1],
      qmdata.solver.currSolution[3*qmdata.origV.rows()+4*i+2],
      qmdata.solver.currSolution[3*qmdata.origV.rows()+4*i+3];
    }
  
  }
}


#endif
