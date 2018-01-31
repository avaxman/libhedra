// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_COMPLEX_MOEBIUS_DEFORM_H
#define HEDRA_COMPLEX_MOEBIUS_DEFORM_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <hedra/Moebius2DEdgeDeviationTraits.h>
#include <hedra/EigenSolverWrapper.h>
#include <hedra/LMSolver.h>
//#include <hedra/complex_cross_ratio.h>

using namespace Eigen;
using namespace std;

typedef std::complex<double> Complex;

void Complex2Coords(const VectorXcd& CV, MatrixXd& V)
{
  V.resize(CV.rows(), 3);
  for (int i=0;i<CV.rows();i++)
    V.row(i)<<CV(i).real(), CV(i).imag(),0.0;
}


void Coords2Complex(const MatrixXd& V, VectorXcd& CV)
{
  CV.resize(V.rows());
  for (int i=0;i<V.rows();i++)
    CV(i)=std::complex<double>(V(i,0), V(i,1));
}


namespace hedra{
  
  
  struct ComplexMoebiusData{
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
    VectorXcd origVc;
    VectorXcd deformVc;
    
    //Edge deviation method variables
    VectorXcd deformY;  //candidate reciprocals of the Mobius transformation.
    VectorXcd deformE;  //errors of vertex reciprocals vs. actual transformed vertices
    
    //for inversion control and other usages
    SparseMatrix<Complex> d0;
    
    //of positional handles
    VectorXi constIndices;
    VectorXcd complexConstPoses;
    
    //optimization operators
    hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > > deformLinearSolver;
    hedra::optimization::Moebius2DEdgeDeviationTraits deformTraits;
    hedra::optimization::LMSolver<hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > >,hedra::optimization::Moebius2DEdgeDeviationTraits> deformSolver;
    
  };
  
  IGL_INLINE void complex_moebius_precompute(const Eigen::VectorXi& h,
                                             const bool isExactDC,
                                             const bool isExactIAP,
                                             const double rigidityFactor,
                                             struct ComplexMoebiusData& mdata){
    
    mdata.constIndices=h;
    mdata.isExactDC=isExactDC;
    mdata.isExactIAP=isExactIAP;
    mdata.rigidityFactor=rigidityFactor;
    mdata.deformTraits.init(mdata.origVc, mdata.D, mdata.F, mdata.extEV,isExactDC, isExactIAP, mdata.constIndices,rigidityFactor);
    mdata.deformSolver.init(&mdata.deformLinearSolver, &mdata.deformTraits);
    
  }
  


  IGL_INLINE void complex_moebius_setup(const Eigen::MatrixXd& V,
                                        const Eigen::VectorXi& D,
                                        const Eigen::MatrixXi& F,
                                        const Eigen::MatrixXi& TF,
                                        const Eigen::MatrixXi& EV,
                                        const Eigen::MatrixXi& EF,
                                        const Eigen::MatrixXi& EFi,
                                        const Eigen::MatrixXi& FE,
                                        const Eigen::MatrixXd& FEs,
                                        const Eigen::VectorXi& innerEdges,
                                        struct ComplexMoebiusData& mdata)
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
    Coords2Complex(mdata.origV,mdata.origVc);
    mdata.deformVc=mdata.origVc;
    
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
    mdata.deformY=VectorXcd::Ones(mdata.deformVc.rows());
    mdata.deformE=VectorXcd::Ones(mdata.extEV.rows());
    
  }
  

  
  IGL_INLINE void complex_moebius_deform(struct ComplexMoebiusData& mdata,
                                         const Eigen::MatrixXd& qh,
                                         const int numIterations,
                                         Eigen::MatrixXd& q){
    
    Coords2Complex(qh, mdata.complexConstPoses);
    
    //feeding initial solution as the previous one
      mdata.deformTraits.complexConstPoses=mdata.complexConstPoses;
    mdata.deformTraits.currPositions=mdata.deformVc;
    mdata.deformTraits.currY=mdata.deformY;
    mdata.deformTraits.currE=mdata.deformE;
    
    mdata.deformSolver.solve(true);
    mdata.deformVc=mdata.deformTraits.finalPositions;
    mdata.deformY=mdata.deformTraits.finalY;
    mdata.deformE=mdata.deformTraits.finalE;
    
    Complex2Coords(mdata.deformVc, mdata.deformV);
    q=mdata.deformV;
  }
}


#endif
