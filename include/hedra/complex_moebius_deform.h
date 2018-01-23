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
#include <hedra/complex_cross_ratio.h>

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
    MatrixXi T;     //triangulated mesh (only for comparisons and visualization)
    VectorXi TF;
    MatrixXi edgeT;  //the edge mesh for MC/IAP error etc. visualization

    MatrixXi quadVertexIndices;  //a list of all four consecutive vertices, corresponding to edge-mesh, that on which the cross ratio is defined

    Eigen::MatrixXd origV, deformV;
    Eigen::MatrixXd origEdgeV, deformEdgeV
    Eigen::MatrixXi F, T, edgeT, edgeTE;
    Eigen::VectorXi D, TF;
    Eigen::MatrixXi EV, EF, FE, EFi;
    Eigen::MatrixXd FEs;
    Eigen::VectorXi innerEdges;
    
    bool isExactMC;
    bool isExactIAP
    double rigidityFactor;
    
    //Raw 3D positions
    MatrixXd origV;    //original vertex positions
    MatrixXd deformV;  //current deformed mesh
    
    //Complex Positions
    VectorXcd OrigVc;
    VectorXcd DeformVc;
    
    //Edge deviation method variables
    VectorXcd deformY;  //candidate reciprocals of the Mobius transformation.
    VectorXcd deformE;  //errors of vertex reciprocals vs. actual transformed vertices
    
    //for inversion control and other usages
    SparseMatrix<Complex> d0;
    
    //of positional handles
    VectorXi constIndices;
    VectorXcd complexConstPoses;
    
    //optimization operators
    hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > > DeformLinearSolver;
    hedra::optimization::Moebius2DEdgeDeviationTraits DeformTraits;
    hedra::optimization::LMSolver<hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > >,hedra::optimization::Moebius2DEdgeDeviationTraits> DeformSolver;
    
  }
  
  IGL_INLINE void complex_moebius_precompute(const double rigidityFactor,
                                             const bool isExactMC,
                                        const bool isExactIAP,
                                             const Eigen::VectorXi& h,
                                             struct ComplexMoebiusData& mdata){
    
    mdata.constIndices=h;
    
    mdata.deformTraits.init(mdata..origVc, mdata.D, mdata.F, mdata.extEV, isExactMC, isExactIAP, h);
    mdata.deformTraits.rigidityFactor=rigidityFactor;
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
                                        const Eigen::MatrixXd& FEs
                                        const Eigen::vectorXi& InnerEdges,
                                        struct ComplexMoebiusData& mdata)
  {
    mdata.F=F;
    mdata.D=D;
    origV=V;
    deformV=V;
    mdata.TF=TF;
    mdata.EV=EV;
    mdata.EF=EF;
    mdata.EFi=EFi;
    mdata.FE=FE;
    mdata.FEs=FEs;
    mdata.InnerEdges=InnerEdges;

    //complex representations
    Coords2Complex(origV,mdata.origVc);
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
    
    
    mdata.quadVertexIndices.resize(innerEdges.size(),4);
    for (int i=0;i<innerEdges.rows();i++){
      int f=EF(innerEdges(i),0);
      int g=EF(innerEdges(i),1);
      int vi=EV(innerEdges(i),0);
      int vk=EV(innerEdges(i),1);
      int vj=F(g,(EFi(innerEdges(i),1)+2)%D(g));
      int vl=F(f,(EFi(innerEdges(i),0)+2)%D(f));
      int vf=mdata.origV.rows()+EF(innerEdges(i),0);
      int vg=mdata.origV.rows()+EF(innerEdges(i),1);
      
      mdata.quadVertexIndices.row(i)<<vi,vj,vk,vl;
    }
    
  }
  
  
  IGL_INLINE void complex_moebius_deform(struct ComplexMoebiusData& mdata,
                                         const Eigen::MatrixXd& qh,
                                         const int numIterations,
                                         Eigen::MatrixXd& q){
    
    Coords2Complex(qh, mdata.complexConstPoses);
    //mdata.deformTraits.smoothFactor=100.0;
    //mdata.deformTraits.posFactor=10.0;
    
    //feeding initial solution as the previous one
    
    DeformSolver.solve(true);
    DeformVc=DeformTraits.finalPositions;
    mdata.DeformY=DeformTraits.finalY;
    mdata.DeformE=DeformTraits.finalE;
    
    Complex2Coords(DeformVc, DeformV);
    ComputeCR(DeformVc, D, F,QuadVertexIndices, DeformECR, DeformFCR);
    
    EdgeDeformV<<DeformV, GetCenters(DeformV, D, F);
    
    
  }
}


#endif
