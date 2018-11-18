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
#include <hedra/quaternionic_operations.h>
#include <hedra/quat_cross_ratio.h>
#include <hedra/quat_normals.h>

namespace hedra
{
  
  struct MoebiusRegularData{
    
    double fullCRFactor;
    double fullFNFactor;
    double unitFactor;
    double HFactor;
    double vecCRFactor;
    double lengthCRFactor;
    double lengthFNFactor;
    double phaseCRFactor;
    double phaseFNFactor;
    
    Eigen::MatrixXi quadVertexIndices;  //a list of all four consecutive vertices, corresponding to edge-mesh, that on which the cross ratio is defined
    Eigen::MatrixXi quadFaceIndices;  //in case of non-triangular faces that are optimized to be cocircular
    Eigen::MatrixXi faceTriads;  //a list of face-consecutive triads on which corner normals are defined
    
    Eigen::VectorXi D;  //face degrees
    Eigen::MatrixXi F;
    Eigen::MatrixXi EV,EF,FE, EFi;
    Eigen::MatrixXi extEV;  //extended with diagonals
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
    Eigen::MatrixXd VDeform, QDeform;    //original vertex positions
    
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
                                        const Eigen::MatrixXi& T,
                                        const Eigen::MatrixXi& EV,
                                        const Eigen::MatrixXi& FE,
                                        const Eigen::MatrixXi& EF,
                                        const Eigen::MatrixXi& EFi,
                                        const Eigen::MatrixXd& FEs,
                                        const Eigen::VectorXi& innerEdges,
                                        const Eigen::VectorXi& constIndices,
                                        MoebiusRegularData& MRData){
    
    using namespace Eigen;
    using namespace std;
    
    MRData.F=F;
    MRData.D=D;
    
    MRData.EV =EV;
    MRData.FE=FE;
    MRData.EF=EF;
    MRData.EFi=EFi;
    MRData.FEs=FEs;
    MRData.innerEdges=innerEdges;
    MRData.constIndices=constIndices;
    
    MRData.VOrig=VOrig;
    MRData.VDeform=VOrig;
    
    MRData.convErrors.resize(1);
    MRData.convErrors(0)=0.0;
    
    //quaternionic representations
    Coords2Quat(MRData.VOrig, MRData.QOrig);
    MRData.QDeform=MRData.QOrig;
    
    //creating full edge list
    vector<pair<int,int>>  Diagonals;
    for (int i=0;i<F.rows();i++)
      for (int j=0;j<D(i);j++)
        for (int k=j+1;k<D(i);k++)
          Diagonals.push_back(pair<int,int>(F(i,j),F(i,k)));
    
    MRData.extEV.resize(EV.rows()+Diagonals.size(),2);
    MRData.extEV.block(0,0,EV.rows(),2)=EV;
    for (int i=0;i<Diagonals.size();i++)
      MRData.extEV.row(EV.rows()+i)<<Diagonals[i].first, Diagonals[i].second;
    
    MRData.quadVertexIndices.resize(2*innerEdges.size(),6);
    MRData.quadFaceIndices.resize(D.sum()-3*D.size(),5);  //for pure triangular meshes - zero
    MRData.faceTriads.resize(D.sum(),4);
    MRData.cornerF.resize(F.rows(), D.maxCoeff());
    MRData.cornerF.setConstant(-1);
    int currTriad=0;
    for (int i=0;i<D.rows();i++){
      for (int j=0;j<D(i);j++){
        MRData.faceTriads.row(currTriad)<<F(i,j), F(i,(j+1)%D(i)), F(i,(j+2)%D(i)),i;
        MRData.cornerF(i,j)=currTriad++;
      }
    }
    
    int currFaceQuad=0;
    for (int i=0;i<D.rows();i++)
      for (int j=0;j<D(i)-3;j++)
        MRData.quadFaceIndices.row(currFaceQuad++)<<F(i,j), F(i,j+1),F(i,j+2), F(i,j+3),i;
    

    for (int i=0;i<innerEdges.rows();i++){
      int f=EF(innerEdges(i),0);
      int g=EF(innerEdges(i),1);
      
      //from the side i->k
      int vis=EV(innerEdges(i),0);
      int vks=EV(innerEdges(i),1);
      int vjs=F(g,(EFi(innerEdges(i),1)+2)%D(g));
      int vls=F(f,(EFi(innerEdges(i),0)+D(f)-1)%D(f));
      
      //from the side k->i
      int vit=EV(innerEdges(i),1);
      int vkt=EV(innerEdges(i),0);
      int vjt=F(f,(EFi(innerEdges(i),0)+2)%D(f));
      int vlt=F(g,(EFi(innerEdges(i),1)+D(g)-1)%D(g));
      
      
      int vf=VOrig.rows()+EF(innerEdges(i),0);
      int vg=VOrig.rows()+EF(innerEdges(i),1);
      
      MRData.quadVertexIndices.row(2*i)  <<vis, vjs, vks, vls, g, f;
      MRData.quadVertexIndices.row(2*i+1)<<vit, vjt, vkt, vlt, f, g;
    }
    
    
    MRData.vertexValences.resize(VOrig.rows());
    MRData.vertexValences.setZero();
    for (int i=0;i<D.rows();i++)
      for (int j=0;j<D(i);j++)
        MRData.vertexValences(F(i,j))++;
    
    //computing one-rings - currently not with any order (do I need it?)
    MRData.oneRings.resize(VOrig.rows(), MRDATA.vertexValences.maxCoeff());
    MRData.oneRings.setConstant(-1);
    VectorXi ringIndices(VOrig.rows()); ringIndices.setZero();
    for (int i=0;i<MRData.quadVertexIndices.rows();i++){
      MRData.oneRings(MRData.quadVertexIndices(i,0),ringIndices(MRDATA.quadVertexIndices(i,0)))=i;
      ringIndices(MRData.quadVertexIndices(i,0))++;
    }
  
    MRDATA.oneRingVertices.resize(VOrig.rows(), MRDATA.vertexValences.maxCoeff()+1);
    ringIndices.setZero();
    for (int i=0;i<EV.rows();i++){
      //cout<<"RingIndices(E2V(i,0)): "<<RingIndices(E2V(i,0))<<endl;
      //cout<<"E2V(i,0): "<<E2V(i,0)<<endl;
      MRData.oneRingVertices(EV(i,0),ringIndices(EV(i,0)))=EV(i,1);
      MRData.oneRingVertices(EV(i,1),ringIndices(EV(i,1)))=EV(i,0);
      ringIndices(EV(i,0))++;
      ringIndices(EV(i,1))++;
    }
    
    vector<vector<int> > boundaryList;
    igl::boundary_loop(T, boundaryList);
    
    MRDATA.boundaryMask.resize(VOrig.rows()); MRDATA.boundaryMask.setZero();
    for (int i=0;i<boundaryList.size();i++)
      for (int j=0;j<boundaryList[i].size();j++)
        MRDATA.boundaryMask(boundaryList[i][j])=1;
    
    MRData.vertexValences+=boundaryMask;  //VValences is about vertices
    
    MRData.constMask=VectorXi::Zero(VOrig.rows());
    
    MRData.boundaryVertices.resize(MRData.boundaryMask.sum());
    int currBoundVertex=0;
    for (int i=0;i<OrigV.rows();i++)
      if (boundaryMask[i]==1)
        MRData.boundaryVertices[currBoundVertex++]=i;
    
  
    /***************Estimating original CR and FN values*********************/
    MRDATA.VCR.resize(VOrig.rows(),3);
    MRDATA.FN.resize(F.rows(),3);
    hedra::quat_cross_ratio(MRData.QOrig,QuadVertexIndices, MRData.origECR);
    //ComputeCR(OrigVq, QuadVertexIndices, OrigECR);
    hedra::quat_normals(MRData.QOrig, MRData.FaceTriads, MRData.OrigCFN);
    //ComputeFN(OrigVq, FaceTriads, OrigCFN);
    
    
    MRData.estimate_combinatorial_intrinsics(OrigVq, D, VValences, QuadVertexIndices, OneRings, BoundaryMask, PatternCRLengths, PatternCRAngles, true);
    
    MRData.estimate_common_ratio_vectors(VValences, OrigECR, OneRings, BoundaryMask, VCR);
    MRData.estimate_common_ratio_vectors(D, OrigCFN, CornerF, VectorXi::Zero(D.size()), FN);
    
    //completing "ear" vertices VCR
    for (int i=0;i<D.rows();i++){
      for (int j=0;j<D(i);j++){
        if (MRData.vertexValences(F(i,j))!=1)
          continue;
        
        RowVector3d vi=VOrig.row(F(i,(j+D(i)-1)%D(i)));
        RowVector3d vj=VOrig.row(F(i,j));
        RowVector3d vk=VOrig.row(F(i,(j+1)%D(i)));
        
        MRData.VCR.row(F(i,j))=((vk-vj).cross(vi-vj)).normalized();
        //cout<<"Valences 2 VCR: "<<VCR.row(F(i,j))<<endl;
      }
      
    }
    
    //estimating intrinsics in every face is trivial
    MRData.patternFNLengths.resize(MRData.faceTriads.rows()); MRData.patternFNLengths.setOnes();
    MRData.patternFNAngles.resize(MRData.faceTriads.rows());
    for (int i=0;i<MRData.cornerF.rows();i++){
      double angle=igl::M_PI*((double)D(i)-2.0)/(double)D(i);
      for (int j=0;j<D(i);j++)
        MRData.patternFNAngles(MRData.cornerF(i,j))=igl::M_PI-angle;
      
    }
    
    currFaceQuad=0;
    MRData.patternFaceCRLengths.resize(MRData.quadFaceIndices.rows()); MRData.patternFaceCRLengths.setOnes();
    MRData.patternFaceCRLengths.resize(MRData.quadFaceIndices.rows());
    MRData.patternFaceCRAngles=MRData.patternFaceCRLengths;
    for (int i=0;i<F.rows();i++){
      double angle=igl::M_PI*((double)D(i)-2.0)/(double)D(i);
      double oppositeLength=1+2*sin(angle-igl::M_PI/2);
      for (int j=0;j<D(i)-3;j++){
        MRData.patternFaceCRLengths(currFaceQuad)=1.0/pppositeLength;
        MRData.patternFaceCRAngles(currFaceQuad++)=igl::M_PI;
      }
      
    }
    
    MRData.deformECR=MRData.origECR;
    MRData.deformCFN=MRData.origCFN;
    
    //prescribed lengths are the originals initially (until externally modified)
    MRData.prescribedLengths.resize(EV.rows());
    for (int i=0;i<EV.rows();i++)
      MRData.prescribedLengths(i)=(VOrig.row(EV(i,0))-VOrig.row(EV(i,1))).norm();
    
    /****************Computing initial energies**************************/
    
    MRData.origMR.resize(VOrig.rows());
    MRData.origW=MRData.origMR;
    MRData.origER.resize(F.rows());
    MRData.compute_ratio_diff_energy(VValences, OrigECR, OneRings, BoundaryMask, PatternCRLengths, PatternCRAngles, OrigMR, false);
    
    MRData.compute_ratio_diff_energy(VValences, OrigECR, OneRings, BoundaryMask, PatternCRLengths, PatternCRAngles, OrigW, true);
    
    MRData.compute_ratio_diff_energy(D, OrigCFN, CornerF, VectorXi::Zero(D.size()), PatternFNLengths, PatternFNAngles, OrigER, false);
    MRData.deformMR=OrigMR;
    MRData.deformER=OrigER;
    MRData.deformW=OrigW;
    
    //ComputeMeanCurvature(VValences, QuadVertexIndices, OrigVq, H);
    
    CSolver.Init(MRData.QOrig, D, F, EV, MRData.quadVertexIndices, MRData.quadFaceIndices, MRData.faceTriads);
    
    MRData.constIndices = constIndices;
    CSolver.set_constant_handles(constIndices);
    return true;
  }
  
  
  IGL_INLINE bool compute_moebius_regular(const MoebiusRegularData& MRData,
                                          const double MRCoeff,
                                          const double ERCoeff,
                                          const Eigen::MatrixXd& constPoses,
                                          Eigen::MatrixXd& VRegular)
  {
    
    /*MatrixXd ConstPoses(InConstPoses.size(),3);
    for (int i=0;i<InConstPoses.size();i++)
      ConstPoses.row(i)=InConstPoses[i];*/
    

    //composing initial solution
    
    for (int i=0;i<MRData.QOrig.rows();i++)
      for (int j=0;j<3;j++)
        MRData.CSolver.currSolution[3*i+j]=MRData.VDeform(i,j);
    
    for (int i=0;i<MRData.QOrig.rows();i++)
      for (int j=0;j<3;j++)
        MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*i+j]=MRData.VCR(i,j);
    
    for (int i=0;i<F.rows();i++)
      for (int j=0;j<3;j++)
        MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*MRData.QOrig.rows()+3*i+j]=MRData.FN(i,j);
    
    
    for (int i=0;i<MRData.constIndices.size();i++)
      for (int j=0;j<3;j++)
        MRData.CSolver.currSolution[3*ConstPosIndices(i)+j]=ConstPoses(i,j);
    
    
    MRData.CSolver.CRLengths<<MRData.patternCRLengths;
    MRData.CSolver.CRAngles<<MRData.patternCRAngles;
    MRData.CSolver.FNLengths<<MRData.patternFNLengths;
    MRData.CSolver.FNAngles<<MRData.patternFNAngles;
    MRData.CSolver.faceCRLengths<<MRData.patternFaceCRLengths;
    MRData.CSolver.faceCRAngles<<MRData.patternFaceCRAngles;
  
    MRData.CSolver.Solve(MRFactor, ERFactor);
    
    for (int i=0;i<MRData.QOrig.rows();i++)
      MRData.VDeform.row(i)<<MRData.CSolver.currSolution[3*i],MRData.CSolver.currSolution[3*i+1],MRData.CSolver.currSolution[3*i+2];
    
    for (int i=0;i<MRData.QOrig.rows();i++)
      MRData.VCR.row(i)<<MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*i],MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*i+1],MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*i+2];
    
    for (int i=0;i<F.rows();i++)
      MRData.FN.row(i)<<MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*MRData.QOrig.rows()+3*i],MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*MRData.QOrig.rows()+3*i+1],MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*MRData.QOrig.rows()+3*i+2];
    
    Coords2Quat(ConstPoses, QuatConstPoses);
    Coords2Quat(DeformV, DeformVq);
    //ComputeCR(DeformVq, QuadVertexIndices, DeformECR);
    //ComputeFN(DeformVq, FaceTriads, DeformCFN);
    
    hedra::quat_cross_ratio(MRData.QDeform,QuadVertexIndices, MRData.deformECR);
    hedra::quat_normals(MRData.QDeform, MRData.FaceTriads, MRData.deformCFN);
    
    MRData.compute_ratio_diff_energy(VValences, DeformECR, OneRings, BoundaryMask, PatternCRLengths, PatternCRAngles, DeformMR, false);
    MRData.compute_ratio_diff_energy(VValences, DeformECR, OneRings, BoundaryMask, PatternCRLengths, PatternCRAngles, DeformW, true);
    MRData.compute_ratio_diff_energy(D, DeformCFN, CornerF, VectorXi::Zero(D.size()), PatternFNLengths, PatternFNAngles, DeformER, false);
    
    return true;
    
  }
  
}
#endif
