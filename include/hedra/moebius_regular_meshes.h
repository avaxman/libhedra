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
#include <hedra/regularity.h>
#include <hedra/dcel.h>
#include <hedra/planarity.h>
#include <hedra/willmore_energy.h>

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
    
    Eigen::VectorXi HV, HF, HE, VH, nextH, prevH, twinH;
    Eigen::MatrixXi EH, FH;
    
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
    
 
    //assuming the angle is [0, pi] always.
    void factorize_quaternion(const Eigen::RowVector4d& q, double& length, double& angle, Eigen::RowVector3d& vec)
    {
      length=q.norm();
      if (length==0.0){
        std::cout<<"Zero quaternion! "<<std::endl;
        vec<<0.0,0.0,0.0;
        angle=M_PI/2;
        return;
      }
      double cosangle=q(0)/length;
      cosangle=(cosangle>1.0 ? 1.0 : cosangle);
      cosangle=(cosangle<-1.0 ? -1.0 : cosangle);
      angle=acos(cosangle);
      if (abs(sin(angle))>10e-6)
        vec=q.tail(3)/(length*sin(angle));
      else
        vec<<0.0,0.0,0.0;
      if ((isnan(vec(0)))||(isnan(vec(1)))||(isnan(vec(2)))){
        std::cout<<"nan in vec!!: "<<q<<std::endl;
        exit(0);
      }
    }
    
    /*double cot(double x){
      if (abs(sin(x))>10e-5)
        return cos(x)/sin(x);
      else return 0.0;
    }*/
    
    /*void compute_mean_curvature(const Eigen::MatrixXi& VValences,
                                const Eigen::MatrixXi& QuadVertexIndices,
                                const Eigen::MatrixXd& Vq,
                                Eigen::VectorXd& H)
    {
      using namespace Eigen;
      H.resize(Vq.rows());
      H.setZero();
      VectorXd VertCotWeights(Vq.rows());
      VertCotWeights.setZero();
      for (int i=0;i<QuadVertexIndices.rows();i++){
        RowVector4d qi=Vq.row(QuadVertexIndices(i,0));
        RowVector4d qj=Vq.row(QuadVertexIndices(i,1));
        RowVector4d qk=Vq.row(QuadVertexIndices(i,2));
        RowVector4d ql=Vq.row(QuadVertexIndices(i,3));
        
        
        RowVector4d Nj=QMult(qj-qi, QInv(qk-qj));
        RowVector4d Nl=QMult(ql-qk, QInv(qi-ql));
        
        
        RowVector4d qCommutator=QMult(Nj,Nl)-QMult(Nl,Nj);
        double lambda=qCommutator.norm()/(qi-qk).norm();
        lambda*=(qCommutator.dot(qi-qk)>0 ? 1 : -1);
        
        double Cosj=Nj(0)/Nj.norm();
        double Cosl=Nl(0)/Nl.norm();
        
        lambda*=(Nj.tail(3).dot(QMult(Nj,Nl).tail(3))> 0.0 ? -1.0 : 1.0);
        Cosj=(Cosj>1.0 ? 1.0 : Cosj);
        Cosj=(Cosj<-1.0 ? -1.0 : Cosj);
        Cosl=(Cosl>1.0 ? 1.0 : Cosl);
        Cosl=(Cosl<-1.0 ? -1.0 : Cosl);
        double CotWeight=cot(M_PI-acos(Cosj))+cot(M_PI-acos(Cosl));
        VertCotWeights(QuadVertexIndices(i,0))+=CotWeight;
        if(isnan(VertCotWeights(QuadVertexIndices(i,0))))
          exit(0);
        
        H(QuadVertexIndices(i,0))+=CotWeight*lambda/QMult(Nj,Nl).tail(3).norm();
        
      }
      
      H=H.cwiseQuotient(VertCotWeights);
      for (int i=0;i<H.size();i++)
        if (isnan(H(i))) H(i)=0;  //TODO: fix isolated vertices!
    }*/
    
    
    
    
    //Compute the energy of quaternionic differences of ratios around a loop indicated by OneRings (unless it's a boundary).
    //the difference is factored with the proper lengths and angles per ratio, so we are measuring the difference from a prescribed 1-ring.
    void compute_ratio_diff_energy(const Eigen::VectorXi& VValences,
                                   const Eigen::MatrixXd& Ratios,
                                   const Eigen::MatrixXi& OneRings,
                                   const Eigen::VectorXi& BoundaryMask,
                                   const Eigen::VectorXd& Lengths,
                                   const Eigen::VectorXd& Angles,
                                   Eigen::VectorXd& W)
    {
      using namespace Eigen;
      
      W.setZero();
      for (int i=0;i<OneRings.rows();i++){
        int NumFlaps=VValences(i)-2*BoundaryMask(i);
        //double avgAngle=((double)numFlaps-2.0)*M_PI/(double)numFlaps;
        for (int j=0;j<NumFlaps;j++){
          double length1, length2, angle1, angle2;
          RowVector3d vec1, vec2;
          //cout<<"Ratios.row(OneRings(i,j)): "<<Ratios.row(OneRings(i,j))<<endl;
          factorize_quaternion(Ratios.row(OneRings(i,j)),length1, angle1, vec1);
          //cout<<"length1, angle1, vec1: "<<length1<<","<<angle1<<","<<vec1<<endl;
          //cout<<"Ratios.row(OneRings(i,j+1)): "<<Ratios.row(OneRings(i,j+1))<<endl;
          factorize_quaternion(Ratios.row(OneRings(i,(j+1)%NumFlaps)), length2, angle2, vec2);
          //cout<<"length2, angle2, vec2: "<<length2<<","<<angle2<<","<<vec2<<endl;
          
          //cout<<"lengths(oneRings(i,j)): "<<Lengths(OneRings(i,j))<<endl;
          //cout<<"Angles(oneRings(i,j)): "<<Angles(OneRings(i,j))<<endl;
          //cout<<"lengths(oneRings(i,j+1)): "<<Lengths(OneRings(i,j+1))<<endl;
          //cout<<"Angles(oneRings(i,j+1)): "<<Angles(OneRings(i,j+1))<<endl;
          
          RowVector4d CompRatio1; CompRatio1<<length1/Lengths(OneRings(i,j))*cos(angle1-Angles(OneRings(i,j))), length1/Lengths(OneRings(i,j))*sin(angle1-Angles(OneRings(i,j)))*vec1;
          RowVector4d CompRatio2; CompRatio2<<length2/Lengths(OneRings(i,(j+1)%NumFlaps))*cos(angle2-Angles(OneRings(i,(j+1)%NumFlaps))), length2/Lengths(OneRings(i,(j+1)%NumFlaps))*sin(angle2-Angles(OneRings(i,(j+1)%NumFlaps)))*vec2;
          
          //cout<<"CompRatio1: "<<CompRatio1<<endl;
          //cout<<"CompRatio2: "<<CompRatio2<<endl;
          
          
          //if (isDirection)
          //  W(i)+=avgAngle-angle//(vec1.cross(vec2)).norm()/(NumFlaps-1);
          //else
          W(i)+=(CompRatio1-CompRatio2).squaredNorm();
          
          //cout<<"Curr W("<<i<<"): "<<W(i)<<endl;
          
          /*RowVector4d Inv1; Inv1<<cos(M_PI/2-Angles(OneRings(i,j))), sin(M_PI/2-Angles(OneRings(i,j)))*CornerCR.row(OneRings(i,j)).tail(3).normalized();
           RowVector4d Inv2; Inv2<<cos(M_PI/2-Angles(OneRings(i,j+1))), sin(M_PI/2-Angles(OneRings(i,j+1)))*CornerCR.row(OneRings(i,j+1)).tail(3).normalized();
           Inv1/=Lengths(OneRings(i,j));
           Inv2/=Lengths(OneRings(i,j+1));
           if (!isDirection)
           W(i)+=(QMult1(CornerCR.row(OneRings(i,j+1)),Inv2)-QMult1(CornerCR.row(OneRings(i,j)),Inv1)).norm()/(NumFlaps-1);
           else
           W(i)=(CornerCR.row(OneRings(i,j)).tail(3).normalized().array().square()-CornerCR.row(OneRings(i,j+1)).tail(3).normalized().array().square()).matrix().norm()/(NumFlaps-1);*/
          
        }
        W(i)=sqrt(W(i));
      }
      
      //cout<<"Total W: "<<W<<endl;
    }
    
    //averages the ratio vectors to find the common one. Factors out the given lengths and angles
    
    void estimate_common_ratio_vectors(const Eigen::VectorXi& VValences,
                                       const Eigen::MatrixXd& Ratios,
                                       const Eigen::MatrixXi& OneRings,
                                       const Eigen::VectorXi& BoundaryMask,
                                       Eigen::MatrixXd& CommonRatios)
    {
      using namespace Eigen;
      CommonRatios.resize(OneRings.rows(),3); CommonRatios.setZero();
      for (int i=0;i<OneRings.rows();i++){
        int NumFlaps=VValences(i)-2*BoundaryMask(i);
        for (int j=0;j<NumFlaps;j++){
          double length, angle;
          RowVector3d vec;
          factorize_quaternion(Ratios.row(OneRings(i,j)),length, angle, vec);
          CommonRatios.row(i)+=vec;
          if ((isnan(CommonRatios(i,0)))||(isnan(CommonRatios(i,1)))||(isnan(CommonRatios(i,2)))){
            std::cout<<"isnan CommonRatios: "<<vec<<CommonRatios.row(i)<<std::endl;
            exit(0);
          }
        }
        if (CommonRatios.row(i).lpNorm<Infinity>()<10e-6)
          CommonRatios(i,0)=1.0;  //just random
      }
      CommonRatios.rowwise().normalize();
    }
    
    
    //getting the radius of the unique circumcircle for the polygon with these edge lengths
    double get_radius_from_lengths(Eigen::VectorXd& lengths,
                                   double AngleSum)
    {
      double Precision=10e-5;
      double MinRadius=lengths.maxCoeff()/2;
      double MaxRadius=MinRadius*100;
     
      double MidAngleSum=1000.0;
      while (abs(MidAngleSum/(AngleSum/2.0)-1.0)>Precision){
        double MidRadius=(MaxRadius+MinRadius)/2.0;
        MidAngleSum=asin(lengths.array()/(2.0*MidRadius)).matrix().sum();
        
        if (MidAngleSum<(AngleSum/2.0))
          MaxRadius=MidRadius;
        else
          MinRadius=MidRadius;
        
      }
      return (MinRadius+MaxRadius)/2.0;
    }
    
    void estimate_combinatorial_intrinsics(const Eigen::MatrixXd& Vq,
                                           const Eigen::VectorXi& D,
                                           const Eigen::VectorXi& VValences,
                                           const Eigen::MatrixXi& QuadVertexIndices,
                                           const Eigen::MatrixXi& OneRings,
                                           const Eigen::VectorXi& BoundaryMask,
                                           Eigen::VectorXd& Lengths,
                                           Eigen::VectorXd& Angles,
                                           bool Smooth)
    {
      using namespace Eigen;
      using namespace std;
      Lengths.resize(QuadVertexIndices.rows()); Lengths.setZero();
      Angles=Lengths;
      for (int i=0;i<OneRings.rows();i++){
        int NumFlaps=VValences(i)-2*BoundaryMask(i);
        if (NumFlaps==0)  //an "ear" of the mesh, sitting on a single face. no given intrinsics
          continue;
        
        set<int> SetFaces;
        VectorXd NewSectorAngles;
        double NaturalSum=0;
        for (int j=0;j<NumFlaps;j++){
          SetFaces.insert(QuadVertexIndices(OneRings(i,j),4));
          SetFaces.insert(QuadVertexIndices(OneRings(i,j),5));
        }
        
        VectorXi Faces(SetFaces.size());
        int SetCounter=0;
        for (set<int>::iterator si=SetFaces.begin();si!=SetFaces.end();si++)
          Faces(SetCounter++)=*si;
        
        //cout<<"Faces Degrees: ";
        VectorXd RingLengths(Faces.size());
        for (int k=0;k<Faces.size();k++){
          //cout<<D(Faces(k))<<","<<endl;
          double RegAngle=(D(Faces(k))-2.0)*M_PI/(double)D(Faces(k));
          RingLengths(k)=2*sin(RegAngle/2.0);
        }
        
        
        double Radius=get_radius_from_lengths(RingLengths, (2-BoundaryMask(i))*M_PI);
        
        
        for (int j=0;j<NumFlaps;j++){
          
          double RegAngle1=(D(QuadVertexIndices(OneRings(i,j),4))-2)*M_PI/(double)D(QuadVertexIndices(OneRings(i,j),4));
          double RegAngle2=(D(QuadVertexIndices(OneRings(i,j),5))-2)*M_PI/(double)D(QuadVertexIndices(OneRings(i,j),5));
         
          double lf=2*sin(RegAngle1/2.0);
          double lg=2*sin(RegAngle2/2.0);
          
          double SectorAnglef=2*asin(lf/(2*Radius));
          double SectorAngleg=2*asin(lg/(2*Radius));
          
          Lengths(OneRings(i,j))= lg/lf; //Lengths(OneRings(i,j)+1)=
         
          Angles(OneRings(i,j))=M_PI-(SectorAnglef+SectorAngleg)/2;
         
          RowVector4d qi=Vq.row(QuadVertexIndices(OneRings(i,j),0));
          RowVector4d qj=Vq.row(QuadVertexIndices(OneRings(i,j),1));
          RowVector4d qk=Vq.row(QuadVertexIndices(OneRings(i,j),2));
          RowVector4d ql=Vq.row(QuadVertexIndices(OneRings(i,j),3));
          
          RowVector4d LocalCR=QMult(QMult(qj-qi, QInv(qk-qj)),QMult(ql-qk, QInv(qi-ql)));
          
        }
        
        //testing
        double SumAngles=0;
        double ProdLengths=1.0;
        for (int j=0;j<NumFlaps;j++){
          ProdLengths*=Lengths(OneRings(i,j));
          SumAngles+=Angles(OneRings(i,j));
        }
      }
    }

    
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
    MRData.oneRings.resize(VOrig.rows(), MRData.vertexValences.maxCoeff());
    MRData.oneRings.setConstant(-1);
    VectorXi ringIndices(VOrig.rows()); ringIndices.setZero();
    for (int i=0;i<MRData.quadVertexIndices.rows();i++){
      MRData.oneRings(MRData.quadVertexIndices(i,0),ringIndices(MRData.quadVertexIndices(i,0)))=i;
      ringIndices(MRData.quadVertexIndices(i,0))++;
    }
  
    MRData.oneRingVertices.resize(VOrig.rows(), MRData.vertexValences.maxCoeff()+1);
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
    
    MRData.boundaryMask.resize(VOrig.rows()); MRData.boundaryMask.setZero();
    for (int i=0;i<boundaryList.size();i++)
      for (int j=0;j<boundaryList[i].size();j++)
        MRData.boundaryMask(boundaryList[i][j])=1;
    
    MRData.vertexValences+=MRData.boundaryMask;  //VValences is about vertices
    
    MRData.constMask=VectorXi::Zero(VOrig.rows());
    
    MRData.boundaryVertices.resize(MRData.boundaryMask.sum());
    int currBoundVertex=0;
    for (int i=0;i<VOrig.rows();i++)
      if (MRData.boundaryMask(i)==1)
        MRData.boundaryVertices[currBoundVertex++]=i;
    
  
    /***************Estimating original CR and FN values*********************/
    MRData.VCR.resize(VOrig.rows(),3);
    MRData.FN.resize(F.rows(),3);
    hedra::quat_cross_ratio(MRData.VOrig,MRData.quadVertexIndices, MRData.origECR);
    //ComputeCR(OrigVq, QuadVertexIndices, OrigECR);
    hedra::quat_normals(MRData.QOrig, MRData.faceTriads, MRData.origCFN);
    //ComputeFN(OrigVq, FaceTriads, OrigCFN);
    
    
    MRData.estimate_combinatorial_intrinsics(MRData.QOrig, D, MRData.vertexValences, MRData.quadVertexIndices, MRData.oneRings, MRData.boundaryMask, MRData.patternCRLengths, MRData.patternCRAngles, true);
    
    MRData.estimate_common_ratio_vectors(MRData.vertexValences, MRData.origECR, MRData.oneRings, MRData.boundaryMask, MRData.VCR);
    MRData.estimate_common_ratio_vectors(D, MRData.origCFN, MRData.cornerF, VectorXi::Zero(D.size()), MRData.FN);
    
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
      double angle=igl::PI*((double)D(i)-2.0)/(double)D(i);
      for (int j=0;j<D(i);j++)
        MRData.patternFNAngles(MRData.cornerF(i,j))=igl::PI-angle;
      
    }
    
    currFaceQuad=0;
    MRData.patternFaceCRLengths.resize(MRData.quadFaceIndices.rows()); MRData.patternFaceCRLengths.setOnes();
    MRData.patternFaceCRLengths.resize(MRData.quadFaceIndices.rows());
    MRData.patternFaceCRAngles=MRData.patternFaceCRLengths;
    for (int i=0;i<F.rows();i++){
      double angle=igl::PI*((double)D(i)-2.0)/(double)D(i);
      double oppositeLength=1+2*sin(angle-igl::PI/2);
      for (int j=0;j<D(i)-3;j++){
        MRData.patternFaceCRLengths(currFaceQuad)=1.0/oppositeLength;
        MRData.patternFaceCRAngles(currFaceQuad++)=igl::PI;
      }
      
    }
    
    MRData.deformECR=MRData.origECR;
    MRData.deformCFN=MRData.origCFN;
    
    //prescribed lengths are the originals initially (until externally modified)
    MRData.prescribedLengths.resize(EV.rows());
    for (int i=0;i<EV.rows();i++)
      MRData.prescribedLengths(i)=(VOrig.row(EV(i,0))-VOrig.row(EV(i,1))).norm();
    
    
    hedra::dcel(MRData.D, MRData.F,MRData.EV,MRData.EF,MRData.EFi,MRData.innerEdges,MRData.VH,MRData.EH,MRData.FH,MRData.HV,MRData.HE,MRData.HF,MRData.nextH, MRData.prevH,MRData.twinH);
    
    /****************Computing initial energies**************************/
    
    MRData.origMR.resize(VOrig.rows());
    MRData.origW=MRData.origMR;
    MRData.origER.resize(F.rows());
    MRData.compute_ratio_diff_energy(MRData.vertexValences, MRData.origECR, MRData.oneRings, MRData.boundaryMask, MRData.patternCRLengths, MRData.patternCRAngles, MRData.origMR);
    hedra::willmore_energy(MRData.VOrig, MRData.VH, MRData.HV, MRData.HE, MRData.HF, MRData.twinH, MRData.nextH, MRData.prevH, MRData.origW);
    //MRData.compute_ratio_diff_energy(MRData.vertexValences, MRData.origECR, MRData.oneRings, MRData.boundaryMask, MRData.patternCRLengths, MRData.patternCRAngles, MRData.origW);
    
    //MRData.compute_ratio_diff_energy(D, MRData.origCFN, MRData.cornerF, VectorXi::Zero(D.size()), MRData.patternFNLengths, MRData.patternFNAngles, MRData.origER);
    
    hedra::regularity(VOrig,MRData.D,MRData.F,MRData.origER);
    MRData.deformMR=MRData.origMR;
    MRData.deformER=MRData.origER;
    MRData.deformW=MRData.origW;
    
    MRData.CSolver.CRLengths=MRData.patternCRLengths;
    MRData.CSolver.CRAngles=MRData.patternCRAngles;
    MRData.CSolver.FNLengths=MRData.patternFNLengths;
    MRData.CSolver.FNAngles=MRData.patternFNAngles;
    MRData.CSolver.faceCRLengths=MRData.patternFaceCRLengths;
    MRData.CSolver.faceCRAngles=MRData.patternFaceCRAngles;
    
    //ComputeMeanCurvature(VValences, QuadVertexIndices, OrigVq, H);
    
    MRData.CSolver.init(MRData.QOrig, D, F, EV, MRData.quadVertexIndices, MRData.quadFaceIndices, MRData.faceTriads);
    
    MRData.constIndices = constIndices;
    MRData.CSolver.set_constant_handles(constIndices);
    return true;
  }
  
  
  IGL_INLINE bool compute_moebius_regular(MoebiusRegularData& MRData,
                                          const double MRCoeff,
                                          const double ERCoeff,
                                          const Eigen::MatrixXd& constPoses,
                                          const bool outputProgress,
                                          Eigen::MatrixXd& VRegular)
  {
    
    //composing initial solution
    
    for (int i=0;i<MRData.QOrig.rows();i++)
      for (int j=0;j<3;j++)
        MRData.CSolver.currSolution[3*i+j]=MRData.VDeform(i,j);
    
    for (int i=0;i<MRData.QOrig.rows();i++)
      for (int j=0;j<3;j++)
        MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*i+j]=MRData.VCR(i,j);
    
    for (int i=0;i<MRData.F.rows();i++)
      for (int j=0;j<3;j++)
        MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*MRData.QOrig.rows()+3*i+j]=MRData.FN(i,j);
    
    
    for (int i=0;i<MRData.constIndices.size();i++)
      for (int j=0;j<3;j++)
        MRData.CSolver.currSolution[3*MRData.constIndices(i)+j]=constPoses(i,j);
    
    
    MRData.CSolver.solve(MRCoeff, ERCoeff, outputProgress);
    
    for (int i=0;i<MRData.QOrig.rows();i++)
      MRData.VDeform.row(i)<<MRData.CSolver.currSolution[3*i],MRData.CSolver.currSolution[3*i+1],MRData.CSolver.currSolution[3*i+2];
    
    VRegular = MRData.VDeform;
    
    for (int i=0;i<MRData.QOrig.rows();i++)
      MRData.VCR.row(i)<<MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*i],MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*i+1],MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*i+2];
    
    for (int i=0;i<MRData.F.rows();i++)
      MRData.FN.row(i)<<MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*MRData.QOrig.rows()+3*i],MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*MRData.QOrig.rows()+3*i+1],MRData.CSolver.currSolution[3*MRData.QOrig.rows()+3*MRData.QOrig.rows()+3*i+2];
    
    Coords2Quat(constPoses, MRData.quatConstPoses);
    Coords2Quat(MRData.VDeform, MRData.QDeform);
    
    hedra::quat_cross_ratio(MRData.VDeform,MRData.quadVertexIndices, MRData.deformECR);
    hedra::quat_normals(MRData.QDeform, MRData.faceTriads, MRData.deformCFN);
    
    MRData.compute_ratio_diff_energy(MRData.vertexValences, MRData.deformECR, MRData.oneRings, MRData.boundaryMask, MRData.patternCRLengths, MRData.patternCRAngles, MRData.deformMR);
    //MRData.compute_ratio_diff_energy(MRData.vertexValences, MRData.deformECR, MRData.oneRings, MRData.boundaryMask, MRData.patternCRLengths, MRData.patternCRAngles, MRData.deformW, true);
    hedra::willmore_energy(MRData.VDeform, MRData.VH, MRData.HV, MRData.HE, MRData.HF, MRData.twinH, MRData.nextH, MRData.prevH, MRData.deformW);
    
    //hedra::moebius_regularity(VRegular, MRData.F, MRData)
    hedra::regularity(VRegular,MRData.D,MRData.F,MRData.deformER);
   //MRData.compute_ratio_diff_energy(MRData.D, MRData.deformCFN, MRData.cornerF, Eigen::VectorXi::Zero(MRData.D.size()), MRData.patternFNLengths, MRData.patternFNAngles, MRData.deformER);
    
    return true;
    
  }
  
}
#endif
