// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_MOEBIUS_TRAITS_H
#define HEDRA_MOEBIUS_TRAITS_H
#include <igl/igl_inline.h>
#include "quaternionic_derivatives.h"
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>
#include <set>


namespace hedra { namespace optimization {
    

    //This traits class implements the "ConstraintTraits" concept, for the energy and constraints of the quaternionic system of deformation in [Vaxman et. al 2015] (see Tutorial).
    class MoebiusTraits{
    public:
        
        //class requirements
        Eigen::VectorXi JERows, JECols;  //rows and column indices for the jacobian matrix
        Eigen::VectorXd JEVals;         //values for the jacobian matrix.
        Eigen::VectorXd EVec;          //energy vector
        int xSize;                  //size of the solution
        //for the constraints
        Eigen::VectorXi JCRows, JCCols;  //rows and column indices for the jacobian matrix
        Eigen::VectorXd JCVals;         //values for the jacobian matrix.
        Eigen::VectorXd CVec;          //energy vector

        Eigen::MatrixXi D, F;
        Eigen::VectorXi ConstIndices;
        Eigen::MatrixXi CornerPairs;
        Eigen::VectorXi CornerOffset;
        Eigen::SparseMatrix<double> d0;
        Eigen::MatrixXd OrigVq;
        Eigen::MatrixXd ConstPoses;
        
        Eigen::MatrixXi EdgeCornerPairs;  //for rigidity energy
        Eigen::MatrixXi EdgeCornerVertices;  //the vertices themselves, for compatibility
        
        Eigen::VectorXd prevSolution;   //from any previous running of the algorithm (in "raw" optimization format, with the corner variables etc.)
        
        Eigen::MatrixXd fullSolution;   //contains the locations of the vertices after the current run.
        
        int NumCorners;
        int NumEdgePairs;
        
        bool isExactMC;

        double SmoothFactor;
        double CloseFactor;
        double PosFactor;
        double RigidRatio;
        
        RowVector4d UnitQuat;
        
        //intermediate variables
        VectorXd CurrX;
        VectorXd CurrLocations;
        VectorXd AMAPVec;
        VectorXd RigidVec;
        VectorXd CompVec;
        VectorXd CloseVec;
        VectorXd PosVec;
        VectorXd MCVec;
 
        int AMAPTriOffset, AMAPRowOffset;
        int RigidTriOffset, RigidRowOffset;
        int CloseTriOffset, CloseRowOffset;
        int CompTriOffset, CompRowOffset;
        int PosTriOffset, PosRowOffset;
        int MCTriOffset, MCRowOffset;

        void init(const Eigen::MatrixXd& _VOrig,
                  const Eigen::MatrixXi& _D ,
                  const Eigen::MatrixXi& _F ,
                  const Eigen::VectorXi& _ConstIndices,
                  const bool _isExactMC,
                  const bool resetPrevSolution=true)
        {
            
            using namespace Eigen;
            
            //imaginary quaternionic representation
            OrigVq.resize(VOrig.rows(),4);
            OrigVq.setZero();
            OrigVq.block(0,1,VOrig.rows(),3)=VOrig;
            F=_F; D=_D;
            ConstIndices=_ConstIndices;
            
            UnitQuat<<1.0,0.0,0.0,0.0;

            //MobConstMat=inMobConstMat;
            isExactMC=_isExactMC;

            //computing relevant topology
            NumCorners=D.sum();
            CornerOffset.resize(F.rows());
            CornerOffset(0)=0;
            
            for (int i=1;i<D.rows();i++)
                CornerOffset(i)=CornerOffset(i-1)+D(i-1);
            
            MatrixXi AdjCorners(OrigV.rows(),12);
            VectorXi Valences(OrigV.rows()); Valences.setZero();
            for (int i=0;i<D.rows();i++){
                for (int j=0;j<D(i);j++){
                    AdjCorners(F(i,j),Valences(F(i,j)))=CornerOffset(i)+j;
                    Valences(F(i,j))++;
                }
            }

            vector<pair<int,int> > CornerPairList;
            for (int i=0;i<OrigV.rows();i++)
                for (int j=0;j<Valences(i)-1;j++)
                    CornerPairList.push_back(pair<int,int>(AdjCorners(i,j),AdjCorners(i,(j+1)%Valences(i))));
            
            CornerPairs.resize(CornerPairList.size(),2);
            for (int i=0;i<CornerPairList.size();i++)
                CornerPairs.row(i)<<CornerPairList[i].first, CornerPairList[i].second;
            
            xSize=4*NumCorners+3*OrigVq.rows();
            CurrX.resize(4*NumCorners);
            CurrLocations.resize(3*OrigVq.rows());
            
            if (resetPrevSolution){
                prevSolution.conservativeResize(xSize);
                for (int i=0;i<DeformX.rows();i++)
                    prevSolution.segment(4*i,4)=UnitQuat;    //corner variables are trivial
                
                for (int i=0;i<DeformVq.rows();i++)
                    prevSolution.segment(4*NumCorners+3*i,3)=_VOrig.row(i);
            }
            
            VectorXi FaceNums=D;
            FaceNums=(FaceNums.cwiseAbs2()-FaceNums);
            NumEdgePairs=FaceNums.sum()/2;
            
            EdgeCornerPairs.resize(NumEdgePairs,2);
            EdgeCornerVertices.resize(NumEdgePairs,2);
            int CurrPair=0;
            for (int i=0;i<D.rows();i++)
                for (int j=0;j<D(i);j++)
                    for (int k=j+1;k<D(i);k++){
                        EdgeCornerPairs.row(CurrPair)<<CornerOffset(i)+j, CornerOffset(i)+k;
                        EdgeCornerVertices.row(CurrPair++)<<F(i,j), F(i,k);
                    }
            
            
            AMAPVec.resize(4*CornerPairs.rows());
            RigidVec.resize(4*EdgeCornerPairs.rows());
            CompVec.resize(4*EdgeCornerPairs.rows());
            CloseVec.resize(SolutionSize);
            PosVec.resize(3*ConstIndices.size());
            if (isExactMC){
                MCVec.resize(CornerPairs.rows());
            else
                MCVec.resize(0);
            
            EVec.resize(AMAPVec.size()+RigidVec.size()+CloseVec.size());
            CVec.resize(CompVec.size()+PosVec.size()+MCVec.size());
            
                
            //TODO: I think I can get rid of this because of Augmented Lagrangian
            CloseFactor=10e-6;
     
            //Constructing Gradient Pattern
            
            JERows.resize(2*4*CornerPairs.rows()+2*4*EdgeCornerPairs.rows()+xSize);
            JECols.resize(JERows.size());
            JEVals.resize(JERows.size());
            
            if (!isExactMC){
                JCRows.resize(38*EdgeCornerPairs.rows()+3*ConstIndices.size());
                JCCols.resize(GradRows.size());
                JCVals.resize(GradRows.size());
            } else{
                JCRows.resize(38*EdgeCornerPairs.rows()+3*ConstIndices.size()+2*4*CornerPairs.rows());
                JCCols.resize(GradRows.size());
                JCVals.resize(GradRows.size());
            }
            
            /*******************************AMAP Energy********************************************/
            AMAPTriOffset=0;
            AMAPRowOffset=0;
            for (int i=0;i<CornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    JERows(AMAPTriOffset+2*(4*i+j))=AMAPRowOffset+4*i+j;
                    JECols(AMAPTriOffset+2*(4*i+j))=4*CornerPairs(i,0)+j;
                    JERows(AMAPTriOffset+2*(4*i+j)+1)=AMAPRowOffset+4*i+j;
                    JECols(AMAPTriOffset+2*(4*i+j)+1)=4*CornerPairs(i,1)+j;
                }
            }
            
            /*******************************Rigidity Energy*****************************************/
            RigidTriOffset=AMAPTriOffset+2*4*CornerPairs.rows();
            RigidRowOffset=AMAPRowOffset+4*CornerPairs.rows();
            for (int i=0;i<EdgeCornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    JERows(RigidTriOffset+2*(4*i+j))=RigidRowOffset+4*i+j;
                    JECols(RigidTriOffset+2*(4*i+j))=4*EdgeCornerPairs(i,0)+j;
                    JERows(RigidTriOffset+2*(4*i+j)+1)=RigidRowOffset+4*i+j;
                    JECols(RigidTriOffset+2*(4*i+j)+1)=4*EdgeCornerPairs(i,1)+j;
                }
            }

            
            /****************************Closeness Energy*******************/
            CloseTriOffset=RigidTriOffset+2*4*EdgeCornerPairs.rows();
            CloseRowOffset=RigidRowOffset+4*EdgeCornerPairs.rows();
            for (int i=0;i<SolutionSize;i++){
                JERows(CloseTriOffset+i)=CloseRowOffset+i;
                JECols(CloseTriOffset+i)=i;
                JEVals(CloseTriOffset+i)=CloseFactor;
            }
            
            /****************************Compatibility Constraints*****************/
            CompTriOffset=0;
            CompRowOffset=0;
            int CompTriCounter=CompTriOffset;
            Vector4i XiTriPoses; XiTriPoses<<0,8,18,28;
            Vector4i XjTriPoses; XjTriPoses<<4,12,22,32;
            for (int i=0;i<EdgeCornerPairs.rows();i++){
                int ColCorneri=4*EdgeCornerPairs(i,0);
                int ColCornerj=4*EdgeCornerPairs(i,1);
                int CurrRowOffset=4*i;
                //derivative of Xi
                quatDerivativeIndices(JCRows, JCCols, CompTriCounter, XiTriPoses, CompRowOffset+CurrRowOffset, ColCorneri);
                
                //Derivative of Xj
                quatDerivativeIndices(JCRows, JCCols, CompTriCounter, XjTriPoses,CompRowOffset+CurrRowOffset, ColCornerj);
                
                for (int k=0;k<3;k++){
                    //wj derivative
                    GradRows(CompTriCounter+16+10*k)=CompRowOffset+CurrRowOffset+1+k;
                    GradCols(CompTriCounter+16+10*k)=4*NumCorners+3*EdgeCornerVertices(i,1)+k;
                    GradValues(CompTriCounter+16+10*k)=-1.0;
                    
                    //wi derivative
                    GradRows(CompTriCounter+17+10*k)=CompRowOffset+CurrRowOffset+1+k;
                    GradCols(CompTriCounter+17+10*k)=4*NumCorners+3*EdgeCornerVertices(i,0)+k;
                    GradValues(CompTriCounter+17+10*k)=1.0;
                }
                
                CompTriCounter+=38;
            }
            
            /****************************Positional Constraints*******************/
            PosTriOffset=CompTriOffset+38*EdgeCornerPairs.rows();
            PosRowOffset=CompRowOffset+4*EdgeCornerPairs.rows();
            
            for (int i=0;i<ConstIndices.size();i++){
                for (int k=0;k<3;k++){
                    JCRows(PosTriOffset+3*i+k)=PosRowOffset+3*i+k;
                    JCCols(PosTriOffset+3*i+k)=4*NumCorners+3*ConstIndices(i)+k;
                }
            }
            
            
            /****************************Metric-Conformal Constraints*************/
            if (isExactMC){
                MCTriOffset=PosTriOffset+3*ConstIndices.size();
                MCRowOffset=PosRowOffset+3*ConstIndices.size();
                for (int i=0;i<CornerPairs.rows();i++){
                    for (int j=0;j<4;j++){
                        JCRows(MCTriOffset+2*(4*i+j))=MCRowOffset+i;
                        JCCols(MCTriOffset+2*(4*i+j))=4*CornerPairs(i,0)+j;
                        JCRows(MCTriOffset+2*(4*i+j)+1)=MCRowOffset+i;
                        JCCols(MCTriOffset+2*(4*i+j)+1)=4*CornerPairs(i,1)+j;
                    }
                }
            }
        }
        
        void update_energy(const Eigen::VectorXd& x){
            using namespace Eigen;
            CurrX<<x.head(4*NumCorners);
            CurrLocations<<x.tail(3*OrigVq.rows());
            
            for (int i=0;i<CornerPairs.rows();i++)
                AMAPVec.segment(4*i,4)=(CurrX.segment(4*CornerPairs(i,1),4)-CurrX.segment(4*CornerPairs(i,0),4));
            
            AMAPVec.array()*=SmoothFactor;
            
            for (int i=0;i<EdgeCornerPairs.rows();i++)
                RigidVec.segment(4*i,4)=(CurrX.segment(4*EdgeCornerPairs(i,1),4)-CurrX.segment(4*EdgeCornerPairs(i,0),4));

            
            RigidVec.array()*=SmoothFactor*RigidRatio;
            
            CloseVec<<CloseFactor*(CurrSolution-InitSolution);
        
            EVec<<AMAPVec, RigidVec, CloseVec

        }
        
        void update_constraints(const Eigen::VectorXd& x){
            using namespace Eigen;
            CurrX<<x.head(4*NumCorners);
            CurrLocations<<x.tail(3*OrigVq.rows());
            
            for (int i=0;i<EdgeCornerPairs.rows();i++){
                RowVector4d Xi=CurrX.segment(4*EdgeCornerPairs(i,0),4).transpose();
                RowVector4d Xj=CurrX.segment(4*EdgeCornerPairs(i,1),4).transpose();
                RowVector4d qi=OrigVq.row(EdgeCornerVertices(i,0));
                RowVector4d qj=OrigVq.row(EdgeCornerVertices(i,1));
                RowVector3d wi=CurrLocations.segment(3*EdgeCornerVertices(i,0),3).transpose();
                RowVector3d wj=CurrLocations.segment(3*EdgeCornerVertices(i,1),3).transpose();
                RowVector4d CurrEdgeVector=QMult1(QMult1(QConj1(Xi),qj-qi),Xj);
                CurrEdgeVector.tail(3)-=(wj-wi);
                CompVec.segment(4*i,4)<<CurrEdgeVector.transpose();
            }
            
            
            for (int i=0;i<ConstIndices.size();i++)
                PosVec.segment(3*i,3)<<PosFactor*(CurrLocations.segment(3*ConstIndices(i),3)-ConstPoses.row(i).transpose());
            
            if (!isExactMC){
                CVec<<CompVec, PosVec;
            } else {
                for (int i=0;i<CornerPairs.rows();i++)
                    MCVec(i)=CurrX.segment(4*CornerPairs(i,1),4).squaredNorm()-CurrX.segment(4*CornerPairs(i,0),4).squaredNorm();
                CVec<<CompVec, PosVec, MCVec;
            }
        }
        
        void update_jacobian(const Eigen::VectorXd& x){
            using namespace Eigen;
            CurrX<<CurrSolution.head(4*NumCorners);
            CurrLocations<<CurrSolution.tail(3*OrigVq.rows());
            
            /*******************************AMAP Energy********************************************/
            for (int i=0;i<CornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    JEVals(AMAPTriOffset+2*(4*i+j))=-SmoothFactor;
                    JEVals(AMAPTriOffset+2*(4*i+j)+1)=SmoothFactor;
                }
            }
            
            /*******************************Rigidity Energy*****************************************/
            for (int i=0;i<EdgeCornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    JEVals(RigidTriOffset+2*(4*i+j))=-SmoothFactor*RigidRatio;
                    JEVals(RigidTriOffset+2*(4*i+j)+1)=SmoothFactor*RigidRatio;
                }
            }
            
            //closeness energy is constant
            
            /****************************Compatibility Constraints*****************/
            int CompTriCounter=CompTriOffset;
            Vector4i XiTriPoses; XiTriPoses<<0,8,18,28;
            Vector4i XjTriPoses; XjTriPoses<<4,12,22,32;
            for (int i=0;i<EdgeCornerPairs.rows();i++){
                RowVector4d Xi=CurrX.segment(4*EdgeCornerPairs(i,0),4).transpose();
                RowVector4d Xj=CurrX.segment(4*EdgeCornerPairs(i,1),4).transpose();
                RowVector4d RightPart=QMult1(OrigVq.row(EdgeCornerVertices(i,1))-OrigVq.row(EdgeCornerVertices(i,0)),Xj);
                RowVector4d LeftPart=QMult1(QConj1(Xi),OrigVq.row(EdgeCornerVertices(i,1))-OrigVq.row(EdgeCornerVertices(i,0)));
                
      
                //derivative of Xi
                quatDerivativeValues(JCVals, CompTriCounter, XiTriPoses, UnitQuat, RightPart, true, false);
                
                //Derivative of Xj
                quatDerivativeValues(JCVals, CompTriCounter, XjTriPoses, LeftPart, UnitQuat, false, false);
                
                //the other compatibility values are constant
                CompTriCounter+=38;
                
                
            }
        
            
            
            /****************************Positional Constraints*******************/
            for (int i=0;i<ConstIndices.size();i++)
                for (int k=0;k<3;k++)
                    JCVals(PosTriOffset+3*i+k)=PosFactor;
            
            /****************************Metric-Conformal Constraints*************/
            if (isExactMC){
                for (int i=0;i<CornerPairs.rows();i++){
                    RowVector4d Xi=CurrX.segment(4*CornerPairs(i,0),4).transpose();
                    RowVector4d Xj=CurrX.segment(4*CornerPairs(i,1),4).transpose();
                    for (int k=0;k<4;k++){
                        JCVals(MCTriOffset+2*(4*i+k))=-2*Xi(k);
                        JCVals(MCTriOffset+2*(4*i+k)+1)=2*Xj(k);
                    }
                }
            }

        }
        
        //provide the initial solution to the solver
        void initial_solution(Eigen::VectorXd& x0){
            x0=prevSolution;
        }
        
        void pre_iteration(const Eigen::VectorXd& prevx){}
        bool post_iteration(const Eigen::VectorXd& x){return false;  /*never stop after an iteration*/}
        
        bool post_optimization(const Eigen::VectorXd& x){
            
            fullSolution.conservativeResize(OrigVq.rows(),3);
            for (int i=0;i<DeformVq.rows();i++)
                fullSolution.row(i)=x.segment(4*NumCorners+3*i,3);
            
            prevSolution=fullSolution;
            
            return true;  //this traits doesn't have any more stop requirements
        }
        
        MoebiusTraits(){}
        ~MoebiusTraits(){}

        
        
        
    };
}



#endif
