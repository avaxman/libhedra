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
        
        //concept requirements
        Eigen::VectorXi JERows, JECols;  //rows and column indices for the jacobian matrix
        Eigen::VectorXd JEVals;         //values for the jacobian matrix.
        Eigen::VectorXd EVec;          //energy vector
        int xSize;                  //size of the solution
        //for the constraints
        Eigen::VectorXi JCRows, JCCols;  //rows and column indices for the jacobian matrix
        Eigen::VectorXd JCVals;         //values for the jacobian matrix.
        Eigen::VectorXd CVec;          //energy vector

        
        Eigen::MatrixXd VOrigq;
        Eigen::MatrixXi D, F;
        Eigen::VectorXi constIndices;
        Eigen::MatrixXi cornerPairs;
        Eigen::VectorXi cornerOffset;

        Eigen::MatrixXd constPoses;
        
        Eigen::MatrixXi edgeCornerPairs;  //for rigidity energy
        Eigen::MatrixXi edgeCornerVertices;  //the vertices themselves, for compatibility
        
        Eigen::VectorXd prevSolution;   //from any previous running of the algorithm (in "raw" optimization format, with the corner variables etc.)
        
        Eigen::MatrixXd fullSolution;   //contains the locations of the vertices after the current run.
        
        int numCorners;
        int numEdgePairs;
        
        bool isExactMC;

        double smoothFactor;
        //double closeFactor;
        double posFactor;
        double rigidRatio;
        
        RowVector4d unitQuat;
        
        //intermediate variables
        VectorXd currX;
        VectorXd currLocations;
        VectorXd AMAPVec;
        VectorXd rigidVec;
        VectorXd compVec;
        //VectorXd closeVec;
        VectorXd posVec;
        VectorXd MCVec;
 
        int AMAPTriOffset, AMAPRowOffset;
        int rigidTriOffset, rigidRowOffset;
        //int closeTriOffset, CloseRowOffset;
        int compTriOffset, compRowOffset;
        int posTriOffset, posRowOffset;
        int MCTriOffset, MCRowOffset;

        void init(const Eigen::MatrixXd& _VOrig,
                  const Eigen::MatrixXi& _D ,
                  const Eigen::MatrixXi& _F ,
                  const Eigen::VectorXi& _constIndices,
                  const bool _isExactMC,
                  const bool resetPrevSolution=true)
        {
            
            using namespace Eigen;
            
            //imaginary quaternionic representation
            VOrigq.resize(VOrig.rows(),4);
            VOrigq.setZero();
            VOrigq.block(0,1,VOrig.rows(),3)=VOrig;
            F=_F; D=_D;
            constIndices=_constIndices;
            
            unitQuat<<1.0,0.0,0.0,0.0;

            //MobConstMat=inMobConstMat;
            isExactMC=_isExactMC;

            //computing relevant topology
            numCorners=D.sum();
            cornerOffset.resize(F.rows());
            cornerOffset(0)=0;
            
            for (int i=1;i<D.rows();i++)
                cornerOffset(i)=cornerOffset(i-1)+D(i-1);
            
            MatrixXi adjCorners(VOrig.rows(),12);
            VectorXi valences(VOrig.rows()); Valences.setZero();
            for (int i=0;i<D.rows();i++){
                for (int j=0;j<D(i);j++){
                    adjCorners(F(i,j),valences(F(i,j)))=cornerOffset(i)+j;
                    valences(F(i,j))++;
                }
            }

            vector<pair<int,int> > cornerPairList;
            for (int i=0;i<VOrig.rows();i++)
                for (int j=0;j<valences(i)-1;j++)
                    cornerPairList.push_back(pair<int,int>(adjCorners(i,j),adjCorners(i,(j+1)%valences(i))));
            
            cornerPairs.resize(cornerPairList.size(),2);
            for (int i=0;i<cornerPairList.size();i++)
                cornerPairs.row(i)<<cornerPairList[i].first, cornerPairList[i].second;
            
            xSize=4*numCorners+3*VOrigq.rows();
            currX.resize(4*numCorners);
            currLocations.reize(3*VOrigq.rows());
            
            if (resetPrevSolution){
                prevSolution.conservativeResize(xSize);
                for (int i=0;i<Vorigq.rows();i++)
                    prevSolution.segment(4*i,4)=unitQuat;    //corner variables are trivial
                
                for (int i=0;i<Vorigq.rows();i++)
                    prevSolution.segment(4*NumCorners+3*i,3)=_VOrig.row(i);
            }
            
            VectorXi faceNums=D;
            faceNums=(faceNums.cwiseAbs2()-faceNums);
            numEdgePairs=faceNums.sum()/2;
            
            edgeCornerPairs.resize(numEdgePairs,2);
            edgeCornerVertices.resize(numEdgePairs,2);
            int currPair=0;
            for (int i=0;i<D.rows();i++)
                for (int j=0;j<D(i);j++)
                    for (int k=j+1;k<D(i);k++){
                        edgeCornerPairs.row(currPair)<<cornerOffset(i)+j, cornerOffset(i)+k;
                        edgeCornerVertices.row(currPair++)<<F(i,j), F(i,k);
                    }
            
            
            AMAPVec.resize(4*cornerPairs.rows());
            rigidVec.resize(4*edgeCornerPairs.rows());
            compVec.resize(4*edgeCornerPairs.rows());
            //CloseVec.resize(xSize);
            posVec.resize(3*constIndices.size());
            if (isExactMC)
                MCVec.resize(cornerPairs.rows());
            else
                MCVec.resize(0);
            
            EVec.resize(AMAPVec.size()+rigidVec.size());//+CloseVec.size());
            CVec.resize(compVec.size()+posVec.size()+MCVec.size());
            
                
            //TODO: I think I can get rid of this because of Augmented Lagrangian
            //CloseFactor=10e-6;
     
            //Constructing Gradient Pattern
            
            JERows.resize(2*4*cornerPairs.rows()+2*4*edgeCornerPairs.rows());//+xSize);
            JECols.resize(JERows.size());
            JEVals.resize(JERows.size());
            
            if (!isExactMC){
                JCRows.resize(38*edgeCornerPairs.rows()+3*constIndices.size());
                JCCols.resize(JCRows.size());
                JCVals.resize(JCRows.size());
            } else{
                JCRows.resize(38*edgeCornerPairs.rows()+3*constIndices.size()+2*4*cornerPairs.rows());
                JCCols.resize(JCRows.size());
                JCVals.resize(JCRows.size());
            }
            
            /*******************************AMAP Energy********************************************/
            AMAPTriOffset=0;
            AMAPRowOffset=0;
            for (int i=0;i<cornerPairs.rows();i++){
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
            for (int i=0;i<edgeCornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    JERows(RigidTriOffset+2*(4*i+j))=rigidRowOffset+4*i+j;
                    JECols(RigidTriOffset+2*(4*i+j))=4*edgeCornerPairs(i,0)+j;
                    JERows(RigidTriOffset+2*(4*i+j)+1)=rigidRowOffset+4*i+j;
                    JECols(RigidTriOffset+2*(4*i+j)+1)=4*edgeCornerPairs(i,1)+j;
                }
            }

            
            /****************************Closeness Energy*******************/
            /*CloseTriOffset=rigidTriOffset+2*4*edgeCornerPairs.rows();
            CloseRowOffsetr=rigidRowOffset+4*edgeCornerPairs.rows();
            for (int i=0;i<SolutionSize;i++){
                JERows(CloseTriOffset+i)=CloseRowOffset+i;
                JECols(CloseTriOffset+i)=i;
                JEVals(CloseTriOffset+i)=CloseFactor;
            }*/
            
            /****************************Compatibility Constraints*****************/
            compTriOffset=0;
            compRowOffset=0;
            int compTriCounter=compTriOffset;
            Vector4i XiTriPoses; XiTriPoses<<0,8,18,28;
            Vector4i XjTriPoses; XjTriPoses<<4,12,22,32;
            for (int i=0;i<edgeCornerPairs.rows();i++){
                int colCorneri=4*edgeCornerPairs(i,0);
                int colCornerj=4*edgeCornerPairs(i,1);
                int currRowOffset=4*i;
                //derivative of Xi
                quatDerivativeIndices(JCRows, JCCols, compTriCounter, XiTriPoses, compRowOffset+currRowOffset, colCorneri);
                
                //Derivative of Xj
                quatDerivativeIndices(JCRows, JCCols, compTriCounter, XjTriPoses,compRowOffset+currRowOffset, colCornerj);
                
                for (int k=0;k<3;k++){
                    //wj derivative
                    JCRows(compTriCounter+16+10*k)=compRowOffset+currRowOffset+1+k;
                    JCCols(compTriCounter+16+10*k)=4*numCorners+3*edgeCornerVertices(i,1)+k;
                    JCVals(compTriCounter+16+10*k)=-1.0;
                    
                    //wi derivative
                    JCRows(compTriCounter+17+10*k)=compRowOffset+currRowOffset+1+k;
                    JCCols(compTriCounter+17+10*k)=4*numCorners+3*edgeCornerVertices(i,0)+k;
                    JCValues(compTriCounter+17+10*k)=1.0;
                }
                
                compTriCounter+=38;
            }
            
            /****************************Positional Constraints*******************/
            posTriOffset=compTriOffset+38*edgeCornerPairs.rows();
            posRowOffset=compRowOffset+4*edgeCornerPairs.rows();
            
            for (int i=0;i<constIndices.size();i++){
                for (int k=0;k<3;k++){
                    JCRows(posTriOffset+3*i+k)=posRowOffset+3*i+k;
                    JCCols(posTriOffset+3*i+k)=4*numCorners+3*constIndices(i)+k;
                }
            }
            
            
            /****************************Metric-Conformal Constraints*************/
            if (isExactMC){
                MCTriOffset=posTriOffset+3*constIndices.size();
                MCRowOffset=posRowOffset+3*constIndices.size();
                for (int i=0;i<cornerPairs.rows();i++){
                    for (int j=0;j<4;j++){
                        JCRows(MCTriOffset+2*(4*i+j))=MCRowOffset+i;
                        JCCols(MCTriOffset+2*(4*i+j))=4*cornerPairs(i,0)+j;
                        JCRows(MCTriOffset+2*(4*i+j)+1)=MCRowOffset+i;
                        JCCols(MCTriOffset+2*(4*i+j)+1)=4*cornerPairs(i,1)+j;
                    }
                }
            }
        }
        
        void update_energy(const Eigen::VectorXd& x){
            using namespace Eigen;
            currX<<x.head(4*numCorners);
            currLocations<<x.tail(3*VOrigq.rows());
            
            for (int i=0;i<cornerPairs.rows();i++)
                AMAPVec.segment(4*i,4)=(currX.segment(4*cornerPairs(i,1),4)-currX.segment(4*cornerPairs(i,0),4));
            
            AMAPVec.array()*=smoothFactor;
            
            for (int i=0;i<edgeCornerPairs.rows();i++)
                rigidVec.segment(4*i,4)=(currX.segment(4*EdgeCornerPairs(i,1),4)-currX.segment(4*edgeCornerPairs(i,0),4));

            
            rigidVec.array()*=smoothFactor*rigidRatio;
            
            //CloseVec<<CloseFactor*(CurrSolution-InitSolution);
        
            EVec<<AMAPVec, rigidVec;//, CloseVec

        }
        
        void update_constraints(const Eigen::VectorXd& x){
            using namespace Eigen;
            currX<<x.head(4*numCorners);
            currLocations<<x.tail(3*VOrigq.rows());
            
            for (int i=0;i<edgeCornerPairs.rows();i++){
                RowVector4d Xi=currX.segment(4*edgeCornerPairs(i,0),4).transpose();
                RowVector4d Xj=currX.segment(4*edgeCornerPairs(i,1),4).transpose();
                RowVector4d qi=VOrigq.row(edgeCornerVertices(i,0));
                RowVector4d qj=VOrigq.row(edgeCornerVertices(i,1));
                RowVector3d wi=currLocations.segment(3*edgeCornerVertices(i,0),3).transpose();
                RowVector3d wj=currLocations.segment(3*edgeCornerVertices(i,1),3).transpose();
                RowVector4d currEdgeVector=QMult1(QMult1(QConj1(Xi),qj-qi),Xj);
                currEdgeVector.tail(3)-=(wj-wi);
                compVec.segment(4*i,4)<<currEdgeVector.transpose();
            }
            
            
            for (int i=0;i<constIndices.size();i++)
                posVec.segment(3*i,3)<<posFactor*(currLocations.segment(3*constIndices(i),3)-constPoses.row(i).transpose());
            
            if (!isExactMC){
                CVec<<compVec, posVec;
            } else {
                for (int i=0;i<cornerPairs.rows();i++)
                    MCVec(i)=currX.segment(4*cornerPairs(i,1),4).squaredNorm()-currX.segment(4*cornerPairs(i,0),4).squaredNorm();
                CVec<<compVec, posVec, MCVec;
            }
        }
        
        void update_jacobian(const Eigen::VectorXd& x){
            using namespace Eigen;
            currX<<currSolution.head(4*numCorners);
            currLocations<<currSolution.tail(3*VOrigq.rows());
            
            /*******************************AMAP Energy********************************************/
            for (int i=0;i<cornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    JEVals(AMAPTriOffset+2*(4*i+j))=-smoothFactor;
                    JEVals(AMAPTriOffset+2*(4*i+j)+1)=smoothFactor;
                }
            }
            
            /*******************************Rigidity Energy*****************************************/
            for (int i=0;i<edgeCornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    JEVals(rigidTriOffset+2*(4*i+j))=-smoothFactor*rigidRatio;
                    JEVals(rigidTriOffset+2*(4*i+j)+1)=smoothFactor*rigidRatio;
                }
            }
            
            //closeness energy is constant
            
            /****************************Compatibility Constraints*****************/
            int compTriCounter=compTriOffset;
            Vector4i XiTriPoses; XiTriPoses<<0,8,18,28;
            Vector4i XjTriPoses; XjTriPoses<<4,12,22,32;
            for (int i=0;i<edgeCornerPairs.rows();i++){
                RowVector4d Xi=CurrX.segment(4*edgeCornerPairs(i,0),4).transpose();
                RowVector4d Xj=CurrX.segment(4*edgeCornerPairs(i,1),4).transpose();
                RowVector4d rightPart=QMult1(VOrigq.row(edgeCornerVertices(i,1))-VOrigq.row(edgeCornerVertices(i,0)),Xj);
                RowVector4d leftPart=QMult1(QConj1(Xi),VOrigq.row(edgeCornerVertices(i,1))-VOrigq.row(edgeCornerVertices(i,0)));
                
      
                //derivative of Xi
                quatDerivativeValues(JCVals, compTriCounter, XiTriPoses, unitQuat, rightPart, true, false);
                
                //Derivative of Xj
                quatDerivativeValues(JCVals, compTriCounter, XjTriPoses, leftPart, unitQuat, false, false);
                
                //the other compatibility values are constant
                compTriCounter+=38;
            }
        
            
            
            /****************************Positional Constraints*******************/
            for (int i=0;i<constIndices.size();i++)
                for (int k=0;k<3;k++)
                    JCVals(posTriOffset+3*i+k)=posFactor;
            
            /****************************Metric-Conformal Constraints*************/
            if (isExactMC){
                for (int i=0;i<cornerPairs.rows();i++){
                    RowVector4d Xi=currX.segment(4*cornerPairs(i,0),4).transpose();
                    RowVector4d Xj=currX.segment(4*cornerPairs(i,1),4).transpose();
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
