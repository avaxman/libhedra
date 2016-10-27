// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_MOEBIUS_3D_CORNER_VARS_TRAITS_H
#define HEDRA_MOEBIUS_3D_CORNER_VARS_TRAITS_H
#include <igl/igl_inline.h>
#include <hedra/quaternionic_derivatives.h>
#include <hedra/quaternionic_operations.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>
#include <set>


namespace hedra { namespace optimization {
    

    //This traits class implements the "unonstraintTraits" concept, for the energy and constraints of the quaternionic system of deformation in [Vaxman et. al 2015] (see Tutorial).
    //TODO: fully integrate this as a constraint class
    class Moebius3DCornerVarsTraits{
    public:
        
        //concept requirements
        Eigen::VectorXi JRows, JCols;  //rows and column indices for the jacobian matrix
        Eigen::VectorXd JVals;         //values for the jacobian matrix.
        Eigen::VectorXd EVec;          //energy vector
        int xSize;                  //size of the solution
        
        Eigen::MatrixXd origVq;
        Eigen::MatrixXi D, F;
        Eigen::VectorXi constIndices;
        Eigen::MatrixXi cornerPairs;
        Eigen::VectorXi cornerOffset;

        Eigen::MatrixXd constPoses;
        
        Eigen::MatrixXi edgeCornerPairs;  //for rigidity energy
        Eigen::MatrixXi edgeCornerVertices;  //the vertices themselves, for compatibility
        
        Eigen::VectorXd initSolution;   //from any previous running of the algorithm (in "raw" optimization format, with the corner variables etc.)
        
        Eigen::MatrixXd finalPositions;   //contains the positions of the vertices after the current run.
        Eigen::MatrixXd finalX;             //the final corner variables
        
        int numCorners;
        int numEdgePairs;
        
        bool isExactMC;

        double smoothFactor;
        double rigidRatio;
        double closeFactor;
        double posFactor;
        double constTolerance;
        double prevError;
        
        Eigen::RowVector4d unitQuat;
        
        //intermediate variables
        Eigen::VectorXd currX;
        Eigen::VectorXd currLocations;
        Eigen::VectorXd AMAPVec;
        Eigen::VectorXd rigidVec;
        Eigen::VectorXd compVec;
        Eigen::VectorXd closeVec;
        Eigen::VectorXd posVec;
        Eigen::VectorXd MCVec;
        
        Eigen::VectorXd constVec;
 
        int AMAPTriOffset, AMAPRowOffset;
        int rigidTriOffset, rigidRowOffset;
        int closeTriOffset, closeRowOffset;
        int compTriOffset, compRowOffset;
        int posTriOffset, posRowOffset;
        int MCTriOffset, MCRowOffset;

        //if constIndices is empty, the initial solution is the original mesh
        void init(const Eigen::MatrixXd& _origV,
                  const Eigen::MatrixXi& _D ,
                  const Eigen::MatrixXi& _F ,
                  const bool _isExactMC,
                  const Eigen::VectorXi& _constIndices=Eigen::VectorXi::Zero(0))
        {
            
            using namespace Eigen;
            
            //imaginary quaternionic representation
            Coords2Quat(_origV, origVq);
            F=_F; D=_D;
            constIndices=_constIndices;
            constPoses.conservativeResize(constIndices.size(),3);
            
            unitQuat<<1.0,0.0,0.0,0.0;

            //MobConstMat=inMobConstMat;
            isExactMC=_isExactMC;

            //computing relevant topology
            numCorners=D.sum();
            cornerOffset.resize(F.rows());
            cornerOffset(0)=0;
            
            for (int i=1;i<D.rows();i++)
                cornerOffset(i)=cornerOffset(i-1)+D(i-1);
            
            MatrixXi adjCorners(origVq.rows(),12);
            VectorXi valences(origVq.rows()); valences.setZero();
            for (int i=0;i<D.rows();i++){
                for (int j=0;j<D(i);j++){
                    adjCorners(F(i,j),valences(F(i,j)))=cornerOffset(i)+j;
                    valences(F(i,j))++;
                }
            }
            
          
            std::vector<std::pair<int,int> > cornerPairList;
            for (int i=0;i<origVq.rows();i++)
                for (int j=0;j<valences(i)-1;j++)
                    cornerPairList.push_back(std::pair<int,int>(adjCorners(i,j),adjCorners(i,(j+1)%valences(i))));
            
            cornerPairs.resize(cornerPairList.size(),2);
            for (int i=0;i<cornerPairList.size();i++)
                cornerPairs.row(i)<<cornerPairList[i].first, cornerPairList[i].second;
            
            xSize=4*numCorners+3*origVq.rows();
            currX.resize(4*numCorners);
            currLocations.resize(3*origVq.rows());
            
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
            closeVec.resize(xSize);
            posVec.resize(3*constIndices.size());
            if (isExactMC)
                MCVec.resize(cornerPairs.rows());
            else
                MCVec.resize(0);
                
            EVec.resize(AMAPVec.size()+rigidVec.size()+closeVec.size()+compVec.size()+posVec.size()+MCVec.size());
            constVec.resize(compVec.size()+posVec.size()+MCVec.size());
            std::cout<<"EVec size: "<<EVec.size()<<std::endl;
            std::cout<<"constVec size: "<<constVec.size()<<std::endl;
            
            closeFactor=1e-6;
            constTolerance=1e-7;
            
            //Constructing Gradient Pattern
            
            if (!isExactMC)
                JRows.resize(2*4*cornerPairs.rows()+2*4*edgeCornerPairs.rows()+xSize+38*edgeCornerPairs.rows()+3*constIndices.size());
            else
                JRows.resize(2*4*cornerPairs.rows()+2*4*edgeCornerPairs.rows()+xSize+38*edgeCornerPairs.rows()+3*constIndices.size()+2*4*cornerPairs.rows());
            
            
            JCols.resize(JRows.size());
            JVals.resize(JRows.size());
            
            
            /*******************************AMAP Energy********************************************/
            AMAPTriOffset=0;
            AMAPRowOffset=0;
            for (int i=0;i<cornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    JRows(AMAPTriOffset+2*(4*i+j))=AMAPRowOffset+4*i+j;
                    JCols(AMAPTriOffset+2*(4*i+j))=4*cornerPairs(i,0)+j;
                    JRows(AMAPTriOffset+2*(4*i+j)+1)=AMAPRowOffset+4*i+j;
                    JCols(AMAPTriOffset+2*(4*i+j)+1)=4*cornerPairs(i,1)+j;
                }
            }
            
            /*******************************Rigidity Energy*****************************************/
            rigidTriOffset=AMAPTriOffset+2*4*cornerPairs.rows();
            rigidRowOffset=AMAPRowOffset+4*cornerPairs.rows();
            for (int i=0;i<edgeCornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    JRows(rigidTriOffset+2*(4*i+j))=rigidRowOffset+4*i+j;
                    JCols(rigidTriOffset+2*(4*i+j))=4*edgeCornerPairs(i,0)+j;
                    JRows(rigidTriOffset+2*(4*i+j)+1)=rigidRowOffset+4*i+j;
                    JCols(rigidTriOffset+2*(4*i+j)+1)=4*edgeCornerPairs(i,1)+j;
                }
            }
                
            /****************************Closeness Energy*******************/
            closeTriOffset=rigidTriOffset+2*4*edgeCornerPairs.rows();
            closeRowOffset=rigidRowOffset+4*edgeCornerPairs.rows();
            for (int i=0;i<xSize;i++){
                JRows(closeTriOffset+i)=closeRowOffset+i;
                JCols(closeTriOffset+i)=i;
                JVals(closeTriOffset+i)=closeFactor;   //TODO: move that since this might not be constant
            }
                

            
            /****************************Compatibility Constraints*****************/
                
            compTriOffset=closeTriOffset+xSize;
            compRowOffset=closeRowOffset+xSize;;
            int compTriCounter=compTriOffset;
            Vector4i XiTriPoses; XiTriPoses<<0,8,18,28;
            Vector4i XjTriPoses; XjTriPoses<<4,12,22,32;
            for (int i=0;i<edgeCornerPairs.rows();i++){
                int colCorneri=4*edgeCornerPairs(i,0);
                int colCornerj=4*edgeCornerPairs(i,1);
                int currRowOffset=4*i;
                double pairLength=(origVq.row(edgeCornerVertices(i,1))-origVq.row(edgeCornerVertices(i,0))).norm();
                //derivative of Xi
                quatDerivativeIndices(JRows, JCols, compTriCounter, XiTriPoses, compRowOffset+currRowOffset, colCorneri);
                
                //Derivative of Xj
                quatDerivativeIndices(JRows, JCols, compTriCounter, XjTriPoses,compRowOffset+currRowOffset, colCornerj);
                
                for (int k=0;k<3;k++){
                    //wj derivative
                    JRows(compTriCounter+16+10*k)=compRowOffset+currRowOffset+1+k;
                    JCols(compTriCounter+16+10*k)=4*numCorners+3*edgeCornerVertices(i,1)+k;
                    JVals(compTriCounter+16+10*k)=-1.0;///pairLength;
                    
                    //wi derivative
                    JRows(compTriCounter+17+10*k)=compRowOffset+currRowOffset+1+k;
                    JCols(compTriCounter+17+10*k)=4*numCorners+3*edgeCornerVertices(i,0)+k;
                    JVals(compTriCounter+17+10*k)=1.0;///pairLength;
                }
                
                compTriCounter+=38;
            }
            
            /****************************Positional Constraints*******************/
            posTriOffset=compTriOffset+38*edgeCornerPairs.rows();
            posRowOffset=compRowOffset+4*edgeCornerPairs.rows();
            
            for (int i=0;i<constIndices.size();i++){
                for (int k=0;k<3;k++){
                    JRows(posTriOffset+3*i+k)=posRowOffset+3*i+k;
                    JCols(posTriOffset+3*i+k)=4*numCorners+3*constIndices(i)+k;
                }
            }
            
            
            /****************************Metric-Conformal Constraints*************/
            if (isExactMC){
                MCTriOffset=posTriOffset+3*constIndices.size();
                MCRowOffset=posRowOffset+3*constIndices.size();
                for (int i=0;i<cornerPairs.rows();i++){
                    for (int j=0;j<4;j++){
                        JRows(MCTriOffset+2*(4*i+j))=MCRowOffset+i;
                        JCols(MCTriOffset+2*(4*i+j))=4*cornerPairs(i,0)+j;
                        JRows(MCTriOffset+2*(4*i+j)+1)=MCRowOffset+i;
                        JCols(MCTriOffset+2*(4*i+j)+1)=4*cornerPairs(i,1)+j;
                    }
                }
            }
            
        
            
            //calibrating prevSolution
            initSolution.conservativeResize(xSize);
            if (constIndices.size()==0){
                for (int i=0;i<numCorners;i++)
                    initSolution.segment(4*i,4)=unitQuat;    //corner variables are trivial
                
                for (int i=0;i<origVq.rows();i++){
                    initSolution.segment(4*numCorners+3*i,3)=_origV.row(i);
                }
                finalPositions=_origV;
                finalX.resize(numCorners,4);
                finalX.setZero();
                finalX.col(0)=VectorXd::Constant(numCorners, 1.0);
            } else {
                update_constraints(initSolution);
                prevError=constVec.lpNorm<Infinity>();
            }
        }
        
        void update_energy(const Eigen::VectorXd& x){
            using namespace Eigen;
            currX<<x.head(4*numCorners);
            currLocations<<x.tail(3*origVq.rows());

            
            for (int i=0;i<cornerPairs.rows();i++)
                AMAPVec.segment(4*i,4)=(currX.segment(4*cornerPairs(i,1),4)-currX.segment(4*cornerPairs(i,0),4));
            
            for (int i=0;i<edgeCornerPairs.rows();i++)
                rigidVec.segment(4*i,4)=(currX.segment(4*edgeCornerPairs(i,1),4)-currX.segment(4*edgeCornerPairs(i,0),4));
            
            closeVec<<x-initSolution;
            
            update_constraints(x);
            
            EVec<<smoothFactor*AMAPVec, smoothFactor*rigidRatio*rigidVec, closeFactor*closeVec, constVec;

        }
        
        void update_constraints(const Eigen::VectorXd& x){
            using namespace Eigen;
            currX<<x.head(4*numCorners);
            currLocations<<x.tail(3*origVq.rows());
            
            for (int i=0;i<edgeCornerPairs.rows();i++){
                RowVector4d Xi=currX.segment(4*edgeCornerPairs(i,0),4).transpose();
                RowVector4d Xj=currX.segment(4*edgeCornerPairs(i,1),4).transpose();
                RowVector4d qi=origVq.row(edgeCornerVertices(i,0));
                RowVector4d qj=origVq.row(edgeCornerVertices(i,1));
                RowVector3d wi=currLocations.segment(3*edgeCornerVertices(i,0),3).transpose();
                RowVector3d wj=currLocations.segment(3*edgeCornerVertices(i,1),3).transpose();
                double pairLength=(qj-qi).norm();
                RowVector4d currEdgeVector=QMult(QMult(QConj(Xi),qj-qi),Xj);///pairLength;
                currEdgeVector.tail(3)-=(wj-wi);///pairLength;
                compVec.segment(4*i,4)<<currEdgeVector.transpose();
            }
            
            
            for (int i=0;i<constIndices.size();i++)
                posVec.segment(3*i,3)<<(currLocations.segment(3*constIndices(i),3)-constPoses.row(i).transpose());
            
            if (!isExactMC){
                if (posVec.size()!=0)
                    constVec<<compVec, posFactor*posVec;
                else{
                    constVec<<compVec;
                }
            } else {
                for (int i=0;i<cornerPairs.rows();i++)
                    MCVec(i)=currX.segment(4*cornerPairs(i,1),4).squaredNorm()-currX.segment(4*cornerPairs(i,0),4).squaredNorm();
                constVec<<compVec, posFactor*posVec, MCVec;
            }
        }
        
        void update_jacobian(const Eigen::VectorXd& x){
            using namespace Eigen;
            currX<<x.head(4*numCorners);
            currLocations<<x.tail(3*origVq.rows());
            
            /*******************************AMAP Energy********************************************/
            for (int i=0;i<cornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    JVals(AMAPTriOffset+2*(4*i+j))=-smoothFactor;
                    JVals(AMAPTriOffset+2*(4*i+j)+1)=smoothFactor;
                }
            }
            
            /*******************************Rigidity Energy*****************************************/
            for (int i=0;i<edgeCornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    JVals(rigidTriOffset+2*(4*i+j))=-smoothFactor*rigidRatio;
                    JVals(rigidTriOffset+2*(4*i+j)+1)=smoothFactor*rigidRatio;
                }
            }
            
            //closeness energy is constant
            
            
            /****************************Compatibility Constraints*****************/
            int compTriCounter=compTriOffset;
            Vector4i XiTriPoses; XiTriPoses<<0,8,18,28;
            Vector4i XjTriPoses; XjTriPoses<<4,12,22,32;
            for (int i=0;i<edgeCornerPairs.rows();i++){
                RowVector4d Xi=currX.segment(4*edgeCornerPairs(i,0),4).transpose();
                RowVector4d Xj=currX.segment(4*edgeCornerPairs(i,1),4).transpose();
                RowVector4d qi=origVq.row(edgeCornerVertices(i,0));
                RowVector4d qj=origVq.row(edgeCornerVertices(i,1));
                double pairLength=(qj-qi).norm();
                RowVector4d rightPart=QMult(origVq.row(edgeCornerVertices(i,1))-origVq.row(edgeCornerVertices(i,0)),Xj);///pairLength;
                RowVector4d leftPart=QMult(QConj(Xi),origVq.row(edgeCornerVertices(i,1))-origVq.row(edgeCornerVertices(i,0)));///pairLength;
                
      
                //derivative of Xi
                quatDerivativeValues(JVals, compTriCounter, XiTriPoses, unitQuat, rightPart, true, false);
                
                //Derivative of Xj
                quatDerivativeValues(JVals, compTriCounter, XjTriPoses, leftPart, unitQuat, false, false);
                
                //the other compatibility values are constant
                compTriCounter+=38;
            }
        
            
            
            /***************************Positional Constraints*******************/
            for (int i=0;i<constIndices.size();i++)
                for (int k=0;k<3;k++)
                    JVals(posTriOffset+3*i+k)=posFactor;
            
            /****************************Metric-Conformal Constraints*************/
            if (isExactMC){
                for (int i=0;i<cornerPairs.rows();i++){
                    RowVector4d Xi=currX.segment(4*cornerPairs(i,0),4).transpose();
                    RowVector4d Xj=currX.segment(4*cornerPairs(i,1),4).transpose();
                    for (int k=0;k<4;k++){
                        JVals(MCTriOffset+2*(4*i+k))=-2*Xi(k);
                        JVals(MCTriOffset+2*(4*i+k)+1)=2*Xj(k);
                    }
                }
            }
            
            for (int i=0;i<JVals.size();i++)
                if (isnan(JVals(i)))
                    std::cout<<"nan in JVals("<<i<<")"<<std::endl;

        }
        
        //provide the initial solution to the solver
        void initial_solution(Eigen::VectorXd& x0){
            x0=initSolution;
        }
        
        void pre_iteration(const Eigen::VectorXd& prevx){
            initSolution=prevx;
            update_constraints(prevx);
            prevError=constVec.lpNorm<Eigen::Infinity>();
        }
        bool post_iteration(const Eigen::VectorXd& x){
            //when error is halved, the smoothness is reduced by slowest, and when error change is zero, smoothness is halved.
            initSolution=x;
            update_constraints(x);
            double rate=constVec.lpNorm<Eigen::Infinity>()/prevError;
            double reduceRate=std::min(rate/2.0,1.0);
            
            smoothFactor*=0.75;//-0.7*(1.0-reduceRate);
            std::cout<<"smoothFactor: "<<smoothFactor<<std::endl;
           
            
            return (constVec.lpNorm<Eigen::Infinity>()<constTolerance);
        }
        
        bool post_optimization(const Eigen::VectorXd& x){
            
            initSolution=x;
            finalPositions.conservativeResize(origVq.rows(),3);
            for (int i=0;i<origVq.rows();i++)
                finalPositions.row(i)<<x.segment(4*numCorners+3*i,3).transpose();
            
            finalX.conservativeResize(4*numCorners,4);
            for (int i=0;i<numCorners;i++)
                finalX.row(i)<<x.segment(4*i,4).transpose();
            
            update_energy(x);
            double finalTotalError=EVec.lpNorm<Eigen::Infinity>();
            double finalConstError=constVec.lpNorm<Eigen::Infinity>();
            std::cout<<"Final Const Error:"<<finalTotalError<<std::endl;
            std::cout<<"Final Total Error:"<<finalConstError<<std::endl;
            return true;
        }
        
        Moebius3DCornerVarsTraits(){}
        ~Moebius3DCornerVarsTraits(){}

    };
} }



#endif
