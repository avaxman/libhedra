// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_CORNER_VARS_MOEBIUS_TRAITS_H
#define HEDRA_CORNER_VARS_MOEBIUS_TRAITS_H
#include <igl/igl_inline.h>
#include "quaternionic_derivatives.h"
#include "QuaternionOps.h"
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>
#include <set>


namespace hedra { namespace optimization {
    

    //This traits class implements the "ConstraintTraits" concept, for the energy and constraints of the quaternionic system of deformation in [Vaxman et. al 2015] (see Tutorial).
    class MoebiusCornerVarsTraits{
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
        
        Eigen::VectorXi fullJERows, fullJECols;
        Eigen::VectorXd fullJEVals;
        Eigen::VectorXi fullJCRows, fullJCCols;
        Eigen::VectorXd fullJCVals;
        Eigen::VectorXi a2x;            //from the entire set of variables "a" to the free variables in the optimization "x".
        Eigen::VectorXi colMap;         //raw map of a2x into the columns from fullJ*Cols into J*Cols
        
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
        double rigidRatio;
        
        Eigen::RowVector4d unitQuat;
        double averageEdgeLength;
        
        //intermediate variables
        Eigen::VectorXd currX;
        Eigen::VectorXd currLocations;
        Eigen::VectorXd AMAPVec;
        Eigen::VectorXd rigidVec;
        Eigen::VectorXd compVec;
        //VectorXd closeVec;
        //Eigen::VectorXd posVec;
        Eigen::VectorXd MCVec;
 
        int AMAPTriOffset, AMAPRowOffset;
        int rigidTriOffset, rigidRowOffset;
        int compTriOffset, compRowOffset;
        int posTriOffset, posRowOffset;
        int MCTriOffset, MCRowOffset;

        
        //if constIndices is empty, the initial solution is the original mesh
        void init(const Eigen::MatrixXd& _VOrig,
                  const Eigen::MatrixXi& _D ,
                  const Eigen::MatrixXi& _F ,
                  const bool _isExactMC,
                  const Eigen::VectorXi& _constIndices=Eigen::VectorXi::Zero(0))
        {
            
            using namespace Eigen;
            
            //imaginary quaternionic representation
            Coords2Quat(_VOrig, VOrigq);
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
            
            MatrixXi adjCorners(VOrigq.rows(),12);
            VectorXi valences(VOrigq.rows()); valences.setZero();
            for (int i=0;i<D.rows();i++){
                for (int j=0;j<D(i);j++){
                    adjCorners(F(i,j),valences(F(i,j)))=cornerOffset(i)+j;
                    valences(F(i,j))++;
                }
            }
            
            averageEdgeLength=0.0;
            for (int i=0;i<D.rows();i++)
                for (int j=0;j<D(i);j++)
                    averageEdgeLength+=(VOrigq.row(F(i,j))-VOrigq.row(F(i,(j+1)%D(i)))).norm();
            
            averageEdgeLength/=(double)D.sum();

            std::vector<std::pair<int,int> > cornerPairList;
            for (int i=0;i<VOrigq.rows();i++)
                for (int j=0;j<valences(i)-1;j++)
                    cornerPairList.push_back(std::pair<int,int>(adjCorners(i,j),adjCorners(i,(j+1)%valences(i))));
            
            cornerPairs.resize(cornerPairList.size(),2);
            for (int i=0;i<cornerPairList.size();i++)
                cornerPairs.row(i)<<cornerPairList[i].first, cornerPairList[i].second;
            
            xSize=4*numCorners+3*(VOrigq.rows()-constIndices.size());
            currX.resize(4*numCorners);
            currLocations.resize(3*VOrigq.rows());
            
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
            //posVec.resize(3*constIndices.size());
            if (isExactMC)
                MCVec.resize(cornerPairs.rows());
            else
                MCVec.resize(0);
            
            EVec.resize(AMAPVec.size()+rigidVec.size());//+CloseVec.size());
            CVec.resize(compVec.size()+MCVec.size()); //posVec.size()+
            std::cout<<"EVec size: "<<EVec.size()<<std::endl;
            std::cout<<"CVec size: "<<CVec.size()<<std::endl;
            
            //Constructing Gradient Pattern
            
            fullJERows.resize(2*4*cornerPairs.rows()+2*4*edgeCornerPairs.rows());//+xSize);
            fullJECols.resize(fullJERows.size());
            fullJEVals.resize(fullJERows.size());
            
            if (!isExactMC){
                fullJCRows.resize(38*edgeCornerPairs.rows());//+3*constIndices.size());
                fullJCCols.resize(fullJCRows.size());
                fullJCVals.resize(fullJCRows.size());
            } else{
                fullJCRows.resize(38*edgeCornerPairs.rows()+2*4*cornerPairs.rows());//+3*constIndices.size());
                fullJCCols.resize(fullJCRows.size());
                fullJCVals.resize(fullJCRows.size());
            }
            
            /*******************************AMAP Energy********************************************/
            AMAPTriOffset=0;
            AMAPRowOffset=0;
            for (int i=0;i<cornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    fullJERows(AMAPTriOffset+2*(4*i+j))=AMAPRowOffset+4*i+j;
                    fullJECols(AMAPTriOffset+2*(4*i+j))=4*cornerPairs(i,0)+j;
                    fullJERows(AMAPTriOffset+2*(4*i+j)+1)=AMAPRowOffset+4*i+j;
                    fullJECols(AMAPTriOffset+2*(4*i+j)+1)=4*cornerPairs(i,1)+j;
                }
            }
            
            /*******************************Rigidity Energy*****************************************/
            rigidTriOffset=AMAPTriOffset+2*4*cornerPairs.rows();
            rigidRowOffset=AMAPRowOffset+4*cornerPairs.rows();
            for (int i=0;i<edgeCornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    fullJERows(rigidTriOffset+2*(4*i+j))=rigidRowOffset+4*i+j;
                    fullJECols(rigidTriOffset+2*(4*i+j))=4*edgeCornerPairs(i,0)+j;
                    fullJERows(rigidTriOffset+2*(4*i+j)+1)=rigidRowOffset+4*i+j;
                    fullJECols(rigidTriOffset+2*(4*i+j)+1)=4*edgeCornerPairs(i,1)+j;
                }
            }

            
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
                double pairLength=(VOrigq.row(edgeCornerVertices(i,1))-VOrigq.row(edgeCornerVertices(i,0))).norm();
                //derivative of Xi
                quatDerivativeIndices(fullJCRows, fullJCCols, compTriCounter, XiTriPoses, compRowOffset+currRowOffset, colCorneri);
                
                //Derivative of Xj
                quatDerivativeIndices(fullJCRows, fullJCCols, compTriCounter, XjTriPoses,compRowOffset+currRowOffset, colCornerj);
                
                for (int k=0;k<3;k++){
                    //wj derivative
                    fullJCRows(compTriCounter+16+10*k)=compRowOffset+currRowOffset+1+k;
                    fullJCCols(compTriCounter+16+10*k)=4*numCorners+3*edgeCornerVertices(i,1)+k;
                    fullJCVals(compTriCounter+16+10*k)=-1.0/pairLength;
                    
                    //wi derivative
                    fullJCRows(compTriCounter+17+10*k)=compRowOffset+currRowOffset+1+k;
                    fullJCCols(compTriCounter+17+10*k)=4*numCorners+3*edgeCornerVertices(i,0)+k;
                    fullJCVals(compTriCounter+17+10*k)=1.0/pairLength;
                }
                
                compTriCounter+=38;
            }
            
            /****************************Positional Constraints*******************/
            /*posTriOffset=compTriOffset+38*edgeCornerPairs.rows();
            posRowOffset=compRowOffset+4*edgeCornerPairs.rows();
            
            for (int i=0;i<constIndices.size();i++){
                for (int k=0;k<3;k++){
                    fullJCRows(posTriOffset+3*i+k)=posRowOffset+3*i+k;
                    fullJCCols(posTriOffset+3*i+k)=4*numCorners+3*constIndices(i)+k;
                }
            }*/
            
            
            /****************************Metric-Conformal Constraints*************/
            if (isExactMC){
                MCTriOffset=compTriOffset+38*edgeCornerPairs.rows(); //posTriOffset+3*constIndices.size();
                MCRowOffset=compRowOffset+4*edgeCornerPairs.rows(); //posRowOffset+3*constIndices.size();
                for (int i=0;i<cornerPairs.rows();i++){
                    for (int j=0;j<4;j++){
                        fullJCRows(MCTriOffset+2*(4*i+j))=MCRowOffset+i;
                        fullJCCols(MCTriOffset+2*(4*i+j))=4*cornerPairs(i,0)+j;
                        fullJCRows(MCTriOffset+2*(4*i+j)+1)=MCRowOffset+i;
                        fullJCCols(MCTriOffset+2*(4*i+j)+1)=4*cornerPairs(i,1)+j;
                    }
                }
            }
            
            a2x=Eigen::VectorXi::Zero(VOrigq.rows());
            int CurrIndex=0;
            for (int i=0;i<constIndices.size();i++)
                a2x(constIndices(i))=-1;
            
            for (int i=0;i<VOrigq.rows();i++)
                if (a2x(i)!=-1)
                    a2x(i)=CurrIndex++;
            
            colMap.resize(4*numCorners+3*VOrigq.rows());
            for (int i=0;i<4*numCorners;i++)
                colMap(i)=i;
            for (int i=0;i<VOrigq.rows();i++)
                if (a2x(i)!=-1)
                    colMap.segment(4*numCorners+3*i,3)<<4*numCorners+3*a2x(i),4*numCorners+3*a2x(i)+1,4*numCorners+3*a2x(i)+2;
                else
                    colMap.segment(4*numCorners+3*i,3)<<-1,-1,-1;
            
            //setting up the Jacobian rows and columns
            int actualGradCounter=0;
            for (int i=0;i<fullJECols.size();i++){
                if (colMap(fullJECols(i))!=-1)  //not a removed variable
                    actualGradCounter++;
            }
            
            JERows.resize(actualGradCounter);
            JECols.resize(actualGradCounter);
            JEVals.resize(actualGradCounter);
            
            actualGradCounter=0;
            for (int i=0;i<fullJCCols.size();i++){
                if (colMap(fullJCCols(i))!=-1)  //not a removed variable
                    actualGradCounter++;
            }
            
            JCRows.resize(actualGradCounter);
            JCCols.resize(actualGradCounter);
            JCVals.resize(actualGradCounter);

            actualGradCounter=0;
            for (int i=0;i<fullJECols.size();i++){
                if (colMap(fullJECols(i))!=-1){  //not a removed variable
                    JERows(actualGradCounter)=fullJERows(i);
                    JECols(actualGradCounter)=colMap(fullJECols(i));
                    JEVals(actualGradCounter++)=fullJEVals(i);
                }
            }
            
            actualGradCounter=0;
            for (int i=0;i<fullJCCols.size();i++){
                if (colMap(fullJCCols(i))!=-1){  //not a removed variable
                    JCRows(actualGradCounter)=fullJCRows(i);
                    JCCols(actualGradCounter)=colMap(fullJCCols(i));
                    JCVals(actualGradCounter++)=fullJCVals(i);
                }
            }
            
            //recalibrating prevSolution
            prevSolution.conservativeResize(xSize);
            if (constIndices.size()==0){
                for (int i=0;i<numCorners;i++)
                    prevSolution.segment(4*i,4)=unitQuat;    //corner variables are trivial
                
                for (int i=0;i<VOrigq.rows();i++){
                    prevSolution.segment(4*numCorners+3*i,3)=_VOrig.row(i);
                    fullSolution=_VOrig;
                }
            } else {
                //the corner variables stay the same
                for (int i=0;i<a2x.size();i++){
                    if (a2x(i)!=-1)
                        prevSolution.segment(4*numCorners+3*a2x(i),3)=fullSolution.row(i);
                }
            }

        }
        
        void update_energy(const Eigen::VectorXd& x){
            using namespace Eigen;
            currX<<x.head(4*numCorners);
            for (int i=0;i<a2x.size();i++)
                if (a2x(i)!=-1)
                    currLocations.segment(3*i,3)<<x.segment(4*numCorners+3*a2x(i),3);
            
            for (int i=0;i<constIndices.size();i++)
                currLocations.segment(3*constIndices(i),3)=constPoses.row(i).transpose();
            
            for (int i=0;i<cornerPairs.rows();i++)
                AMAPVec.segment(4*i,4)=(currX.segment(4*cornerPairs(i,1),4)-currX.segment(4*cornerPairs(i,0),4));
            
            AMAPVec.array()*=smoothFactor;
            
            for (int i=0;i<edgeCornerPairs.rows();i++)
                rigidVec.segment(4*i,4)=(currX.segment(4*edgeCornerPairs(i,1),4)-currX.segment(4*edgeCornerPairs(i,0),4));

            
            rigidVec.array()*=smoothFactor*rigidRatio;
            
            EVec<<AMAPVec, rigidVec;//, CloseVec

        }
        
        void update_constraints(const Eigen::VectorXd& x){
            using namespace Eigen;
            currX<<x.head(4*numCorners);
            for (int i=0;i<a2x.size();i++)
                if (a2x(i)!=-1)
                    currLocations.segment(3*i,3)<<x.segment(4*numCorners+3*a2x(i),3);
            
            for (int i=0;i<constIndices.size();i++)
                currLocations.segment(3*constIndices(i),3)=constPoses.row(i).transpose();
    
            for (int i=0;i<edgeCornerPairs.rows();i++){
                RowVector4d Xi=currX.segment(4*edgeCornerPairs(i,0),4).transpose();
                RowVector4d Xj=currX.segment(4*edgeCornerPairs(i,1),4).transpose();
                RowVector4d qi=VOrigq.row(edgeCornerVertices(i,0));
                RowVector4d qj=VOrigq.row(edgeCornerVertices(i,1));
                RowVector3d wi=currLocations.segment(3*edgeCornerVertices(i,0),3).transpose();
                RowVector3d wj=currLocations.segment(3*edgeCornerVertices(i,1),3).transpose();
                double pairLength=(qj-qi).norm();
                RowVector4d currEdgeVector=QMult(QMult(QConj(Xi),qj-qi),Xj)/pairLength;
                currEdgeVector.tail(3)-=(wj-wi)/pairLength;
                compVec.segment(4*i,4)<<currEdgeVector.transpose();
            }
            
            
            //for (int i=0;i<constIndices.size();i++)
            //    posVec.segment(3*i,3)<<(currLocations.segment(3*constIndices(i),3)-constPoses.row(i).transpose())/averageEdgeLength;
            
            if (!isExactMC){
                CVec<<compVec;//11, posVec;
            } else {
                for (int i=0;i<cornerPairs.rows();i++)
                    MCVec(i)=currX.segment(4*cornerPairs(i,1),4).squaredNorm()-currX.segment(4*cornerPairs(i,0),4).squaredNorm();
                CVec<<compVec, /*posVec,*/ MCVec;
            }
        }
        
        void update_jacobian(const Eigen::VectorXd& x){
            using namespace Eigen;
            currX<<x.head(4*numCorners);
            for (int i=0;i<a2x.size();i++)
                if (a2x(i)!=-1)
                    currLocations.segment(3*i,3)<<x.segment(4*numCorners+3*a2x(i),3);
            
            for (int i=0;i<constIndices.size();i++)
                currLocations.segment(3*constIndices(i),3)=constPoses.row(i).transpose();
            
            /*******************************AMAP Energy********************************************/
            for (int i=0;i<cornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    fullJEVals(AMAPTriOffset+2*(4*i+j))=-smoothFactor;
                    fullJEVals(AMAPTriOffset+2*(4*i+j)+1)=smoothFactor;
                }
            }
            
            /*******************************Rigidity Energy*****************************************/
            for (int i=0;i<edgeCornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    fullJEVals(rigidTriOffset+2*(4*i+j))=-smoothFactor*rigidRatio;
                    fullJEVals(rigidTriOffset+2*(4*i+j)+1)=smoothFactor*rigidRatio;
                }
            }
            
            
            /****************************Compatibility Constraints*****************/
            int compTriCounter=compTriOffset;
            Vector4i XiTriPoses; XiTriPoses<<0,8,18,28;
            Vector4i XjTriPoses; XjTriPoses<<4,12,22,32;
            for (int i=0;i<edgeCornerPairs.rows();i++){
                RowVector4d Xi=currX.segment(4*edgeCornerPairs(i,0),4).transpose();
                RowVector4d Xj=currX.segment(4*edgeCornerPairs(i,1),4).transpose();
                RowVector4d qi=VOrigq.row(edgeCornerVertices(i,0));
                RowVector4d qj=VOrigq.row(edgeCornerVertices(i,1));
                double pairLength=(qj-qi).norm();
                RowVector4d rightPart=QMult(VOrigq.row(edgeCornerVertices(i,1))-VOrigq.row(edgeCornerVertices(i,0)),Xj)/pairLength;
                RowVector4d leftPart=QMult(QConj(Xi),VOrigq.row(edgeCornerVertices(i,1))-VOrigq.row(edgeCornerVertices(i,0)))/pairLength;
                
      
                //derivative of Xi
                quatDerivativeValues(fullJCVals, compTriCounter, XiTriPoses, unitQuat, rightPart, true, false);
                
                //Derivative of Xj
                quatDerivativeValues(fullJCVals, compTriCounter, XjTriPoses, leftPart, unitQuat, false, false);
                
                //the other compatibility values are constant
                compTriCounter+=38;
            }
        
            
            
            /***************************Positional Constraints*******************/
            /*for (int i=0;i<constIndices.size();i++)
                for (int k=0;k<3;k++)
                    fullJCVals(posTriOffset+3*i+k)=1.0/averageEdgeLength;*/
            
            /****************************Metric-Conformal Constraints*************/
            if (isExactMC){
                for (int i=0;i<cornerPairs.rows();i++){
                    RowVector4d Xi=currX.segment(4*cornerPairs(i,0),4).transpose();
                    RowVector4d Xj=currX.segment(4*cornerPairs(i,1),4).transpose();
                    for (int k=0;k<4;k++){
                        fullJCVals(MCTriOffset+2*(4*i+k))=-2*Xi(k);
                        fullJCVals(MCTriOffset+2*(4*i+k)+1)=2*Xj(k);
                    }
                }
            }
            
            int actualGradCounter=0;
            for (int i=0;i<fullJECols.size();i++)
                if (colMap(fullJECols(i))!=-1)  //not a removed variable
                    JEVals(actualGradCounter++)=fullJEVals(i);
            
            for (int i=0;i<JEVals.size();i++)
                if (isnan(JEVals(i)))
                    std::cout<<"nan in JEVals("<<i<<")"<<std::endl;
            
            actualGradCounter=0;
            for (int i=0;i<fullJCCols.size();i++)
                if (colMap(fullJCCols(i))!=-1)  //not a removed variable
                    JCVals(actualGradCounter++)=fullJCVals(i);
            
            for (int i=0;i<JCVals.size();i++)
                if (isnan(JCVals(i)))
                    std::cout<<"nan in JEVals("<<i<<")"<<std::endl;

        }
        
        //provide the initial solution to the solver
        void initial_solution(Eigen::VectorXd& x0){
            x0=prevSolution;
        }
        
        void pre_iteration(const Eigen::VectorXd& prevx){}
        bool post_iteration(const Eigen::VectorXd& x){return false;  /*never stop after an iteration*/}
        
        bool post_optimization(const Eigen::VectorXd& x){
            
            fullSolution.conservativeResize(VOrigq.rows(),3);
            for (int i=0;i<a2x.size();i++)
                if (a2x(i)!=-1)
                    fullSolution.row(i)<<x.segment(4*numCorners+3*a2x(i),3).transpose();
            
            for (int i=0;i<constIndices.size();i++)
                fullSolution.row(constIndices(i))=constPoses.row(i);
            
            double prevError; update_constraints(prevSolution);
            prevError=CVec.lpNorm<Eigen::Infinity>();
            double currError; update_constraints(x);
            currError=CVec.lpNorm<Eigen::Infinity>();
            std::cout<<"prevError, currError: "<<prevError<<","<<currError<<std::endl;
            prevSolution=x;
            return true;  //this traits doesn't have any more stop requirements
        }
        
        MoebiusCornerVarsTraits(){}
        ~MoebiusCornerVarsTraits(){}

        
        
        
    };
} }



#endif
