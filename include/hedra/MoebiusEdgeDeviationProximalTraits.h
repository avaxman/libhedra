// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_EDGE_DEVIATION_PROXIMAL_MOEBIUS_TRAITS_H
#define HEDRA_EDGE_DEVIATION_PROXIMAL_MOEBIUS_TRAITS_H
#include <igl/igl_inline.h>
#include "quaternionic_derivatives.h"
#include "QuaternionOps.h"
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>
#include <set>


namespace hedra { namespace optimization {
    
    
    //This traits class implements the "ConstrainedTraits" concept, for the energy of the quaternionic system of deformation in [Vaxman et. al 2015] (see Tutorial).
    class MoebiusEdgeDeviationProximalTraits{
    public:
        
        //concept requirements
        Eigen::VectorXi JERows, JECols;  //rows and column indices for the jacobian matrix
        Eigen::VectorXd JEVals;         //values for the jacobian matrix.
        Eigen::VectorXd EVec;          //energy vector
        int xSize;                  //size of the solution
        
        Eigen::VectorXi JCRows, JCCols;  //rows and column indices for the jacobian matrix
        Eigen::VectorXd JCVals;         //values for the jacobian matrix.
        Eigen::VectorXd CVec;          //energy vector
        
        /*Eigen::VectorXi fullJRows, fullJCols;
        Eigen::VectorXd fullJVals;
        Eigen::VectorXi a2x;            //from the entire set of variables "a" to the free variables in the optimization "x".
        Eigen::VectorXi colMap;         //raw map of a2x into the columns from fullJ*Cols into J*Cols*/
        
        Eigen::MatrixXd VOrigq;
        Eigen::MatrixXi D, F, EV;
        Eigen::VectorXi constIndices;
        Eigen::MatrixXi cornerPairs;
        Eigen::VectorXi cornerOffset;
        
        Eigen::MatrixXd constPoses;
        
        Eigen::VectorXd prevSolution;   //from any previous running of the algorithm (in "raw" optimization format, with the corner variables etc.)
        
        Eigen::MatrixXd fullSolution;   //contains the locations of the vertices after the current run.
        
        double smoothFactor;
        double rigidRatio;
        //double posFactor;
        double closeFactor;
        double prevPosError;
        int bigIter;
        
        Eigen::RowVector4d unitQuat;
        
        //intermediate variables
        Eigen::VectorXd currY;
        Eigen::VectorXd prevY;
        Eigen::VectorXd currLocations;
        Eigen::VectorXd AMAPVec;
        Eigen::VectorXd rigidVec;
        Eigen::VectorXd posVec;
        Eigen::VectorXd closeVec;
        
        int AMAPTriOffset, AMAPRowOffset;
        int rigidTriOffset, rigidRowOffset;
        int posTriOffset, posRowOffset;
        int closeTriOffset, closeRowOffset;
        
        //if constIndices is empty, the initial solution is the original mesh
        void init(const Eigen::MatrixXd& _VOrig,
                  const Eigen::MatrixXi& _D ,
                  const Eigen::MatrixXi& _F ,
                  const Eigen::MatrixXi& _EV,
                  const Eigen::VectorXi& _constIndices=Eigen::VectorXi::Zero(0))
        {
            
            using namespace Eigen;
            
            //imaginary quaternionic representation
            Coords2Quat(_VOrig, VOrigq);
            F=_F; D=_D;
            constIndices=_constIndices;
            constPoses.conservativeResize(constIndices.size(),3);
     
            //enriching list of edges with diagonals
            std::vector<std::pair<int,int> > EVList;
            for (int i=0;i<F.rows();i++)
                for (int j=0;j<D(i);j++)
                    for (int k=j+2;k<D(i)-1;k++){
                        EVList.push_back(std::pair<int, int>(F(i,j), F(i,k)));
                        //std::cout<<"j,k, F(i,j), F(i,k)"<<j<<","<<k<<","<<F(i,j)<<","<<F(i,k)<<std::endl;
                    }
            
            EV.resize(_EV.rows()+EVList.size(),2);
            EV.block(0,0, _EV.rows(),2)=_EV;
            for (int i=0;i<EVList.size();i++)
                EV.row(_EV.rows()+i)<<EVList[i].first, EVList[i].second;
            
            unitQuat<<1.0,0.0,0.0,0.0;            
            xSize=4*VOrigq.rows()+3*(VOrigq.rows());//-constIndices.size());
            currY.resize(4*VOrigq.rows());
            prevY.resize(4*VOrigq.rows());
            currLocations.resize(3*VOrigq.rows());
            
            AMAPVec.resize(4*EV.rows());
            rigidVec.resize(4*EV.rows());
            closeVec.resize(xSize);
            posVec.resize(3*constIndices.rows());
            
            EVec.resize(AMAPVec.size()+rigidVec.size()+closeVec.size());
            CVec.resize(posVec.size());
            
            //Constructing Gradient Pattern
            
            JERows.resize(38*EV.rows()+2*4*EV.rows()+xSize);
            JECols.resize(JERows.size());
            JEVals.resize(JERows.size());
            
            JCRows.resize(3*constIndices.size());
            JCCols.resize(JCRows.size());
            JCVals.resize(JCRows.size());
            

            /*******************************AMAP Energy********************************************/
            AMAPTriOffset=0;
            AMAPRowOffset=0;
            int AMAPTriCounter=AMAPTriOffset;
            Vector4i YiTriPoses; YiTriPoses<<0,8,18,28;
            Vector4i YjTriPoses; YjTriPoses<<4,12,22,32;
            for (int i=0;i<EV.rows();i++){
                int coli=4*EV(i,0);
                int colj=4*EV(i,1);
                int currRowOffset=4*i;
                double pairLength=(VOrigq.row(EV(i,1))-VOrigq.row(EV(i,0))).norm();
                //derivative of Xi
                quatDerivativeIndices(JERows, JECols, AMAPTriCounter, YiTriPoses, AMAPRowOffset+currRowOffset, coli);
                
                //Derivative of Xj
                quatDerivativeIndices(JERows, JECols, AMAPTriCounter, YjTriPoses,AMAPRowOffset+currRowOffset, colj);
                
                for (int k=0;k<3;k++){
                    //wj derivative
                    JERows(AMAPTriCounter+16+10*k)=AMAPRowOffset+currRowOffset+1+k;
                    JECols(AMAPTriCounter+16+10*k)=4*VOrigq.rows()+3*EV(i,1)+k;
                    //JVals(AMAPTriCounter+16+10*k)=-1.0;///pairLength;
                    
                    //wi derivative
                    JERows(AMAPTriCounter+17+10*k)=AMAPRowOffset+currRowOffset+1+k;
                    JECols(AMAPTriCounter+17+10*k)=4*VOrigq.rows()+3*EV(i,0)+k;
                    //JVals(AMAPTriCounter+17+10*k)=1.0;///pairLength;
                }
                
                AMAPTriCounter+=38;
            }

            
            /*******************************Rigidity Energy*****************************************/
            rigidTriOffset=AMAPTriOffset+38*EV.rows();
            rigidRowOffset=AMAPRowOffset+4*EV.rows();
            for (int i=0;i<EV.rows();i++){
                for (int j=0;j<4;j++){
                    JERows(rigidTriOffset+2*(4*i+j))=rigidRowOffset+4*i+j;
                    JECols(rigidTriOffset+2*(4*i+j))=4*EV(i,0)+j;
                    JERows(rigidTriOffset+2*(4*i+j)+1)=rigidRowOffset+4*i+j;
                    JECols(rigidTriOffset+2*(4*i+j)+1)=4*EV(i,1)+j;
                }
            }
            
            /******************************closeness to previous solution***************************/
            closeTriOffset=rigidTriOffset+2*4*EV.rows();
            closeRowOffset=rigidRowOffset+4*EV.rows();
            
            for (int i=0;i<xSize;i++){
                JERows(closeTriOffset+i)=closeRowOffset+i;
                JECols(closeTriOffset+i)=i;
            }
            
            /****************************positional constraints*******************************/
            posTriOffset=0;
            posRowOffset=0;
            
            for (int i=0;i<constIndices.size();i++){
                for (int k=0;k<3;k++){
                    JCRows(posTriOffset+3*i+k)=posRowOffset+3*i+k;
                    JCCols(posTriOffset+3*i+k)=4*VOrigq.rows()+3*constIndices(i)+k;
                }
            }
            
        
            /*a2x=Eigen::VectorXi::Zero(VOrigq.rows());
            int CurrIndex=0;
            for (int i=0;i<constIndices.size();i++)
                a2x(constIndices(i))=-1;
            
            for (int i=0;i<VOrigq.rows();i++)
                if (a2x(i)!=-1)
                    a2x(i)=CurrIndex++;
            
            colMap.resize(4*VOrigq.rows()+3*VOrigq.rows());
            for (int i=0;i<4*VOrigq.rows();i++)
                colMap(i)=i;
            for (int i=0;i<VOrigq.rows();i++)
                if (a2x(i)!=-1)
                    colMap.segment(4*VOrigq.rows()+3*i,3)<<4*VOrigq.rows()+3*a2x(i),4*VOrigq.rows()+3*a2x(i)+1,4*VOrigq.rows()+3*a2x(i)+2;
                else
                    colMap.segment(4*VOrigq.rows()+3*i,3)<<-1,-1,-1;
            
            //setting up the Jacobian rows and columns
            int actualGradCounter=0;
            for (int i=0;i<fullJCols.size();i++){
                if (colMap(fullJCols(i))!=-1)  //not a removed variable
                    actualGradCounter++;
            }
            
            JRows.resize(actualGradCounter);
            JCols.resize(actualGradCounter);
            JVals.resize(actualGradCounter);
            
            actualGradCounter=0;
            for (int i=0;i<fullJCols.size();i++){
                if (colMap(fullJCols(i))!=-1){  //not a removed variable
                    JRows(actualGradCounter)=fullJRows(i);
                    JCols(actualGradCounter)=colMap(fullJCols(i));
                    JVals(actualGradCounter++)=fullJVals(i);
                }
            }*/
            
            //recalibrating prevSolution
            prevSolution.conservativeResize(xSize);
            if (constIndices.size()==0){
                for (int i=0;i<VOrigq.rows();i++)
                    prevSolution.segment(4*i,4)=unitQuat;    //corner variables are trivial
                
                for (int i=0;i<VOrigq.rows();i++){
                    prevSolution.segment(4*VOrigq.rows()+3*i,3)=_VOrig.row(i);
                }
                fullSolution=_VOrig;
                //std::cout<<"Initialization without constraints: "<<std::endl;
            } /*else {
                //the corner variables stay the same
                
                for (int i=0;i<a2x.size();i++){
                    if (a2x(i)!=-1)
                        prevSolution.segment(4*VOrigq.rows()+3*a2x(i),3)=fullSolution.row(i).transpose();
                }
                //std::cout<<"prevSolution: "<<prevSolution<<std::endl;
                //std::cout<<"Initialization with constraints:"<<std::endl;
            }*/
            
        }
        
        void update_energy(const Eigen::VectorXd& x){
            using namespace Eigen;
            currY<<x.head(4*VOrigq.rows());
            prevY<<prevSolution.head(4*VOrigq.rows());
            currLocations<<x.tail(3*VOrigq.rows());
            /*for (int i=0;i<a2x.size();i++)
                if (a2x(i)!=-1)
                    currLocations.segment(3*i,3)<<x.segment(4*VOrigq.rows()+3*a2x(i),3);
            
            for (int i=0;i<constIndices.size();i++)
                currLocations.segment(3*constIndices(i),3)=constPoses.row(i).transpose();*/
            
            for (int i=0;i<EV.rows();i++){
                RowVector4d Yi=prevY.segment(4*EV(i,0),4).transpose();
                RowVector4d Yj=currY.segment(4*EV(i,1),4).transpose();
                RowVector4d qi=VOrigq.row(EV(i,0));
                RowVector4d qj=VOrigq.row(EV(i,1));
                RowVector3d wi=currLocations.segment(3*EV(i,0),3).transpose();
                RowVector3d wj=currLocations.segment(3*EV(i,1),3).transpose();
                double pairLength=(qj-qi).norm();
                RowVector4d currEdgeVector=QMult(QMult(QConj(Yi),qj-qi),Yj);///pairLength;
                currEdgeVector.tail(3)-=(wj-wi);///pairLength;
                AMAPVec.segment(4*i,4)<<currEdgeVector.transpose();
            }

            //AMAPVec.array()*=smoothFactor;
            //std::cout<<"AMAP.maxCoeffs(): "<<AMAPVec.lpNorm<Infinity>()<<std::endl;
            
            for (int i=0;i<EV.rows();i++)
                rigidVec.segment(4*i,4)=(currY.segment(4*EV(i,1),4)-currY.segment(4*EV(i,0),4));
            

            
            //rigidVec.array()*=smoothFactor*rigidRatio;
            //std::cout<<"rigidVec.maxCoeff(): "<<rigidVec.lpNorm<Infinity>()<<std::endl;
            
            //
            
            
            closeVec<<(x-prevSolution);
            
            std::cout<<"rigidVec max: "<<rigidVec.lpNorm<Infinity>()<<std::endl;
            std::cout<<"AMAPVec max: "<<AMAPVec.lpNorm<Infinity>()<<std::endl;
            std::cout<<"closeVec max: "<<closeVec.lpNorm<Infinity>()<<std::endl;
            
            EVec<<smoothFactor*AMAPVec, smoothFactor*rigidRatio*rigidVec, closeFactor*closeVec;
            
        }
        
        void update_constraints(const Eigen::VectorXd& x){
            
            for (int i=0;i<constIndices.size();i++)
                posVec.segment(3*i,3)<<currLocations.segment(3*constIndices(i),3)-constPoses.row(i).transpose();
            
            if (posVec.size()!=0)
                CVec<<posVec;
            
            std::cout<<"posVec max: "<<posVec.lpNorm<Eigen::Infinity>()<<std::endl;
        }
        
       
        void update_jacobian(const Eigen::VectorXd& x){
            using namespace Eigen;
            currY<<x.head(4*VOrigq.rows());
            prevY<<prevSolution.head(4*VOrigq.rows());
            currLocations<<x.tail(3*VOrigq.rows());
            /*for (int i=0;i<a2x.size();i++)
                if (a2x(i)!=-1)
                    currLocations.segment(3*i,3)<<x.segment(4*VOrigq.rows()+3*a2x(i),3);
            
            for (int i=0;i<constIndices.size();i++)
                currLocations.segment(3*constIndices(i),3)=constPoses.row(i).transpose();*/
            
            /****************************AMAP Energy*****************/
            int AMAPTriCounter=AMAPTriOffset;
            Vector4i YiTriPoses; YiTriPoses<<0,8,18,28;
            Vector4i YjTriPoses; YjTriPoses<<4,12,22,32;
            for (int i=0;i<EV.rows();i++){
                RowVector4d Yi=prevY.segment(4*EV(i,0),4).transpose();
                RowVector4d Yj=currY.segment(4*EV(i,1),4).transpose();
                RowVector4d qi=VOrigq.row(EV(i,0));
                RowVector4d qj=VOrigq.row(EV(i,1));
                double pairLength=(qj-qi).norm();
                //std::cout<<"pairLength: "<<pairLength<<std::endl;
                RowVector4d rightPart; rightPart.setZero(); //smoothFactor*QMult(VOrigq.row(EV(i,1))-VOrigq.row(EV(i,0)),Yj)/pairLength;   //this is constant in the iteration
                RowVector4d leftPart=smoothFactor*QMult(QConj(Yi),VOrigq.row(EV(i,1))-VOrigq.row(EV(i,0)));///pairLength;
                
                //derivative of Yi
                quatDerivativeValues(JEVals, AMAPTriCounter, YiTriPoses, unitQuat, rightPart, true, false);
                
                //Derivative of Yj
                quatDerivativeValues(JEVals, AMAPTriCounter, YjTriPoses, leftPart, unitQuat, false, false);
                
                for (int k=0;k<3;k++){
                    JEVals(AMAPTriCounter+16+10*k)=-smoothFactor;///pairLength;
                    JEVals(AMAPTriCounter+17+10*k)=smoothFactor;///pairLength;
                }
                AMAPTriCounter+=38;
            }
            
            
            
            /*******************************Rigidity Energy*****************************************/
            for (int i=0;i<EV.rows();i++){
                for (int j=0;j<4;j++){
                    JEVals(rigidTriOffset+2*(4*i+j))=-smoothFactor*rigidRatio;
                    JEVals(rigidTriOffset+2*(4*i+j)+1)=smoothFactor*rigidRatio;
                }
            }
            
            /*int actualGradCounter=0;
            for (int i=0;i<fullJCols.size();i++)
                if (colMap(fullJCols(i))!=-1)  //not a removed variable
                    JVals(actualGradCounter++)=fullJVals(i);*/
            

            

            /***************************Closeness energy*******************/
            for (int i=0;i<xSize;i++){
                JEVals(closeTriOffset+i)=closeFactor;
            }
            
            /***************************Positional Constraints*******************/
            for (int i=0;i<constIndices.size();i++)
                for (int k=0;k<3;k++)
                    JCVals(posTriOffset+3*i+k)=1.0;
            
            for (int i=0;i<JEVals.size();i++)
                if (isnan(JEVals(i)))
                    std::cout<<"nan in JEVals("<<i<<")"<<std::endl;
            
            for (int i=0;i<JCVals.size();i++)
                if (isnan(JCVals(i)))
                    std::cout<<"nan in JCVals("<<i<<")"<<std::endl;
      
        }
        
        //provide the initial solution to the solver
        void initial_solution(Eigen::VectorXd& x0){
            x0=prevSolution;
        }
        
        void pre_iteration(const Eigen::VectorXd& prevx){prevSolution=prevx; std::cout<<"updating prevY"<<std::endl;}
        bool post_iteration(const Eigen::VectorXd& x){return false;  /*never stop after an iteration*/}
        
        bool post_optimization(const Eigen::VectorXd& x){
            
            fullSolution.conservativeResize(VOrigq.rows(),3);
            
            for (int i=0;i<VOrigq.rows();i++)
                fullSolution.row(i)<<x.segment(4*VOrigq.rows()+3*i,3).transpose();
            
            std::cout<<"updating fullSolution"<<std::endl;
            /*for (int i=0;i<a2x.size();i++)
                if (a2x(i)!=-1)
                    fullSolution.row(i)<<x.segment(4*VOrigq.rows()+3*a2x(i),3).transpose();
            
            for (int i=0;i<constIndices.size();i++)
                fullSolution.row(constIndices(i))=constPoses.row(i);*/
            
            //if the smoothness inner iteration did not converge, let it
            return true;
        }
        
        MoebiusEdgeDeviationProximalTraits(){}
        ~MoebiusEdgeDeviationProximalTraits(){}
        
        
        
        
    };
} }



#endif
