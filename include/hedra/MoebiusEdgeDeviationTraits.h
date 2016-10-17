// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_EDGE_DEVIATION_MOEBIUS_TRAITS_H
#define HEDRA_EDGE_DEVIATION_MOEBIUS_TRAITS_H
#include <igl/igl_inline.h>
#include "quaternionic_derivatives.h"
#include "QuaternionOps.h"
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>
#include <set>


namespace hedra { namespace optimization {
    
    
    //This traits class implements the "UnconstrainedTraits" concept, for the energy of the quaternionic system of deformation in [Vaxman et. al 2015] (see Tutorial).
    class MoebiusEdgeDeviationTraits{
    public:
        
        //concept requirements
        Eigen::VectorXi JRows, JCols;  //rows and column indices for the jacobian matrix
        Eigen::VectorXd JVals;         //values for the jacobian matrix.
        Eigen::VectorXd EVec;          //energy vector
        int xSize;                  //size of the solution
        
        Eigen::VectorXi fullJRows, fullJCols;
        Eigen::VectorXd fullJVals;
        Eigen::VectorXi a2x;            //from the entire set of variables "a" to the free variables in the optimization "x".
        Eigen::VectorXi colMap;         //raw map of a2x into the columns from fullJ*Cols into J*Cols
        
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
        
        Eigen::RowVector4d unitQuat;
        
        //intermediate variables
        Eigen::VectorXd currY;
        Eigen::VectorXd currLocations;
        Eigen::VectorXd AMAPVec;
        Eigen::VectorXd rigidVec;
        
        int AMAPTriOffset, AMAPRowOffset;
        int rigidTriOffset, rigidRowOffset;
        
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
            F=_F; D=_D; EV=_EV;
            constIndices=_constIndices;
            constPoses.conservativeResize(constIndices.size(),3);
            
            unitQuat<<1.0,0.0,0.0,0.0;
            
            xSize=4*VOrigq.rows()+3*(VOrigq.rows()-constIndices.size());
            currY.resize(4*VOrigq.rows());
            currLocations.resize(3*VOrigq.rows());
            
            AMAPVec.resize(4*EV.rows());
            rigidVec.resize(4*EV.rows());
            
            EVec.resize(AMAPVec.size()+rigidVec.size());
            
            //Constructing Gradient Pattern
            
            fullJRows.resize(38*EV.rows()+2*4*EV.rows());
            fullJCols.resize(fullJRows.size());
            fullJVals.resize(fullJRows.size());
            

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
                quatDerivativeIndices(fullJRows, fullJCols, AMAPTriCounter, YiTriPoses, AMAPRowOffset+currRowOffset, coli);
                
                //Derivative of Xj
                quatDerivativeIndices(fullJRows, fullJCols, AMAPTriCounter, YjTriPoses,AMAPRowOffset+currRowOffset, colj);
                
                for (int k=0;k<3;k++){
                    //wj derivative
                    fullJRows(AMAPTriCounter+16+10*k)=AMAPRowOffset+currRowOffset+1+k;
                    fullJCols(AMAPTriCounter+16+10*k)=4*VOrigq.rows()+3*EV(i,1)+k;
                    fullJVals(AMAPTriCounter+16+10*k)=-1.0/pairLength;
                    
                    //wi derivative
                    fullJRows(AMAPTriCounter+17+10*k)=AMAPRowOffset+currRowOffset+1+k;
                    fullJCols(AMAPTriCounter+17+10*k)=4*VOrigq.rows()+3*EV(i,0)+k;
                    fullJVals(AMAPTriCounter+17+10*k)=1.0/pairLength;
                }
                
                AMAPTriCounter+=38;
            }

            
            /*******************************Rigidity Energy*****************************************/
            rigidTriOffset=AMAPTriOffset+38*EV.rows();
            rigidRowOffset=AMAPRowOffset+4*EV.rows();
            for (int i=0;i<EV.rows();i++){
                for (int j=0;j<4;j++){
                    fullJRows(rigidTriOffset+2*(4*i+j))=rigidRowOffset+4*i+j;
                    fullJCols(rigidTriOffset+2*(4*i+j))=4*EV(i,0)+j;
                    fullJRows(rigidTriOffset+2*(4*i+j)+1)=rigidRowOffset+4*i+j;
                    fullJCols(rigidTriOffset+2*(4*i+j)+1)=4*EV(i,1)+j;
                }
            }
            
        
            a2x=Eigen::VectorXi::Zero(VOrigq.rows());
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
            }
            
            //recalibrating prevSolution
            prevSolution.conservativeResize(xSize);
            if (constIndices.size()==0){
                for (int i=0;i<VOrigq.rows();i++)
                    prevSolution.segment(4*i,4)=unitQuat;    //corner variables are trivial
                
                for (int i=0;i<VOrigq.rows();i++){
                    prevSolution.segment(4*VOrigq.rows()+3*i,3)=_VOrig.row(i);
                    fullSolution=_VOrig;
                }
            } else {
                //the corner variables stay the same
                for (int i=0;i<a2x.size();i++){
                    if (a2x(i)!=-1)
                        prevSolution.segment(4*VOrigq.rows()+3*a2x(i),3)=fullSolution.row(i);
                }
            }
            
        }
        
        void update_energy(const Eigen::VectorXd& x){
            using namespace Eigen;
            currY<<x.head(4*VOrigq.rows());
            for (int i=0;i<a2x.size();i++)
                if (a2x(i)!=-1)
                    currLocations.segment(3*i,3)<<x.segment(4*VOrigq.rows()+3*a2x(i),3);
            
            for (int i=0;i<constIndices.size();i++)
                currLocations.segment(3*constIndices(i),3)=constPoses.row(i).transpose();
            
            for (int i=0;i<EV.rows();i++){
                RowVector4d Xi=currY.segment(4*EV(i,0),4).transpose();
                RowVector4d Xj=currY.segment(4*EV(i,1),4).transpose();
                RowVector4d qi=VOrigq.row(EV(i,0));
                RowVector4d qj=VOrigq.row(EV(i,1));
                RowVector3d wi=currLocations.segment(3*EV(i,0),3).transpose();
                RowVector3d wj=currLocations.segment(3*EV(i,1),3).transpose();
                double pairLength=(qj-qi).norm();
                RowVector4d currEdgeVector=QMult(QMult(QConj(Xi),qj-qi),Xj)/pairLength;
                currEdgeVector.tail(3)-=(wj-wi)/pairLength;
                AMAPVec.segment(4*i,4)<<currEdgeVector.transpose();
            }

            AMAPVec.array()*=smoothFactor;
            
            for (int i=0;i<EV.rows();i++)
                rigidVec.segment(4*i,4)=(currY.segment(4*EV(i,1),4)-currY.segment(4*EV(i,0),4));
            
            
            rigidVec.array()*=smoothFactor*rigidRatio;
            
            EVec<<AMAPVec, rigidVec;
            
        }
        
       
        void update_jacobian(const Eigen::VectorXd& x){
            using namespace Eigen;
            currY<<x.head(4*VOrigq.rows());
            for (int i=0;i<a2x.size();i++)
                if (a2x(i)!=-1)
                    currLocations.segment(3*i,3)<<x.segment(4*VOrigq.rows()+3*a2x(i),3);
            
            for (int i=0;i<constIndices.size();i++)
                currLocations.segment(3*constIndices(i),3)=constPoses.row(i).transpose();
            
            /****************************AMAP Energy*****************/
            int AMAPTriCounter=AMAPTriOffset;
            Vector4i YiTriPoses; YiTriPoses<<0,8,18,28;
            Vector4i YjTriPoses; YjTriPoses<<4,12,22,32;
            for (int i=0;i<EV.rows();i++){
                RowVector4d Yi=currY.segment(4*EV(i,0),4).transpose();
                RowVector4d Yj=currY.segment(4*EV(i,1),4).transpose();
                RowVector4d qi=VOrigq.row(EV(i,0));
                RowVector4d qj=VOrigq.row(EV(i,1));
                double pairLength=(qj-qi).norm();
                RowVector4d rightPart=QMult(VOrigq.row(EV(i,1))-VOrigq.row(EV(i,0)),Yj)/pairLength;
                RowVector4d leftPart=QMult(QConj(Yi),VOrigq.row(EV(i,1))-VOrigq.row(EV(i,0)))/pairLength;
                
                
                //derivative of Xi
                quatDerivativeValues(fullJVals, AMAPTriCounter, YiTriPoses, unitQuat, rightPart, true, false);
                
                //Derivative of Xj
                quatDerivativeValues(fullJVals, AMAPTriCounter, YjTriPoses, leftPart, unitQuat, false, false);
                
                //the other compatibility values are constant
                AMAPTriCounter+=38;
            }
            
            
            
            /*******************************Rigidity Energy*****************************************/
            for (int i=0;i<EV.rows();i++){
                for (int j=0;j<4;j++){
                    fullJVals(rigidTriOffset+2*(4*i+j))=-smoothFactor*rigidRatio;
                    fullJVals(rigidTriOffset+2*(4*i+j)+1)=smoothFactor*rigidRatio;
                }
            }
            
            int actualGradCounter=0;
            for (int i=0;i<fullJCols.size();i++)
                if (colMap(fullJCols(i))!=-1)  //not a removed variable
                    JVals(actualGradCounter++)=fullJVals(i);
            
            for (int i=0;i<JVals.size();i++)
                if (isnan(JVals(i)))
                    std::cout<<"nan in JVals("<<i<<")"<<std::endl;
            
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
                    fullSolution.row(i)<<x.segment(4*VOrigq.rows()+3*a2x(i),3).transpose();
            
            for (int i=0;i<constIndices.size();i++)
                fullSolution.row(constIndices(i))=constPoses.row(i);
            
            prevSolution=x;
            
            return true;  //this traits doesn't have any more stop requirements
        }
        
        MoebiusEdgeDeviationTraits(){}
        ~MoebiusEdgeDeviationTraits(){}
        
        
        
        
    };
} }



#endif
