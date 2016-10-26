// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_MOEBIUS_2D_EDGE_DEVIATION_TRAITS_H
#define HEDRA_MOEBIUS_2D_EDGE_DEVIATION_TRAITS_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>
#include <set>
#include <iostream>


namespace hedra { namespace optimization {
    

    //This traits class implements the "unonstraintTraits" concept, for the energy and constraints of the complex system of deformation in [Vaxman et. al 2015] (see Tutorial).
    //TODO: fully integrate this as a constraint class
    
    class Moebius2DEdgeDeviationTraits{
    public:
        
        
        //concept requirements
        Eigen::VectorXi JRows, JCols;  //rows and column indices for the jacobian matrix
        Eigen::VectorXd JVals;         //values for the jacobian matrix.
        Eigen::VectorXd EVec;          //energy vector
        int xSize;                  //size of the solution

        typedef std::complex<double> Complex;
        
        Eigen::VectorXcd origVc;
        Eigen::MatrixXi EV, D, F;
        Eigen::VectorXi constIndices;//, VarIndices;
        Eigen::VectorXcd complexConstPoses;
        Eigen::VectorXcd initSolution;
        Eigen::VectorXcd finalPositions;  //final vertex positions
        Eigen::VectorXcd finalY;  //final vertex reciprocals
        Eigen::VectorXcd finalE;  //final edge deivations
        
        double smoothFactor;
        double closeFactor;
        double constTolerance;
        double rigidRatio;
        double posFactor;
        double prevError;
        
        bool isExactMC;
        bool isExactIAP;
        
        //Complex variables
        Eigen::VectorXcd currSolution;
        Eigen::VectorXcd currY;     //vertex reciprocals
        Eigen::VectorXcd currLocations;
        Eigen::VectorXcd currE;    //edge deviations
        Eigen::VectorXcd currEdges;
        Eigen::VectorXcd AMAPVec;
        Eigen::VectorXcd mobVec;
        Eigen::VectorXcd posVec;
        Eigen::VectorXcd rigidVec;
        Eigen::VectorXcd closeVec;
        Eigen::VectorXcd deviationVec;
        Eigen::VectorXcd origFCR;
        
        //real variables
        Eigen::VectorXd MCVec;
        Eigen::VectorXd IAPVec;
        
        
        Eigen::SparseMatrix<Complex> d0;
        
        Eigen::VectorXi complexJRows;
        Eigen::VectorXi complexJCols;
        Eigen::VectorXcd complexJVals;
        
        //into the complex values
        int AMAPTriOffset, AMAPRowOffset;
        int rigidTriOffset, rigidRowOffset;
        int posTriOffset, posRowOffset;
        int mobTriOffset, mobRowOffset;
        int closeTriOffset, closeRowOffset;
        int deviationTriOffset, deviationRowOffset;
        int complexRowOffset, complexTriOffset;
        
        //into the real values
        int MCTriOffset, MCRowOffset;
        int IAPTriOffset, IAPRowOffset;
        
        Eigen::VectorXcd constVec;
        

        void init(const Eigen::VectorXcd& _origVc,
                  const Eigen::MatrixXi& _D,
                  const Eigen::MatrixXi& _F,
                  const Eigen::MatrixXi& _EV,
                  bool _isExactMC,
                  bool _isExactIAP,
                  const Eigen::VectorXi& _constIndices=Eigen::VectorXi::Zero(0)){
            
            using namespace Eigen;
            using namespace std;
            
            F=_F; D=_D; EV=_EV;
            constIndices=_constIndices;
            origVc=_origVc;
            isExactMC=_isExactMC;
            isExactIAP=_isExactIAP;
            
            
            xSize=2*(origVc.rows()+origVc.rows());
            
            if (isExactMC || isExactIAP)
                xSize+=2*EV.rows();
            
            closeFactor=10e-6;
            constTolerance=10e-9;
            
            //creating difference operator
            d0.resize(EV.rows(), origVc.rows());
            d0.reserve(2*EV.rows());
            vector<Triplet<Complex> > d0Tris(2*EV.rows());
            for (int i=0;i<EV.rows();i++){
                d0Tris[2*i]=Triplet<Complex>(i,EV(i,0), -1.0);
                d0Tris[2*i+1]=Triplet<Complex>(i,EV(i,1), 1.0);
            }
            
            d0.setFromTriplets(d0Tris.begin(), d0Tris.end());
            
            //Allocating intermediate and output vectors
            AMAPVec.resize(EV.rows());
            rigidVec.resize(EV.rows());
            posVec.resize(constIndices.size());
            closeVec.resize(xSize);
            mobVec.resize(D.sum()-3*D.rows());
            if (isExactMC || isExactIAP)
                deviationVec.resize(EV.rows());
            if (isExactMC){
                MCVec.resize(EV.rows());
                IAPVec.resize(0);
            } else if (isExactIAP){
                IAPVec.resize(EV.rows());
                MCVec.resize(0);
            } else {
                IAPVec.resize(0);
                MCVec.resize(0);
            }
            currSolution.resize(xSize/2);
            currY.resize(origVc.rows());
            currLocations.resize(origVc.rows());
            currE.resize(EV.rows());
            currEdges.resize(EV.rows());
            
            
            //original face-based cross ratios
            origFCR.resize(D.sum()-3*D.size());
            int mobCounter=0;
            for (int i=0;i<D.rows();i++){
                for (int j=0;j<D(i)-3;j++){
                    Complex zi=origVc(F(i,j));
                    Complex zj=origVc(F(i,j+1));
                    Complex zk=origVc(F(i,j+2));
                    Complex zl=origVc(F(i,j+3));
                    origFCR(mobCounter++)=(zj-zi)/(zk-zj)*(zl-zk)/(zi-zl);
                }
            }
            
            if (!isExactMC && !isExactIAP)
                EVec.resize(2*(AMAPVec.size()+rigidVec.size()+closeVec.size()+posVec.size()+mobVec.size()));
            else
                EVec.resize(2*(AMAPVec.size()+rigidVec.size()+closeVec.size()+posVec.size()+mobVec.size()+deviationVec.size())+MCVec.size()+IAPVec.size());
            
            
            if (!isExactMC && !isExactIAP)
                constVec.resize(posVec.size()+mobVec.size());
            else
                constVec.resize(posVec.size()+mobVec.size()+deviationVec.size()+MCVec.size()+IAPVec.size());
                                       
            
            /**************************************************Creating Gradient Pattern*******************************************/
            
            if (!isExactMC && !isExactIAP){
                complexJRows.resize(4*EV.rows()+2*EV.rows()+xSize/2+constIndices.size()+4*origFCR.size());
                complexJCols.resize(complexJRows.size());
                complexJVals.resize(complexJRows.size());
                
                JRows.resize(4*complexJRows.size());
                JCols.resize(4*complexJRows.size());
                JVals.resize(4*complexJRows.size());
            } else {
                complexJRows.resize(4*EV.rows()+2*EV.rows()+xSize/2+constIndices.size()+4*origFCR.size()+5*EV.rows());
                complexJCols.resize(complexJRows.size());
                complexJVals.resize(complexJRows.size());
                
                if (isExactMC){
                    JRows.resize(4*complexJRows.size()+2*EV.rows());
                    JCols.resize(4*complexJRows.size()+2*EV.rows());
                    JVals.resize(4*complexJRows.size()+2*EV.rows());
                } else{
                    JRows.resize(4*complexJRows.size()+EV.rows());
                    JCols.resize(4*complexJRows.size()+EV.rows());
                    JVals.resize(4*complexJRows.size()+EV.rows());
                }
            }
            

            /*********************AMAPEnergy**********************/
            //Xj*Xi*zij part
            AMAPTriOffset=0;
            AMAPRowOffset=0;
            for (int i=0;i<EV.rows();i++){
                complexJRows(AMAPTriOffset+4*i)=AMAPRowOffset+i;
                complexJCols(AMAPTriOffset+4*i)=EV(i,0);
                
                complexJRows(AMAPTriOffset+4*i+1)=AMAPRowOffset+i;
                complexJCols(AMAPTriOffset+4*i+1)=EV(i,1);
                
                complexJRows(AMAPTriOffset+4*i+2)=AMAPRowOffset+i;
                complexJCols(AMAPTriOffset+4*i+2)=origVc.rows()+EV(i,0);
                
                complexJRows(AMAPTriOffset+4*i+3)=AMAPRowOffset+i;
                complexJCols(AMAPTriOffset+4*i+3)=origVc.rows()+EV(i,1);
                
                
            }
            
            /*******************Rigidity Energy*******************/
            rigidTriOffset=AMAPTriOffset+4*EV.rows();
            rigidRowOffset=AMAPRowOffset+EV.rows();
            
            for (int i=0;i<EV.rows();i++){
                complexJRows(rigidTriOffset+2*i)=rigidRowOffset+i;
                complexJCols(rigidTriOffset+2*i)=EV(i,0);
                
                complexJRows(rigidTriOffset+2*i+1)=rigidRowOffset+i;
                complexJCols(rigidTriOffset+2*i+1)=EV(i,1);

            }
            
            /************************Closeness Energy******************/
            closeTriOffset=rigidTriOffset+2*EV.rows();
            closeRowOffset=rigidRowOffset+EV.rows();
            for (int i=0;i<xSize/2;i++){
                complexJRows(closeTriOffset+i)=closeRowOffset+i;
                complexJCols(closeTriOffset+i)=i;
                complexJVals(closeTriOffset+i)=closeFactor;
            }
            
            /******************Positional Soft Constraints************/
            
            posTriOffset=closeTriOffset+xSize/2;
            posRowOffset=closeRowOffset+xSize/2;
            for (int i=0;i<constIndices.size();i++){
                complexJRows(posTriOffset+i)=posRowOffset+i;
                complexJCols(posTriOffset+i)=origVc.rows()+constIndices(i);
            }
            
            /*****************Mobius-Equivalent Constraints********/
            
            mobTriOffset=posTriOffset+constIndices.size();
            mobRowOffset=posRowOffset+constIndices.size();
            int mobTriCounter=0;
            for (int i=0;i<D.rows();i++){
                for (int j=0;j<D(i)-3;j++){
                    for (int k=0;k<4;k++){
                        complexJRows(mobTriOffset+4*mobTriCounter+k)=mobRowOffset+mobTriCounter;
                        complexJCols(mobTriOffset+4*mobTriCounter+k)=origVc.rows()+F(i,j+k);
                    }
                    mobTriCounter++;
                }
            }
            
            complexRowOffset=mobRowOffset+origFCR.size();
            complexTriOffset=mobTriOffset+4*origFCR.size();

            if (isExactMC || isExactIAP){  //deviation constraints are included
                
                deviationTriOffset=complexTriOffset;
                deviationRowOffset=complexRowOffset;
                for (int i=0;i<EV.rows();i++){
                    complexJRows(deviationTriOffset+5*i)=deviationRowOffset+i;
                    complexJCols(deviationTriOffset+5*i)=EV(i,0);
                    
                    complexJRows(deviationTriOffset+5*i+1)=deviationRowOffset+i;
                    complexJCols(deviationTriOffset+5*i+1)=EV(i,1);
                    
                    complexJRows(deviationTriOffset+5*i+2)=deviationRowOffset+i;
                    complexJCols(deviationTriOffset+5*i+2)=origVc.rows()+EV(i,0);
                    
                    complexJRows(deviationTriOffset+5*i+3)=deviationRowOffset+i;
                    complexJCols(deviationTriOffset+5*i+3)=origVc.rows()+EV(i,1);
                    
                    complexJRows(deviationTriOffset+5*i+4)=deviationRowOffset+i;
                    complexJCols(deviationTriOffset+5*i+4)=2*origVc.rows()+i;
                }
                
                complexRowOffset+=EV.rows();
                complexTriOffset+=5*EV.rows();
            }
            
            //creating the real-valued pattern and adding MC\IAP constraints
            //[Real -imag; imag real]
            for (int i=0;i<complexJRows.size();i++){
                
                //real upper left
                JRows(2*i)=complexJRows(i);
                JCols(2*i)=complexJCols(i);
                
                //-imag upper right
                JRows(2*i+1)=complexJRows(i);
                JCols(2*i+1)=xSize/2+complexJCols(i);
                
                //imag lower left
                JRows(2*i+2*complexJRows.size())=complexRowOffset+complexJRows(i);
                JCols(2*i+2*complexJRows.size())=complexJCols(i);
                
                //real lower right
                JRows(2*i+2*complexJRows.size()+1)=complexRowOffset+complexJRows(i);
                JCols(2*i+2*complexJRows.size()+1)=xSize/2+complexJCols(i);
            }
            
            
            
            /*************************Metric-Conformal and Intersection-Angle Preserving Constraints********************************/
            if (isExactMC){
                MCRowOffset=complexRowOffset*2;
                MCTriOffset=4*complexJRows.size();
                
                
                //actual values are updated in the gradient function
                for (int i=0;i<EV.rows();i++){
                    JRows(MCTriOffset+2*i)=MCRowOffset+i;
                    JCols(MCTriOffset+2*i)=2*origVc.rows()+i;  //real part of e_ij
                    
                    JRows(MCTriOffset+2*i+1)=MCRowOffset+i;
                    JCols(MCTriOffset+2*i+1)=xSize/2+2*origVc.rows()+i;  //imaginary part of e_ij
                }
            } else if (isExactIAP){
                IAPRowOffset=2*complexRowOffset;
                IAPTriOffset=4*complexJRows.size();
                for (int i=0;i<EV.rows();i++){
                    JRows(IAPTriOffset+i)=IAPRowOffset+i;
                    JCols(IAPTriOffset+i)=xSize/2+2*origVc.rows()+i;  //imaginary part of e_ij
                    JVals(IAPTriOffset+i)=1.0;
                }
            }
            
            //calibrating initSolution
            initSolution.conservativeResize(xSize);
            if (constIndices.size()==0){
                initSolution.head(origVc.rows())=VectorXcd::Zero(origVc.rows());
                initSolution.segment(origVc.rows(),origVc.rows())=origVc;
                finalPositions=_origVc;
                finalY=VectorXcd::Constant(origVc.rows(), Complex(1.0,0.0));
                if (isExactMC || isExactIAP)
                    finalE=VectorXcd::Constant(EV.rows(), Complex(1.0,0.0));
                    
                prevError=0.0;
            } else {
                //update_constraints(initSolution);
                //prevError=constVec.lpNorm<Infinity>();
            }
            
        }
        
        void update_energy(const Eigen::VectorXd& x){
            
            using namespace Eigen;
            using namespace std;
            
            currSolution.array().real()<<x.head(x.size()/2);
            currSolution.array().imag()<<x.tail(x.size()/2);
            
            currY<<currSolution.head(origVc.rows());
            currLocations<<currSolution.segment(origVc.rows(),origVc.rows());
            
            for (int i=0;i<EV.rows();i++)
                AMAPVec(i)=(currY(EV(i,0))*(origVc(EV(i,1))-origVc(EV(i,0)))*currY(EV(i,1))-(currLocations(EV(i,1))-currLocations(EV(i,0))));
            
            rigidVec<<d0*currY;
            closeVec<<(currSolution-initSolution);
        
            update_constraints(x);
            
            if (!isExactMC && !isExactIAP)
                EVec<<smoothFactor*AMAPVec.real(), smoothFactor*rigidRatio*rigidVec.real(), closeFactor*closeVec.real(), posFactor*posVec.real(), mobVec.real(),
                          smoothFactor*AMAPVec.imag(), smoothFactor*rigidRatio*rigidVec.imag(), closeFactor*closeVec.imag(), posFactor*posVec.imag(), mobVec.imag();
            
            else
                EVec<<smoothFactor*AMAPVec.real(), smoothFactor*rigidRatio*rigidVec.real(), closeFactor*closeVec.real(), posFactor*posVec.real(), mobVec.real(), deviationVec.real(),
                          smoothFactor*AMAPVec.imag(), smoothFactor*rigidRatio*rigidVec.imag(), closeFactor*closeVec.imag(), posFactor*posVec.imag(), mobVec.imag(), deviationVec.imag(),MCVec, IAPVec;
            
        }
        
        void update_constraints(const Eigen::VectorXd& x)
        {
            currSolution.array().real()<<x.head(x.size()/2);
            currSolution.array().imag()<<x.tail(x.size()/2);
            
            currY<<currSolution.head(origVc.rows());
            currLocations<<currSolution.segment(origVc.rows(),origVc.rows());
            
            for (int i=0;i<constIndices.size();i++)
                posVec(i)=currLocations(constIndices(i))-complexConstPoses(i);
            
            int mobCounter=0;
            for (int i=0;i<D.rows();i++){
                for (int j=0;j<D(i)-3;j++){
                    Complex wi=currLocations(F(i,j));
                    Complex wj=currLocations(F(i,j+1));
                    Complex wk=currLocations(F(i,j+2));
                    Complex wl=currLocations(F(i,j+3));
                    mobVec(mobCounter)=(wj-wi)*(wl-wk)-origFCR(mobCounter)*(wi-wl)*(wk-wj);
                    mobCounter++;
                }
            }
            
            if (!isExactMC && !isExactIAP){
                constVec<<posVec, mobVec;
            } else{
                //Deviation constraints
                currE<<currSolution.segment(origVc.rows()+origVc.rows(),EV.rows());
                currEdges<<d0*currLocations;
                for (int i=0;i<EV.rows();i++)
                    deviationVec(i)=currY(EV(i,0))*(origVc(EV(i,1))-origVc(EV(i,0)))*currY(EV(i,1))-currE(i)*currEdges(i);
                
                if (isExactMC)
                    MCVec<<currE.array().abs2()-1;
                if (isExactIAP)
                    IAPVec<<currE.array().imag();
                
                constVec<<posVec, mobVec, deviationVec, MCVec.cast<Complex>(), IAPVec.cast<Complex>();
            }
            
        }

        
        void update_jacobian(const Eigen::VectorXd& x){
            
            using namespace Eigen;
            using namespace std;
            
            currSolution.array().real()<<x.head(x.size()/2);
            currSolution.array().imag()<<x.tail(x.size()/2);
            
            currY<<currSolution.head(origVc.rows());
            currLocations<<currSolution.segment(origVc.rows(),origVc.rows());

            /*********************AMAPEnergy**********************/
            for (int i=0;i<EV.rows();i++){
                //Yi Derivative
                complexJVals(4*i)=smoothFactor*currY(EV(i,1))*(origVc(EV(i,1))-origVc(EV(i,0)));
                //Yj Derivative
                complexJVals(4*i+1)=smoothFactor*currY(EV(i,0))*(origVc(EV(i,1))-origVc(EV(i,0)));
                //wi derivative
                complexJVals(4*i+2)=smoothFactor;
                //wj derivative
                complexJVals(4*i+3)=-smoothFactor;
            }
            
            /*******************Rigidity Energy******************/
            
            for (int i=0;i<EV.rows();i++){
                //Yi derivative
                complexJVals(rigidTriOffset+2*i)=-smoothFactor*rigidRatio;
                //Yj derivative
                complexJVals(rigidTriOffset+2*i+1)=smoothFactor*rigidRatio;
            }
            
            //Closeness gradient is constant
            
            /********************Positional soft Constraints**************************/
            for (int i=0;i<constIndices.size();i++){
                complexJVals(posTriOffset+i)=posFactor;
            }
            
            /*****************Mobius-Equivalent Constraints********/
            
            //the actual values are updated in the gradient
            int mobCounter=0;
            for (int i=0;i<D.rows();i++){
                for (int j=0;j<D(i)-3;j++){
                    Complex wi=currLocations(F(i,j));
                    Complex wj=currLocations(F(i,j+1));
                    Complex wk=currLocations(F(i,j+2));
                    Complex wl=currLocations(F(i,j+3));
                    //MobVec(Counter++)=(Zj-Zi)*(Zl-Zk)-FCR(MobCounter)*(Zi-Zl)*(Zk-Zj);
                    
                    //wi derivative
                    complexJVals(mobTriOffset+4*mobCounter)=-(wl-wk)-origFCR(mobCounter)*(wk-wj);
                    //wj derivative
                    complexJVals(mobTriOffset+4*mobCounter+1)=(wl-wk)+origFCR(mobCounter)*(wi-wl);
                    //wk derivative
                    complexJVals(mobTriOffset+4*mobCounter+2)=-(wj-wi)-origFCR(mobCounter)*(wi-wl);
                    //wl derivative
                    complexJVals(mobTriOffset+4*mobCounter+3)=(wj-wi)+origFCR(mobCounter)*(wk-wj);
                    
                    mobCounter++;
                    
                }
            }
            
            
            if (isExactMC || isExactIAP){
                /*****************Deviation constraints***************/
                currE<<currSolution.segment(origVc.rows()+origVc.rows(),EV.rows());
                for (int i=0;i<EV.rows();i++){
                    //Yi derivative
                    complexJVals(deviationTriOffset+5*i)=currY(EV(i,1))*(origVc(EV(i,1))-origVc(EV(i,0)));
                    //Yj derivative
                    complexJVals(deviationTriOffset+5*i+1)=currY(EV(i,0))*(origVc(EV(i,1))-origVc(EV(i,0)));
        
                    //wi derivative
                    complexJVals(deviationTriOffset+5*i+2)=currE(i);
                    
                    //wj derivative
                    complexJVals(deviationTriOffset+5*i+3)=-currE(i);
                    
                    //E_ij derivatve
                    complexJVals(deviationTriOffset+5*i+4)=-(currLocations(EV(i,1))-currLocations(EV(i,0)));
                }
            }
            
            //updating real values from complex ones
            for (int i=0;i<complexJRows.size();i++){
                JVals(2*i)=   complexJVals(i).real();
                JVals(2*i+1)=-complexJVals(i).imag();
                JVals(2*i+2*complexJRows.size())= complexJVals(i).imag();
                JVals(2*i+1+2*complexJRows.size())= complexJVals(i).real();
            }
            
            
            
            /*************************Metric-Conformal and Intersection-Angle Preserving Constraints********************************/
            if (isExactMC){
                for (int i=0;i<EV.rows();i++){
                    JVals(MCTriOffset+2*i)=2.0*currE(i).real();
                    JVals(MCTriOffset+2*i+1)=2.0*currE(i).imag();
                    
                }
            }
            
            //IAP constraint jacobian are constant
            
        }
        
        void initial_solution(Eigen::VectorXd& x0){
            x0.conservativeResize(xSize);
            x0<<initSolution.real(), initSolution.imag();
        }
        
        void pre_iteration(const Eigen::VectorXd& prevx){
            initSolution.array().real()=prevx.head(prevx.size()/2);
            initSolution.array().imag()=prevx.tail(prevx.size()/2);
            update_constraints(prevx);
            prevError=constVec.lpNorm<Eigen::Infinity>();
        }
        bool post_iteration(const Eigen::VectorXd& x){
            //when error is halved, the smoothness is reduced by slowest, and when error change is zero, smoothness is halved.
            initSolution.array().real()=x.head(x.size()/2);
            initSolution.array().imag()=x.tail(x.size()/2);
            update_constraints(x);
            double rate=constVec.lpNorm<Eigen::Infinity>()/prevError;
            double reduceRate=std::min(rate/2.0,1.0);
            
            smoothFactor*=0.9;//-0.7*(1.0-reduceRate);
            std::cout<<"smoothFactor: "<<smoothFactor<<std::endl;
            return (constVec.lpNorm<Eigen::Infinity>()<constTolerance);
        }
        
        bool post_optimization(const Eigen::VectorXd& x){
            initSolution.array().real()=x.head(x.size()/2);
            initSolution.array().imag()=x.tail(x.size()/2);
   
            finalPositions=initSolution.segment(origVc.rows(),origVc.rows());
            finalY=initSolution.head(origVc.rows());
            if (isExactMC || isExactIAP)
                finalE=initSolution.tail(EV.rows());
          
            update_energy(x);
            double finalTotalError=EVec.lpNorm<Eigen::Infinity>();
            double finalConstError=constVec.lpNorm<Eigen::Infinity>();
            std::cout<<"Final Const Error:"<<finalTotalError<<std::endl;
            std::cout<<"Final Total Error:"<<finalConstError<<std::endl;
            return true;
        }
        
        Moebius2DEdgeDeviationTraits(){}
        ~Moebius2DEdgeDeviationTraits(){}
    };
    
}}



#endif
