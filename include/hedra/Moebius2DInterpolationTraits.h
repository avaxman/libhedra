// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_MOEBIUS_2D_INTERPOLATION_TRAITS_H
#define HEDRA_MOEBIUS_2D_INTERPOLATION_TRAITS_H
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
    
    class Moebius2DInterpolationTraits{
    public:
        
        //concept requirements
        Eigen::VectorXi JRows, JCols;  //rows and column indices for the jacobian matrix
        Eigen::VectorXd JVals;         //values for the jacobian matrix.
        Eigen::VectorXd EVec;          //energy vector
        int xSize;                  //size of the solution
        
        typedef std::complex<double> Complex;

        Eigen::VectorXcd origVc;
        Eigen::MatrixXi D, F;
        Eigen::VectorXi innerEdges;
        Eigen::MatrixXi faceCornerIndices;  //rows of f1, f2, zi, zj (oriented)
        Eigen::VectorXcd presJumps;
        Eigen::Vector2cd firstcd;
        Eigen::VectorXcd initSolution;
        Eigen::VectorXcd currSolution;
        Eigen::VectorXcd finalcd;
        bool isExactMC;
        
        double smoothFactor;
        double closeFactor;
        double constTolerance;
        double prevError;

        Eigen::VectorXcd presVec;
        Eigen::VectorXcd compVec;
        Eigen::VectorXcd closeVec;
        Eigen::Vector2cd firstcdVec;
        Eigen::VectorXcd constVec;
        Eigen::VectorXd MCVec;
        
        Eigen::VectorXi complexJRows, complexJCols;
        Eigen::VectorXcd complexJVals;
        
        int presTriOffset, presRowOffset;
        int closeTriOffset, closeRowOffset;
        int compTriOffset, compRowOffset;
        int firstcdTriOffset, firstcdRowOffset;
        int MCTriOffset,MCRowOffset;
        int complexTriOffset, complexRowOffset;
        
        
        void init(const Eigen::VectorXcd& _origVc,
                  const Eigen::MatrixXi& _D,
                  const Eigen::MatrixXi& _F,
                  const Eigen::MatrixXi& EV,
                  const Eigen::MatrixXi& EF,
                  const Eigen::VectorXi& innerEdges,
                  const bool& _isExactMC){
            
            using namespace Eigen;
            using namespace std;
            
            F=_F; D=_D;
            origVc=_origVc;
            isExactMC=_isExactMC;

            closeFactor=10e-6;
            constTolerance=10e-7;
            
            xSize=4*F.rows();
            
            firstcd<<Complex(0.0), Complex(1.0);
            
            faceCornerIndices.resize(innerEdges.size(),4);
            for (int i=0;i<innerEdges.size();i++)
                faceCornerIndices.row(i)<<EF.row(innerEdges(i)), EV(innerEdges(i),1), EV(innerEdges(i),0);
            
            if (isExactMC)
                MCVec.resize(2*faceCornerIndices.rows());
            else
                MCVec.resize(0);
            
            closeVec.resize(xSize/2);
            presVec.resize(2*faceCornerIndices.rows());
            compVec.resize(faceCornerIndices.rows());
            constVec.resize(compVec.size()+firstcdVec.size()+MCVec.size());
            currSolution.resize(xSize/2);
            
            if (!isExactMC)
                EVec.resize(2*(closeVec.size()+presVec.size()+compVec.size()+firstcdVec.size()));
            else
                EVec.resize(2*(closeVec.size()+presVec.size()+compVec.size()+firstcdVec.size())+MCVec.size());
            
            
            /************************Creating Gradient Pattern*******************/
            
            complexJRows.resize(8*faceCornerIndices.rows()+xSize/2+4*faceCornerIndices.rows()+2);
            complexJCols.resize(complexJRows.size());
            complexJVals.resize(complexJRows.size());
            
             if (isExactMC){
                 JRows.resize(4*complexJRows.size()+16*faceCornerIndices.rows());
                 JCols.resize(4*complexJRows.size()+16*faceCornerIndices.rows());
                 JVals.resize(4*complexJRows.size()+16*faceCornerIndices.rows());
             } else{
                 JRows.resize(4*complexJRows.size());
                JCols.resize(4*complexJRows.size());
                JVals.resize(4*complexJRows.size());
             }
          
            /******************Prescription Energy**********************/
            presTriOffset=0;
            presRowOffset=0;
            for (int i=0;i<faceCornerIndices.rows();i++){
                complexJRows(presTriOffset+8*i)=presRowOffset+2*i;
                complexJCols(presTriOffset+8*i)=2*faceCornerIndices(i,0);
                
                complexJRows(presTriOffset+8*i+1)=presRowOffset+2*i;
                complexJCols(presTriOffset+8*i+1)=2*faceCornerIndices(i,0)+1;
                
                complexJRows(presTriOffset+8*i+2)=presRowOffset+2*i;
                complexJCols(presTriOffset+8*i+2)=2*faceCornerIndices(i,1);
                
                complexJRows(presTriOffset+8*i+3)=presRowOffset+2*i;
                complexJCols(presTriOffset+8*i+3)=2*faceCornerIndices(i,1)+1;
                
                complexJRows(presTriOffset+8*i+4)=presRowOffset+2*i+1;
                complexJCols(presTriOffset+8*i+4)=2*faceCornerIndices(i,0);
                
                complexJRows(presTriOffset+8*i+5)=presRowOffset+2*i+1;
                complexJCols(presTriOffset+8*i+5)=2*faceCornerIndices(i,0)+1;
                
                complexJRows(presTriOffset+8*i+6)=presRowOffset+2*i+1;
                complexJCols(presTriOffset+8*i+6)=2*faceCornerIndices(i,1);
                
                complexJRows(presTriOffset+8*i+7)=presRowOffset+2*i+1;
                complexJCols(presTriOffset+8*i+7)=2*faceCornerIndices(i,1)+1;
                
                
            }
            
            /*********************Closeness Energy**************************/
            closeTriOffset=presTriOffset+8*faceCornerIndices.rows();
            closeRowOffset=presRowOffset+2*faceCornerIndices.rows();
            for (int i=0;i<xSize/2;i++){
                complexJRows(closeTriOffset+i)=closeRowOffset+i;
                complexJCols(closeTriOffset+i)=i;
                complexJVals(closeTriOffset+i)=closeFactor;
            }
            
            /**********************Compatibility Constraint*****************************/
            compTriOffset=closeTriOffset+xSize/2;
            compRowOffset=closeRowOffset+xSize/2;
            for (int i=0;i<faceCornerIndices.rows();i++){
                complexJRows(compTriOffset+4*i)=compRowOffset+i;
                complexJCols(compTriOffset+4*i)=2*faceCornerIndices(i,0);
                
                complexJRows(compTriOffset+4*i+1)=compRowOffset+i;
                complexJCols(compTriOffset+4*i+1)=2*faceCornerIndices(i,0)+1;
                
                complexJRows(compTriOffset+4*i+2)=compRowOffset+i;
                complexJCols(compTriOffset+4*i+2)=2*faceCornerIndices(i,1);
                
                complexJRows(compTriOffset+4*i+3)=compRowOffset+i;
                complexJCols(compTriOffset+4*i+3)=2*faceCornerIndices(i,1)+1;
                
            }
            
            /************************Firstcd Constraint***********************/
            firstcdTriOffset=compTriOffset+4*faceCornerIndices.rows();
            firstcdRowOffset=compRowOffset+faceCornerIndices.rows();
            
            complexJRows(firstcdTriOffset)=firstcdRowOffset;
            complexJCols(firstcdTriOffset)=0;
            complexJVals(firstcdTriOffset)=1.0;
            
            complexJRows(firstcdTriOffset+1)=firstcdRowOffset+1;
            complexJCols(firstcdTriOffset+1)=1;
            complexJVals(firstcdTriOffset+1)=1.0;
            
            
            complexRowOffset=firstcdRowOffset+2;
            complexTriOffset=firstcdTriOffset+2;
            
    
            //creating the real-valued pattern and adding MC\IAP constraints
            //[Real -imag; imag real]
            for (int i=0;i<complexJRows.size();i++){
                
                //real upper left
                JRows(2*i)=complexJRows(i);
                JCols(2*i)=complexJCols(i);
                JVals(2*i)=complexJVals(i).real();
                
                //-imag upper right
                JRows(2*i+1)=complexJRows(i);
                JCols(2*i+1)=xSize/2+complexJCols(i);
                JVals(2*i+1)=-complexJVals(i).imag();
                
                //imag lower left
                JRows(2*i+2*complexJRows.size())=complexRowOffset+complexJRows(i);
                JCols(2*i+2*complexJRows.size())=complexJCols(i);
                JVals(2*i+2*complexJRows.size())=complexJVals(i).imag();
                
                //real lower right
                JRows(2*i+2*complexJRows.size()+1)=complexRowOffset+complexJRows(i);
                JCols(2*i+2*complexJRows.size()+1)=xSize/2+complexJCols(i);
                JVals(2*i+2*complexJRows.size()+1)=complexJVals(i).real();
            }
            
            if (isExactMC){
                MCRowOffset=complexRowOffset*2;
                MCTriOffset=4*complexJRows.size();
                
                for (int i=0;i<faceCornerIndices.rows();i++){
                    for (int j=0;j<2;j++){
                        JRows(MCTriOffset+16*i+8*j)=MCRowOffset+2*i+j;
                        JCols(MCTriOffset+16*i+8*j)=2*faceCornerIndices(i,0);
                        
                        JRows(MCTriOffset+16*i+1+8*j)=MCRowOffset+2*i+j;
                        JCols(MCTriOffset+16*i+1+8*j)=2*F.rows()+2*faceCornerIndices(i,0);
                        
                        JRows(MCTriOffset+16*i+2+8*j)=MCRowOffset+2*i+j;
                        JCols(MCTriOffset+16*i+2+8*j)=2*faceCornerIndices(i,0)+1;
                        
                        JRows(MCTriOffset+16*i+3+8*j)=MCRowOffset+2*i+j;
                        JCols(MCTriOffset+16*i+3+8*j)=2*F.rows()+2*faceCornerIndices(i,0)+1;
                        
                        JRows(MCTriOffset+16*i+4+8*j)=MCRowOffset+2*i+j;
                        JCols(MCTriOffset+16*i+4+8*j)=2*faceCornerIndices(i,1);
                        
                        JRows(MCTriOffset+16*i+5+8*j)=MCRowOffset+2*i+j;
                        JCols(MCTriOffset+16*i+5+8*j)=2*F.rows()+2*faceCornerIndices(i,1);
                        
                        JRows(MCTriOffset+16*i+6+8*j)=MCRowOffset+2*i+j;
                        JCols(MCTriOffset+16*i+6+8*j)=2*faceCornerIndices(i,1)+1;
                        
                        JRows(MCTriOffset+16*i+7+8*j)=MCRowOffset+2*i+j;
                        JCols(MCTriOffset+16*i+7+8*j)=2*F.rows()+2*faceCornerIndices(i,1)+1;
                    }
                }
                
            }
        }
        
        void update_energy(const Eigen::VectorXd& x){
            
            using namespace Eigen;
            using namespace std;
            
            currSolution.array().real()<<x.head(x.size()/2);
            currSolution.array().imag()<<x.tail(x.size()/2);
            
            for (int i=0;i<faceCornerIndices.rows();i++){
                Complex c1=currSolution(2*faceCornerIndices(i,0));
                Complex d1=currSolution(2*faceCornerIndices(i,0)+1);
                Complex c2=currSolution(2*faceCornerIndices(i,1));
                Complex d2=currSolution(2*faceCornerIndices(i,1)+1);
                Complex zi=origVc(faceCornerIndices(i,2));
                Complex zj=origVc(faceCornerIndices(i,3));
                Complex g=presJumps(i);
                presVec(2*i)=smoothFactor*((c1*zi+d1)*g-(c2*zi+d2));
                presVec(2*i+1)=smoothFactor*((c2*zj+d2)*g-(c1*zj+d1));
                
                compVec(i)=(c1*zi+d1)*(c1*zj+d1)-(c2*zi+d2)*(c2*zj+d2);
            }
            
            closeVec<<closeFactor*(currSolution-initSolution);
            
            update_constraints(x);
            
            if (!isExactMC)
                EVec<<presVec.real(), closeVec.real(), compVec.real(), firstcdVec.real(),
                presVec.imag(), closeVec.imag(), compVec.imag(), firstcdVec.imag();
            else
                EVec<<presVec.real(), closeVec.real(), compVec.real(), firstcdVec.real(),
                presVec.imag(), closeVec.imag(), compVec.imag(), firstcdVec.imag(),  MCVec;
            
        }
        void update_constraints(const Eigen::VectorXd& x)
        {
            using namespace Eigen;
            using namespace std;
            VectorXcd currSolution(xSize/2);
            currSolution.array().real()=x.head(x.size()/2);
            currSolution.array().imag()=x.tail(x.size()/2);
            
            for (int i=0;i<faceCornerIndices.rows();i++){
                Complex c1=currSolution(2*faceCornerIndices(i,0));
                Complex d1=currSolution(2*faceCornerIndices(i,0)+1);
                Complex c2=currSolution(2*faceCornerIndices(i,1));
                Complex d2=currSolution(2*faceCornerIndices(i,1)+1);
                Complex zi=origVc(faceCornerIndices(i,2));
                Complex zj=origVc(faceCornerIndices(i,3));
                
                compVec(i)=(c1*zi+d1)*(c1*zj+d1)-(c2*zj+d2)*(c2*zi+d2);
            }
            
            firstcdVec=currSolution.head(2)-firstcd;
            
            if (isExactMC){
                for (int i=0;i<faceCornerIndices.rows();i++){
                    Complex c1=currSolution(2*faceCornerIndices(i,0));
                    Complex d1=currSolution(2*faceCornerIndices(i,0)+1);
                    Complex c2=currSolution(2*faceCornerIndices(i,1));
                    Complex d2=currSolution(2*faceCornerIndices(i,1)+1);
                    Complex zi=origVc(faceCornerIndices(i,2));
                    Complex zj=origVc(faceCornerIndices(i,3));
                    double gabs2=presJumps(i).real()*presJumps(i).real()+presJumps(i).imag()*presJumps(i).imag();
                    
                    MCVec(2*i)=norm(c1*zi+d1)*gabs2-norm(c2*zi+d2);
                    MCVec(2*i+1)=norm(c1*zj+d1)-norm(c2*zj+d2) * gabs2;
                }
                constVec<<compVec, firstcdVec, MCVec.cast<Complex>();
            } else
                constVec<<compVec, firstcdVec;
        }

        
        void update_jacobian(const Eigen::VectorXd& x){
            
            using namespace Eigen;
            using namespace std;
            
            currSolution.array().real()<<x.head(x.size()/2);
            currSolution.array().imag()<<x.tail(x.size()/2);
            
            
            //Energy and compatibility
            for (int i=0;i<faceCornerIndices.rows();i++){
                
                Complex c1=currSolution(2*faceCornerIndices(i,0));
                Complex d1=currSolution(2*faceCornerIndices(i,0)+1);
                Complex c2=currSolution(2*faceCornerIndices(i,1));
                Complex d2=currSolution(2*faceCornerIndices(i,1)+1);
                Complex zi=origVc(faceCornerIndices(i,2));
                Complex zj=origVc(faceCornerIndices(i,3));
                Complex g=presJumps(i);
                
                complexJVals(presTriOffset+8*i)=smoothFactor*zi*g;
                complexJVals(presTriOffset+8*i+1)=smoothFactor*g;
                complexJVals(presTriOffset+8*i+2)=-smoothFactor*zi;
                complexJVals(presTriOffset+8*i+3)=-smoothFactor*1.0;
                
                complexJVals(presTriOffset+8*i+4)=-smoothFactor*zj;
                complexJVals(presTriOffset+8*i+5)=-smoothFactor*1.0;
                complexJVals(presTriOffset+8*i+6)=smoothFactor*g*zj;
                complexJVals(presTriOffset+8*i+7)=smoothFactor*g;
                
                
                complexJVals(compTriOffset+4*i)=zi*(c1*zj+d1)+zj*(c1*zi+d1);
                complexJVals(compTriOffset+4*i+1)=(c1*zi+d1)+(c1*zj+d1);
                complexJVals(compTriOffset+4*i+2)=-(zj*(c2*zi+d2)+zi*(c2*zj+d2));
                complexJVals(compTriOffset+4*i+3)=-((c2*zi+d2)+(c2*zj+d2));
                
            }
            
            
            //Firstcd and closeness are constant
            
            //updating real values from complex ones
            for (int i=0;i<complexJRows.size();i++){
                JVals(2*i)=   complexJVals(i).real();
                JVals(2*i+1)=-complexJVals(i).imag();
                JVals(2*i+2*complexJRows.size())= complexJVals(i).imag();
                JVals(2*i+1+2*complexJRows.size())= complexJVals(i).real();
            }
            
            if (isExactMC){
                //[ 2*cx*zx^2 + 2*dx*zx + 2*cx*zy^2 + 2*dy*zy,
                //2*cy*zx^2 + 2*dy*zx + 2*cy*zy^2 - 2*dx*zy,
                //2*dx + 2*cx*zx - 2*cy*zy,
                //2*dy + 2*cx*zy + 2*cy*zx]
                for (int i=0;i<faceCornerIndices.rows();i++){
                    Complex c1=currSolution(2*faceCornerIndices(i,0));
                    Complex d1=currSolution(2*faceCornerIndices(i,0)+1);
                    Complex c2=currSolution(2*faceCornerIndices(i,1));
                    Complex d2=currSolution(2*faceCornerIndices(i,1)+1);
                    Complex zi=origVc(faceCornerIndices(i,2));
                    Complex zj=origVc(faceCornerIndices(i,3));
                    double gabs2=presJumps(i).real()*presJumps(i).real()+presJumps(i).imag()*presJumps(i).imag();
                    
                    //derivatives of |(c1*zi+d1)|^2*|g|^2-|(c2*zi+d2)|^2
                    //c1.real() derivative
                    JVals(MCTriOffset+16*i)=(2*c1.real()*zi.real()*zi.real()  + 2*d1.real()*zi.real() + 2*c1.real()*zi.imag()*zi.imag() + 2*d1.imag()*zi.imag())*gabs2;
                    //c1.imag() derivative
                    JVals(MCTriOffset+16*i+1)=(2*c1.imag()*zi.real()*zi.real() + 2*d1.imag()*zi.real() + 2*c1.imag()*zi.imag()*zi.imag() - 2*d1.real()*zi.imag())*gabs2,
                    //d1.real() derivative
                    JVals(MCTriOffset+16*i+2)=(2*d1.real() + 2*c1.real()*zi.real() - 2*c1.imag()*zi.imag())*gabs2;
                    //d1.real() derivative
                    JVals(MCTriOffset+16*i+3)=(2*d1.imag() + 2*c1.real()*zi.imag() + 2*c1.imag()*zi.real())*gabs2;
                    
                    //c2.real() derivative
                    JVals(MCTriOffset+16*i+4)=-(2*c2.real()*zi.real()*zi.real() + 2*d2.real()*zi.real() + 2*c2.real()*zi.imag()*zi.imag() + 2*d2.imag()*zi.imag());
                    //c2.imag() derivative
                    JVals(MCTriOffset+16*i+5)=-(2*c2.imag()*zi.real()*zi.real()+ 2*d2.imag()*zi.real() + 2*c2.imag()*zi.imag()*zi.imag() - 2*d2.real()*zi.imag());
                    //d2.real() derivative
                    JVals(MCTriOffset+16*i+6)=-(2*d2.real() + 2*c2.real()*zi.real() - 2*c2.imag()*zi.imag());
                    //d2.real() derivative
                    JVals(MCTriOffset+16*i+7)=-(2*d2.imag() + 2*c2.real()*zi.imag() + 2*c2.imag()*zi.real());
                    
                    //derivatives of |(c1*zj+d1)|^2-|(c2*zj+d2)|^2 * |g|^2
                    //c1.real() derivative
                    JVals(MCTriOffset+16*i+8)=2*c1.real()*zj.real()*zj.real() + 2*d1.real()*zj.real() + 2*c1.real()*zj.imag()*zj.imag() + 2*d1.imag()*zj.imag();
                    //c1.imag() derivative
                    JVals(MCTriOffset+16*i+9)=2*c1.imag()*zj.real()*zj.real() + 2*d1.imag()*zj.real() + 2*c1.imag()*zj.imag()*zj.imag() - 2*d1.real()*zj.imag(),
                    //d1.real() derivative
                    JVals(MCTriOffset+16*i+10)=2*d1.real() + 2*c1.real()*zj.real() - 2*c1.imag()*zj.imag();
                    //d1.real() derivative
                    JVals(MCTriOffset+16*i+11)=2*d1.imag() + 2*c1.real()*zj.imag() + 2*c1.imag()*zj.real();
                    
                    //c2.real() derivative
                    JVals(MCTriOffset+16*i+12)=-(2*c2.real()*zj.real()*zj.real() + 2*d2.real()*zj.real() + 2*c2.real()*zj.imag()*zj.imag() + 2*d2.imag()*zj.imag())*gabs2;
                    //c2.imag() derivative
                    JVals(MCTriOffset+16*i+13)=-(2*c2.imag()*zj.real()*zj.real() + 2*d2.imag()*zj.real() + 2*c2.imag()*zj.imag()*zj.imag() - 2*d2.real()*zj.imag())*gabs2;
                    //d2.real() derivative
                    JVals(MCTriOffset+16*i+14)=-(2*d2.real() + 2*c2.real()*zj.real() - 2*c2.imag()*zj.imag())*gabs2;
                    //d2.real() derivative
                    JVals(MCTriOffset+16*i+15)=-(2*d2.imag() + 2*c2.real()*zj.imag() + 2*c2.imag()*zj.real())*gabs2;
                }
            }
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
            
            finalcd=initSolution;
            update_energy(x);
            double finalTotalError=EVec.lpNorm<Eigen::Infinity>();
            double finalConstError=constVec.lpNorm<Eigen::Infinity>();
            std::cout<<"Final Const Error:"<<finalTotalError<<std::endl;
            std::cout<<"Final Total Error:"<<finalConstError<<std::endl;
            return true;
        }
        
        Moebius2DInterpolationTraits(){}
        ~Moebius2DInterpolationTraits(){}

    };
}}



#endif
