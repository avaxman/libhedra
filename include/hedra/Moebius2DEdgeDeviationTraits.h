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
  
  
  //This traits class implements the "constraintTraits" concept, for the energy and constraints of the complex system of deformation in [Vaxman et. al 2015] (see Tutorial).

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
    Eigen::VectorXcd finalPositions;  //final vertex positions
    Eigen::VectorXcd finalY;  //final vertex reciprocals
    Eigen::VectorXcd finalE;  //final edge deivations
    
    double smoothFactor;
    //double closeFactor;
    double constTolerance;
    double rigidityFactor;
    double posFactor;
    double prevError;
    
    bool isExactDC;
    bool isExactIAP;
    
    //Complex variables
    Eigen::VectorXcd currSolution;
    Eigen::VectorXcd currY;     //vertex reciprocals
    Eigen::VectorXcd currPositions;
    Eigen::VectorXcd currE;    //edge deviations
    Eigen::VectorXcd currEdges;
    Eigen::VectorXcd AMAPVec;
    Eigen::VectorXcd mobVec;
    Eigen::VectorXcd posVec;
    Eigen::VectorXcd rigidVec;
    //Eigen::VectorXcd closeVec;
    Eigen::VectorXcd deviationVec;
    Eigen::VectorXcd origFCR;
    
    //real variables
    Eigen::VectorXd DCVec;
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
    //int closeTriOffset, closeRowOffset;
    int deviationTriOffset, deviationRowOffset;
    int complexRowOffset, complexTriOffset;
    
    //into the real values
    int DCTriOffset, DCRowOffset;
    int IAPTriOffset, IAPRowOffset;
    
    Eigen::VectorXcd constVec;
    
    void init(const Eigen::VectorXcd& _origVc,
              const Eigen::MatrixXi& _D,
              const Eigen::MatrixXi& _F,
              const Eigen::MatrixXi& _EV,
              bool _isExactDC,
              bool _isExactIAP,
              const Eigen::VectorXi& _constIndices,
              const double _rigidityFactor){
      
      using namespace Eigen;
      using namespace std;
      
      F=_F; D=_D; EV=_EV;
      constIndices=_constIndices;
      origVc=_origVc;
      isExactDC=_isExactDC;
      isExactIAP=_isExactIAP;
      rigidityFactor=_rigidityFactor;
      
      xSize=2*(origVc.rows()+origVc.rows());
      
      if (isExactDC || isExactIAP)
        xSize+=2*EV.rows();
      
      //closeFactor=1e-6;
      constTolerance=1e-7;
      
      //creating difference operator
      d0.conservativeResize(EV.rows(), origVc.rows());
      d0.reserve(2*EV.rows());
      vector<Triplet<Complex> > d0Tris(2*EV.rows());
      for (int i=0;i<EV.rows();i++){
        d0Tris[2*i]=Triplet<Complex>(i,EV(i,0), -1.0);
        d0Tris[2*i+1]=Triplet<Complex>(i,EV(i,1), 1.0);
      }
      
      d0.setFromTriplets(d0Tris.begin(), d0Tris.end());
      
      //Allocating intermediate and output vectors
      AMAPVec.conservativeResize(EV.rows());
      rigidVec.conservativeResize(EV.rows());
      posVec.conservativeResize(constIndices.size());
      //closeVec.conservativeResize(xSize/2);
      mobVec.conservativeResize(D.sum()-3*D.rows());
      if (isExactDC || isExactIAP)
        deviationVec.conservativeResize(EV.rows());
      if (isExactDC){
        DCVec.conservativeResize(EV.rows());
        IAPVec.conservativeResize(0);
      } else if (isExactIAP){
        IAPVec.conservativeResize(EV.rows());
        DCVec.conservativeResize(0);
      } else {
        IAPVec.conservativeResize(0);
        DCVec.conservativeResize(0);
      }
      currSolution.conservativeResize(xSize/2);
      currY.conservativeResize(origVc.rows());
      currPositions.conservativeResize(origVc.rows());
      currE.conservativeResize(EV.rows());
      currEdges.conservativeResize(EV.rows());
      
      
      //original face-based cross ratios
      origFCR.conservativeResize(D.sum()-3*D.size());
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
      
      if (!isExactDC && !isExactIAP)
        EVec.conservativeResize(2*(AMAPVec.size()+rigidVec.size()+/*closeVec.size()+*/posVec.size()+mobVec.size()));
      else
        EVec.conservativeResize(2*(AMAPVec.size()+rigidVec.size()+/*closeVec.size()+*/posVec.size()+mobVec.size()+deviationVec.size())+DCVec.size()+IAPVec.size());
      
      
      if (!isExactDC && !isExactIAP)
        constVec.conservativeResize(posVec.size()+mobVec.size());
      else
        constVec.conservativeResize(posVec.size()+mobVec.size()+deviationVec.size()+DCVec.size()+IAPVec.size());
      
      
      /**************************************************Creating Gradient Pattern*******************************************/
      
      if (!isExactDC && !isExactIAP){
        complexJRows.conservativeResize(4*EV.rows()+2*EV.rows()/*+xSize/2*/+constIndices.size()+4*origFCR.size());
        complexJCols.conservativeResize(complexJRows.size());
        complexJVals.conservativeResize(complexJRows.size());
        
        JRows.conservativeResize(4*complexJRows.size());
        JCols.conservativeResize(4*complexJRows.size());
        JVals.conservativeResize(4*complexJRows.size());
      } else {
        complexJRows.conservativeResize(4*EV.rows()+2*EV.rows()/*+xSize/2*/+constIndices.size()+4*origFCR.size()+5*EV.rows());
        complexJCols.conservativeResize(complexJRows.size());
        complexJVals.conservativeResize(complexJRows.size());
        
        if (isExactDC){
          JRows.conservativeResize(4*complexJRows.size()+2*EV.rows());
          JCols.conservativeResize(4*complexJRows.size()+2*EV.rows());
          JVals.conservativeResize(4*complexJRows.size()+2*EV.rows());
        } else{
          JRows.conservativeResize(4*complexJRows.size()+EV.rows());
          JCols.conservativeResize(4*complexJRows.size()+EV.rows());
          JVals.conservativeResize(4*complexJRows.size()+EV.rows());
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
      /*closeTriOffset=rigidTriOffset+2*EV.rows();
      closeRowOffset=rigidRowOffset+EV.rows();
      for (int i=0;i<xSize/2;i++){
        complexJRows(closeTriOffset+i)=closeRowOffset+i;
        complexJCols(closeTriOffset+i)=i;
        complexJVals(closeTriOffset+i)=closeFactor;
      }*/
      
      /******************Positional Soft Constraints************/
      
      posTriOffset=rigidTriOffset+2*EV.rows();//closeTriOffset+xSize/2;
      posRowOffset=rigidRowOffset+EV.rows();//closeRowOffset+xSize/2;
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
      
      if (isExactDC || isExactIAP){  //deviation constraints are included
        
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
      
      //creating the real-valued pattern and adding DC\IAP constraints
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
      if (isExactDC){
        DCRowOffset=complexRowOffset*2;
        DCTriOffset=4*complexJRows.size();
        
        
        //actual values are updated in the gradient function
        for (int i=0;i<EV.rows();i++){
          JRows(DCTriOffset+2*i)=DCRowOffset+i;
          JCols(DCTriOffset+2*i)=2*origVc.rows()+i;  //real part of e_ij
          
          JRows(DCTriOffset+2*i+1)=DCRowOffset+i;
          JCols(DCTriOffset+2*i+1)=xSize/2+2*origVc.rows()+i;  //imaginary part of e_ij
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
      finalPositions.conservativeResize(origVc.rows());
      finalY.conservativeResize(origVc.rows());
      finalE.conservativeResize(EV.rows());
      //if (constIndices.size()==0){
      prevError=0.0;
      //} else {
      //  if (isExactDC || isExactIAP)
         // initSolution.tail(EV.rows())=finalE;
        
        //update_constraints(initSolution);
        //prevError=constVec.lpNorm<Infinity>();
     // }
      
    }
    
    void update_energy(const Eigen::VectorXd& x){
      
      using namespace Eigen;
      using namespace std;
      
      currSolution.array().real()<<x.head(x.size()/2);
      currSolution.array().imag()<<x.tail(x.size()/2);
      
      currY<<currSolution.head(origVc.rows());
      currPositions<<currSolution.segment(origVc.rows(),origVc.rows());
      
      for (int i=0;i<EV.rows();i++)
        AMAPVec(i)=(currY(EV(i,0))*(origVc(EV(i,1))-origVc(EV(i,0)))*currY(EV(i,1))-(currPositions(EV(i,1))-currPositions(EV(i,0))));
      
      rigidVec<<d0*currY;
      //closeVec<<(currSolution-initSolution);
      
      update_constraints(x);
      
      if (mobVec.size()!=0){
        if (!isExactDC && !isExactIAP)
          EVec<<smoothFactor*AMAPVec.real(), smoothFactor*rigidityFactor*rigidVec.real(), /*closeFactor*closeVec.real(),*/ posFactor*posVec.real(), mobVec.real(),
          smoothFactor*AMAPVec.imag(), smoothFactor*rigidityFactor*rigidVec.imag(), /*closeFactor*closeVec.imag(), */posFactor*posVec.imag(), mobVec.imag();
        
        else{
          if (isExactDC){
            EVec<<smoothFactor*AMAPVec.real(), smoothFactor*rigidityFactor*rigidVec.real(), /*closeFactor*closeVec.real(), */posFactor*posVec.real(), mobVec.real(), deviationVec.real(),
            smoothFactor*AMAPVec.imag(), smoothFactor*rigidityFactor*rigidVec.imag(), /*closeFactor*closeVec.imag(),*/ posFactor*posVec.imag(), mobVec.imag(), deviationVec.imag(),DCVec;
            
          }else if (isExactIAP){
            EVec<<smoothFactor*AMAPVec.real(), smoothFactor*rigidityFactor*rigidVec.real(), /*closeFactor*closeVec.real(),*/ posFactor*posVec.real(), mobVec.real(), deviationVec.real(),
            smoothFactor*AMAPVec.imag(), smoothFactor*rigidityFactor*rigidVec.imag(), /*closeFactor*closeVec.imag(),*/ posFactor*posVec.imag(), mobVec.imag(), deviationVec.imag(),IAPVec;
          }
        }
      } else {
        if (!isExactDC && !isExactIAP)
          EVec<<smoothFactor*AMAPVec.real(), smoothFactor*rigidityFactor*rigidVec.real(), /*closeFactor*closeVec.real(),*/ posFactor*posVec.real(),
          smoothFactor*AMAPVec.imag(), smoothFactor*rigidityFactor*rigidVec.imag(), /*closeFactor*closeVec.imag(), */posFactor*posVec.imag();
        
        else{
          if (isExactDC){
            EVec<<smoothFactor*AMAPVec.real(), smoothFactor*rigidityFactor*rigidVec.real(), /*closeFactor*closeVec.real(),*/ posFactor*posVec.real(), deviationVec.real(),
            smoothFactor*AMAPVec.imag(), smoothFactor*rigidityFactor*rigidVec.imag(), /*closeFactor*closeVec.imag(),*/ posFactor*posVec.imag(), deviationVec.imag(),DCVec;
            
          }else if (isExactIAP){
            EVec<<smoothFactor*AMAPVec.real(), smoothFactor*rigidityFactor*rigidVec.real(), /*closeFactor*closeVec.real(), */posFactor*posVec.real(), deviationVec.real(),
            smoothFactor*AMAPVec.imag(), smoothFactor*rigidityFactor*rigidVec.imag(), /*closeFactor*closeVec.imag(), */posFactor*posVec.imag(), deviationVec.imag(),IAPVec;
          }
        }
        
      }
      
      
    }
    
    void update_constraints(const Eigen::VectorXd& x)
    {
      currSolution.array().real()<<x.head(x.size()/2);
      currSolution.array().imag()<<x.tail(x.size()/2);
      
      currY<<currSolution.head(origVc.rows());
      currPositions<<currSolution.segment(origVc.rows(),origVc.rows());
      
      for (int i=0;i<constIndices.size();i++)
        posVec(i)=currPositions(constIndices(i))-complexConstPoses(i);
      
      int mobCounter=0;
      for (int i=0;i<D.rows();i++){
        for (int j=0;j<D(i)-3;j++){
          Complex wi=currPositions(F(i,j));
          Complex wj=currPositions(F(i,j+1));
          Complex wk=currPositions(F(i,j+2));
          Complex wl=currPositions(F(i,j+3));
          mobVec(mobCounter)=(wj-wi)*(wl-wk)-origFCR(mobCounter)*(wi-wl)*(wk-wj);
          mobCounter++;
        }
      }
      
      if (!isExactDC && !isExactIAP){
        if (mobVec.size()!=0)
          constVec<<posVec, mobVec;
        else
          constVec<<posVec;
      } else{
        //Deviation constraints
        currE<<currSolution.segment(origVc.rows()+origVc.rows(),EV.rows());
        currEdges<<d0*currPositions;
        for (int i=0;i<EV.rows();i++)
          deviationVec(i)=currY(EV(i,0))*(origVc(EV(i,1))-origVc(EV(i,0)))*currY(EV(i,1))-currE(i)*currEdges(i);
        
        if (mobVec.size()!=0){
          if (isExactDC){
            DCVec<<currE.array().abs2()-1;
            constVec<<posVec, mobVec, deviationVec, DCVec.cast<Complex>();
          }if (isExactIAP){
            IAPVec<<currE.array().imag();
            constVec<<posVec, mobVec, deviationVec, IAPVec.cast<Complex>();
          }
        } else {
          if (isExactDC){
            DCVec<<currE.array().abs2()-1;
            constVec<<posVec, deviationVec, DCVec.cast<Complex>();
          }if (isExactIAP){
            IAPVec<<currE.array().imag();
            constVec<<posVec, deviationVec, IAPVec.cast<Complex>();
          }
        }
      }
      
    }
    
    
    void update_jacobian(const Eigen::VectorXd& x){
      
      using namespace Eigen;
      using namespace std;
      
      currSolution.array().real()<<x.head(x.size()/2);
      currSolution.array().imag()<<x.tail(x.size()/2);
      
      currY<<currSolution.head(origVc.rows());
      currPositions<<currSolution.segment(origVc.rows(),origVc.rows());
      
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
        complexJVals(rigidTriOffset+2*i)=-smoothFactor*rigidityFactor;
        //Yj derivative
        complexJVals(rigidTriOffset+2*i+1)=smoothFactor*rigidityFactor;
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
          Complex wi=currPositions(F(i,j));
          Complex wj=currPositions(F(i,j+1));
          Complex wk=currPositions(F(i,j+2));
          Complex wl=currPositions(F(i,j+3));
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
      
      
      if (isExactDC || isExactIAP){
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
          complexJVals(deviationTriOffset+5*i+4)=-(currPositions(EV(i,1))-currPositions(EV(i,0)));
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
      if (isExactDC){
        for (int i=0;i<EV.rows();i++){
          JVals(DCTriOffset+2*i)=2.0*currE(i).real();
          JVals(DCTriOffset+2*i+1)=2.0*currE(i).imag();
          
        }
      }
      
      //IAP constraint jacobian are constant
      
    }
    
    void initial_solution(Eigen::VectorXd& x0){
      Eigen::VectorXcd initSolution(xSize/2);
      initSolution.head(origVc.rows())=currY;
      initSolution.segment(origVc.rows(),origVc.rows())=currPositions;
      if (isExactDC || isExactIAP)
        initSolution.tail(EV.rows())=currE;
      x0.conservativeResize(xSize);
      x0<<initSolution.real(), initSolution.imag();
      smoothFactor=100.0;
      posFactor=10.0;
    }
    
    void pre_iteration(const Eigen::VectorXd& prevx){
      currSolution.array().real()=prevx.head(prevx.size()/2);
      currSolution.array().imag()=prevx.tail(prevx.size()/2);
      update_constraints(prevx);
      prevError=constVec.lpNorm<Eigen::Infinity>();
    }
    bool post_iteration(const Eigen::VectorXd& x){
      //when error is halved, the smoothness is reduced by slowest, and when error change is zero, smoothness is halved.
      currSolution.array().real()<<x.head(x.size()/2);
      currSolution.array().imag()<<x.tail(x.size()/2);
      update_constraints(x);
      //double rate=constVec.lpNorm<Eigen::Infinity>()/prevError;
      //double reduceRate=std::min(rate/2.0,1.0);
      
      smoothFactor*=0.9;//-0.7*(1.0-reduceRate);
      std::cout<<"smoothFactor: "<<smoothFactor<<std::endl;
      return (constVec.lpNorm<Eigen::Infinity>()<constTolerance);
    }
    
    bool post_optimization(const Eigen::VectorXd& x){
      currSolution.array().real()=x.head(x.size()/2);
      currSolution.array().imag()=x.tail(x.size()/2);
      
      finalPositions=currSolution.segment(origVc.rows(),origVc.rows());
      finalY=currSolution.head(origVc.rows());
      if (isExactDC || isExactIAP)
        finalE=currSolution.tail(EV.rows());
      
      update_energy(x);
      double finalTotalError=EVec.lpNorm<Eigen::Infinity>();
      double finalConstError=constVec.lpNorm<Eigen::Infinity>();
      std::cout<<"Final Total Error:"<<finalTotalError<<std::endl;
      std::cout<<"Final Const Error:"<<finalConstError<<std::endl;
      return true;
    }
    
    Moebius2DEdgeDeviationTraits(){}
    ~Moebius2DEdgeDeviationTraits(){}
  };
  
}}



#endif
