// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_GAUSS_NEWTON_SOLVER_H
#define HEDRA_GAUSS_NEWTON_SOLVER_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra::optimization
{
    
    template<class LinearSolver, class SolverTraits>
    class GNSolver{
    private:
        Eigen::MatrixXd x;      //current solution; always updated
        Eigen::MatrixXd prevx;  //the solution of the previous iteration
        Eigen::MatrixXd x0;     //the initial solution to the system
        VectorXd d;             //the direction taken.
        VectorXd currEnergy;    //the current value of the energy
        VectorXd prevEnergy;    //the previous value of the energy
        
        Eigen::VectorXd lambda;  //Lagrange multipliers
        Eigen::VectorXi Ie, Je;  //(row,col) pairs for energy matrix
        Eigen::VectorXd Se;      //values for energy matrix
        
        Eigen::VectorXi Ic, Jc;  //(row,col) pairs for constraint matrix
        Eigen::VectorXd Sc;      //values for constraint matrix
    
        LinearSolver LS;
        SolverTraits ST;
        double penalty;
        int maxIterations;
        
        VectorXi MatRows, MatCols;
        VectorXd MatValues;
        VectorXd Rhs;
        MatrixXi Single2Double;
        

    public:
        
        ALSolver(){};
        
        
        void init(){
        
            //analysing pattern
            
            MatrixPattern(st->GradRows, st->GradCols,MatRows,MatCols,Single2Double);
            MatValues.resize(MatRows.size());
            
            lLS
            
            ps.set_type(2);
            ps.set_pattern(MatRows,MatCols,VectorXd::Ones(MatRows.size()));
            ps.analyze_pattern();
            
            Direction.resize(SolverTraits->SolutionSize);
            CurrSolution.resize(SolverTraits->SolutionSize);
            PrevSolution.resize(SolverTraits->SolutionSize);
            CurrEnergy.resize(SolverTraits->TotalVec.size());
            PrevEnergy.resize(SolverTraits->TotalVec.size());
            
            int MaxCol=MatCols.maxCoeff()+1;
            Rhs.resize(MaxCol);
            
            //TestMatrixOperations();
        }
        
        
        VectorXd Solve(const VectorXd& InitSolution, int MaxIterations){
            
            vector<double> ConvErrorsList;
            
            PrevSolution<<InitSolution;
            int CurrIter=0;
            bool isValid;
            
            SolverTraits->UpdateConstraints(InitSolution);
            SolverTraits->Lambda=SolverTraits->mu*SolverTraits->ConstVec;
            SolverTraits->InitSolution=InitSolution;
            
            SolverTraits->UpdateEnergy(InitSolution);
            cout<<"Initial Squared Energy: "<<SolverTraits->TotalVec.squaredNorm()<<endl;
            //cout<<"SolverTraits->TotalVec"<<SolverTraits->TotalVec<<endl;
            /*if (SolverTraits->ConstVec.size()!=0)
             cout<<"Initial Consts: "<<SolverTraits->ConstVec.template lpNorm<Infinity>()<<endl;*/
            
            double CurrEnergyMaxError, PrevEnergyMaxError;
            
            bool EnergyStop=false;
            bool SolutionStop=false;
            bool GradStop=false;
            bool Stop=false;
            
            double DiagAddition=0; //10e-2;
            double Inith=1.0;
            
            if (GradientDescent)
                DiagAddition=1;
            
            do{
                cout<<"Energy Minimization Iteration"<<endl;
                cout<<"*****************************"<<endl;
                do{
                    SolverTraits->UpdateEnergy(PrevSolution);
                    SolverTraits->UpdateGradient(PrevSolution);
                    MatrixValues(MatRows, MatCols, SolverTraits->GradValues, Single2Double, MatValues, DiagAddition);
                    MultiplyAdjointVector(SolverTraits->GradRows, SolverTraits->GradCols, SolverTraits->GradValues, -SolverTraits->TotalVec, Rhs);
                    
                    //cout<<"GradValues: "<<SolverTraits->GradValues<<endl;
                    //cout<<"TotalVec: "<<SolverTraits->TotalVec<<endl;
                    //cout<<"Rhs: "<<Rhs<<endl;
                    
                    //solving to get the direction
                    ps.update_a(MatValues);
                    if(!ps.factorize()) {
                        // decomposition failed
                        cout<<"Solver Failed to factorize! "<<endl;
                    }
                    
                    ps.solve(Rhs,Direction);
                    
                    //doing a line search
                    double h=Inith;
                    PrevEnergy<<SolverTraits->TotalVec;
                    PrevEnergyMaxError=PrevEnergy.squaredNorm();
                    //cout<<"PrevEnergyMaxError: "<<PrevEnergyMaxError<<endl;
                    double CurrEnergyMaxError;
                    double StepSize;
                    do{
                        CurrSolution<<PrevSolution+h*Direction;
                        SolverTraits->UpdateEnergy(CurrSolution);
                        CurrEnergy<<SolverTraits->TotalVec;
                        CurrEnergyMaxError=CurrEnergy.squaredNorm();
                        //cout<<"CurrEnergySquaredError: "<<CurrEnergyMaxError<<endl;
                        h*=0.75;
                        StepSize=h*Direction.squaredNorm();
                        //cout<<"h:"<<h<<endl;
                    }while ((CurrEnergyMaxError>PrevEnergyMaxError)&&(h>10e-9));
                    
                    cout<<"Step Energy: "<<CurrEnergyMaxError<<endl;
                    cout<<"StepSize"<<StepSize<<endl;
                    /*if ((CurrEnergyMaxError>=PrevEnergyMaxError)&&(DiagAddition<10e5)){
                     cout<<"DiagAddition: "<<DiagAddition<<endl;
                     DiagAddition*=10.0;
                     continue;
                     }
                     
                     if (!GradientDescent){
                     DiagAddition/=10.0;
                     cout<<"DiagAddition: "<<DiagAddition<<endl;
                     }*/
                    
                    CurrIter++;
                    ConvErrorsList.push_back(CurrEnergyMaxError);
                    
                    EnergyStop=(PrevEnergyMaxError-CurrEnergyMaxError) < SqrtEpsilon*(1+CurrEnergyMaxError);
                    cout<<"Energy Difference Norm: "<<(PrevEnergyMaxError-CurrEnergyMaxError)<<endl;
                    SolutionStop=(PrevSolution-CurrSolution).squaredNorm() < SqrtEpsilon*(1+CurrEnergyMaxError);
                    cout<<"Solution difference norm: "<<(PrevSolution-CurrSolution).squaredNorm()<<endl;
                    //cout<<"Rhs Error: "<<Rhs.norm()<<endl;
                    GradStop=Rhs.squaredNorm() < Sqrt3Epsilon*(1+CurrEnergyMaxError);
                    
                    cout<<"EnergyStop, SolutionStop, GradStop: "<<EnergyStop<<","<<SolutionStop<<","<<GradStop<<endl;
                    
                    Stop=(EnergyStop);// && SolutionStop);
                    
                    PrevSolution<<CurrSolution;
                    SolverTraits->InitSolution=PrevSolution;
                    
                }while ((CurrIter<MaxIterations)&&(!Stop));
                
                cout<<"End of Inner Iteration"<<endl;
                cout<<"**********************"<<endl;
                //cout<<"Validity: "<<(SolverTraits->ConstVec).template lpNorm<Infinity>()<<endl;
                //at this point, a local minimum was reached. If the constraints are well, it is time to terminate
                if (SolverTraits->isValid())
                    break;
                
                //Otherwise, reformulating and solving anew from that same point
                if (SolverTraits->ConstVec.size()!=0)
                    cout<<"Validity of constraints: "<<(SolverTraits->ConstVec).template lpNorm<Infinity>()<<endl;
                SolverTraits->UpdateEnergy(CurrSolution);
                SolverTraits->UpdateGradient(CurrSolution);
                SolverTraits->Reformulate(CurrIter, MaxIterations, CurrSolution);
            }while((CurrIter<MaxIterations));
            
            SolverTraits->UpdateEnergy(CurrSolution);
            double FinalTotalError=(SolverTraits->TotalVec).template lpNorm<Infinity>();
            if (SolverTraits->ConstVec.size()!=0){
                double FinalConstError=(SolverTraits->ConstVec).template lpNorm<Infinity>();
                cout<<"Final Const Error:"<<FinalConstError<<endl;
            }
            cout<<"Final Total Error:"<<FinalTotalError<<endl;
            cout<<"Number of Iterations and max iterations:"<<CurrIter<<","<<MaxIterations<<endl;
            ConvErrors.resize(ConvErrorsList.size());
            for (int i=0;i<ConvErrorsList.size();i++)
                ConvErrors(i)=ConvErrorsList[i];
            
            return CurrSolution;
            
        }

        
        
    
    };

}


#endif


