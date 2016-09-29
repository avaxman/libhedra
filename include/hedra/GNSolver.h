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
    public:
        Eigen::MatrixXd x;      //current solution; always updated
        Eigen::MatrixXd prevx;  //the solution of the previous iteration
        Eigen::MatrixXd x0;     //the initial solution to the system
        VectorXd d;             //the direction taken.
        VectorXd currEnergy;    //the current value of the energy
        VectorXd prevEnergy;    //the previous value of the energy
    
        Eigen::VectorXi HRows, HCols;  //(row,col) pairs for H=J^T*J matrix
        Eigen::VectorXd HVals;      //values for H matrix
        
        LinearSolver LS;
        SolverTraits ST;
        double penalty;
        int maxIterations;
        

        //Input: pattern of matrix M by (iI,iJ) representation
        //Output: pattern of matrix M^T*M by (oI, oJ) representation
        //        map between values in the input to values in the output (Single2Double). The map is aggregating values from future iS to oS
        //prerequisite: iI are sorted by rows (not necessary columns)
        void MatrixPattern(const Eigen::VectorXi& iI,
                           const Eigen::VectorXi& iJ,
                           Eigen::VectorXi& oI,
                           Eigen::VectorXi& oJ,
                           Eigen::MatrixXi& S2D)
        {
            int CurrTri=0;
            using namespace Eigen;
            vector<int> oIlist;
            vector<int> oJlist;
            vector<pair<int, int> > S2Dlist;
            do{
                int CurrRow=iI(CurrTri);
                int NumCurrTris=0;
                while ((CurrTri+NumCurrTris<iI.size())&&(iI(CurrTri+NumCurrTris)==CurrRow))
                    NumCurrTris++;
                
                for (int i=CurrTri;i<CurrTri+NumCurrTris;i++){
                    for (int j=CurrTri;j<CurrTri+NumCurrTris;j++){
                        if (iJ(j)>=iJ(i)){
                            oIlist.push_back(iJ(i));
                            oJlist.push_back(iJ(j));
                            S2Dlist.push_back(pair<int,int>(i,j));
                        }
                    }
                }
                CurrTri+=NumCurrTris;
            }while (CurrTri!=iI.size());
            
            oI.resize(oIlist.size());
            oJ.resize(oJlist.size());
            s2d.resize(S2Dlist.size(),2);
            
            for (int i=0;i<oIlist.size();i++){
                oI(i)=oIlist[i];
                oJ(i)=oJlist[i];
                S2D.row(i)<<S2Dlist[i].first, S2Dlist[i].second;
            }
        }
        
        //returns the values of M^T*M by multiplication and aggregating from Single2double list.
        //prerequisite - oS is allocated
        void MatrixValues(const Eigen::VectorXi& oI,
                          const Eigen::VectorXi& oJ,
                          const Eigen::VectorXd& iS,
                          const Eigen::MatrixXi& S2D,
                          Eigen::VectorXd& oS)
        {
            for (int i=0;i<s2d.rows();i++)
                oS(i)=iS(S2D(i,0))*iS(S2D(i,1));
        }
        
        //returns M^t*ivec by (I,J,S) representation
        void MultiplyAdjointVector(const Eigen::VectorXi& iI,
                                   const Eigen::VectorXi& iJ,
                                   const Eigen::VectorXd& iS,
                                   const Eigen::VectorXd& iVec,
                                   Eigen::VectorXd& oVec)
        {
            oVec.setZero();
            for (int i=0;i<iI.size();i++)
                oVec(iJ(i))+=iS(i)*iVec(iI(i));
        }
        

    public:
        
        ALSolver(){};
        
        
        void init(){
        
            //analysing pattern
            LS.init();
            ST.init();
            
            MatrixPattern(ST->JRows, ST->JCols,HRows,HCols,S2D);
            HVals.resize(HRows.size());
            
            LS.analyze(HRows,HCols);
            
            d.resize(ST->xSize);
            x.resize(ST->xSize);
            prevx.resize(ST->solutionSize);
            currEnergy.resize(ST->EVec.size());
            prevEnergy.resize(ST->EVec.size());
            
        }
        
        
        bool Solve(const bool verbose) {
            
            using namespace Eigen;
            using namespace std;
            prevx0<<x0;
            
            int CurrIter=0;
          
            ST.UpdateEnergyJacobian(x0);
            
            if (verbose)
                std::cout<<"Initial Energy: "<<ST->EVec.squaredNorm()<<endl;
            
            
            double CurrEnergyMaxError, PrevEnergyMaxError;
            VectorXd rhs(ST.xSize);
            
            do{
                ST.UpdateEnergyJacobian(prevx0);
                MatrixValues(HRows, HCols, ST.JValues, S2D, HValues);
                MultiplyAdjointVector(ST.JRows, ST.JCols, ST.JVals, -ST.EVec, rhs);
                
                //solving to get the direction
                LS.update_a(MatValues);
                if(!ps.factorize()) {
                    // decomposition failed
                    cout<<"Solver Failed to factorize! "<<endl;
                    return false;
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
            
            SolverTraits->UpdateEnergy(CurrSolution);
            double FinalTotalError=(SolverTraits->TotalVec).template lpNorm<Infinity>();
            if (SolverTraits->ConstVec.size()!=0){
                double FinalConstError=(SolverTraits->ConstVec).template lpNorm<Infinity>();
                cout<<"Final Const Error:"<<FinalConstError<<endl;
            }
            cout<<"Final Total Error:"<<FinalTotalError<<endl;
            cout<<"Number of Iterations and max iterations:"<<CurrIter<<","<<MaxIterations<<endl;
            
        }

        
        
    
    };

}


#endif