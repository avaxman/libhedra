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
#include <iostream>

namespace hedra {
    namespace optimization
    {
        
        template<class LinearSolver, class SolverTraits>
        class GNSolver{
        public:
            Eigen::VectorXd x;      //current solution; always updated
            Eigen::VectorXd prevx;  //the solution of the previous iteration
            Eigen::VectorXd x0;     //the initial solution to the system
            Eigen::VectorXd d;             //the direction taken.
            Eigen::VectorXd currEnergy;    //the current value of the energy
            Eigen::VectorXd prevEnergy;    //the previous value of the energy
            
            Eigen::VectorXi HRows, HCols;  //(row,col) pairs for H=J^T*J matrix
            Eigen::VectorXd HVals;      //values for H matrix
            Eigen::MatrixXi S2D;        //single J to J^J indices
            
            LinearSolver* LS;
            SolverTraits* ST;
            int maxIterations;
            double hTolerance;
            double xTolerance;
            double fooTolerance;
            
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
                std::vector<int> oIlist;
                std::vector<int> oJlist;
                std::vector<std::pair<int, int> > S2Dlist;
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
                                S2Dlist.push_back(std::pair<int,int>(i,j));
                            }
                        }
                    }
                    CurrTri+=NumCurrTris;
                    //std::cout<<"CurrTri, NumCurrTris: "<<CurrTri<<","<<NumCurrTris<<std::endl;
                }while (CurrTri!=iI.size());
                
                oI.resize(oIlist.size());
                oJ.resize(oJlist.size());
                S2D.resize(S2Dlist.size(),2);
                
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
                for (int i=0;i<S2D.rows();i++)
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
            
            GNSolver(){};
            
            void init(LinearSolver* _LS,
                      SolverTraits* _ST,
                      int _maxIterations=100,
                      double _xTolerance=10e-6,
                      double _hTolerance=10e-9,
                      double _fooTolerance=10e7){
                
                LS=_LS;
                ST=_ST;
                maxIterations=_maxIterations;
                hTolerance=_hTolerance;
                xTolerance=_xTolerance;
                fooTolerance=_fooTolerance;
                //analysing pattern
                MatrixPattern(ST->JRows, ST->JCols,HRows,HCols,S2D);
                HVals.resize(HRows.size());
                
                LS->analyze(HRows,HCols,true);
                
                d.resize(ST->xSize);
                x.resize(ST->xSize);
                x0.resize(ST->xSize);
                prevx.resize(ST->xSize);
                currEnergy.resize(ST->EVec.size());
                prevEnergy.resize(ST->EVec.size());
            }
            
            
            bool solve(const bool verbose) {
                
                using namespace Eigen;
                using namespace std;
                ST->initial_solution(x0);
                prevx<<x0;
                int currIter=0;
                bool stop=false;
                double currError, prevError;
                VectorXd rhs(ST->xSize);
                MatrixXd direction;
                if (verbose)
                    cout<<"******Beginning Optimization******"<<endl;
                
                do{
                    currIter=0;
                    stop=false;
                    do{
                        ST->pre_iteration(prevx);
                        ST->update_energy(prevx);
                        ST->update_jacobian(prevx);
                        if (verbose)
                            cout<<"Initial Energy for Iteration "<<currIter<<": "<<ST->EVec.template lpNorm<Infinity>()<<endl;
                        MatrixValues(HRows, HCols, ST->JVals, S2D, HVals);
                        MultiplyAdjointVector(ST->JRows, ST->JCols, ST->JVals, -ST->EVec, rhs);
                        
                        //solving to get the GN direction
                        if(!LS->factorize(HVals, true)) {
                            // decomposition failed
                            cout<<"Solver Failed to factorize! "<<endl;
                            return false;
                        }
                        
                        LS->solve(rhs,direction);
                        cout<<"direction max"<<direction.template lpNorm<Infinity>()<<endl;
                        
                        //doing a line search by decreasing by half until the energy goes down
                        //TODO: more effective line search
                        prevEnergy<<ST->EVec;
                        prevError=prevEnergy.template lpNorm<Infinity>();
                        cout<<"prevError: "<<prevError<<endl;
                        double h=1.0;
                        double t=0.0; //10e-4*direction.dot(rhs);
                        do{
                            x<<prevx+h*direction;
                            ST->update_energy(x);
                            currEnergy<<ST->EVec;
                            currError=currEnergy.template lpNorm<Infinity>();
                            //cout<<"currError: "<<currError<<endl;
                            //cout<<"h="<<h<<endl;
                            if (prevError-currError>=0.0)//h*t)
                                break;
                            
                            h*=0.5;
                            
                        }while (h>hTolerance);
                        
                        if (verbose){
                            //cout<<"currError: "<<currError<<endl;
                        }
                        
                        currIter++;
                        double xDiff=(x-prevx).template lpNorm<Infinity>();
                        double firstOrderOptimality=rhs.lpNorm<Infinity>();
                        stop=(firstOrderOptimality<fooTolerance)&&(xDiff<xTolerance);
                        if (verbose){
                            //cout<<"xDiff: "<<xDiff<<endl;
                            //cout<<"firstOrderOptimality: "<<firstOrderOptimality<<endl;
                            //cout<<"stop: "<<stop<<endl;
                        }
                        prevx<<x;
                        //The SolverTraits can order the optimization to stop by giving "true" of to continue by giving "false"
                        bool stopFromTraits=ST->post_iteration(x);
                        stop = /*stop ||*/ stopFromTraits;
                        if (stopFromTraits)
                            cout<<"ST->Post_iteration() gave a stop"<<endl;
                    }while ((currIter<=maxIterations)&&(!stop));
                }while (!ST->post_optimization(x));
                return stop;
            }
        };
        
    }
}


#endif
