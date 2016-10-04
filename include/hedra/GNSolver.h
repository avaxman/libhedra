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
        int maxIterations;
        double htolerance;
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
        
        GNSolver(){};
        
        void init(double _hTolerance=10e-5, double _xTolerance=10e-7, double _fooTolerance=10e7):hTolerance(_hTolerance), xTolerance(_xTolerance), fooTolerance(_fooTolerance) {
        
            //analysing pattern
            LS.init();
            ST.init();
            
            MatrixPattern(ST->JRows, ST->JCols,HRows,HCols,S2D);
            HVals.resize(HRows.size());
            
            LS.analyze(HRows,HCols);
            
            d.resize(ST->xSize);
            x.resize(ST->xSize);
            x0.resize(ST->xSize);
            prevx.resize(ST->solutionSize);
            currEnergy.resize(ST->EVec.size());
            prevEnergy.resize(ST->EVec.size());
        }
        
        
        bool solve(const bool verbose) {
            
            using namespace Eigen;
            using namespace std;
            prevx<<x0;
            
            int currIter=0;
            bool stop=false;
            double currMaxError, prevMaxError;
            VectorXd rhs(ST.xSize);
            
            do{
                ST.init_iteration(prevx);
                ST.update_energy_jacobian(prevx);
                MatrixValues(HRows, HCols, ST.JValues, S2D, HValues);
                MultiplyAdjointVector(ST.JRows, ST.JCols, ST.JVals, -ST.EVec, rhs);
                
                //solving to get the GN direction
                if(!LS.factorize(MatValues)) {
                    // decomposition failed
                    cout<<"Solver Failed to factorize! "<<endl;
                    return false;
                }
                
                LS.solve(rhs,direction);
                
                //doing a line search
                //doing a "lion in the desert" sampling as a heuristic to find the maximal h that still reduces the energy
                double hmin=0.0, hmax=10.0;  //TODO: arbitrary values
                prevEnergy<<ST->EVec;
                prevMaxError=prevEnergy.lpNorm<Infinity>();
                double h;
                do{
                    h=(hmin+hmax)/2.0;
                    x<<prevx+h*direction;
                    ST->update_energy_jacobian(x);
                    currEnergy<<ST->EVec;
                    currMaxError=currEnergy.lpNorm<Infinity>();
                    if (currMaxError<prevMaxError)
                        hmin=h;
                    else
                        hmax=h;
                }while ((hmax-hmin)>htolerance);
                
                currIter++;
                bool xDiff=(x-prevx).lpNorm<infinity>();
                bool firstOrderOptimality=rhs.lpNorm<infinity>();
                stop=(firstOrderOptimality<fooTolerance)&&(xDiff<xTolerance);
                prevx<<x;
                PrevSolution<<CurrSolution;
                
            }while ((CurrIter<=MaxIterations)&&(!stop));
            return stop;
        }
    };

}


#endif
