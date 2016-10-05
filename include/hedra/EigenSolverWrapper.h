// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_EIGEN_SOLVER_WRAPPER_H
#define HEDRA_EIGEN_SOLVER_WRAPPER_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra::optimization
{
    //a templated wrapper to all sparse solvers by Eigen. Not doing much and not entirely efficient since the matrix has to be initialized twice, but not too bad.
    
    //TODO: perhaps better to invalidate the analysis stage and do it all in the factorization.
    template<class EigenSparseSolver>
    class EigenSolverWrapper{
    public:
        EigenSparseSolver solver;
        Eigen::SparseMatrix<double> A;
        Eigen::VectorXi rows, cols;
        
        bool analyze(const Eigen::VectorXi& _rows,
                     const Eigen::VectorXi& _cols):rows(_rows), cols(_cols){
            A.resize(rows.maxCoeff(), cols.maxCoeff());
            std::vector<Eigen::Triplets<double> > triplets;
            for (int i=0;i<rows.size();i++)
                triplets.push_Back(Triplet<double> (rows(i), cols(i), 1.0)));  //it's just a pattern
            A.setFromTriplets(triplets);
            solver.analyzePattern(A);
        }
        
        bool factorize(const Eigen::VectorXd& values){
            std::vector<Eigen::Triplets<double> > triplets;
            for (int i=0;i<rows.size();i++)
                triplets.push_Back(Triplet<double> (rows(i), cols(i), values(i))));
            A.setFromTriplets(triplets);
            solver.factorize(A);
            
        }
        
        bool solve(const Eigen::MatrixXd& rhs,
                   Eigen::MatrixXd& x){
            
            x1 = solver.solve(rhs);
        }
    };

}


#endif
