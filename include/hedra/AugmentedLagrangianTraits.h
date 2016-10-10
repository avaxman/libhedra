// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_AUGMENTED_LAGRANGIAN_TRAITS_H
#define HEDRA_AUGMENTED_LAGRANGIAN_TRAITS_H
#include <igl/igl_inline.h>
#include <igl/harmonic.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>
#include <set>


namespace hedra { namespace optimization {
        
        //this class is a traits class for optimization of discrete shells deformation by given positional constraints. It is an implementation on [Froehlich and Botsch 2012] for general polyhedral meshes, using a triangulation of them
    
        //the solution vector is assumed to be arranged as xyzxyzxyz... where each triplet is a coordinate of the free vertices.
    
    template<class ConstraintTraits>
        class AugmentedLagrangianTraits{
        public:
            
            //Requisites of the Gauss-Newton traits class
            Eigen::VectorXi JRows, JCols;  //rows and column indices for the jacobian matrix
            Eigen::VectorXd JVals;         //values for the jacobian matrix.
            int xSize;                  //size of the solution
            Eigen::VectorXd EVec;          //energy vector
            
            Eigen::VectorXd lambda;         //Lagrange multiplers
            double penalty;                 //"miu" penalty
            ConstraintTraits* CT;
            double constTolerance;          //tolerance to max(constraint)
            int maxBigIterations;            //max iterations of lambda correction + GN solve.
            int currBigIteration;            //current big iteration
            
            void init(ConstraintTraits* _CT, int _maxBigIterations=10, double constTolerance=10e-6){
                
                CT=_CT;
                xSize=CT->xSize;
                maxBigIterations=_maxBigIterations;
                lambda.resize(CT->CVec.size());
                EVec.resize(CT->EVec.size()+CT->CVec.size());
                JRows.resize(CT->JERows.size()+CT->JCRows.size());
                JCols.resize(CT->JECols.size()+CT->JCCols.size());
                JVals.resize(CT->JEVals.size()+CT->JCVals.size());
                
                penalty=1.0;  //TODO: something more sophisticated
                lambda.setZero();
                
                JRows<<CT->JERows, CT->JCRows.array()+CT->EVec.size();
                JCols<<CT->JECols, CT->JCCols;
                
                std::cout<<"max JRows: "<<JRows.maxCoeff()<<std::endl;
                std::cout<<"max JERows: "<<CT->JERows.maxCoeff()<<std::endl;
                std::cout<<"max JCRows: "<<CT->JCRows.maxCoeff()<<std::endl;
                
                std::cout<<"EVec size: "<<EVec.size()<<std::endl;
                
            }
            
            void initial_solution(Eigen::VectorXd& x0){
                CT->initial_solution(x0);
                currBigIteration=0;
            }
            
            void pre_iteration(const Eigen::VectorXd& prevx){
                CT->pre_iteration(prevx);
                
            }
            bool post_iteration(const Eigen::VectorXd& x){
                CT->post_iteration(x);
                return false;  //do not stop
            }
            
            
            //updating the energy vector and the jacobian values for a given current solution
            void update_energy(const Eigen::VectorXd& x){
                
                CT->update_energy(x);
                CT->update_constraints(x);
                EVec<<CT->EVec, sqrt(penalty/2.0)*(CT->CVec-lambda/penalty);
            }
            
            void update_jacobian(const Eigen::VectorXd& x){
                
                CT->update_jacobian(x);
                JVals<<CT->JEVals, sqrt(penalty/2.0)*CT->JCVals;
                
            }
            
            bool post_optimization(const Eigen::VectorXd& x){
                //updating the lagrangian function
                currBigIteration++;
                CT->update_constraints(x);
                lambda=lambda-penalty*CT->CVec;
                std::cout<<"Constraint Error: "<<CT->CVec.template lpNorm<Eigen::Infinity>()<<std::endl;
                
                if ((CT->CVec.template lpNorm<Eigen::Infinity>()<constTolerance)||(currBigIteration>=maxBigIterations))
                    return CT->post_optimization(x);  //Only stopping if the ConstraintTraits wants to stop
                else
                    return false;  ///do another optimization process, since we have not reached the constraints

            }
            
            AugmentedLagrangianTraits(){}
            ~AugmentedLagrangianTraits(){}
        };
        
        
    } }


#endif
