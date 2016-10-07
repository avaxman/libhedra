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
            
            void init(ConstraintTraits* _CT){
                
                CT=_CT:
                lambda.resize(CT->CVec.size());
                EVec.resize(CT->EVec.size()+CT->CVec.size());
                JRows.resize(CT->JERows.size()+CT->JCRows.size());
                JCols.resize(CT->JECols.size()+CT->JCCols.size());
                JVals.resize(CT->JEVals.size()+CT->JCVals.size());
                
                penalty=1.0;  //TODO: something more sophisticated
                lambda.setZero();
                
                JRows<<CT->JERows, CT->JCRows;
                JCols<<CT->JECols, CT->JCCols;
                
            }
            
            void initial_solution(Eigen::VectorXd& x0){
                CT->initial_solution(x0);
            }
            
            void pre_iteration(const Eigen::VectorXd& prevx){
                CT->pre_iteration(prevx);
                
            }
            void post_iteration(const Eigen::VectorXd& x){
                CT->post_iteration(prevx);
                
                //updating the lagrangian function
                CT->update_constraints();
                lambda=lambda-penalty*CT->CVec;
             
                //TODO: maybe update penalty
            }
            
            
            //updating the energy vector and the jacobian values for a given current solution
            void update_energy_jacobian(const Eigen::VectorXd& x){
                //TODO
                
            }
            
            AugmentedLagrangianTraits(){}
            ~AugmentedLagrangianTraits(){}
        };
        
        
    } }


#endif
