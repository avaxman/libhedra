// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_DISCRETE_SHELLS_TRAITS_H
#define HEDRA_DISCRETE_SHELLS_TRAITS_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>
#include <set>


namespace hedra { namespace optimization {
        
        //this class is a traits class for optimization of discrete shells deformation by given positional constraints. It is an implementation on [Froehlich and Botsch 2012] for general polyhedral meshes, when taking into account all diagonals inside a face
        
        //prerequisite before operating the optimizer with a pointer to an instance of this calss: filling out V, F, h, b
        
        //the solution vector is assumed to be arranged as xyzxyzxyz... where each triplet is a coordinate of the free vertices.
        
        class DiscreteShellsTraits{
        public:
            
            //Requisites of the traits class
            Eigen::VectorXi JRows, JCols;  //rows and column indices for the jacobian matrix
            Eigen::VectorXd JVals;         //values for the jacobian matrix.
            int xSize;                  //size of the solution
            Eigen::VectorXd EVec;          //energy vector
            
            //These are for the the full Jacobian matrix without removing the handles
            Eigen::VectorXi fullJRows, fullJCols;
            Eigen::VectorXd fullJVals;
            Eigen::MatrixXi edgeIndices;    //all pairs of edges for edge preservation
            Eigen::VectorXi h;              //list of handles
            Eigen::MatrixXd qh;             //#h by 3 positions
            Eigen::VectorXi a2x;            //from the entire set of variables "a" to the free variables in the optimization "x".
            Eigen::VectorXi colMap;         //raw map of a2x into the columns from fullJCols into JCols
            Eigen::VectorXd origLengths;    //the original edge and diagonal lengths.
            Eigen::MatrixXd VOrig;          //original positions
            
            Eigen::VectorXd x0;                 //the initial solution to the optimization
            Eigen::MatrixXd fullSolution;       //The final solution of the last optimization
            
            void init(const Eigen::MatrixXd& V,
                      const Eigen::VectorXi& D,
                      const Eigen::MatrixXi& F,
                      const Eigen::VectorXi& _h,
                      const Eigen::MatrixXi& EF,
                      const Eigen::MatrixXi& EFi,
                      const Eigen::VectorXi& innerEdges){
                
                using namespace std;
                using namespace Eigen;
                
                std::set<pair<int, int> > edgeIndicesList;
                
                VOrig=V;
                h=_h;
                //inside each face
                for (int i=0;i<D.size();i++)
                    for (int j=0;j<D(i);j++)
                        for (int k=j+1;k<D(i);k++)
                            edgeIndicesList.insert(pair<int, int> (F(i,j)>F(i,k) ? F(i,k) : F(i,j), F(i,j)>F(i,k) ? F(i,j) : F(i,k)));
                
                //across edges
                for (int i=0;i<innerEdges.rows();i++){
                    int f=EF(innerEdges(i),0);
                    int g=EF(innerEdges(i),1);
                    
                    //from the side i->k
                    int vjs=F(g,(EFi(innerEdges(i),1)+2)%D(g));
                    int vls=F(f,(EFi(innerEdges(i),0)+D(f)-1)%D(f));
                    
                    //from the side k->i
                    int vjt=F(f,(EFi(innerEdges(i),0)+2)%D(f));
                    int vlt=F(g,(EFi(innerEdges(i),1)+D(g)-1)%D(g));
                    
                    edgeIndicesList.insert(pair<int, int> (vjs>vls ? vls : vjs, vjs>vls ? vjs : vls));
                    edgeIndicesList.insert(pair<int, int> (vjt>vlt ? vlt : vjt, vjt>vlt ? vjt : vlt));
                    
                }
                edgeIndices.resize(edgeIndicesList.size(),2);
                int currIndex=0;
                for (set<pair<int, int> >::iterator si=edgeIndicesList.begin(); si!=edgeIndicesList.end();si++)
                    edgeIndices.row(currIndex++)<<si->first, si->second;
                
                EVec.resize(edgeIndices.rows());
                origLengths.resize(edgeIndices.rows());
                
                for (int i=0;i<edgeIndices.rows();i++)
                    origLengths(i)=(VOrig.row(edgeIndices(i,1))-VOrig.row(edgeIndices(i,0))).norm();
                
                
                //creating the Jacobian pattern
                xSize=3*(VOrig.rows()-h.size());
                fullJRows.resize(6*edgeIndices.rows());
                fullJCols.resize(6*edgeIndices.rows());
                fullJVals.resize(6*edgeIndices.rows());
                
                for (int i=0;i<edgeIndices.rows();i++){
                    fullJRows.segment(6*i,6).setConstant(i);
                    fullJCols.segment(6*i,3)<<3*edgeIndices(i,0),3*edgeIndices(i,0)+1,3*edgeIndices(i,0)+2;
                    fullJCols.segment(6*i+3,3)<<3*edgeIndices(i,1),3*edgeIndices(i,1)+1,3*edgeIndices(i,1)+2;
                }
                
                //cout<<"fullJRows: "<<fullJRows<<endl;
                //cout<<"fullJCols: "<<fullJCols<<endl;
                
                a2x.resize(VOrig.rows());
                int CurrIndex=0;
                for (int i=0;i<h.size();i++)
                    a2x(h(i))=-1;
                
                for (int i=0;i<V.rows();i++)
                    if (a2x(i)!=-1)
                        a2x(i)=CurrIndex++;
                
                colMap.resize(3*VOrig.rows());
                for (int i=0;i<VOrig.rows();i++)
                    if (a2x(i)!=-1)
                        colMap.segment(3*i,3)<<3*a2x(i),3*a2x(i)+1,3*a2x(i)+2;
                    else
                        colMap.segment(3*i,3)<<-1,-1,-1;
                
                //setting up the Jacobian rows and columns
                int actualGradCounter=0;
                for (int i=0;i<fullJCols.size();i++){
                    if (colMap(fullJCols(i))!=-1)  //not a removed variable
                        actualGradCounter++;
                }
                
                JRows.resize(actualGradCounter);
                JCols.resize(actualGradCounter);
                JVals.resize(actualGradCounter);
                
                actualGradCounter=0;
                for (int i=0;i<fullJCols.size();i++){
                    if (colMap(fullJCols(i))!=-1){  //not a removed variable
                        JRows(actualGradCounter)=fullJRows(i);
                        JCols(actualGradCounter++)=colMap(fullJCols(i));
                    }
                }
                
                cout<<"JRows maxCoeff: "<<JRows.maxCoeff()<<endl;
                cout<<"JCols.maxCoeff(): "<<JCols.maxCoeff()<<endl;
                cout<<"xSize: "<<xSize<<endl;
                
                
                //the "last found" solution is just the original mesh
                x0.resize(xSize);
                for (int i=0;i<VOrig.rows();i++)
                    if (a2x(i)!=-1)
                        x0.segment(3*a2x(i),3)<<VOrig.row(i).transpose();
            }
            
            //provide the initial solution to the solver
            void initial_solution(Eigen::VectorXd& _x0){
                _x0=x0;
            }
            
            void pre_iteration(const Eigen::VectorXd& prevx){}
            void post_iteration(const Eigen::VectorXd& x){}
            
            
            //updating the energy vector and the jacobian values for a given current solution
            void update_energy_jacobian(const Eigen::VectorXd& x){
                
                using namespace std;
                using namespace Eigen;
                
                MatrixXd fullx(xSize+h.size(),3);
                for (int i=0;i<a2x.size();i++)
                    if (a2x(i)!=-1)
                        fullx.row(i)<<x.segment(3*a2x(i),3).transpose();
                
                for (int i=0;i<h.size();i++)
                    fullx.row(h(i))=qh.row(i);
                
                for (int i=0;i<edgeIndices.rows();i++)
                    EVec(i)=((fullx.row(edgeIndices(i,1))-fullx.row(edgeIndices(i,0))).norm()-origLengths(i))/(origLengths(i)*origLengths(i));
                
                //Jacobian
                for (int i=0;i<edgeIndices.rows();i++){
                    RowVector3d normedEdgeVector=(fullx.row(edgeIndices(i,1))-fullx.row(edgeIndices(i,0))).normalized();
                    fullJVals.segment(6*i,3)<<-normedEdgeVector.transpose()/(origLengths(i)*origLengths(i));
                    fullJVals.segment(6*i+3,3)<<normedEdgeVector.transpose()/(origLengths(i)*origLengths(i));
                }
                
                int actualGradCounter=0;
                for (int i=0;i<fullJCols.size();i++)
                    if (colMap(fullJCols(i))!=-1)  //not a removed variable
                        JVals(actualGradCounter++)=fullJVals(i);
            }
            
            void post_optimization(const Eigen::VectorXd& x){
                x0=x;  //the last solution will be used in consequent optimizations
                fullSolution.conservativeResize(xSize+h.size(),3);
                for (int i=0;i<a2x.size();i++)
                    if (a2x(i)!=-1)
                        fullSolution.row(i)<<x.segment(3*a2x(i),3).transpose();
                
                for (int i=0;i<h.size();i++)
                    fullSolution.row(h(i))=qh.row(i);
            }
            
            DiscreteShellsTraits(){}
            ~DiscreteShellsTraits(){}
        };
        
        
    } }


#endif
