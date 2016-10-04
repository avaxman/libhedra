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
#include <hedra/dihedral_angles>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra::optimization
{
    
    //this class is a traits class for optimization of discrete shells deformation by given positional constraints. It is an implementation on [Froehlich and Botsch 2012] for general polyhedral meshes, when taking into account all diagonals inside a face
    
    //prerequisite before operating the optimizer with a pointer to an instance of this calss: filling out V, F, h, b
    
    //the solution vector is assumed to be arranged as xyzxyzxyz... where each triplet is a coordinate of the free vertices.
    
    class DiscreteShellsTraits{
        
        Eigen::VectorXi JRows, JCols;  //rows and column indices for the jacobian matrix
        Eigen::VectorXd JVals;         //values for the jacobian matrix.
        double xSize;                  //size of the solution
        Eigen::VectorXd EVec;          //energy vector
        
        Eigen::VectorXd V, D, F;        //libhedra mesh representation
        Eigen::MatrixXi edgeIndices;  //all pairs of edges for edge preservation
        Eigen::VectorXd EF, EFi, innerEdges;
        Eigen::VectorXi h;
        Eigen::VectorXi a2x;  //from the entire set of variables "a" to the free variables in the optimization "x".
        
        Eigen::VectorXd origLengths;
        
        void init(const Eigen::VectorXd& x0){
            
            using namespace std;
            using namespace Eigen;
            
            set<pair<int, int> > edgeIndicesList;
            
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
            for (set<pair<int, int> >::iterator si=edgeIndicesList.begin(); si!=edgeIndicesList.end();si++)
                edgeIndices.row(i)<<si->first, si->second;
            
            EVec.resize(edgeIndices.rows());
            origLengths.resize(edgeIndices.rows());
            
            for (int i=0;i<edgeIndices.rows();i++)
                origLengths(i)=(x.row(edgeIndices(i,0))-x.row(edgeIndices(i,1))).norm();
            
            int CurrIndex=0;
            for (int i=0;i<h.size();i++)
                a2x(h(i))=-1;
            
            for (int i=0;i<V.rows();i++)
                if (a2x(i)!=-1)
                    a2x(i)=CurrIndex++;
            
            for (int i=0;i<F.rows();i++)
                a2x(V.rows()+i)=CurrIndex++;

        }
        //updating the energy vector and the jacobian values for a given current solution
        
        void update_energy_jacobian(const Eigen::VectorXd& x){
            
            using namespace std;
            using namespace Eigen;
            
            MatrixXd fullx(xSize+h.size(),3);
            for (int i=0;i<a2x.size();i++)
                if (a2x(i)!=-1)
                    fullx.row(i)<<x0.segment(3*a2x(i),3).transpose();
            
            for (int i=0;i<h.size();i++)
                fullx.row(h(i))=qh.row(i);
            
            
            for (int i=0;i<edgeIndices.rows();i++)
                EVec(i)=lengthCoeff*((fullx.row(EV(i,0))-fullx.row(EV(i,1))).norm()-origLengths(i));
            
            Eigen::VectorXd currDihedralAngles;
            dihedral_angles(x, quadIndices,currDihedralAngles);
            
            EVec.segment(EV.rows(),quadIndices.rows())=bendCoeff*(currDihedralAngles-origDihedralAngles);
            
            //Jacobian
            for (int i=0;i<EV.rows();i++)
                
                
        }
        
        DiscreteShellsTraits(){}
        ~DiscreteShellsTraits(){}
    };
    
    
}


#endif
