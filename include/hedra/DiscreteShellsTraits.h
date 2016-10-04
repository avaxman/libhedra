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
    
    class DiscreteShellsTraits{
        
        Eigen::VectorXi JRows, JCols;  //rows and column indices for the jacobian matrix
        Eigen::VectorXd JVals;         //values for the jacobian matrix.
        double xSize;                  //size of the solution
        Eigen::VectorXd EVec;          //energy vector
        
        Eigen::VectorXd V, D, F;        //libhedra mesh representation
        Eigen::MatrixXi quadIndices;  //all quadruplets of vertices that are measured for bending, if form (i,j,k,l) where the bending is of diagonal (i,k) between triangles (i,j,k) and (k,l,i)
        Eigen::VectorXd EV, EF, EFi, innerEdges;
        Eigen::VectorXi h;         //
        
        double lengthCoeff, bendCoeff;  //coefficients the weight bending and length energies
        Eigen::VectorXd origLengths, origDihedralAngles;  //of the original mesh
        
        void init(const Eigen::VectorXd& x0){
            
            using namespace std;
            using namespace Eigen;
            
            vector<RowVector4i> quadIndicesList;
            
            //for each non-triangular face
            RowVector4i indices;
            for (int i=0;i<D.size();i++){
                for (int j=0;j<D(i)-3;j++){
                    //both directions of diagonals
                     indices<<F(i,j), F(i,(j+1)%D(i)), F(i,(j+2)%D(i)), F(i,(j+3)%D(i));
                    quadIndicesList.push_back(indices);
                    indices<<F(i,j), F(i,(j+2)%D(i)), F(i,(j+1)%D(i)),  F(i,(j+3)%D(i));
                    quadIndicesList.push_back(indices);
                }
            }
            //for each edge
            for (int i=0;i<innerEdges.rows();i++){
                int f=EF(innerEdges(i),0);
                int g=EF(innerEdges(i),1);
                
                //from the side i->k
                int vis=EV(innerEdges(i),0);
                int vks=EV(innerEdges(i),1);
                int vjs=F(g,(EFi(innerEdges(i),1)+2)%D(g));
                int vls=F(f,(EFi(innerEdges(i),0)+D(f)-1)%D(f));
                
                //from the side k->i
                int vit=EV(innerEdges(i),1);
                int vkt=EV(innerEdges(i),0);
                int vjt=F(f,(EFi(innerEdges(i),0)+2)%D(f));
                int vlt=F(g,(EFi(innerEdges(i),1)+D(g)-1)%D(g));
                
                
                indices<<vis, vjs, vks, vls;
                quadIndicesList.push_back(indices);
                indices<<vit, vjt, vkt, vlt;
                quadIndicesList.push_back(indices);

            }
            quadIndices.resize(quadIndicesList.size(),4);
            for (int i=0;<quadIndicesList.size();i++)
                quadIndices.row(i)<<quadIndicesList[i];
            
            EVec.resize(EV.rows()+quadIndices.rows());
            
            for (int i=0;i<EV.rows();i++)
                origLengths(i)=(x.row(EV(i,0))-x.row(EV(i,1))).norm();
            
            dihedral_angles(x0, quadIndices,origDihedralAngles);

        }
        //updating the energy vector and the jacobian values for a given current solution
        update_energy_jacobian(const Eigen::VectorXd& x){
            
            //Energy
            for (int i=0;i<EV.rows();i++)
                EVec(i)=lengthCoeff*((x.row(EV(i,0))-x.row(EV(i,1))).norm()-origLengths(i));
            
            Eigen::VectorXd currDihedralAngles;
            dihedral_angles(x, quadIndices,currDihedralAngles);
            
            EVec.segment(EV.rows(),quadIndices.rows())=bendCoeff*(currDihedralAngles-origDihedralAngles);
            
            //Jacobian
        }
        
        DiscreteShellsTraits(){}
        ~DiscreteShellsTraits(){}
    };
    
    
}


#endif
