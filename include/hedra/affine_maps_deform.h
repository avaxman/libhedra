// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_AFFINE_MAPS_DEFORM_H
#define HEDRA_AFFINE_MAPS_DEFORM_H

#include <igl/igl_inline.h>
#include <igl/setdiff.h>
#include <Eigen/Core>
#include <Eigen/SparseLU>
#include <Eigen/Geometry>
#include <vector>


//This file implements the deformation algorithm described in:

//Amir Vaxman
//Modeling Polyhedral Meshes with Affine Maps
//Computer Graphics Forum (Proc. SGP) 31(5), 2012

namespace hedra
{
    // Deform a polyhedral mesh with a single affine map per face
    
    enum AffineEnergyTypes{
        ARAP,  //as-rigid-as-possible energy
        ASAP  //as-similar-as-possible energy ("conformal")
    };
    
    struct AffineData{
        Eigen::SparseMatrix<double> AFull, A;  //energy matrices (full and substituted)
        Eigen::SparseMatrix<double> CFull, C;  //constraint matrices
        Eigen::MatrixXd toB, toD;  //constant additions to right-hand sides due to handles
        Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> >  solver;
        AffineEnergyTypes aet;
        double bendFactor;
        Eigen::VectorXi h;  //handle indices
        Eigen::VectorXi x2f;  //variables to free variables
        int FSize, VSize;
    };
    
    
    //precomputation of the necessary matrices
    
    //input:
    //  V  eigen double matrix  #V by 3 original mesh coordinates
    //  D  eigen int vector     #F by 1 - face degrees
    //  F  eigen int matrix     #F by max(D) - vertex indices in face
    //  EF eigen int matrix     #E by 2 - map from edges to adjacent faces
    //  EV eigen int matrix     #E by 2 - map from edges to end vertices
    //  h eigen int vector      #constraint vertex indices (handles)
    //  bendFactor double       #the relative similarty between affine maps on adjacent faces
  
    // Output:
    // adata struct AffineData     the data necessary to solve the linear system.

    //TODO: Currently uniform weights. Make them geometric.
    IGL_INLINE void affine_maps_precompute(const Eigen::MatrixXd& V,
                                           const Eigen::VectorXi& D,
                                           const Eigen::MatrixXi& F,
                                           const Eigen::MatrixXi& EF,
                                           const Eigen::MatrixXi& FE,
                                           const Eigen::MatrixXi& EV,
                                           const Eigen::VectorXi& h,
                                           struct AffineData& adata)
    {
        
        
        //creating independent sets of edges
        
        
        //Assembling the constraint matrix C
        int CRows=0;
        int NumFullVars=3*F.rows()+V.rows();  //every dimension is seperable.
        int NumVars=NumFullVars-h.size();
        
        Eigen::MatrixXd FaceNormals(F.rows(),3);
        for (int i=0;i<F.rows();i++){
            Eigen::Matrix3d v; v<<V.row(F(i,0)), V.row(F(i,1)), V.row(F(i,2));
                                                
            FaceNormals.row(i)=(v.row(2)-v.row(1)).cross(v.row(0)-v.row(1)).normalized();
        }
            
        
        std::vector<Eigen::Triplet<double> > CTripletList;
        for(int i=0;i<EF.rows();i++)
        {
            std::cout<<"EF.row(i): "<<EF.row(i)<<std::endl;
            Eigen::RowVector3d EdgeVector=V.row(EV(i,1))-V.row(EV(i,0));
            for (int j=0;j<2;j++){
                if (EF(i,j)!=-1)
                    for (int k=0;k<3;k++){
                        if (3*EF(i,j)+k==317)
                            std::cout<<"EdgeVector(k)"<<EdgeVector(k)<<std::endl;
                        CTripletList.push_back(Eigen::Triplet<double>(CRows,3*EF(i,j)+k,EdgeVector(k)));
                        CTripletList.push_back(Eigen::Triplet<double>(CRows,3*D.size()+EV(i,0),-1.0));
                        CTripletList.push_back(Eigen::Triplet<double>(CRows,3*D.size()+EV(i,1),1.0));
                    }
                CRows++;
            }
        }
        //the phantom codntitions for the normals
        for (int i=0;i<F.rows();i++){
            for (int j=0;j<3;j++)
                CTripletList.push_back(Eigen::Triplet<double>(CRows, 3*F.rows()+j,FaceNormals(i,j)));
            CRows++;
        }
        
    
        //Assembling the energy matrix E
        int ARows=0;
        //NumVars is the same
        std::vector<Eigen::Triplet<double> > ATripletList;
        
        //prescription to a given matrix per face - just an identity
        //TODO: stack up a stock identity matrix instead of wording it here?
        for(int i=0;i<3*F.rows();i++){
            ATripletList.push_back(Eigen::Triplet<double>(i,i,1.0));
        }
        ARows+=3*F.rows();
        
        //"bending" energy to difference of adjacent matrices
        for(int i=0;i<EF.rows();i++)
        {
            if ((EF(i,0)==-1)||(EF(i,1)==-1))  //boundary edge
                continue;
            
            for (int k=0;k<3;k++){
                ATripletList.push_back(Eigen::Triplet<double>(ARows,3*EF(i,0)+k,-1.0));
                ATripletList.push_back(Eigen::Triplet<double>(ARows,3*EF(i,1)+k,1.0));
            }
            ARows++;
        }
        
        Eigen::MatrixXd xconst=Eigen::MatrixXd::Zero(NumFullVars,3);
        
       
        //constructing Ax-b
        adata.AFull.resize(ARows,NumFullVars);
        adata.AFull.setFromTriplets(ATripletList.begin(), ATripletList.end());
        adata.toB=-adata.AFull*xconst;
        
        
        adata.CFull.resize(CRows,NumFullVars);
        adata.CFull.setFromTriplets(CTripletList.begin(), CTripletList.end());
        adata.toD=-adata.CFull*xconst;
        
        
        adata.FSize=F.rows();
        adata.VSize=V.rows();
        
        //removing the columns of the matrices by filtering the triplets (cheaper than the slice function)
        adata.x2f=Eigen::VectorXi::Zero(3*F.rows()+V.rows(), 1);  //vertex index to free vertex index
        for (int i=0;i<3*F.rows();i++)
            adata.x2f(i)=i;
        int CurrIndex=3*F.rows();
        for (int i=0;i<h.size();i++)
            adata.x2f(h(i)+3*F.rows())=-1;
        
        for (int i=0;i<V.rows();i++)
            if (adata.x2f(i+3*F.rows())!=-1)
                adata.x2f(i+3*F.rows())=CurrIndex++;
        
        std::cout<<"CurrIndex: "<<CurrIndex<<std::endl;
        std::cout<<"adata.x2f: "<<adata.x2f<<std::endl;
        
        std::vector<Eigen::Triplet<double> > temp=ATripletList;
        ATripletList.clear();
        for (int i=0;i<temp.size();i++)
            if (adata.x2f(temp[i].col())!=-1)
                ATripletList.push_back(Eigen::Triplet<double>(temp[i].row(), adata.x2f(temp[i].col()), temp[i].value()));
        
        temp=CTripletList;
        CTripletList.clear();
        for (int i=0;i<temp.size();i++)
            if (adata.x2f(temp[i].col())!=-1){
                CTripletList.push_back(Eigen::Triplet<double>(temp[i].row(), adata.x2f(temp[i].col()), temp[i].value()));
                std::cout<<CTripletList[i].row()<<","<<CTripletList[i].col()<<std::endl;
            }
        
        adata.A.resize(ARows,NumVars);
        adata.A.setFromTriplets(ATripletList.begin(), ATripletList.end());
        
        std::cout<<"NumVars: "<<NumVars<<std::endl;
        adata.C.resize(CRows,NumVars);
        adata.C.setFromTriplets(CTripletList.begin(), CTripletList.end());
        

        
        //igl::min_quad_with_fixed_precompute<double>(adata.E.transpose()*adata.E, h,adata.C,true,adata.mqwfd);
        
        //constructing the Lagrangian Matrix
        
        
        Eigen::SparseMatrix<double, Eigen::ColMajor> BigMat(CRows+NumVars, CRows+NumVars);
        std::vector<Eigen::Triplet<double> > BigMatTris;
        Eigen::SparseMatrix<double> AtA=adata.A.transpose()*adata.A;
        
        //block EtE (there is no better way to do this?)
        for (int k=0; k < adata.A.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(adata.A,k); it; ++it)
                BigMatTris.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
        

        
        //block C and Ct
        for (int k=0; k < adata.C.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(adata.C,k); it; ++it)
            {
                BigMatTris.push_back(Eigen::Triplet<double>(it.row()+NumVars, it.col(), it.value()));
                BigMatTris.push_back(Eigen::Triplet<double>(it.col(), it.row()+NumVars, it.value()));
            }
        
        std::cout<<"3*F.rows()"<<3*F.rows()<<std::endl;
        std::cout<<"BigMat(:,317)"<<BigMat.col(317)<<std::endl;
        std::cout<<"C(:,317)"<<adata.C.col(317)<<std::endl;
        BigMat.setFromTriplets(BigMatTris.begin(), BigMatTris.end());
        adata.solver.analyzePattern(BigMat);
        adata.solver.factorize(BigMat);
    }
    
    
    //Computing the deformation.
    //Prerequisite: affine_maps_precompute is called, and
    //qh input values are matching those in h.
    
    //input:
    // adata struct AffineData     the data necessary to solve the linear system.
    // qh eigen double matrix           h by 3 new handle positions
    // q0 eigen double matrix           v by 3 initial solution
    
    //output:
    // A eigen double matrix            3*F by 3 affine maps (stacked 3x3 per face)
    // q eigen double matrix            V by 3 new vertex positions (note: include handles)
    
    
    //currently: solving only one global system (thus, initial solution is not used).
    IGL_INLINE void affine_maps_deform(struct AffineData& adata,
                                       const Eigen::MatrixXd& qh,
                                       const Eigen::MatrixXd& q0,
                                       Eigen::MatrixXd A,
                                       Eigen::MatrixXd q)
    {
    
        
        //currently trying to prescribe the identity transformation alone.
        Eigen::MatrixXd Brhs(adata.A.rows(),3);
        for (int i=0;i<adata.FSize;i++)
            Brhs.block(3*i,0,3,3)=Eigen::Matrix3d::Identity();
        
        Eigen::MatrixXd B=adata.A.transpose()*(Brhs+adata.toB);
        
        Eigen::MatrixXd D=adata.toD;
        
        Eigen::MatrixXd rhs(B.rows()+D.rows(),3); rhs<<B,D;
        Eigen::MatrixXd RawResult = adata.solver.solve(rhs);
        
        Eigen::MatrixXd RawFullResult(3*adata.FSize+adata.VSize,3);
        
        for (int i=0;i<adata.x2f.size();i++)
            RawFullResult.row(i)=RawResult.row(adata.x2f(i));
        
        for (int i=0;i<adata.h.size();i++)
            RawFullResult.row(adata.h(i))=qh.row(i);
            
        
        A=RawFullResult.block(0,0,3*adata.FSize,3);
        q=RawFullResult.block(3*adata.FSize,0,adata.VSize,3);
    }
}


#endif
