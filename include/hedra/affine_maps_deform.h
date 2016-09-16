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
#include <Eigen/SVD>
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
    //  EFi eigen int matrix    #E by 2 - map from edges to relative index in faces
    //  EV eigen int matrix     #E by 2 - map from edges to end vertices
    //  h eigen int vector      #constraint vertex indices (handles)
    //  bendFactor double       #the relative similarty between affine maps on adjacent faces
  
    // Output:
    // adata struct AffineData     the data necessary to solve the linear system.

    //TODO: Currently uniform weights. Make them geometric.
    IGL_INLINE void affine_maps_precompute(const Eigen::MatrixXd& V,
                                           const Eigen::VectorXi& D,
                                           const Eigen::MatrixXi& F,
                                           const Eigen::MatrixXi& EV,
                                           const Eigen::MatrixXi& EF,
                                           const Eigen::MatrixXi& EFi,
                                           const Eigen::MatrixXi& FE,
                                           const Eigen::MatrixXi& EV,
                                           const Eigen::VectorXi& h,
                                           struct AffineData& adata)
    
    {
        
        
        //creating independent sets of edges
        
        using namespace Eigen;
        using namespace std;
        //Assembling the constraint matrix C
        int CRows=0;
        int NumFullVars=V.rows()+F.rows();  //every dimension is seperable.
        int NumVars=NumFullVars-h.size();
        MatrixXd OrigNormals(D.rows(),3);
        
        for (int i=0;i<D.rows();i++){
            RowVector3d FaceNormal; FaceNormal<<0.0,0.0,0.0;
            for (int j=0;j<D(i);j++){
                RowVectorXd vn=V.row(F(i,(j+D(i)-1)%D(i)));
                RowVectorXd v0=V.row(F(i,j));
                RowVectorXd v1=V.row(F(i,(j+1)%D(i)));
                
                FaceNormal=FaceNormal+((v1-v0).cross(v0-vn)).normalized();
            }
            
            OrigNormals.row(i)=FaceNormal.normalized();
        }
        
        
        /********************Assembling the full constraint matrix********************************************/
        vector<Triplet<double> > CTriplets;
        for(int i=0;i<F.rows();i++){
            for (int j=0;j<D(i)-3;j++){  //in case of triangle, nothing happens
                int vi[4];
                for (int k=0;k<4;k++)
                    vi[k]=F(i,(j+k)%D(i));
                
                Matrix3d Coeffs1; Coeffs1<<V.row(vi[2])-V.row(vi[1]),V.row(vi[1])-V.row(vi[0]), OrigNormals.row(i);
                Matrix3d Coeffs2; Coeffs2<<V.row(vi[3])-V.row(vi[2]),V.row(vi[2])-V.row(vi[1]), OrigNormals.row(i);
                
                Coeffs1=Coeffs1.inverse();
                Coeffs2=Coeffs2.inverse();
                
                for (int k=0;k<3;k++){
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[2], Coeffs1(k,0)));
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[1], -Coeffs1(k,0)));
                    
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[1], Coeffs1(k,1)));
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[0], -Coeffs1(k,1)));
                    
                    CTriplets.push_back(Triplet<double>(CRows+k, V.rows()+i, Coeffs1(k,2)));
                    
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[3], -Coeffs2(k,0)));
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[2], Coeffs2(k,0)));
                    
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[2], -Coeffs2(k,1)));
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[1], Coeffs2(k,1)));
                    
                    CTriplets.push_back(Triplet<double>(CRows+k, V.rows()+i, -Coeffs2(k,2)));
                }
                
                CRows+=3;
            }
        }
        
        
        /**************Assembling full energy matrix*************/
        
        //prescription matrix
    
        int ARows=0;
        vector<Triplet<double> > ATriplets;
        for(int i=0;i<F.rows();i++){
            for (int j=0;j<D(i);j++){
                int vi[3];
                for (int k=0;k<3;k++)
                    vi[k]=F(i,(j+k)%D(i));
                
                Matrix3d Coeffs; Coeffs<<V.row(vi[2])-V.row(vi[1]),V.row(vi[1])-V.row(vi[0]), OrigNormals.row(i);
               
                Coeffs1=Coeffs1.inverse();
  
                for (int k=0;k<3;k++){
                    ATriplets.push_back(Triplet<double>(ARows+k, vi[2], Coeffs1(k,0)));
                    ATriplets.push_back(Triplet<double>(ARows+k, vi[1], -Coeffs1(k,0)));
                    
                    ATriplets.push_back(Triplet<double>(ARows+k, vi[1], Coeffs1(k,1)));
                    ATriplets.push_back(Triplet<double>(ARows+k, vi[0], -Coeffs1(k,1)));
                    
                    ATriplets.push_back(Triplet<double>(ARows+k, V.rows()+i, Coeffs1(k,2)));
                }
                
                ARows+=3;
            }
        }

        
        //"bending" energy to difference of adjacent matrices
    
        for(int i=0;i<EF.rows();i++)
        {
            if ((EF(i,0)==-1)||(EF(i,1)==-1))  //boundary edge
                continue;
            
            
            int vi[3], vj[3];
            for (int k=0;k<3;k++){
                vi[k]=F(EF(i,0),(EFi(i,0)+k)%D(EF(i,0)));
                vi[k]=F(EF(i,1),(EFi(i,1)+k)%D(EF(i,1)));
            }
            
            Matrix3d Coeffsi; Coeffsi<<V.row(vi[2])-V.row(vi[1]),V.row(vi[1])-V.row(vi[0]), OrigNormals.row(EF(i,0));
            Matrix3d Coeffsj; Coeffsj<<V.row(vj[2])-V.row(vj[1]),V.row(vj[1])-V.row(vj[0]), OrigNormals.row(EF(i,1));
            
            Coeffsi=Coeffsi.inverse();
            Coeffsj=Coeffsj.inverse();
            
            for (int k=0;k<3;k++){
                ATriplets.push_back(Triplet<double>(ARows+k, vi[2], Coeffsi(k,0)));
                ATriplets.push_back(Triplet<double>(ARows+k, vi[1], -Coeffsi(k,0)));
                
                ATriplets.push_back(Triplet<double>(ARows+k, vi[1], Coeffsi(k,1)));
                ATriplets.push_back(Triplet<double>(ARows+k, vi[0], -Coeffsi(k,1)));
                
                ATriplets.push_back(Triplet<double>(ARows+k, V.rows()+EF(i,0), Coeffsi(k,2)));
                
                ATriplets.push_back(Triplet<double>(ARows+k, vj[2], Coeffsj(k,0)));
                ATriplets.push_back(Triplet<double>(ARows+k, vj[1], -Coeffsj(k,0)));
                
                ATriplets.push_back(Triplet<double>(ARows+k, vj[1], Coeffsj(k,1)));
                ATriplets.push_back(Triplet<double>(ARows+k, vj[0], -Coeffsj(k,1)));
                
                ATriplets.push_back(Triplet<double>(ARows+k, V.rows()+EF(i,1), Coeffsj(k,2)));
            }
            
            ARows+=3;
        }
        
        Eigen::MatrixXd xconst=Eigen::MatrixXd::Zero(NumFullVars,3);
        
        //constructing Ax-b
        adata.AFull.resize(ARows,NumFullVars);
        adata.AFull.setFromTriplets(ATriplets.begin(), ATriplets.end());
        
        adata.CFull.resize(CRows,NumFullVars);
        adata.CFull.setFromTriplets(CTriplets.begin(), CTripletsend());
        
        adata.FSize=F.rows();
        adata.VSize=V.rows();
        
        //removing the columns of the matrices by filtering the triplets (cheaper than the slice function)
        adata.x2f=Eigen::VectorXi::Zero(V.rows()+F.rows, 1);  //vertex index to free vertex index

        int CurrIndex=0;
        for (int i=0;i<h.size();i++)
            adata.x2f(h(i))=-1;
        
        for (int i=0;i<V.rows();i++)
            if (adata.x2f(i)!=-1)
                adata.x2f(i)=CurrIndex++;
        
        for (int i=0;i<F.rows();i++)
            adata(V.rows()+i)=CurrIndex++;
        
        std::cout<<"CurrIndex: "<<CurrIndex<<std::endl;
        //std::cout<<"adata.x2f: "<<adata.x2f<<std::endl;
        
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
        
        BigMat.setFromTriplets(BigMatTris.begin(), BigMatTris.end());
        adata.solver.analyzePattern(BigMat);
        adata.solver.factorize(BigMat);
    }
    
    
    //Computing a valid transformation that is as close as possible to a prescribed one.
    //This is the core part of the deformation and the interpoaltion algorithm
    //Prerequisite: affine_maps_precompute is called, and
    //qh input values are matching those in h.
    
    //input:
    // adata struct AffineData          the data necessary to solve the linear system.
    // qh eigen double matrix           h by 3 new handle positions
    
    //output:
    // A eigen double matrix            prescribed 3*F by 3 affine maps (stacked 3x3 per face)
    // q eigen double matrix            V by 3 new vertex positions (note: include handles)
    
    
    //currently: solving only one global system (thus, initial solution is not used).
    IGL_INLINE void affine_maps_prescribe(struct AffineData& adata,
                                          const Eigen::MatrixXd& q0,
                                          const Eigen::MatrixXd& A,
                                          Eigen::MatrixXd& q)
    {
    
        
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
