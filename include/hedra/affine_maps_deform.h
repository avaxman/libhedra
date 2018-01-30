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
#include <Eigen/SparseQR>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <vector>
#include <igl/matlab_format.h>


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
        Eigen::SparseMatrix<double> W;          //weight matrix for energy
        Eigen::SparseQR<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> >  solver;
        AffineEnergyTypes aet;
        double sqrtBendFactor;
        Eigen::VectorXi h;  //handle indices
        Eigen::VectorXi x2f;  //variables to free variables
        Eigen::MatrixXd VOrig;
        Eigen::VectorXi D;
        Eigen::MatrixXi F;
        Eigen::MatrixXd OrigNormals;
    };
    
    
    //currently only ARAP. Also, assuming A is allocated properly
    void getIdealAffineTransformation(const struct AffineData& adata,
                                      const Eigen::MatrixXd& q,
                                      Eigen::MatrixXd& A)
    {
        
        for(int i=0;i<adata.F.rows();i++){
            Eigen::MatrixXd P(adata.D(i),3);
            Eigen::MatrixXd Q(adata.D(i),3);
            for (int j=0;j<adata.D(i);j++){
                P.row(j)=adata.VOrig.row(adata.F(i,(j+1)%adata.D(i)))-adata.VOrig.row(adata.F(i,j));
                Q.row(j)=q.row(adata.F(i,(j+1)%adata.D(i)))-q.row(adata.F(i,j));
            }
            
            Eigen::Matrix3d S=P.transpose()*Q;
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(S, Eigen::ComputeFullU | Eigen::ComputeFullV);
            A.block(3*i,0,3,3)=svd.matrixU()*svd.matrixV().transpose();
        }
    }
    
    
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
                                           const Eigen::VectorXi& h,
                                           const double bendFactor,
                                           struct AffineData& adata)
    
    {
        
        
        //creating independent sets of edges
        
        using namespace Eigen;
        using namespace std;
        //Assembling the constraint matrix C
        int CRows=0;
        int NumFullVars=V.rows()+F.rows();  //every dimension is seperable.
        int NumVars=NumFullVars-h.size();
        adata.OrigNormals.resize(D.rows(),3);
        
        for (int i=0;i<D.rows();i++){
            RowVector3d FaceNormal; FaceNormal<<0.0,0.0,0.0;
            for (int j=0;j<D(i);j++){
                RowVector3d vn=V.row(F(i,(j+D(i)-1)%D(i)));
                RowVector3d v0=V.row(F(i,j));
                RowVector3d v1=V.row(F(i,(j+1)%D(i)));
                
                FaceNormal=FaceNormal+((v1-v0).cross(vn-v0)).normalized();
            }
            
            adata.OrigNormals.row(i)=FaceNormal.normalized();
        }
        
        
        /********************Assembling the full constraint matrix********************************************/
        vector<Triplet<double> > CTriplets;
        for(int i=0;i<F.rows();i++){
            for (int j=0;j<D(i)-3;j++){  //in case of triangle, nothing happens
                int vi[4];
                for (int k=0;k<4;k++)
                    vi[k]=F(i,(j+k)%D(i));
                
                Matrix3d Coeffs1; Coeffs1<<V.row(vi[2])-V.row(vi[1]),V.row(vi[1])-V.row(vi[0]), adata.OrigNormals.row(i);
                Matrix3d Coeffs2; Coeffs2<<V.row(vi[3])-V.row(vi[2]),V.row(vi[2])-V.row(vi[1]), adata.OrigNormals.row(i);
                
                Matrix3d iCoeffs1=Coeffs1.inverse();
                Matrix3d iCoeffs2=Coeffs2.inverse();
                
                for (int k=0;k<3;k++){
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[2], iCoeffs1(k,0)));
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[1], -iCoeffs1(k,0)));
                    
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[1], iCoeffs1(k,1)));
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[0], -iCoeffs1(k,1)));
                    
                    CTriplets.push_back(Triplet<double>(CRows+k, V.rows()+i, iCoeffs1(k,2)));
                    
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[3], -iCoeffs2(k,0)));
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[2], iCoeffs2(k,0)));
                    
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[2], -iCoeffs2(k,1)));
                    CTriplets.push_back(Triplet<double>(CRows+k, vi[1], iCoeffs2(k,1)));
                    
                    CTriplets.push_back(Triplet<double>(CRows+k, V.rows()+i, -iCoeffs2(k,2)));
                }
                
                CRows+=3;
            }
        }
        
        
        /**************Assembling full energy matrix*************/
        
        //prescription matrix
    
        int ARows=0;
        vector<Triplet<double> > ATriplets;
        vector<Triplet<double> > WTriplets;
        for(int i=0;i<F.rows();i++){
            for (int j=0;j<D(i);j++){
                int vi[3];
                for (int k=0;k<3;k++)
                    vi[k]=F(i,(j+k)%D(i));
                
                Matrix3d Coeffs; Coeffs<<V.row(vi[2])-V.row(vi[1]),V.row(vi[1])-V.row(vi[0]), adata.OrigNormals.row(i);
               
                Matrix3d iCoeffs=Coeffs.inverse();
                
                //cout<<"iCoeff  "<<i<<":"<<iCoeffs<<endl;
  
                for (int k=0;k<3;k++){
                    ATriplets.push_back(Triplet<double>(ARows+k, vi[2], iCoeffs(k,0)));
                    ATriplets.push_back(Triplet<double>(ARows+k, vi[1], -iCoeffs(k,0)));
                    
                    ATriplets.push_back(Triplet<double>(ARows+k, vi[1], iCoeffs(k,1)));
                    ATriplets.push_back(Triplet<double>(ARows+k, vi[0], -iCoeffs(k,1)));
                    
                    ATriplets.push_back(Triplet<double>(ARows+k, V.rows()+i, iCoeffs(k,2)));
                    
                    WTriplets.push_back(Triplet<double>(ARows+k, ARows+k, 1.0));
                }
                
                ARows+=3;
            }
        }

        
        //"bending" energy to difference of adjacent matrices
    
        adata.sqrtBendFactor=sqrt(bendFactor);
        for(int i=0;i<EF.rows();i++)
        {
            if ((EF(i,0)==-1)||(EF(i,1)==-1))  //boundary edge
                continue;
            
            
            int vi[3], vj[3];
            for (int k=0;k<3;k++){
                vi[k]=F(EF(i,0),(EFi(i,0)+k)%D(EF(i,0)));
                vj[k]=F(EF(i,1),(EFi(i,1)+k)%D(EF(i,1)));
            }
            
            //std::cout<<"vi: "<<vi[0]<<","<<vi[1]<<","<<vi[2]<<std::endl;
            //std::cout<<"vj: "<<vj[0]<<","<<vj[1]<<","<<vj[2]<<std::endl;
            Matrix3d Coeffsi; Coeffsi<<V.row(vi[2])-V.row(vi[1]),V.row(vi[1])-V.row(vi[0]), adata.OrigNormals.row(EF(i,0));
            Matrix3d Coeffsj; Coeffsj<<V.row(vj[2])-V.row(vj[1]),V.row(vj[1])-V.row(vj[0]), adata.OrigNormals.row(EF(i,1));
            
            Matrix3d iCoeffsi=Coeffsi.inverse();
            Matrix3d iCoeffsj=Coeffsj.inverse();
            
            //cout<<"iCoeffi  "<<i<<":"<<iCoeffsi<<endl;
            //cout<<"iCoeffj  "<<i<<":"<<iCoeffsj<<endl;
            
            for (int k=0;k<3;k++){
                ATriplets.push_back(Triplet<double>(ARows+k, vi[2], iCoeffsi(k,0)));
                ATriplets.push_back(Triplet<double>(ARows+k, vi[1], -iCoeffsi(k,0)));
                
                ATriplets.push_back(Triplet<double>(ARows+k, vi[1], iCoeffsi(k,1)));
                ATriplets.push_back(Triplet<double>(ARows+k, vi[0], -iCoeffsi(k,1)));
                
                ATriplets.push_back(Triplet<double>(ARows+k, V.rows()+EF(i,0), iCoeffsi(k,2)));
                
                ATriplets.push_back(Triplet<double>(ARows+k, vj[2], -iCoeffsj(k,0)));
                ATriplets.push_back(Triplet<double>(ARows+k, vj[1], iCoeffsj(k,0)));
                
                ATriplets.push_back(Triplet<double>(ARows+k, vj[1], -iCoeffsj(k,1)));
                ATriplets.push_back(Triplet<double>(ARows+k, vj[0], iCoeffsj(k,1)));
                
                ATriplets.push_back(Triplet<double>(ARows+k, V.rows()+EF(i,1), -iCoeffsj(k,2)));
                
                WTriplets.push_back(Triplet<double>(ARows+k, ARows+k, adata.sqrtBendFactor));
            }
            
            ARows+=3;
        }
        
        //constructing Ax-b
        adata.AFull.resize(ARows,NumFullVars);
        adata.AFull.setFromTriplets(ATriplets.begin(), ATriplets.end());
        
        adata.CFull.resize(CRows,NumFullVars);
        adata.CFull.setFromTriplets(CTriplets.begin(), CTriplets.end());
        
        adata.W.resize(ARows, ARows);
        adata.W.setFromTriplets(WTriplets.begin(), WTriplets.end());
        
        adata.F=F;
        adata.VOrig=V;
        adata.D=D;
        adata.h=h;
        
        //sanity checks for the full matrices
        /*MatrixXd T=MatrixXd::Random(4,3);
        Eigen::MatrixXd B(adata.AFull.rows(),3);
        int CurrIndex=0;
        for (int i=0;i<F.rows();i++)
            for (int j=0;j<D(i);j++)
                B.block(3*(CurrIndex++),0,3,3)=T.block(0,0,3,3);
        
        Eigen::MatrixXd x(adata.AFull.cols(),3);
        
        for (int i=0;i<V.rows();i++)
            x.row(i)=V.row(i)*T.block(0,0,3,3)+T.row(3);
        
        for (int i=0;i<F.rows();i++)
            x.row(V.rows()+i)=adata.OrigNormals.row(i)*T.block(0,0,3,3);
        
        cout<<"A*x-b"<<adata.W*(adata.AFull*x-B)<<endl;
        cout<<"C*x"<<adata.CFull*x<<endl;
        exit(0);*/
        
        //removing the columns of the matrices by filtering the triplets (cheaper than the slice function)
        adata.x2f=Eigen::VectorXi::Zero(V.rows()+F.rows(), 1);  //vertex index to free vertex index

        int CurrIndex=0;
        for (int i=0;i<h.size();i++)
            adata.x2f(h(i))=-1;
        
        for (int i=0;i<V.rows();i++)
            if (adata.x2f(i)!=-1)
                adata.x2f(i)=CurrIndex++;
        
        for (int i=0;i<F.rows();i++)
            adata.x2f(V.rows()+i)=CurrIndex++;
        
        //std::cout<<"CurrIndex: "<<CurrIndex<<std::endl;
        //std::cout<<"adata.x2f: "<<adata.x2f<<std::endl;
        
        std::vector<Eigen::Triplet<double> > temp=ATriplets;
        ATriplets.clear();
        for (int i=0;i<temp.size();i++)
            if (adata.x2f(temp[i].col())!=-1)
                ATriplets.push_back(Eigen::Triplet<double>(temp[i].row(), adata.x2f(temp[i].col()), temp[i].value()));
        
        temp=CTriplets;
        CTriplets.clear();
        for (int i=0;i<temp.size();i++)
            if (adata.x2f(temp[i].col())!=-1){
                CTriplets.push_back(Eigen::Triplet<double>(temp[i].row(), adata.x2f(temp[i].col()), temp[i].value()));
                //std::cout<<temp[i].row()<<","<<temp[i].col()<<std::endl;
            }
        
        adata.A.resize(ARows,NumVars);
        adata.A.setFromTriplets(ATriplets.begin(), ATriplets.end());
        
        //std::cout<<igl::matlab_format(adata.A,"A")<<std::endl;
        
        //std::cout<<"NumVars: "<<NumVars<<std::endl;
        adata.C.resize(CRows,NumVars);
        adata.C.setFromTriplets(CTriplets.begin(), CTriplets.end());
        
        Eigen::SparseMatrix<double, Eigen::ColMajor> BigMat(CRows+NumVars, CRows+NumVars);
        std::vector<Eigen::Triplet<double> > BigMatTris;
        Eigen::SparseMatrix<double> AtWA=adata.A.transpose()*adata.W*adata.A;
        
        //block At*W*A (there is no better way to do this?)
        for (int k=0; k < AtWA.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(AtWA,k); it; ++it)
                BigMatTris.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
        
        
        //std::cout<<igl::matlab_format(AtA,"AtA")<<std::endl;
        //std::cout<<igl::matlab_format(AtA,"C")<<std::endl;
        

        
        //block C and Ct
        for (int k=0; k < adata.C.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(adata.C,k); it; ++it)
            {
                BigMatTris.push_back(Eigen::Triplet<double>(it.row()+NumVars, it.col(), it.value()));
                BigMatTris.push_back(Eigen::Triplet<double>(it.col(), it.row()+NumVars, it.value()));
            }
        
       
        BigMat.setFromTriplets(BigMatTris.begin(), BigMatTris.end());
         //std::cout<<igl::matlab_format(BigMat,"BigMat")<<std::endl;
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
    

    IGL_INLINE void affine_maps_prescribe(struct AffineData& adata,
                                          const Eigen::MatrixXd& qh,
                                          const Eigen::MatrixXd& A,
                                          Eigen::MatrixXd& q)
    {
    
        
        Eigen::MatrixXd Brhs=Eigen::MatrixXd::Zero(adata.A.rows(),3);
        int ARows=0;
        for(int i=0;i<adata.F.rows();i++){
            for (int j=0;j<adata.D(i);j++){
                Brhs.block(ARows,0,3,3)=A.block(3*i,0,3,3);
                ARows+=3;
            }
        }
        
        Eigen::MatrixXd xconst=Eigen::MatrixXd::Zero(adata.AFull.cols(),3);
        for (int i=0;i<qh.rows();i++)
            xconst.row(adata.h(i))=qh.row(i);
        
        Eigen::MatrixXd toB=-adata.AFull*xconst;
        Eigen::MatrixXd toD=-adata.CFull*xconst;


        Eigen::MatrixXd B=adata.A.transpose()*adata.W*(Brhs+toB);
        
        Eigen::MatrixXd D=toD;
        
        Eigen::MatrixXd rhs(B.rows()+D.rows(),3); rhs<<B,D;
        Eigen::MatrixXd RawResult = adata.solver.solve(rhs);
        
        Eigen::MatrixXd RawFullResult(adata.VOrig.rows()+adata.F.rows(),3);
        
        for (int i=0;i<adata.x2f.size();i++)
            if (adata.x2f(i)!=-1)
                RawFullResult.row(i)=RawResult.row(adata.x2f(i));
        
        for (int i=0;i<adata.h.size();i++)
            RawFullResult.row(adata.h(i))=qh.row(i);
        
        q=RawFullResult.block(0,0,adata.VOrig.rows(),3);

    }
    
    
    //Computing an ARAP/ASAP deformation within the affine-map space.
    
    //input:
    // adata struct AffineData          the data necessary to solve the linear system.
    // qh eigen double matrix           h by 3 new handle positions
    // numIterations                    of one local-global cycle
    
    //output:
    // q eigen double matrix            V by 3 new vertex positions (note: include handles)
    IGL_INLINE void affine_maps_deform(struct AffineData& adata,
                                       const Eigen::MatrixXd& qh,
                                       const int numIterations,
                                       Eigen::MatrixXd& q)
    {
        q.conservativeResize(adata.VOrig.rows(), adata.VOrig.cols());
        Eigen::MatrixXd A(3*adata.F.rows(),3);
        for (int i=0;i<numIterations;i++){
            getIdealAffineTransformation(adata, q, A);
            affine_maps_prescribe(adata,qh,A, q);
        }
    }
}


#endif
