// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_QC_ERROR_H
#define HEDRA_QC_ERROR_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <vector>
#include <cmath> 


namespace hedra
{
    // Computes the quasiconformal error of every polygonal face from source mesh Vs to combinatorially equivalent mesh Vt. This is the norm of a vector where each element is the QC error of a consecutive three vertices (which is just one element in the triangle-mesh case
    
    // Inputs:
    //  V           eigen double matrix     #V by 3 - mesh coordinates
    //  D           eigen int matrix        #F by 1 - face degrees
    //  F           eigen int matrix        #F by max(D)
    // Outputs:
    //  qcError   eigen double matix      #F by 1
    IGL_INLINE bool qc_error(const Eigen::MatrixXd& Vs,
                             const Eigen::MatrixXd& Vt,
                             const Eigen::VectorXi& D,
                             const Eigen::MatrixXi& F,
                             Eigen::VectorXd& qcError)
    {
        using namespace Eigen;
        qcError.conservativeResize(D.size());
        
        for (int i=0;i<D.rows();i++){
            Eigen::VectorXd triQCErrors(D(i));
            for (int j=0;j<D(i);j++){
                RowVector3d vs0=Vs.row(F(i,j));
                RowVector3d vs1=Vs.row(F(i,(j+1)%D(i)));
                RowVector3d vs2=Vs.row(F(i,(j+2)%D(i)));
                
                RowVector3d vt0=Vt.row(F(i,j));
                RowVector3d vt1=Vt.row(F(i,(j+1)%D(i)));
                RowVector3d vt2=Vt.row(F(i,(j+2)%D(i)));
                
                Vector3d vs01=vs1-vs0;
                Vector3d vs12=vs2-vs1;
                
                Vector3d vt01=vt1-vt0;
                Vector3d vt12=vt2-vt1;
                
                Vector3d ns=(vs12.cross(vs01)).normalized();
                Vector3d nt=(vt12.cross(vt01)).normalized();
                
                Matrix3d sMat; sMat<<vs01, vs12, ns;
                Matrix3d tMat; tMat<<vt01, vt12, nt;
                
                Matrix3d A=sMat*tMat.inverse();
                
                JacobiSVD<Matrix3d> svd(A);
                Vector3d sings=(svd.singularValues());
                
                Vector2d trueSings;
                trueSings.setOnes();
                int Count=0;
                for (int i=0;i<3;i++)
                    if (abs(sings(i)-1.0)>10e-5)
                        trueSings(Count++)=abs(sings(i));
                
                //cout<<"Sings: "<<Sings<<endl;
                //cout<<"True Sings: "<<TrueSings<<endl;
                double smax=std::max(trueSings(0),trueSings(1));
                double smin=std::min(trueSings(0), trueSings(1));
                double bigK=smax/smin;
                
                triQCErrors(j)=abs(bigK-1.0);
            }
            qcError(i)=sqrt(triQCErrors.squaredNorm()/(double)D(i));
        }
        return true;
    }
}




#endif


