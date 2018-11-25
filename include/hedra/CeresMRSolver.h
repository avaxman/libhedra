// This file is part of libhedra, a library for polygonal mesh processing
//
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef HEDRA_CERES_MR_SOLVER_H
#define HEDRA_CERES_MR_SOLVER_H

#include "ceres/ceres.h"
#include "glog/logging.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;


//templated versions needed for ceres

template <typename T>
inline Eigen::Matrix<T, 1, 4> QConjT(const Eigen::Matrix<T, 1, 4>& q)
{
  Eigen::Matrix<T, 1, 4> newq(q.rows(), q.cols());
  newq<<q(0), -q.tail(3);
  return newq;
}

template <typename T>
inline Eigen::Matrix<T, 1, 4> QMultT(const Eigen::Matrix<T, 1, 4>& q1, const Eigen::Matrix<T, 1, 4>& q2)
{
  Eigen::Matrix<T, 1, 4> newq;
  T r1=q1(0);
  T r2=q2(0);
  Eigen::Matrix<T, 1, 3> v1=q1.tail(3);
  Eigen::Matrix<T, 1, 3> v2=q2.tail(3);
  newq<<r1*r2-v1.dot(v2), r1*v2+r2*v1+v1.cross(v2);
  return newq;
}

template <typename T>
inline Eigen::Matrix<T, 1, 4> QInvT(const Eigen::Matrix<T, 1, 4>& q)
{
  return(QConjT(q)/q.squaredNorm());
}


struct FullCRError{
public:
  FullCRError(double* _pLength, double* _pAngle, double* _factor):pLength(_pLength), pAngle(_pAngle), factor(_factor){}
  ~FullCRError(){};
  
  template <typename T>
  bool operator()(const T* const pi,
                  const T* const pj,
                  const T* const pk,
                  const T* const pl,
                  const T* const _crvec,
                  T* residuals) const {
    
    Eigen::Matrix< T, 1, 4 > wi; wi<<T(0),pi[0],pi[1],pi[2];
    Eigen::Matrix< T, 1, 4 > wj; wj<<T(0),pj[0],pj[1],pj[2];
    Eigen::Matrix< T, 1, 4 > wk; wk<<T(0),pk[0],pk[1],pk[2];
    Eigen::Matrix< T, 1, 4 > wl; wl<<T(0),pl[0],pl[1],pl[2];
    Eigen::Matrix< T, 1, 3 > crvec; crvec<<_crvec[0],_crvec[1],_crvec[2];
    //Eigen::Matrix< T, 1, 3 > unitcrvec=crvec.normalized();
    Eigen::Matrix< T, 1, 4 > cr; cr<<(T)(*pLength)*cos(*pAngle),(T)((*pLength)*sin(*pAngle))*crvec;
    
    
    Eigen::Matrix< T, 1, 4 > CRError=(QMultT<T>(QMultT<T>(wj-wi,QInvT<T>(wk-wj)), QMultT<T>(wl-wk, QInvT<T>(wi-wl)))-cr);
    for (int i=0;i<4;i++)
      residuals[i]=(T)((*factor)*CRError(i));
    
    return true;
  }
  double* pLength;
  double* pAngle;
  double* factor;
};

struct LengthCRError{
public:
  LengthCRError(double* _pLength,  double* _factor):pLength(_pLength), factor(_factor){}
  ~LengthCRError(){};
  
  template <typename T>
  bool operator()(const T* const pi,
                  const T* const pj,
                  const T* const pk,
                  const T* const pl,
                  T* residuals) const {
    
    Eigen::Matrix< T, 1, 4 > wi; wi<<T(0),pi[0],pi[1],pi[2];
    Eigen::Matrix< T, 1, 4 > wj; wj<<T(0),pj[0],pj[1],pj[2];
    Eigen::Matrix< T, 1, 4 > wk; wk<<T(0),pk[0],pk[1],pk[2];
    Eigen::Matrix< T, 1, 4 > wl; wl<<T(0),pl[0],pl[1],pl[2];
    Eigen::Matrix< T, 1, 4 > cr; cr<<(T)(-(*pLength)), T(0.0), T(0.0), T(0.0);
    
    Eigen::Matrix< T, 1, 4 > CRError=(QMultT<T>(QMultT<T>(wj-wi,QInvT<T>(wk-wj)), QMultT<T>(wl-wk, QInvT<T>(wi-wl)))-cr);
    for (int i=0;i<4;i++)
      residuals[i]=(T)((*factor)*CRError(i));
    return true;
  }
  double* pLength;
  double* factor;
};


struct FullFNError{
public:
  FullFNError(double* _pLength, double* _pAngle, double* _factor):pLength(_pLength), pAngle(_pAngle), factor(_factor){}
  FullFNError(){};
  
  template <typename T>
  bool operator()(const T* const pi,
                  const T* const pj,
                  const T* const pk,
                  const T* const _fnvec,
                  T* residuals) const {
    Eigen::Matrix< T, 1, 4 > wi; wi<<T(0),pi[0],pi[1],pi[2];
    Eigen::Matrix< T, 1, 4 > wj; wj<<T(0),pj[0],pj[1],pj[2];
    Eigen::Matrix< T, 1, 4 > wk; wk<<T(0),pk[0],pk[1],pk[2];
    Eigen::Matrix< T, 1, 3 > fnvec; fnvec<<_fnvec[0],_fnvec[1],_fnvec[2];
    Eigen::Matrix< T, 1, 4 > fn; fn<<(T)(*pLength)*cos(*pAngle),(T)((*pLength)*sin(*pAngle))*fnvec;
    
    Eigen::Matrix< T, 1, 4 > FNError=(QMultT<T>(wj-wi,QInvT<T>(wk-wj))-fn);
    for (int i=0;i<4;i++)
      residuals[i]=(T)((*factor)*FNError(i));
    
    //exit(0);
    
    return true;
    
  }
  double* pLength;
  double* pAngle;
  double* factor;
};


class CeresMRSolver{
public:
  
  CeresMRSolver():currSolution(NULL),problem(NULL){}
  ~CeresMRSolver(){if (problem!=NULL) delete problem; if (currSolution!=NULL) delete[] currSolution;}
  
  Eigen::MatrixXi D, F;
  
  Eigen::MatrixXd QOrig;
  Eigen::MatrixXi EV;
  Eigen::MatrixXi quadVertexIndices;        //rows of wi,wj,wk,wl,CR (CR belongs to wi)
  Eigen::MatrixXi quadFaceIndices;          //rows of wi, wj, wk, wl
  Eigen::MatrixXi faceTriads;               //rows of wi, wj, wk, FN (N belongs to wj)
  
  double* currSolution;   //3*|V| (vertex positions) + 3*|V| (vertex cross-ratio vectors) + 3*|F| (face normals).
  
  //prescribed variables
  Eigen::VectorXd CRLengths;
  Eigen::VectorXd CRAngles;  //angles of cross-ratios, arranged by QuadvertexIndices
  Eigen::VectorXd FNLengths;
  Eigen::VectorXd FNAngles;    //angles of normal ratios, arranged by FaceTriads
  Eigen::VectorXd faceCRLengths;
  Eigen::VectorXd faceCRAngles;
  Eigen::VectorXi constPosIndices;
  
  double CRFactor;
  double FNFactor;
  
  ceres::Problem* problem;
  
  void init(const Eigen::MatrixXd& inQOrig,
            const Eigen::MatrixXi& inD,
            const Eigen::MatrixXi& inF,
            const Eigen::MatrixXi& inEV,
            const Eigen::MatrixXi& inQuadVertexIndices,
            const Eigen::MatrixXi& inQuadFaceIndices,
            const Eigen::MatrixXi& inFaceTriads)
  {
    
    F=inF; D=inD;
    EV=inEV;
    quadVertexIndices=inQuadVertexIndices;
    quadFaceIndices=inQuadFaceIndices;
    faceTriads=inFaceTriads;
    QOrig=inQOrig;
    
    CRLengths.conservativeResize(quadVertexIndices.rows());
    CRAngles.conservativeResize(quadVertexIndices.rows());
    faceCRLengths.conservativeResize(quadFaceIndices.rows());
    faceCRAngles.conservativeResize(quadFaceIndices.rows());
    FNLengths.conservativeResize(faceTriads.rows());
    FNAngles.conservativeResize(faceTriads.rows());
    
    if (currSolution!=NULL)
      delete[] currSolution;
    
    if (problem!=NULL)
      delete[] problem;
    
    problem=new ceres::Problem;
    
    currSolution=new double[3*QOrig.rows()+3*QOrig.rows()+3*F.rows()];
    
    for (int i=0;i<QOrig.rows();i++)
      problem->AddParameterBlock(currSolution+3*i, 3);
    
    for (int i=0;i<QOrig.rows();i++)
      problem->AddParameterBlock(currSolution+3*QOrig.rows()+3*i, 3, new ceres::HomogeneousVectorParameterization(3));
    
    for (int i=0;i<F.rows();i++)
      problem->AddParameterBlock(currSolution+3*QOrig.rows()+3*QOrig.rows()+3*i, 3,new ceres::HomogeneousVectorParameterization(3));
    
    //Vertex CR
    for (int i = 0; i <quadVertexIndices.rows(); ++i) {
      ceres::CostFunction* cost_function=new AutoDiffCostFunction<FullCRError, 4, 3, 3, 3, 3, 3>(new FullCRError(&CRLengths(i), &CRAngles(i), &CRFactor));
      problem->AddResidualBlock(cost_function,
                                NULL, // TODO: update with coefficients somehow,
                                currSolution+3*quadVertexIndices(i,0),
                                currSolution+3*quadVertexIndices(i,1),
                                currSolution+3*quadVertexIndices(i,2),
                                currSolution+3*quadVertexIndices(i,3),
                                currSolution+3*QOrig.rows()+3*quadVertexIndices(i,0));
    }
    
    //Face CR
    for (int i = 0; i<quadFaceIndices.rows(); ++i) {
      ceres::CostFunction* cost_function=new AutoDiffCostFunction<LengthCRError, 4, 3, 3, 3, 3>(new LengthCRError(&faceCRLengths(i), &CRFactor));
      problem->AddResidualBlock(cost_function,
                                NULL, // TODO: update with coefficients somehow,
                                currSolution+3*quadFaceIndices(i,0),
                                currSolution+3*quadFaceIndices(i,1),
                                currSolution+3*quadFaceIndices(i,2),
                                currSolution+3*quadFaceIndices(i,3));
    }
    
    //Face FN
    for (int i = 0; i <faceTriads.rows(); ++i) {
      ceres::CostFunction* cost_function=new AutoDiffCostFunction<FullFNError, 4, 3, 3, 3, 3>(new FullFNError(&FNLengths(i), &FNAngles(i), &FNFactor));
      problem->AddResidualBlock(cost_function,
                                NULL, // TODO: update with coefficients somehow,
                                currSolution+3*faceTriads(i,0),
                                currSolution+3*faceTriads(i,1),
                                currSolution+3*faceTriads(i,2),
                                currSolution+3*QOrig.rows()+3*QOrig.rows()+3*faceTriads(i,3));
    }
  }
  
  void set_constant_handles(const Eigen::VectorXi& _constPosIndices)
  {
    
    constPosIndices=_constPosIndices;
    //clearing constantness //TODO: check if this can be done all at once
    for (int i=0;i<QOrig.rows();i++)
      problem->SetParameterBlockVariable(currSolution+3*i);
    
    for (int i=0;i<constPosIndices.size();i++)
      problem->SetParameterBlockConstant(currSolution+3*constPosIndices(i));
    
  }
  
  //Every position is constant
  void set_constant_positions()
  {
    
    for (int i=0;i<QOrig.rows();i++)
      problem->SetParameterBlockConstant(currSolution+3*i);
    
    //set the ratios as variable
    for (int i=0;i<QOrig.rows();i++)
      problem->SetParameterBlockVariable(currSolution+3*QOrig.rows()+3*i);
    
    for (int i=0;i<F.rows();i++)
      problem->SetParameterBlockVariable(currSolution+3*QOrig.rows()+3*QOrig.rows()+3*i);
    
  }
  
  //previous solution is always the current solution
  //user is responsible to initalize both
  void solve(const double& _CRFactor, const double& _FNFactor, const bool outputProgress){
    
    CRFactor=_CRFactor;
    FNFactor=_FNFactor;
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    options.sparse_linear_algebra_library_type=ceres::SUITE_SPARSE;
    options.minimizer_progress_to_stdout = outputProgress;
    //options.preconditioner_type = ceres::JACOBI;
    //options.use_inner_iterations=true;
    //options.check_gradients=true;
    options.max_num_iterations=250;
    options.num_threads=16;
    //options.num_linear_solver_threads = 16;
    ceres::Solver::Summary summary;
    ceres::Solve(options, problem, &summary);
    if (outputProgress)
      std::cout << summary.FullReport() << "\n";
  }
};


#endif /* CeresSolver_h */
