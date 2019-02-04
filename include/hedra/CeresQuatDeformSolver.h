// This file is part of libhedra, a library for polygonal mesh processing
//
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef HEDRA_CERES_QUAT_DEFORM_SOLVER_H
#define HEDRA_CERES_QUAT_DEFORM_SOLVER_H

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



struct DCError{
public:
  DCError(double* _factor, const double& _x, const double& _y, const double& _z):factor(_factor),x(_x),y(_y),z(_z){}
  DCError(){};
  
  template <typename T>
  bool operator()(const T* const _wi,
                  const T* const _wj,
                  const T* const _Yi,
                  const T* const _Yj,
                  T* residuals) const {
    
    Eigen::Matrix< T, 1, 4 > wi; wi<<T(0),_wi[0],_wi[1],_wi[2];
    Eigen::Matrix< T, 1, 4 > wj; wj<<T(0),_wj[0],_wj[1],_wj[2];
    Eigen::Matrix< T, 1, 4 > Yi; Yi<<_Yi[0], _Yi[1], _Yi[2], _Yi[3];
    Eigen::Matrix< T, 1, 4 > Yj; Yj<<_Yj[0], _Yj[1], _Yj[2], _Yj[3];
    Eigen::Matrix< T, 1, 4 > qij; qij<<T(0),T(x),T(y),T(z);
    
    T DCRes=(wj-wi).squaredNorm() - (QMultT<T>(QMultT<T>(QConjT<T>(Yi),qij),Yj)).squaredNorm();

    residuals[0]=(T)((*factor)*DCRes);
    
    return true;
  }
  double x,y,z;
  double* factor;
};


struct RigidityError{
public:
  RigidityError(double* _factor):factor(_factor){}
  RigidityError(){};
  
  template <typename T>
  bool operator()(const T* const _Yi,
                  const T* const _Yj,
                  T* residuals) const {
    
    Eigen::Matrix< T, 1, 4 > Yi; Yi<<_Yi[0], _Yi[1], _Yi[2], _Yi[3];
    Eigen::Matrix< T, 1, 4 > Yj; Yj<<_Yj[0], _Yj[1], _Yj[2], _Yj[3];
    
    Eigen::Matrix< T, 1, 4 > RigidityRes=Yi-Yj;
    for (int i=0;i<4;i++)
      residuals[i]=(T)((*factor)*RigidityRes(i));
    
    return true;
  }
  double* factor;
};


struct AMAPError{
public:
  AMAPError(double* _factor, const double& _x, const double& _y, const double& _z):factor(_factor),x(_x),y(_y),z(_z){}
  
  template <typename T>
  bool operator()(const T* const _wi,
                  const T* const _wj,
                  const T* const _Yi,
                  const T* const _Yj,
                  T* residuals) const {
    
    Eigen::Matrix< T, 1, 4 > wi; wi<<T(0),_wi[0],_wi[1],_wi[2];
    Eigen::Matrix< T, 1, 4 > wj; wj<<T(0),_wj[0],_wj[1],_wj[2];
    Eigen::Matrix< T, 1, 4 > Yi; Yi<<_Yi[0], _Yi[1], _Yi[2], _Yi[3];
    Eigen::Matrix< T, 1, 4 > Yj; Yj<<_Yj[0], _Yj[1], _Yj[2], _Yj[3];
    Eigen::Matrix< T, 1, 4> qij; qij<<T(0),T(x),T(y),T(z);
    
    Eigen::Matrix< T, 1, 4 > AMAPRes=(wj-wi) - QMultT<T>(QMultT<T>(QConjT<T>(Yi),qij),Yj);
    for (int i=0;i<4;i++)
      residuals[i]=(T)((*factor)*AMAPRes(i));
    
    return true;
  }
  double x,y,z;
  double* factor;
};



class CeresQMDSolver{
public:
  
  CeresQMDSolver():currSolution(NULL),problem(NULL){}
  ~CeresQMDSolver(){if (problem!=NULL) delete problem; if (currSolution!=NULL) delete[] currSolution;}
  
  Eigen::MatrixXi D, F;
  Eigen::MatrixXd VOrig;
  
  double* currSolution;   //3*|V| (vertex positions) + 4*|V| (quaternionic vertex variables)
  
  Eigen::MatrixXi extEV;

  double AMAPFactor;
  double DCFactor;
  double rigidityFactor;
  
  //positional handles
  Eigen::VectorXi constIndices;
  Eigen::MatrixXd constPoses;
  
  ceres::Problem* problem;
  
  void init(const Eigen::MatrixXd& _VOrig,
            const Eigen::MatrixXi& _D,
            const Eigen::MatrixXi& _F,
            const Eigen::MatrixXi& _extEV)
  {
    
    F=_F; D=_D;
    extEV=_extEV;
    VOrig=_VOrig;
    
    if (currSolution!=NULL)
      delete[] currSolution;
    
    if (problem!=NULL)
      delete problem;
    
    problem=new ceres::Problem;
    
    currSolution=new double[3*VOrig.rows()+4*VOrig.rows()];
    
    for (int i=0;i<VOrig.rows();i++)
      problem->AddParameterBlock(currSolution+3*i, 3);
    
    for (int i=0;i<VOrig.rows();i++)
      problem->AddParameterBlock(currSolution+3*VOrig.rows()+4*i, 4);
    

    //AMAP Energy
    for (int i = 0; i <extEV.rows(); ++i) {
      Eigen::RowVector3d vij=VOrig.row(extEV(i,1))-VOrig.row(extEV(i,0));
      ceres::CostFunction* cost_function=new AutoDiffCostFunction<AMAPError, 4, 3, 3, 4,4>(new AMAPError(&AMAPFactor, vij(0), vij(1), vij(2)));
      problem->AddResidualBlock(cost_function,
                                NULL, // TODO: update with coefficients somehow,
                                currSolution+3*extEV(i,0),
                                currSolution+3*extEV(i,1),
                                currSolution+3*VOrig.rows()+4*extEV(i,0),
                                currSolution+3*VOrig.rows()+4*extEV(i,1));
    }
    
    //Rigidity Energy
    for (int i = 0; i<D.rows(); ++i) {
      for (int j = 0; j<D(i); ++j) {
        ceres::CostFunction* cost_function=new AutoDiffCostFunction<RigidityError, 2, 4,4>(new RigidityError(&rigidityFactor));
      problem->AddResidualBlock(cost_function,
                                NULL, // TODO: update with coefficients somehow,
                                currSolution+3*VOrig.rows()+4*F(i,j),
                                currSolution+3*VOrig.rows()+4*F(i,(j+1)%D(i)));
      }
    }
    
    //DC soft constraint
    for (int i = 0; i <extEV.rows(); ++i) {
      Eigen::RowVector3d vij=VOrig.row(extEV(i,1))-VOrig.row(extEV(i,0));
      ceres::CostFunction* cost_function=new AutoDiffCostFunction<DCError, 1, 3, 3, 4,4>(new DCError(&DCFactor,vij(0), vij(1),vij(2)));
      problem->AddResidualBlock(cost_function,
                                NULL, // TODO: update with coefficients somehow,
                                currSolution+3*extEV(i,0),
                                currSolution+3*extEV(i,1),
                                currSolution+3*VOrig.rows()+4*extEV(i,0),
                                currSolution+3*VOrig.rows()+4*extEV(i,1));
    }
  }
  
  void set_constant_handles(const Eigen::VectorXi& _constIndices)
  {
    
    constIndices=_constIndices;
    //clearing constantness //TODO: check if this can be done all at once
    for (int i=0;i<VOrig.rows();i++)
      problem->SetParameterBlockVariable(currSolution+3*i);
    
    for (int i=0;i<constIndices.size();i++)
      problem->SetParameterBlockConstant(currSolution+3*constIndices(i));
    
  }
  
  //Every position is constant
  void set_constant_positions()
  {
    
    for (int i=0;i<VOrig.rows();i++)
      problem->SetParameterBlockConstant(currSolution+3*i);
    
    //set the ratios as variable
    for (int i=0;i<VOrig.rows();i++)
      problem->SetParameterBlockVariable(currSolution+3*VOrig.rows()+3*i);
    
    for (int i=0;i<F.rows();i++)
      problem->SetParameterBlockVariable(currSolution+3*VOrig.rows()+3*VOrig.rows()+3*i);
    
  }
  
  //previous solution is always the current solution
  //user is responsible to initalize both
  void solve(const double& _AMAPFactor, const double& _RigidityFactor,  const double& _DCFactor, const bool outputProgress){
    
    rigidityFactor=_RigidityFactor;
    AMAPFactor=_AMAPFactor;
    DCFactor=_DCFactor;
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
