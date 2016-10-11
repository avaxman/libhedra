#include <hedra/LMSolver.h>
#include <hedra/EigenSolverWrapper.h>
#include <hedra/check_traits.h>
#include <iostream>
#include <Eigen/core>
#include <hedra/check_traits.h>



typedef hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > > LinearSolver;

#define m 1000
#define n 250

class LineaFunctionTraits{
public:
    Eigen::VectorXi JRows, JCols;
    Eigen::VectorXd JVals;
    int xSize;
    Eigen::VectorXd EVec;
    
    Eigen::MatrixXd A;
    
    void init()
    {
        using namespace Eigen;
        xSize=n;
        A.resize(m,n);
        A<<Matrix<double, n, n>::Identity()-MatrixXd::Constant(n,n,2.0/(double)m),
            -MatrixXd::Constant(m-n, n,2.0/(double)m);

        JRows.resize(m*n);
        JCols.resize(m*n);
        JVals.resize(m*n);
        for (int i=0;i<m;i++){
            for (int j=0;j<n;j++){
                JRows(n*i+j)=i;
                JCols(n*i+j)=j;
                JVals(n*i+j)=A(i,j);
            }
        }
        EVec.resize(m);
    }
    
    void initial_solution(Eigen::VectorXd& x0){
        x0=Eigen::VectorXd::Constant(n, 1.0);
    }
    void pre_iteration(const Eigen::VectorXd& prevx){}
    bool post_iteration(const Eigen::VectorXd& x){return false;}
    void update_energy(const Eigen::VectorXd& x){
        EVec<<A*x-Eigen::VectorXd::Constant(m, 1.0);
    }
    void update_jacobian(const Eigen::VectorXd& x){
        //jacobian is constant
    }
    bool post_optimization(const Eigen::VectorXd& x){
        //std::cout<<"x:"<<x<<std::endl;
        return true;
    }
};



LineaFunctionTraits slTraits;
LinearSolver lSolver;
hedra::optimization::LMSolver<LinearSolver,LineaFunctionTraits> lmSolver;


int main(int argc, char *argv[])
{
    
    // Load a mesh in OFF format
    using namespace std;
    using namespace Eigen;
    
    slTraits.init();
    lmSolver.init(&lSolver, &slTraits, 100);
    hedra::optimization::check_traits(slTraits, slTraits.xSize);
    lmSolver.solve(true);
    
    return 0;
}
