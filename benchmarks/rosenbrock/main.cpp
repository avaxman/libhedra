#include <hedra/LMSolver.h>
#include <hedra/EigenSolverWrapper.h>
#include <hedra/check_traits.h>
#include <iostream>
#include <Eigen/core>
#include <hedra/check_traits.h>



typedef hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > > LinearSolver;

#define VALLEY_COEFF 20.0

class RosenbrockTraits{
public:
    Eigen::VectorXi JRows, JCols;
    Eigen::VectorXd JVals;
    int xSize;
    Eigen::Vector2d EVec;
    
    void init()
    {
        xSize=2;
        
        JRows.resize(4);
        JCols.resize(4);
        JVals.resize(4);
        
        //f1(x)
        JRows(0)=0;
        JCols(0)=0;
        JRows(1)=0;
        JCols(1)=1;
        
        //f2(x)
        JRows(2)=1;
        JCols(2)=0;
        JRows(3)=1;
        JCols(3)=1;
        
    }
    
    void initial_solution(Eigen::VectorXd& x0){
        x0<<-1.2,1.0;
    }
    void pre_iteration(const Eigen::VectorXd& prevx){}
    bool post_iteration(const Eigen::VectorXd& x){return false;}
    void update_energy(const Eigen::VectorXd& x){
        EVec<<VALLEY_COEFF*(x(1)-x(0)*x(0)), 1.0-x(0);
    }
    void update_jacobian(const Eigen::VectorXd& x){
        JVals(0)=-2.0*VALLEY_COEFF*x(0);
        JVals(1)=VALLEY_COEFF;
        JVals(2)=-1.0;
        JVals(3)=0.0;
    }
    bool post_optimization(const Eigen::VectorXd& x){
        std::cout<<"x:"<<x<<std::endl;
        return true;
    }
};



RosenbrockTraits slTraits;
LinearSolver lSolver;
hedra::optimization::LMSolver<LinearSolver,RosenbrockTraits> lmSolver;


int main(int argc, char *argv[])
{
    
    // Load a mesh in OFF format
    using namespace std;
    using namespace Eigen;
    
    slTraits.init();
    lmSolver.init(&lSolver, &slTraits, 1000);
    hedra::optimization::check_traits(slTraits);
    lmSolver.solve(true);
    
    return 0;
}
