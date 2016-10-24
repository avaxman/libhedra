#include <hedra/LMSolver.h>
#include <hedra/EigenSolverWrapper.h>
#include <hedra/check_traits.h>
#include <iostream>
#include <Eigen/core>
#include <hedra/AugmentedLagrangianTraits.h>



typedef hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > > LinearSolver;

#define VALLEY_COEFF 5.0

class g11Traits{
public:
    Eigen::VectorXi JERows, JECols;
    Eigen::VectorXd JEVals;
    Eigen::VectorXi JCRows, JCCols;
    Eigen::VectorXd JCVals;
    int xSize;
    Eigen::VectorXd EVec, CVec;
    
    void init()
    {
        xSize=2;
        
        JERows.resize(4);
        JECols.resize(4);
        JEVals.resize(4);
        
        JCRows.resize(2);
        JCCols.resize(2);
        JCVals.resize(2);
        
        //f1(x)
        JERows(0)=0;
        JECols(0)=0;
        JERows(1)=0;
        JECols(1)=1;
        
        //f2(x)
        JERows(2)=1;
        JECols(2)=0;
        JERows(3)=1;
        JECols(3)=1;
        
        JCRows(0)=0;
        JCCols(0)=0;
        
        JCRows(1)=0;
        JCCols(1)=1;
        
        EVec.resize(2);
        CVec.resize(1);
        
    }
    
    void initial_solution(Eigen::VectorXd& x0){
        x0<<100.0,-30.0;
    }
    void pre_iteration(const Eigen::VectorXd& prevx){}
    bool post_iteration(const Eigen::VectorXd& x){return false;}
    void update_energy(const Eigen::VectorXd& x){
        EVec<<x(0),x(1)-1.0;
    }
    void update_constraints(const Eigen::VectorXd& x){
        CVec<<x(1)-x(0)*x(0);
    }
    void update_jacobian(const Eigen::VectorXd& x){
        
        //energy
        JEVals(0)=1.0;
        JEVals(1)=0.0;
        JEVals(2)=0.0;
        JEVals(3)=1.0;
        
        JCVals(0)=-2.0*x(0);
        JCVals(1)=1.0;
    }
    bool post_optimization(const Eigen::VectorXd& x){
        std::cout<<"x:"<<x<<std::endl;
        return true;
    }
};



g11Traits slTraits;
hedra::optimization::AugmentedLagrangianTraits<g11Traits> ctTraits;
LinearSolver lSolver;
hedra::optimization::LMSolver<LinearSolver,hedra::optimization::AugmentedLagrangianTraits<g11Traits> > lmSolver;


int main(int argc, char *argv[])
{
    
    // Load a mesh in OFF format
    using namespace std;
    using namespace Eigen;
    
    slTraits.init();
    ctTraits.init(&slTraits, 100);
    lmSolver.init(&lSolver, &ctTraits, 1000);
    hedra::optimization::check_traits(ctTraits);
    //exit(0);
    lmSolver.solve(true);
    
    return 0;
}
