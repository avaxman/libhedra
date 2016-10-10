#include <igl/viewer/Viewer.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/AugmentedLagrangianTraits.h>
#include <hedra/GNSolver.h>
#include <hedra/EigenSolverWrapper.h>
#include <hedra/OffsetMeshTraits.h>
#include <hedra/check_traits.h>



typedef hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > > LinearSolver;

typedef hedra::optimization::AugmentedLagrangianTraits<hedra::optimization::OffsetMeshTraits> SolverTraits;

hedra::optimization::OffsetMeshTraits osTraits;
SolverTraits slTraits;
LinearSolver lSolver;
hedra::optimization::GNSolver<LinearSolver,SolverTraits> gnSolver;

enum ViewingMode{ORIGINAL, OFFSET, GAUSS_MAP} ViewingMode=ORIGINAL;

Eigen::MatrixXd VOrig, VOffset, VGauss;
Eigen::MatrixXi F, T;
Eigen::VectorXi D, TF;
Eigen::MatrixXi EV, EF, FE, EFi;
Eigen::MatrixXd FEs;
Eigen::VectorXi innerEdges;
Eigen::Vector3d spans;

bool UpdateCurrentView(igl::viewer::Viewer & viewer)
{
    using namespace Eigen;
    using namespace std;
    
    Eigen::MatrixXd V;
    
    switch(ViewingMode){
        case ORIGINAL: V=VOrig; break;
        case OFFSET:    V=VOffset; break;
        case GAUSS_MAP: V=VGauss; break;
    }
    
    viewer.core.show_lines=false;
    Eigen::MatrixXd OrigEdgeColors(EV.rows(),3);
    OrigEdgeColors.col(0)=Eigen::VectorXd::Constant(EV.rows(),0.0);
    OrigEdgeColors.col(1)=Eigen::VectorXd::Constant(EV.rows(),0.0);
    OrigEdgeColors.col(2)=Eigen::VectorXd::Constant(EV.rows(),0.0);
    
    viewer.data.clear();
    viewer.core.align_camera_center(V);
    viewer.data.set_mesh(V,T);
    viewer.data.compute_normals();
    viewer.data.set_edges(V,EV,OrigEdgeColors);
    return true;
}


bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
    switch(key)
    {
        case '1': ViewingMode=ORIGINAL; break;
        case '2': ViewingMode=OFFSET; break;
        case '3': ViewingMode=GAUSS_MAP; break;
    }
    UpdateCurrentView(viewer);
    return false;
}


int main(int argc, char *argv[])
{
    
    // Load a mesh in OFF format
    using namespace std;
    using namespace Eigen;
    
    hedra::polygonal_read_OFF(DATA_PATH "/eyeye_circular.off", VOrig, D, F);
    F.array()-=1;  //this is unfortunately a 1-indexed file.
    hedra::triangulate_mesh(D, F, T, TF);
    hedra::polygonal_edge_topology(D, F, EV, FE, EF,EFi,FEs,innerEdges);
    spans=VOrig.colwise().maxCoeff()-VOrig.colwise().minCoeff();
    
    //solving for offset 1.0
    double d=0.05;
    osTraits.init(VOrig, D, F,  EV,hedra::optimization::OffsetMeshTraits::VERTEX_OFFSET, d);
    slTraits.init(&osTraits, 10);
    gnSolver.init(&lSolver, &slTraits, 150, 10e-6);
    //hedra::optimization::check_traits(slTraits, slTraits.xSize);
    //exit(0);
    gnSolver.solve(true);
    
    VOffset=osTraits.fullSolution;
    VGauss=VOffset-VOrig;
        
    igl::viewer::Viewer viewer;
    viewer.callback_key_down=&key_down;
    viewer.core.background_color<<0.75,0.75,0.75,1.0;
    UpdateCurrentView(viewer);
    viewer.launch();
    
    cout<<"press 1 for original surface"<<endl;
    cout<<"press 2 for offset surface"<<endl;
    cout<<"press 3 for Gauss map (difference between the surfaces)"<<endl;
    

    
}
