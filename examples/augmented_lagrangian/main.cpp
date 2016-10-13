#include <igl/viewer/Viewer.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/AugmentedLagrangianTraits.h>
#include <hedra/LMSolver.h>
#include <hedra/EigenSolverWrapper.h>
#include <hedra/OffsetMeshTraits.h>
#include <hedra/check_traits.h>



typedef hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > > LinearSolver;

typedef hedra::optimization::AugmentedLagrangianTraits<hedra::optimization::OffsetMeshTraits> SolverTraits;

hedra::optimization::OffsetMeshTraits osTraits;
SolverTraits slTraits;
LinearSolver lSolver;
hedra::optimization::LMSolver<LinearSolver,SolverTraits> lmSolver;

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
    
    Eigen::MatrixXd bigV;
    Eigen::MatrixXi bigT;
    Eigen::MatrixXi bigEV;
    Eigen::MatrixXd bigTC;
    
    switch(ViewingMode){
        case ORIGINAL:  {
            bigV=VOrig;
            bigT=T;
            bigEV=EV;
            bigTC.resize(T.rows(),3);
            bigTC.col(0).setConstant(255.0/255.0);
            bigTC.col(1).setConstant(215.0/255.0);
            bigTC.col(2).setConstant(0.0);
            break;
        }
        case OFFSET:{
            bigV.resize(VOrig.rows()+VOffset.rows(),3);
            bigT.resize(2*T.rows(),3);
            bigEV.resize(2*EV.rows(),2);
            bigV<<VOrig,VOffset;
            bigT<<T,T.array()+VOrig.rows();
            bigEV<<EV, EV.array()+VOrig.rows();
            bigTC.resize(2*T.rows(),3);
            bigTC.block(0,0,T.rows(),1).setConstant(255.0/255.0);
            bigTC.block(0,1,T.rows(),1).setConstant(215.0/255.0);
            bigTC.block(0,2,T.rows(),1).setConstant(0.0);
            bigTC.block(T.rows(),0,T.rows(),1).setConstant(201.0/255.0);
            bigTC.block(T.rows(),1,T.rows(),1).setConstant(192.0/255.0);
            bigTC.block(T.rows(),2,T.rows(),1).setConstant(187.0/255.0);
            break;
        }
        case GAUSS_MAP: {
            bigV=VGauss;
            bigT=T;
            bigEV=EV;
            bigTC.resize(T.rows(),3);
            bigTC.resize(T.rows(),3);
            bigTC.col(0).setConstant(255.0/255.0);
            bigTC.col(1).setConstant(215.0/255.0);
            bigTC.col(2).setConstant(0.0);
            break;
        }
    }
    
    viewer.core.show_lines=false;
    Eigen::MatrixXd OrigEdgeColors(bigEV.rows(),3);
    OrigEdgeColors.col(0)=Eigen::VectorXd::Constant(bigEV.rows(),0.0);
    OrigEdgeColors.col(1)=Eigen::VectorXd::Constant(bigEV.rows(),0.0);
    OrigEdgeColors.col(2)=Eigen::VectorXd::Constant(bigEV.rows(),0.0);
    viewer.data.clear();
    viewer.core.align_camera_center(bigV);
    viewer.data.set_mesh(bigV,bigT);
    viewer.data.compute_normals();
    viewer.data.set_colors(bigTC);
    viewer.data.set_edges(bigV,bigEV,OrigEdgeColors);
    return true;
}


bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
    using namespace std;
    switch(key)
    {
        case '1': ViewingMode=ORIGINAL; cout<<"Showing original mesh"<<endl; break;
        case '2': ViewingMode=OFFSET; cout<<"Showing vertex offset mesh"<<endl; break;
        case '3': ViewingMode=GAUSS_MAP; cout<<"Showing Gauss map"<<endl; break;
    }
    UpdateCurrentView(viewer);
    return false;
}


int main(int argc, char *argv[])
{
    
    // Load a mesh in OFF format
    using namespace std;
    using namespace Eigen;
    
    hedra::polygonal_read_OFF(DATA_PATH "/enneper2.off", VOrig, D, F);
    //F.array()-=1;
    hedra::triangulate_mesh(D, F, T, TF);
    hedra::polygonal_edge_topology(D, F, EV, FE, EF,EFi,FEs,innerEdges);
    spans=VOrig.colwise().maxCoeff()-VOrig.colwise().minCoeff();
    
    //solving for offset 1.0
    double d=100.0;
    osTraits.init(VOrig, D, F,  EV,hedra::optimization::OffsetMeshTraits::FACE_OFFSET, d);
    slTraits.init(&osTraits, 5);
    lmSolver.init(&lSolver, &slTraits, 100);
    //hedra::optimization::check_traits(slTraits, slTraits.xSize);
    //exit(0);
    lmSolver.solve(true);
    
    VOffset=osTraits.fullSolution;
    VGauss=(VOffset-VOrig).array()/d;
    cout<<"VGauss.norm().minCoeff(): "<<VGauss.rowwise().norm().minCoeff()<<endl;
    cout<<"VGauss.norm().maxCoeff(): "<<VGauss.rowwise().norm().maxCoeff()<<endl;
    
    cout<<"press 1 for original surface"<<endl;
    cout<<"press 2 for offset surface"<<endl;
    cout<<"press 3 for Gauss map (difference between the surfaces)"<<endl;
    
    igl::viewer::Viewer viewer;
    viewer.callback_key_down=&key_down;
    viewer.core.background_color<<0.75,0.75,0.75,1.0;
    UpdateCurrentView(viewer);
    viewer.launch();
}
