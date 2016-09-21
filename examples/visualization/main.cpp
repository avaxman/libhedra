#include <igl/viewer/Viewer.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/hedra_read_OFF.h>
#include <hedra/hedra_edge_topology.h>
#include <hedra/point_spheres.h>

Eigen::MatrixXd V, FEs;
Eigen::MatrixXi F,EV, EF, FE,T, EFi;
Eigen::VectorXi D,TF, innerEdges;
Eigen::MatrixXd TC;
Eigen::MatrixXd sphereColors;

bool ShowPolygonalEdges=false;
bool ShowVertexSpheres=false;
bool ShowFaceNormals=false;

void UpdateCurrentView(igl::viewer::Viewer& viewer)
{
    viewer.data.clear();
    viewer.data.set_face_based(true);
    
    Eigen::MatrixXd sphereV,sphereTC;
    Eigen::MatrixXi sphereT;
    if (ShowVertexSpheres)
        hedra::point_spheres(V, 0.05, sphereColors, 10, false, sphereV, sphereT, sphereTC);
    
    Eigen::MatrixXd bigV(V.rows()+sphereV.rows(),3);
    Eigen::MatrixXi bigT(T.rows()+sphereT.rows(),3);
    Eigen::MatrixXd bigTC(TC.rows()+sphereTC.rows(),3);
    if (sphereV.rows()!=0){
        std::cout<<"Showing vertex spheres"<<std::endl;
        bigV<<V, sphereV;
        bigT<<T, sphereT+Eigen::MatrixXi::Constant(sphereT.rows(), sphereT.cols(), V.rows());
        bigTC<<TC, sphereTC;
    } else{
        //std::cout<<"Here is good"<<std::endl;
        bigV<<V;
        bigT<<T;
        bigTC<<TC;
    }
    
    viewer.data.clear();
    if (ShowPolygonalEdges){
        viewer.core.show_lines=false;
        std::cout<<"Showing polygonal lines"<<std::endl;
        Eigen::MatrixXd OrigEdgeColors(EV.rows(),3);
        OrigEdgeColors.col(0)=Eigen::VectorXd::Constant(EV.rows(),0.0);
        OrigEdgeColors.col(1)=Eigen::VectorXd::Constant(EV.rows(),0.0);
        OrigEdgeColors.col(2)=Eigen::VectorXd::Constant(EV.rows(),0.0);
        viewer.data.set_edges(V,EV,OrigEdgeColors);
    } else {
        viewer.core.show_lines=true;
    }
    
    //std::cout<<"bigT.size()"<<bigT.size()<<std::endl;
    viewer.data.set_mesh(bigV,bigT);
    //std::cout<<"bigTC.size()"<<bigTC.size()<<std::endl;
    viewer.data.set_colors(bigTC);
}


// This function is called every time a keyboard button is pressed
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
    switch(key){
		case '1': ShowPolygonalEdges=!ShowPolygonalEdges; break;
        case '2': ShowVertexSpheres=!ShowVertexSpheres; break;
        case '3': ShowFaceNormals=!ShowFaceNormals; break;
    }
    
    UpdateCurrentView(viewer);
    return false;
}


int main(int argc, char *argv[])
{
    using namespace std;

    hedra::hedra_read_OFF(DATA_PATH "/rhombitruncated_cubeoctahedron.off", V, D, F);
    
    std::cout<<R"(
    1 Switch polygonal edges/triangulated edges
    2 Switch vertex spheres on/off
    3 Switch face normals on/off
    )";
    
    sphereColors.resize(V.rows(),3);
    for (int i=0;i<V.rows();i++)
        sphereColors.row(i)<<(double)i/(double)V.rows(), 1.0-(double)i/(double)V.rows(), 0.0;
    
   
    hedra::triangulate_mesh(D, F,T,TF);
    hedra::hedra_edge_topology(D,F, EV,FE, EF, EFi, FEs, innerEdges);
    
    TC.resize(T.rows(),3);
    TC.col(0)=Eigen::VectorXd::Constant(T.rows(),1.0);
    TC.col(1)=Eigen::VectorXd::Constant(T.rows(),1.0);
    TC.col(2)=Eigen::VectorXd::Constant(T.rows(),1.0);
    
    
    // Plot the mesh
    igl::viewer::Viewer viewer;
    viewer.callback_key_down = &key_down;
    UpdateCurrentView(viewer);
    viewer.launch();
    
    return 0;
}


