#include <igl/viewer/Viewer.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/hedra_read_OFF.h>
#include <hedra/hedra_edge_topology.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F,EV, EF, FE,T;
Eigen::VectorXi D,TF;

// This function is called every time a keyboard button is pressed
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
    if (key == '1')
    {
        // Clear should be called before drawing the mesh
        std::cout<<"Showing polygonal lines"<<std::endl;
        Eigen::MatrixXd OrigEdgeColors(EV.rows(),3);
        OrigEdgeColors.col(0)=Eigen::VectorXd::Constant(EV.rows(),1.0);
        OrigEdgeColors.col(1)=Eigen::VectorXd::Constant(EV.rows(),1.0);
        OrigEdgeColors.col(2)=Eigen::VectorXd::Constant(EV.rows(),0.0);
        viewer.data.set_edges(V,EV,OrigEdgeColors);
    }
    else if (key == '2')
    {
        std::cout<<"Showing triangulated lines"<<std::endl;
        viewer.data.clear();
        viewer.data.set_mesh(V, T);
    }
    
    return false;
}


int main(int argc, char *argv[])
{
    using namespace std;

    hedra::hedra_read_OFF(DATA_PATH "/rhombitruncated_cubeoctahedron.off", V, D, F);
    
    std::cout<<R"(
    1 Switch to polygonal edges
    2 Switch to triangulated edges
    )";

    hedra::triangulate_mesh(D, F,T,TF);
    hedra::hedra_edge_topology(D,F,EV,FE,EF);

    // Plot the mesh
    igl::viewer::Viewer viewer;
    viewer.callback_key_down = &key_down;
    viewer.data.set_mesh(V, T);
    viewer.data.set_face_based(true);
    viewer.launch();
    
    return 0;
}

