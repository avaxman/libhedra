#include <igl/viewer/Viewer.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/hedra_read_OFF.h>
#include <unistd.h>


int main(int argc, char *argv[])
{
    using namespace std;
    Eigen::MatrixXi F;
    Eigen::VectorXi D;
    Eigen::MatrixXd V;
    char buffer[2000];
    cout<<DATA_PATH<<endl;
    if (!hedra::hedra_read_OFF(DATA_PATH "/rhombitruncated_cubeoctahedron.off", V, D, F))
        cout<<"File not found!"<<endl;
    Eigen::MatrixXi T;
    Eigen::VectorXi TF;
    hedra::triangulate_mesh(D, F,T,TF);
    
    cout<<"V:"<<V<<endl;
    cout<<"F:"<<F<<endl;
    cout<<"D:"<<D<<endl;
    cout<<"T: "<<T<<endl;
    cout<<"TF:"<<TF<<endl;
    
    // Plot the mesh
    igl::viewer::Viewer viewer;
    viewer.data.set_mesh(V, T);
    viewer.data.set_face_based(true);
    viewer.launch();
    
    return 0;
}
