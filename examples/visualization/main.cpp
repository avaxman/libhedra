#include <igl/viewer/Viewer.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/polyhedral_face_normals.h>
#include <hedra/polygonal_face_centers.h>
#include <hedra/point_spheres.h>
#include <hedra/line_cylinders.h>

Eigen::MatrixXd V, FEs;
Eigen::MatrixXi F,EV, EF, FE,T, EFi;
Eigen::VectorXi D,TF, innerEdges;
Eigen::MatrixXd TC;
Eigen::MatrixXd sphereColors, lineColors;
Eigen::MatrixXd faceCenters, faceNormals;

bool ShowPolygonalEdges=true;
bool ShowVertexSpheres=false;
bool ShowFaceNormals=false;

void UpdateCurrentView(igl::viewer::Viewer& viewer)
{
    viewer.data.clear();
    viewer.data.set_face_based(true);
    
    Eigen::MatrixXd sphereV,lineV, sphereTC, lineTC;
    Eigen::MatrixXi sphereT, lineT;
    if (ShowVertexSpheres)
        hedra::point_spheres(V, 0.02, sphereColors, 10, false, sphereV, sphereT, sphereTC);
    
    if (ShowFaceNormals)
        hedra::line_cylinders(faceCenters, faceCenters+faceNormals, 0.02, lineColors, 10, false, lineV, lineT, lineTC);
    
    
    Eigen::MatrixXd bigV(V.rows()+sphereV.rows()+lineV.rows(),3);
    Eigen::MatrixXi bigT(T.rows()+sphereT.rows()+lineT.rows(),3);
    Eigen::MatrixXd bigTC(TC.rows()+sphereTC.rows()+lineTC.rows(),3);
    
    bigV.block(0,0,V.rows(),3)=V;
    bigT.block(0,0,T.rows(),3)=T;
    bigTC.block(0,0,TC.rows(),3)=TC;
    if (sphereV.rows()!=0){
        std::cout<<"Showing vertex spheres"<<std::endl;
        bigV.block(V.rows(),0,sphereV.rows(),3)=sphereV;
        bigT.block(T.rows(),0,sphereT.rows(),3)=sphereT+Eigen::MatrixXi::Constant(sphereT.rows(), sphereT.cols(), V.rows());
        bigTC.block(TC.rows(),0,sphereTC.rows(),3)=sphereTC;
    }
    if (lineV.rows()!=0){
        std::cout<<"Showing face normals"<<std::endl;
        bigV.block(V.rows()+sphereV.rows(),0,lineV.rows(),3)=lineV;
        bigT.block(T.rows()+sphereT.rows(),0,lineT.rows(),3)=lineT+Eigen::MatrixXi::Constant(lineT.rows(), lineT.cols(), V.rows()+sphereV.rows());
        bigTC.block(TC.rows()+sphereTC.rows(),0,lineTC.rows(),3)=lineTC;
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
    
    viewer.data.set_mesh(bigV,bigT);
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

    hedra::polygonal_read_OFF(DATA_PATH "/rhombitruncated_cubeoctahedron_fixed.off", V, D, F);
    
    std::cout<<R"(
    1 Switch polygonal edges/triangulated edges
    2 Switch vertex spheres on/off
    3 Switch face normals on/off
    )";
    
    sphereColors.resize(V.rows(),3);
    for (int i=0;i<V.rows();i++)
        sphereColors.row(i)<<(double)i/(double)V.rows(), 1.0-(double)i/(double)V.rows(), 0.0;
    
    hedra::triangulate_mesh(D, F,T,TF);
    hedra::polygonal_edge_topology(D,F, EV,FE, EF, EFi, FEs, innerEdges);
    hedra::polygonal_face_centers(V,D, F,faceCenters);
    hedra::polyhedral_face_normals(V,D, F,faceNormals);
    
    lineColors.resize(T.rows(),3);
    lineColors.col(0)=Eigen::VectorXd::Constant(T.rows(),0.5);
    lineColors.col(1)=Eigen::VectorXd::Constant(T.rows(),0.5);
    lineColors.col(2)=Eigen::VectorXd::Constant(T.rows(),1.0);
    
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


