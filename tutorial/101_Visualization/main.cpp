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
  
  Eigen::MatrixXd bigV=V;
  Eigen::MatrixXi bigT=T;
  Eigen::MatrixXd bigTC=TC;
  if (ShowVertexSpheres){
    hedra::point_spheres(V, 0.02, sphereColors, 10, false, true, bigV, bigT, bigTC);
    std::cout<<"Showing vertex spheres"<<std::endl;
  }
  
  if (ShowFaceNormals){
    std::cout<<"Showing face normals"<<std::endl;
    hedra::line_cylinders(faceCenters, faceCenters+faceNormals, 0.02, lineColors, 10, false, true, bigV, bigT, bigTC);
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
  
  hedra::polygonal_read_OFF(TUTORIAL_SHARED_PATH "/rhombitruncated_cubeoctahedron_fixed.off", V, D, F);
  
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


