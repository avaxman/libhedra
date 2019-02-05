#include <algorithm>
#include <math.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <igl/edge_topology.h>
#include <hedra/copyleft/cgal/generate_mesh.h>
#include <hedra/polygonal_write_OFF.h>

Eigen::MatrixXd V, TV, newV;
Eigen::MatrixXi F, newF, FTC;
Eigen::VectorXi newD;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  return true;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  
  cout<<"0: show original mesh with parameterization"<<endl<<
  "1: show UV coordinates"<<endl<<
  "2: show arrangement with tiling"<<endl<<
  "3: Show final mesh"<<endl;

  Eigen::MatrixXd TC;
  Eigen::MatrixXd N;
  Eigen::MatrixXd FN;
  Eigen::MatrixXi EF, FE, EV;
  igl::readOBJ(TUTORIAL_SHARED_PATH "/Interior-parametrization.obj", V, TC, N, F, FTC, FN);
  
  //OBJ values are skewed to fit PNG as follows
  //OutputCornerValues=[CornerValues(:,1)/3 CornerValues(:,2)/sqrt(3)];
  TC.col(0).array()*=3;
  TC.col(1).array()*=SQRT3;
  
  igl::edge_topology(V,F,EV,FE,EF);
  
  Eigen::RowVector3d spans=V.colwise().maxCoeff()-V.colwise().minCoeff();
  hedra::copyleft::cgal::generate_mesh(6, V, F, EV, FE, EF, TC, FTC, 0.01*spans.maxCoeff(),newV, newD, newF);
  hedra::polygonal_write_OFF(std::string(TUTORIAL_SHARED_PATH "/Interior-hex.off"),newV,newD,newF);
  
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  
  //viewer.core.background_color<<0.75,0.75,0.75,1.0;
  
  //edges mesh
  viewer.append_mesh();
  
  //control polygon mesh
  viewer.append_mesh();
  
  viewer.selected_data_index=0;
  //update_mesh(viewer);
  viewer.launch();
}
