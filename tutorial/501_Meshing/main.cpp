#include <algorithm>
#include <math.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <hedra/copyleft/cgal/extract_mesh.h>

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
  igl::readOBJ(TUTORIAL_SHARED_PATH "/lilium_param.obj", V, TC, N, F, FTC, FN);
  
  hedra::copyleft::cgal::extract_mesh(V, F, TC, FTC, newV, newD, newF);
  
  
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
