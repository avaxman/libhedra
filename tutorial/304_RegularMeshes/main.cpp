#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject.h>
#include <igl/project.h>
#include <igl/sortrows.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/invert_diag.h>
#include <math.h>
#include <igl/per_vertex_normals.h>
#include <algorithm>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/polygonal_write_OFF.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/moebius_regular_meshes.h>
#include <hedra/visualization_schemes.h>


typedef enum {ORIGINAL_MESH, REGULAR_MESH} MeshModes;
MeshModes meshMode=REGULAR_MESH;

typedef enum {STANDARD, MOEBIUS_REGULARITY, EUCLIDEAN_REGULARITY, WILLMORE_ENERGY} ViewingModes;
ViewingModes viewingMode=STANDARD;

// Mesh
Eigen::MatrixXd VOrig, VRegular, FEs;
Eigen::MatrixXi F, TF;
Eigen::VectorXi D, innerEdges;
Eigen::MatrixXi EV, FE, EF;


Eigen::VectorXi constVertices;
Eigen::MatrixXd constPoses;

double MRCoeff=1.0;
double ERCoeff=0.01;

hedra::MoebiusRegularData MRData;


void update_mesh(igl::opengl::glfw::Viewer& viewer)
{
  
  //mesh
  viewer.data_list[0].clear();
  viewer.data_list[0].set_face_based(true);
  viewer.data_list[0].show_lines=false;
  
  Eigen::MatrixXd C;  //color
  switch (viewingMode){
    case STANDARD: C = hedra::default_mesh_color(); break;
    case MOEBIUS_REGULARITY: C = hedra::scalar2RGB((meshMode==ORIGINAL_MESH ? origMoebRegularity : optMoebRegularity), 0.0,1.0);
      case EUCLIDEAN_REGULARITY: C = hedra::scalar2RGB((meshMode==ORIGINAL_MESH ? origEucRegularity : optEucRegularity), 0.0,1.0);
      case WILLMORE_ENERGY: C = hedra::scalar2RGB((meshMode==ORIGINAL_MESH ? origWillmore : optWillmore), 0.0,1.0);
  }
  
  viewer.data_list[0].set_mesh((meshMode==ORIGINAL_MESH ? VOrig : VRegular), TF);
  viewer.data_list[0].set_colors(C);
  viewer.data_list[0].set_edges(V,EV,hedra::default_edge_color());
  
  hedra::vertex_spheres((meshMode==ORIGINAL_MESH ? VOrig : VRegular), hedra::passive_handle_color(),  VHandles, FHandles, CHandles);
  viewer.data_list[1].set_mesh(VHandles, FHandles);
  viewer.data_list[1].set_colors(CHandles);
  viewer.data_list[1].show_faces = true;
  viewer.data_list[1].show_lines = false;
}

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
  switch(key)
  {
    default:
      return false;
    case '1': {meshMode++; break;}
    case '2': {viewingMode++; break;}
    case '3': {ERCoeff+=0.1; hedra::compute_moebius_regular(MRData, VRegular, MRCoeff, ERCoeff); break;}
    case '4': {ERCoeff-=0.1; hedra::compute_moebius_regular(MRData, VRegular, MRCoeff, ERCoeff); break;}
  }
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  
  cout<<"1  Show original/Moebius regular mesh"<<endl<<
        "2  Toggle measurements "<<endl<<
        "+  Increase Euclidean Regularity "<<endl<<
  "-  Decrease Euclidean Regularity "<<endl;
  hedra::polygonal_read_OFF(TUTORIAL_SHARED_PATH "/intersection.off", V, D, F);
  hedra::polygonal_edge_topology(D, F, EV, FE, EF, EFi, FEs, innerEdges);
  hedra::triangulate_mesh(D, F,T,TF);
  
  //TODO: code to choose entire boundary
  
  
  hedra::setup_moebius_regular(VOrig, D, F, EV, FE, EF, EFi, FEs, innerEdges, constVertices, MRData);
  hedra::compute_moebius_regular(MRData,  MRCoeff, ERCoeff, constPoses, VRegular);
  
  igl::opengl::glfw::Viewer viewer;
  Viewer.callback_mouse_down = &mouse_down;
  Viewer.callback_mouse_up = &mouse_up;
  Viewer.callback_mouse_move = &mouse_move;
  Viewer.callback_key_down = &key_down;
  Viewer.callback_key_up = &key_up;
  Viewer.callback_init = &init;
  Viewer.core().background_color<<0.75,0.75,0.75,1.0;
  update_mesh(viewer);
  viewer.launch();
  
}
