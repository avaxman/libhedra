#include <algorithm>
#include <math.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/avg_edge_length.h>
#include <igl/boundary_loop.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/polygonal_write_OFF.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/visualization_schemes.h>
#include <hedra/catmull_clark.h>
#include <hedra/dual_mesh.h>
#include <hedra/operator_1246.h>
#include <hedra/simplest_subdivision.h>
#include <hedra/dual_truncation.h>
//#include <hedra/quad_kobbelt.h>

int currViewLevel=0;
int subdLevel=0;

// Subdivision Meshes
std::vector<Eigen::MatrixXd> V;
std::vector<Eigen::MatrixXi> F, T;
std::vector<Eigen::VectorXi> D,TF;
std::vector<Eigen::MatrixXi> EV, FE, EF, EFi;

void update_mesh(igl::opengl::glfw::Viewer& viewer)
{
  viewer.data_list[0].clear();
  viewer.data_list[0].set_face_based(true);
  viewer.data_list[0].show_lines=false;
  viewer.data_list[0].set_mesh(V[currViewLevel], T[currViewLevel]);
  viewer.data_list[0].set_edges(V[currViewLevel],EV[currViewLevel],hedra::default_edge_color());
  viewer.data_list[0].set_colors(hedra::default_mesh_color());
}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  Eigen::MatrixXi FNext;
  Eigen::MatrixXd VNext;
  Eigen::VectorXi DNext;
  switch(key)
  {
    default:
      return false;
    case 'A': {
      currViewLevel=std::min(currViewLevel+1,subdLevel);
      std::cout<<"Viewing Level "<<currViewLevel<<std::endl;
      break;}
      
    case 'S': {
      currViewLevel=std::max(currViewLevel-1,0);
      std::cout<<"Viewing Level "<<currViewLevel<<std::endl;
      break;}
      
    case '2': hedra::catmull_clark(V[subdLevel], D[subdLevel], F[subdLevel], hedra::CANONICAL_MOEBIUS_SUBDIVISION, VNext, DNext, FNext); break;
    /*case '3': hedra::kobbelt_quad(V[subdLevel], D[subdLevel], F[subdLevel], VNext, DNext, FNext, hedra::LINEAR_SUBDIVISION); break;
    case '4': hedra::simplest(V[subdLevel], D[subdLevel], F[subdLevel], VNext, DNext, FNext, hedra::LINEAR_SUBDIVISION); break;
    case '5': hedra::dual_truncation(V[subdLevel], D[subdLevel], F[subdLevel], VNext, DNext, FNext, hedra::LINEAR_SUBDIVISION); break;
    case '6': hedra::dual_mesh(V[subdLevel], D[subdLevel], F[subdLevel], VNext, DNext, FNext, hedra::LINEAR_SUBDIVISION); break;
    case '7': hedra::vertex_insertion(V[subdLevel], D[subdLevel], F[subdLevel], VNext, DNext, FNext, hedra::LINEAR_SUBDIVISION); break;
    case '8': hedra::operator_1246(V[subdLevel], D[subdLevel], F[subdLevel], VNext, DNext, FNext, hedra::LINEAR_SUBDIVISION); break;*/
  }
  
  if ((key >='2')&&( key<='8')){
 
    Eigen::MatrixXi TNext;
    Eigen::VectorXi TFNext, innerEdgesNext;
    Eigen::MatrixXi EVNext, FENext, EFNext, EFiNext;
    Eigen::MatrixXd FEsNext;
    
    hedra::polygonal_edge_topology(DNext, FNext, EVNext, FENext, EFNext, EFiNext, FEsNext, innerEdgesNext);
    hedra::triangulate_mesh(DNext,FNext,TNext, TFNext);
    
    V.push_back(VNext);
    F.push_back(FNext);
    D.push_back(DNext);
    T.push_back(TNext);
    TF.push_back(TFNext);
    EV.push_back(EVNext);
    subdLevel++;
    currViewLevel=subdLevel;

  }
  
  update_mesh(viewer);
  return true;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  
  cout<<"1  Show original/Moebius regular mesh"<<endl<<
  "2  Catmull-Clark Subdivision "<<endl<<
  "3  Kobbelt Quad Subdivision "<<endl<<
  "4  Simplest Subdivision"<<endl<<
  "5  Dual Truncation "<<endl<<
  "6  Dual Mesh "<<endl<<
  "7  Vertex Insertion "<<endl<<
  "8  1246 Operator "<<endl<<
  "A  show Next Level"<<endl<<
  "S  show previous level"<<endl;
  
  Eigen::MatrixXd V0;
  Eigen::MatrixXi F0, T0;
  Eigen::VectorXi D0,TF0, innerEdges0;
  Eigen::MatrixXi EV0, FE0, EF0, EFi0;
  Eigen::MatrixXd FEs;
  
  hedra::polygonal_read_OFF(TUTORIAL_SHARED_PATH "/cube.off", V0, D0, F0);
  hedra::polygonal_edge_topology(D0, F0, EV0, FE0, EF0, EFi0, FEs, innerEdges0);
  hedra::triangulate_mesh(D0,F0,T0, TF0);
  
  V.push_back(V0);
  F.push_back(F0);
  D.push_back(D0);
  T.push_back(T0);
  TF.push_back(TF0);
  EV.push_back(EV0);
  /*EF.push_back(EF0);
  EFi.push_back(EFi0);
  innerEdges.push_back(innerEdges0);
  FEs.push_back(FEs0);
  FE.push_back(FE0);*/
  
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  
  viewer.core.background_color<<0.75,0.75,0.75,1.0;
  update_mesh(viewer);
  viewer.launch();
}
