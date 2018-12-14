#include <algorithm>
#include <math.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/avg_edge_length.h>
#include <igl/boundary_loop.h>
#include <hedra/polygonal_edge_lines.h>
#include <hedra/vertex_insertion.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/polygonal_write_OFF.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/visualization_schemes.h>
#include <hedra/catmull_clark.h>
#include <hedra/simplest_subdivision.h>
#include <hedra/dual_mesh.h>
#include <hedra/simplest_subdivision.h>
#include <hedra/dual_truncation.h>
#include <hedra/subdivision_basics.h>
#include <hedra/operator_1264.h>
//#include <hedra/quad_kobbelt.h>

int currViewLevel=0;
int subdLevel=0;
bool showControlPolygon=false;

// Subdivision Meshes
std::vector<Eigen::MatrixXd> V, VEdges, CEdges;
std::vector<Eigen::MatrixXi> F, T, TEdges;
std::vector<Eigen::VectorXi> D,TF;
std::vector<Eigen::MatrixXi> EV, FE, EF, EFi;

void update_mesh(igl::opengl::glfw::Viewer& viewer)
{
  
  //current mesh
  viewer.data_list[0].clear();
  viewer.data_list[0].set_face_based(true);
  viewer.data_list[0].show_lines=false;
  viewer.data_list[0].set_mesh(V[currViewLevel], T[currViewLevel]);
  viewer.data_list[0].set_colors(hedra::default_mesh_color());
  
  //edges
  viewer.data_list[1].clear();
  viewer.data_list[1].set_mesh(VEdges[currViewLevel], TEdges[currViewLevel]);
  viewer.data_list[1].set_face_based(true);
  viewer.data_list[1].set_colors(CEdges[currViewLevel]);
  viewer.data_list[1].show_lines=false;
  
  
  //control polygon
  viewer.data_list[2].clear();
  if (showControlPolygon){
    viewer.data_list[2].set_mesh(VEdges[0], TEdges[0]);
    viewer.data_list[2].set_face_based(true);
    viewer.data_list[2].set_colors(hedra::active_handle_color().replicate(TEdges[0].rows(),1));
    viewer.data_list[2].show_lines=false;
  }
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
      
    case '1': showControlPolygon=!showControlPolygon; break;
    case '2': hedra::catmull_clark(V[subdLevel], D[subdLevel], F[subdLevel], hedra::CANONICAL_MOEBIUS_SUBDIVISION, VNext, DNext, FNext); break;
    case '3': hedra::simplest_subdivision(V[subdLevel], D[subdLevel], F[subdLevel], hedra::CANONICAL_MOEBIUS_SUBDIVISION, VNext, DNext, FNext); break;
    case '4': hedra::vertex_insertion(V[subdLevel], D[subdLevel], F[subdLevel], hedra::CANONICAL_MOEBIUS_SUBDIVISION, VNext, DNext, FNext); break;
    case '5': hedra::dual_truncation(V[subdLevel], D[subdLevel], F[subdLevel], hedra::CANONICAL_MOEBIUS_SUBDIVISION, VNext, DNext, FNext); break;
    case '6': hedra::operator_1264(V[subdLevel], D[subdLevel], F[subdLevel], hedra::CANONICAL_MOEBIUS_SUBDIVISION, VNext, DNext, FNext); break;
    case '7': hedra::dual_mesh(V[subdLevel], D[subdLevel], F[subdLevel], hedra::CANONICAL_MOEBIUS_SUBDIVISION, VNext, DNext, FNext); break;
      
  }
  
  if ((key >='2')&&( key<='8')){
    
    Eigen::MatrixXi TNext, TEdgesNext;
    Eigen::VectorXi TFNext, innerEdgesNext;
    Eigen::MatrixXi EVNext, FENext, EFNext, EFiNext;
    Eigen::MatrixXd FEsNext, VEdgesNext, CEdgesNext;
    
    hedra::polygonal_edge_topology(DNext, FNext, EVNext, FENext, EFNext, EFiNext, FEsNext, innerEdgesNext);
    hedra::triangulate_mesh(DNext,FNext,TNext, TFNext);
    
    hedra::polygonal_edge_lines(VNext, FNext, TNext, EVNext,  VEdgesNext, TEdgesNext, CEdgesNext,0.25);
    
    V.push_back(VNext);
    F.push_back(FNext);
    D.push_back(DNext);
    T.push_back(TNext);
    TF.push_back(TFNext);
    EV.push_back(EVNext);
    VEdges.push_back(VEdgesNext);
    TEdges.push_back(TEdgesNext);
    CEdges.push_back(CEdgesNext);
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
  
  cout<<"1  Show control polygon mesh"<<endl<<
  "2  Catmull-Clark Subdivision "<<endl<<
  "3  Simplest Subdivision"<<endl<<
  "4  Vertex Insertion "<<endl<<
  "5  Dual Truncation "<<endl<<
  "6  1264 Operator "<<endl<<
  "7  Dual Mesh "<<endl<<
  "S  show Next Level"<<endl<<
  "A  show previous level"<<endl;
  
  Eigen::MatrixXd V0, VEdges0;
  Eigen::MatrixXi F0, T0, TEdges0;
  Eigen::VectorXi D0,TF0, innerEdges0;
  Eigen::MatrixXi EV0, FE0, EF0, EFi0;
  Eigen::MatrixXd FEs, CEdges0;
  
  hedra::polygonal_read_OFF(TUTORIAL_SHARED_PATH "/pipes_all_quads.off", V0, D0, F0);
  hedra::polygonal_edge_topology(D0, F0, EV0, FE0, EF0, EFi0, FEs, innerEdges0);
  hedra::triangulate_mesh(D0,F0,T0, TF0);
  hedra::polygonal_edge_lines(V0, F0, T0, EV0,  VEdges0, TEdges0, CEdges0,0.25);
  
  V.push_back(V0);
  F.push_back(F0);
  D.push_back(D0);
  T.push_back(T0);
  TF.push_back(TF0);
  EV.push_back(EV0);
  VEdges.push_back(VEdges0);
  TEdges.push_back(TEdges0);
  CEdges.push_back(CEdges0);
  
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  
  //viewer.core.background_color<<0.75,0.75,0.75,1.0;
  
  //edges mesh
  viewer.append_mesh();
  
  //control polygon mesh
  viewer.append_mesh();
  
  update_mesh(viewer);
  viewer.launch();
}
