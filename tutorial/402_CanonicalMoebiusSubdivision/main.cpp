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
#include <igl/unproject_onto_mesh.h>
//#include <hedra/quad_kobbelt.h>

#define SUBDIVISION_LEVELS 3

int currViewLevel=SUBDIVISION_LEVELS-1;
int subdivisionLevels=3;
bool editingMode = false;
bool choosingHandleMode = false;
double currWinZ;

std::vector<int> handles;
std::vector<Eigen::RowVector3d> handlePoses;
int currHandle;

// Subdivision Meshes
Eigen::MatrixXd V[SUBDIVISION_LEVELS], VEdges[SUBDIVISION_LEVELS], CEdges[SUBDIVISION_LEVELS], FEs[SUBDIVISION_LEVELS];
Eigen::MatrixXi F[SUBDIVISION_LEVELS], T[SUBDIVISION_LEVELS], TEdges[SUBDIVISION_LEVELS];
Eigen::VectorXi D[SUBDIVISION_LEVELS],TF[SUBDIVISION_LEVELS], innerEdges[SUBDIVISION_LEVELS];
Eigen::MatrixXi EV[SUBDIVISION_LEVELS], FE[SUBDIVISION_LEVELS], EF[SUBDIVISION_LEVELS], EFi[SUBDIVISION_LEVELS];

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
  viewer.data_list[2].set_mesh(VEdges[0], TEdges[0]);
  viewer.data_list[2].set_face_based(true);
  viewer.data_list[2].set_colors(hedra::active_handle_color().replicate(TEdges[0].rows(),1));
  viewer.data_list[2].show_lines=false;
}


void update_subdivision()
{
  hedra::polygonal_edge_lines(V[0], F[0], T[0], EV[0],  VEdges[0], TEdges[0], CEdges[0],0.25);
  for (int i=0;i<SUBDIVISION_LEVELS-1;i++){
    hedra::catmull_clark(V[i], D[i], F[i], hedra::CANONICAL_MOEBIUS_SUBDIVISION, V[i+1], D[i+1], F[i+1]);
    hedra::polygonal_edge_topology(D[i+1], F[i+1], EV[i+1], FE[i+1], EF[i+1], EFi[i+1], FEs[i+1], innerEdges[i+1]);
    hedra::triangulate_mesh(D[i+1],F[i+1],T[i+1], TF[i+1]);
    hedra::polygonal_edge_lines(V[i+1], F[i+1], T[i+1], EV[i+1],  VEdges[i+1], TEdges[i+1], CEdges[i+1],0.25);
  }
}

bool mouse_move(igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y)
{
  if (!editingMode)
    return false;
  
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  Eigen::RowVector3f NewPos=igl::unproject<float>(Eigen::Vector3f(x,y,currWinZ),
                                                  viewer.core.view,
                                                  viewer.core.proj,
                                                  viewer.core.viewport);
  
  handlePoses[handlePoses.size()-1]=NewPos.cast<double>();
  
  Eigen::MatrixXd A(3*F[0].rows(),3);
  for (int i=0;i<F[0].rows();i++)
    A.block(3*i,0,3,3)=Eigen::Matrix3d::Identity();
  
  Eigen::MatrixXd bc(handles.size(),V[0].cols());
  for (int i=0;i<handles.size();i++)
    V[0].row(handles[i])=handlePoses[i].transpose();
  update_subdivision();
  update_mesh(viewer);
  return true;
  
}


bool mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
  if (((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Left))
    return false;
  
  editingMode=false;
  
  return true;
}

bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
  if (((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Left))
    return false;
  int vid, fid;
  Eigen::Vector3f bc;
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  if (!choosingHandleMode){
    editingMode=true;
    return false;
  }

  if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core.view,
                              viewer.core.proj, viewer.core.viewport, V[0], F[0], fid, bc))
  {
    //add the closest vertex to the handles
    Eigen::MatrixXf::Index maxRow, maxCol;
    bc.maxCoeff(&maxRow);
    int CurrVertex=F[0](fid, maxRow);
    bool Found=false;
    for (int i=0;i<handles.size();i++)
      if (handles[i]==CurrVertex){
        CurrVertex=handles[i];
        Found=true;
      }
    
    if (!Found){
      handles.push_back(CurrVertex);
      handlePoses.push_back(V[0].row(CurrVertex));
    }
    
    Eigen::Vector3f WinCoords=igl::project<float>(V[0].row(CurrVertex).cast<float>(),
                                                  viewer.core.view,
                                                  viewer.core.proj,
                                                  viewer.core.viewport);
    
    currWinZ=WinCoords(2);
    std::cout<<"Choosing Vertex :"<<CurrVertex<<std::endl;
    
    Eigen::VectorXi b(handles.size());
    for (int i=0;i<handles.size();i++)
      b(i)=handles[i];
    
    update_mesh(viewer);
  }
  return true;
}

bool key_up(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  switch(key)
  {
    case '0': choosingHandleMode=false;
      break;
  }
  return false;
}


bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{

  switch(key)
  {
    case '0': choosingHandleMode=true;
      break;
      
    case '1': {
      currViewLevel=std::min(currViewLevel+1,SUBDIVISION_LEVELS);
      std::cout<<"Viewing Level "<<currViewLevel<<std::endl;
      break;}
      
    case '2': {
      currViewLevel=std::max(currViewLevel-1,0);
      std::cout<<"Viewing Level "<<currViewLevel<<std::endl;
      break;}
    default:
      return false;
  }
  update_mesh(viewer);
  return true;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  
  cout<<"0+Right-click  choose handle"<<endl<<
  "right+drag:  move handle"<<endl<<
  "1  Show previous level"<<endl<<
  "2  Show next level"<<endl;

  hedra::polygonal_read_OFF(TUTORIAL_SHARED_PATH "/pipes_all_quads.off", V[0], D[0], F[0]);
  hedra::polygonal_edge_topology(D[0], F[0], EV[0], FE[0], EF[0], EFi[0], FEs[0], innerEdges[0]);
  hedra::triangulate_mesh(D[0],F[0],T[0], TF[0]);
  hedra::polygonal_edge_lines(V[0], F[0], T[0], EV[0],  VEdges[0], TEdges[0], CEdges[0],0.25);
  
  update_subdivision();
  
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.callback_mouse_move = &mouse_move;
  viewer.callback_mouse_up=&mouse_up;
  
  //viewer.core.background_color<<0.75,0.75,0.75,1.0;
  
  //edges mesh
  viewer.append_mesh();
  
  //control polygon mesh
  viewer.append_mesh();
  
  viewer.selected_data_index=0;
  update_mesh(viewer);
  viewer.launch();
}
