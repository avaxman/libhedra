#include <igl/unproject_onto_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readDMAT.h>
#include <igl/jet.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/affine_maps_deform.h>
#include <hedra/point_spheres.h>
#include <hedra/planarity.h>
#include <hedra/visualization_schemes.h>
#include <hedra/scalar2RGB.h>


std::vector<int> Handles;
std::vector<Eigen::RowVector3d> HandlePoses;
int CurrentHandle;
Eigen::MatrixXd VOrig, V;
Eigen::MatrixXi F, T;
Eigen::VectorXi D, TF;
Eigen::MatrixXi EV, EF, FE, EFi;
Eigen::MatrixXd FEs;
Eigen::VectorXi innerEdges;
Eigen::Vector3d spans;
bool Editing=false;
bool ChoosingHandleMode=false;
double CurrWinZ;
hedra::AffineData affine_data;


bool UpdateCurrentView(igl::opengl::glfw::Viewer & viewer)
{
  using namespace Eigen;
  using namespace std;
  
  MatrixXd bc(Handles.size(),V.cols());
  for (int i=0;i<Handles.size();i++)
    bc.row(i)=HandlePoses[i].transpose();
  
  VectorXd planarity;
  hedra::planarity(V,D,F,planarity);
  
  Eigen::MatrixXd C, TC;
  hedra::scalar2RGB(planarity, 0.0,1.0, C);
  TC.resize(T.rows(),3);
  for (int i=0;i<T.rows();i++)
    TC.row(i)=C.row(TF(i));
  
  double sphereRadius=spans.sum()/200.0;
  MatrixXd sphereGreens(Handles.size(),3);
  sphereGreens.col(0).setZero();
  sphereGreens.col(1).setOnes();
  sphereGreens.col(2).setZero();
  
  MatrixXd VHandles;
  MatrixXi THandles;
  MatrixXd CHandles;
  
  //hedra::point_spheres(bc, sphereRadius, sphereGreens, 10,VHandle, THandle, CHandle);
  
  double radius = 0.25*igl::avg_edge_length(V, T);
  hedra::point_spheres(bc, radius, (hedra::passive_handle_color()).replicate(F.rows(),1), 8, VHandles, FHandles, CHandles);
  viewer.data_list[1].clear();
  viewer.data_list[1].set_mesh(VHandles, FHandles);
  viewer.data_list[1].set_face_based(true);
  viewer.data_list[1].set_colors(CHandles);
  viewer.data_list[1].show_faces = true;
  viewer.data_list[1].show_lines = false;
  
  Eigen::MatrixXd OrigEdgeColors(EV.rows(),3);
  OrigEdgeColors.col(0)=Eigen::VectorXd::Constant(EV.rows(),0.0);
  OrigEdgeColors.col(1)=Eigen::VectorXd::Constant(EV.rows(),0.0);
  OrigEdgeColors.col(2)=Eigen::VectorXd::Constant(EV.rows(),0.0);
  
  viewer.data_list[0].clear();
  viewer.data_list[0].set_mesh(V,T);
  viewer.data_list[0].set_colors(TC);
  //viewer.data().compute_normals();
  viewer.data_list[0].set_edges(V,EV,OrigEdgeColors);
  
  viewer.data_list[1].clear();
  viewer.data_list[1].show_lines=false;
  viewer.data_list[1].set_mesh(VHandle,THandle);
  viewer.data_list[1].set_colors(CHandle);
  
  return true;
}

bool mouse_move(igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y)
{
  if (!Editing)
    return false;
  
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  Eigen::RowVector3f NewPos=igl::unproject<float>(Eigen::Vector3f(x,y,CurrWinZ),
                                                  viewer.core.view,
                                                  viewer.core.proj,
                                                  viewer.core.viewport);
  
  HandlePoses[HandlePoses.size()-1]=NewPos.cast<double>();
  
  Eigen::MatrixXd A(3*F.rows(),3);
  for (int i=0;i<F.rows();i++)
    A.block(3*i,0,3,3)=Eigen::Matrix3d::Identity();
  
  Eigen::MatrixXd bc(Handles.size(),V.cols());
  for (int i=0;i<Handles.size();i++)
    bc.row(i)=HandlePoses[i].transpose();
  hedra::affine_maps_deform(affine_data,bc,5, V);
  UpdateCurrentView(viewer);
  return true;
  
}


bool mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
  if (((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Left))
    return false;
  
  Editing=false;
  
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
  if (!ChoosingHandleMode){
    Editing=true;
    return false;
  }
  if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core.view,
                              viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
  {
    //add the closest vertex to the handles
    Eigen::MatrixXf::Index maxRow, maxCol;
    bc.maxCoeff(&maxRow);
    int CurrVertex=F(fid, maxRow);
    bool Found=false;
    for (int i=0;i<Handles.size();i++)
      if (Handles[i]==CurrVertex){
        CurrVertex=Handles[i];
        Found=true;
      }
    
    if (!Found){
      Handles.push_back(CurrVertex);
      HandlePoses.push_back(V.row(CurrVertex));
    }
    
    Eigen::Vector3f WinCoords=igl::project<float>(V.row(CurrVertex).cast<float>(),
                                                  viewer.core.view,
                                                  viewer.core.proj,
                                                  viewer.core.viewport);
    
    CurrWinZ=WinCoords(2);
    std::cout<<"Choosing Vertex :"<<CurrVertex<<std::endl;
    
    Eigen::VectorXi b(Handles.size());
    for (int i=0;i<Handles.size();i++)
      b(i)=Handles[i];
    hedra::affine_maps_precompute(V,D,F,EV,EF,EFi, FE, b, 3, affine_data);
    
    UpdateCurrentView(viewer);
  }
  return true;
}

bool key_up(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  switch(key)
  {
      
    case '1': ChoosingHandleMode=false;
      break;
  }
  return false;
}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  switch(key)
  {
    case '1': ChoosingHandleMode=true;
      break;
  }
  return false;
}


int main(int argc, char *argv[])
{
  
  // Load a mesh in OFF format
  using namespace std;
  using namespace Eigen;
  
  cout<<"press 1+right button to select new handles"<<endl;
  cout<<"press the right button and drag the edit the mesh"<<endl;
  
  hedra::polygonal_read_OFF(TUTORIAL_SHARED_PATH "/six.off", V, D, F);
  hedra::polygonal_edge_topology(D, F, EV, FE, EF,EFi,FEs,innerEdges);
  hedra::triangulate_mesh(D, F, T, TF);
  spans=V.colwise().maxCoeff()-V.colwise().minCoeff();
  
  VOrig=V;
  igl::opengl::glfw::Viewer viewer;
  
  viewer.append_mesh();
  viewer.selected_data_index=0;
  viewer.callback_mouse_down = &mouse_down;
  viewer.callback_mouse_move = &mouse_move;
  viewer.callback_mouse_up=&mouse_up;
  viewer.callback_key_down=&key_down;
  viewer.callback_key_up=&key_up;
  viewer.core.background_color<<0.75,0.75,0.75,1.0;
  UpdateCurrentView(viewer);
  viewer.launch();
}
