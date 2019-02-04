#include <igl/unproject_onto_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readDMAT.h>
#include <igl/jet.h>
#include <igl/per_vertex_normals.h>
#include <igl/avg_edge_length.h>
#include <hedra/edge_mesh.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/affine_maps_deform.h>
#include <hedra/point_spheres.h>
#include <hedra/concyclity.h>
#include <hedra/dc_error.h>
#include <hedra/ia_error.h>
#include <hedra/qc_error.h>
#include <hedra/concyclity.h>
#include <hedra/scalar2RGB.h>
#include <hedra/quat_moebius_deform.h>
#include <hedra/visualization_schemes.h>
#include <algorithm>

using namespace Eigen;
using namespace std;


igl::opengl::glfw::Viewer viewer;
typedef enum {DC_ERROR, IA_ERROR, FACE_CONCYCLITY, QC_ERROR} ViewingModes;
ViewingModes viewingMode=FACE_CONCYCLITY;

bool showDeformed=true;
std::vector<int> handles;
std::vector<Eigen::RowVector3d> handlePoses;
int currHandle;

Eigen::MatrixXd origV, deformV, currV, VHandles, CHandles;
Eigen::MatrixXd edgeOrigV, edgeDeformV;
Eigen::MatrixXi F, edgeT, currT, polyT, THandles;
Eigen::VectorXi D, polyTF,edgeTE, currTF;
Eigen::MatrixXi EV, EF, FE, EFi, I4;
Eigen::MatrixXd FEs;
Eigen::VectorXi innerEdges;
Eigen::Vector3d spans;
bool isExactMC=false;

bool editing=false;
bool choosingHandleMode=false;
double currWinZ;
hedra::QuatMoebiusData qmdata;

double rigidityFactor=1.0;      //inversion control ratio
double DCSensitivity=0.01;     //ratio of abs(cr)
double IASensitivity=1.0;      //circle IA degrees difference between faces
double QCSensitivity=0.2;      //abs(max_sing/min_sing)
double circSensitivity=5.0;    //circle IA difference

void choose_current_mesh(){
  //choosing mesh
  switch(viewingMode){
    case FACE_CONCYCLITY:
    case QC_ERROR:
      currV=(showDeformed ? deformV : origV);
      currT=polyT;
      currTF=polyTF;
      break;
      
    case DC_ERROR:
    case IA_ERROR:
      currV=(showDeformed ? edgeDeformV : edgeOrigV);
      currT=edgeT;
      currTF=edgeTE;
      break;
  }
}

void UpdateCurrentView()
{
  choose_current_mesh();
  MatrixXd bc(handles.size(),currV.cols());
  for (int i=0;i<handles.size();i++)
    bc.row(i)=handlePoses[i].transpose();
  
  Eigen::MatrixXd C, TC;
  VectorXd scalarMeasure;
  switch(viewingMode){
    case FACE_CONCYCLITY:
      hedra::concyclity(deformV,D,F,scalarMeasure);
      scalarMeasure/=circSensitivity;
      break;
    case QC_ERROR:
      hedra::qc_error(origV, deformV,D,F,scalarMeasure);
      scalarMeasure/=QCSensitivity;
      break;
    case DC_ERROR:
      hedra::dc_error(origV, deformV, D, F, EV, EF, EFi, innerEdges, scalarMeasure);
      scalarMeasure/=DCSensitivity;
      break;
    case IA_ERROR:
      hedra::ia_error(origV, deformV, D, F, EV, EF, EFi, innerEdges, scalarMeasure);
      scalarMeasure/=IASensitivity;
      break;
  }
  
  
  hedra::scalar2RGB(scalarMeasure, 0.0,1.0, C);
  TC.resize(currT.rows(),3);
  for (int i=0;i<currTF.rows();i++)
    TC.row(i)=C.row(currTF(i));

  Eigen::MatrixXd OrigEdgeColors(EV.rows(),3);
  OrigEdgeColors.col(0)=Eigen::VectorXd::Constant(EV.rows(),0.0);
  OrigEdgeColors.col(1)=Eigen::VectorXd::Constant(EV.rows(),0.0);
  OrigEdgeColors.col(2)=Eigen::VectorXd::Constant(EV.rows(),0.0);
  
  viewer.data_list[0].clear();
  viewer.data_list[0].set_mesh(currV,currT);
  viewer.data_list[0].set_colors(TC);
  //viewer.data().compute_normals();
  viewer.data_list[0].set_edges(currV,EV,OrigEdgeColors);
  
  double radius = spans.sum()/200.0;
  hedra::point_spheres(bc, radius, (hedra::passive_handle_color()).replicate(F.rows(),1), 8, VHandles, THandles, CHandles);
  if (!handles.empty()){
    viewer.data_list[1].clear();
    viewer.data_list[1].show_lines=false;
    viewer.data_list[1].set_mesh(VHandles, THandles);
    viewer.data_list[1].set_face_based(true);
    viewer.data_list[1].set_colors(CHandles);
    viewer.data_list[1].show_faces = true;
    viewer.data_list[1].show_lines = false;
  }
}

bool mouse_move(igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y)
{
  if (!editing)
    return false;
  
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  Eigen::RowVector3f NewPos=igl::unproject<float>(Eigen::Vector3f(x,y,currWinZ),
                                                  viewer.core.view,
                                                  viewer.core.proj,
                                                  viewer.core.viewport);
  
  handlePoses[handlePoses.size()-1]=NewPos.cast<double>();
  UpdateCurrentView();
  return true;
  
}


bool mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
  if (((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Left))
    return false;

  if (editing==true){
    Eigen::MatrixXd bc(handles.size(),currV.cols());
    for (int i=0;i<handles.size();i++)
      bc.row(i)=handlePoses[i].transpose();
  
    hedra::quat_moebius_deform(qmdata, 1.0, rigidityFactor, isExactMC, bc,true, deformV);
    hedra::edge_mesh(deformV,D,F,EV, EF, edgeDeformV, edgeT,edgeTE);
    UpdateCurrentView();
  }
  editing=false;
  return true;
}

bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
  if (((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Left))
    return false;
  int vid, fid;
  Eigen::Vector3f bc;
  if (!choosingHandleMode){
    editing=true;
    return false;
  }
  
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core.view,
                              viewer.core.proj, viewer.core.viewport, currV, currT, fid, bc)){
    //add the closest vertex to the handles
    Eigen::MatrixXf::Index maxRow, maxCol;
    bc.maxCoeff(&maxRow);
    int CurrVertex=currT(fid, maxRow);
    bool Found=false;
    for (int i=0;i<handles.size();i++)
      if (handles[i]==CurrVertex){
        CurrVertex=handles[i];
        Found=true;
      }
    
    if (!Found){
      handles.push_back(CurrVertex);
      handlePoses.push_back(currV.row(CurrVertex));
    }
    
    Eigen::Vector3f WinCoords=igl::project<float>(currV.row(CurrVertex).cast<float>(),
                                                  viewer.core.view,
                                                  viewer.core.proj,
                                                  viewer.core.viewport);
    
    currWinZ=WinCoords(2);
    std::cout<<"Choosing Vertex :"<<CurrVertex<<std::endl;
    
    Eigen::VectorXi b(handles.size());
    for (int i=0;i<handles.size();i++)
      b(i)=handles[i];
    
    hedra::quat_moebius_precompute(b, qmdata);
    UpdateCurrentView();
  }
  return true;
}

bool key_up(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  switch(key)
  {
      
    case '1': choosingHandleMode=false;
      break;
  }
  return false;
}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  switch(key)
  {
    case '1':{
      choosingHandleMode=true;
      break;
    }
      
    case '2': {
      viewingMode=(ViewingModes)(((int)viewingMode+1)%4);
      UpdateCurrentView();
      switch(viewingMode){
        case DC_ERROR: cout<<"Showing discrete conformal error"<<endl; break;
        case IA_ERROR: cout<<"Showing intersection-angle error"<<endl; break;
        case FACE_CONCYCLITY: cout<<"Showing face concyclity error"<<endl; break;
        case QC_ERROR: cout<<"Showing quasiconformal error"<<endl; break;
      }
      break;
    }
    case '3': {
      showDeformed=!showDeformed;
      UpdateCurrentView();
      break;
    }
    case '5': {
      switch(viewingMode){
        case DC_ERROR: DCSensitivity+=0.01 ;cout<<"discrete conformal error range: [0,"<<DCSensitivity<<"]"<<endl; break;
        case IA_ERROR: IASensitivity+=1.0 ;cout<<"intersection-angle error range: [0,"<<IASensitivity<<"]"<<endl; break;
        case FACE_CONCYCLITY: circSensitivity+=1.0 ;cout<<"concyclity error range: [0,"<<circSensitivity<<"]"<<endl; break;
        case QC_ERROR: QCSensitivity+=0.05 ;cout<<"quasiconformal error range: [0,"<<QCSensitivity<<"]"<<endl; break;
      }
      UpdateCurrentView();
      break;
    }
    case '4': {
      switch(viewingMode){
        case DC_ERROR:  if (DCSensitivity>0.01) DCSensitivity-=0.01; cout<<"discrete conformal error range: [0,"<<DCSensitivity<<"]"<<endl; break;
        case IA_ERROR: if (IASensitivity>1.0)  IASensitivity-=1.0 ;cout<<"intersection-angle error range: [0,"<<IASensitivity<<"]"<<endl; break;
        case FACE_CONCYCLITY:  if (circSensitivity>1.0) circSensitivity-=1.0 ;cout<<"concyclity error range: [0,"<<circSensitivity<<"]"<<endl; break;
        case QC_ERROR: if (QCSensitivity>0.05) QCSensitivity-=0.05 ;cout<<"quasiconformal error range: [0,"<<QCSensitivity<<"]"<<endl; break;
      }
      UpdateCurrentView();
      break;
    }
      
    case '6':{
      isExactMC=!isExactMC;
      Eigen::MatrixXd bc(handles.size(),currV.cols());
      for (int i=0;i<handles.size();i++)
        bc.row(i)=handlePoses[i].transpose();
      
      hedra::quat_moebius_deform(qmdata, 1.0, rigidityFactor, isExactMC, bc,true, deformV);
      hedra::edge_mesh(deformV,D,F,EV, EF, edgeDeformV, edgeT,edgeTE);
      UpdateCurrentView();
    }
      
  }
  
  return false;
}


int main(int argc, char *argv[])
{
  
  cout<<"press 1+right button to select new handles."<<endl;
  cout<<"press 2 button to toggle between different quality measures."<<endl;
  cout<<"press 3 to toggle between original and deformed mesh."<<endl;
  cout<<"press 4-5 to change the sensitivity of the quality measurement."<<endl;
   cout<<"press 6 to toggle exact Discrete Conformality."<<endl;
  cout<<"press the right button and drag the current handle for deformation."<<endl;
  
  // Load a mesh in OFF format
  hedra::polygonal_read_OFF(TUTORIAL_SHARED_PATH "/moomoo.off", origV, D, F);
  hedra::polygonal_edge_topology(D, F, EV, FE, EF,EFi,FEs,innerEdges);
  hedra::triangulate_mesh(D, F, polyT, polyTF);
  hedra::quat_moebius_setup(origV,D,F,polyTF,EV,EF,EFi,FE,FEs, innerEdges, qmdata);
  hedra::edge_mesh(origV, D, F, EV, EF, edgeOrigV, edgeT, edgeTE);
  
  currV=origV;
  deformV=origV;
  edgeDeformV=edgeOrigV;
  spans=currV.colwise().maxCoeff()-currV.colwise().minCoeff();
  
  viewer.callback_mouse_down = &mouse_down;
  viewer.callback_mouse_move = &mouse_move;
  viewer.callback_mouse_up=&mouse_up;
  viewer.callback_key_down=&key_down;
  viewer.callback_key_up=&key_up;
  viewer.core.background_color<<0.75,0.75,0.75,1.0;
  
  viewer.append_mesh();  //handles mesh
  viewer.selected_data_index=0;
  UpdateCurrentView();
  viewer.launch();
  
  
}
