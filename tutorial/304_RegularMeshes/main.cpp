#include <algorithm>
#include <math.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject.h>
#include <igl/project.h>
#include <igl/avg_edge_length.h>
#include <igl/boundary_loop.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/polygonal_write_OFF.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/moebius_regular_meshes.h>
#include <hedra/visualization_schemes.h>
#include <hedra/point_spheres.h>
#include <hedra/scalar2RGB.h>


typedef enum {ORIGINAL_MESH, REGULAR_MESH} MeshModes;
MeshModes meshMode=ORIGINAL_MESH;

typedef enum {STANDARD, MOEBIUS_REGULARITY, EUCLIDEAN_REGULARITY, WILLMORE_ENERGY} ViewingModes;
ViewingModes viewingMode=STANDARD;

// Mesh
Eigen::MatrixXd VOrig, VRegular, VHandles, FEs, CHandles;
Eigen::MatrixXi F, T, FHandles;
Eigen::VectorXi D, innerEdges, TF;
Eigen::MatrixXi EV, FE, EF, EFi;

Eigen::VectorXi constIndices;
Eigen::MatrixXd constPoses;

double MRCoeff=1.0;
double ERCoeff=0.1;

hedra::MoebiusRegularData MRData;


void update_mesh(igl::opengl::glfw::Viewer& viewer)
{
  
  double maxMR = MRData.origMR.maxCoeff()/2.5;
  double maxER = MRData.origER.maxCoeff()/2.5;
  double maxW = MRData.origW.maxCoeff()/2.5;
  switch(meshMode){
    case ORIGINAL_MESH: std::cout<<"Original Mesh"<<std::endl; break;
    case REGULAR_MESH: std::cout<<"Regular Mesh"<<std::endl; break;
  }
  switch(viewingMode){
    case STANDARD: std::cout<<"Standard coloring"<<std::endl; break;
    case MOEBIUS_REGULARITY: std::cout<<"Moebius Regularity between [0,"<<maxMR<<"]"<<std::endl; break;
    case EUCLIDEAN_REGULARITY: std::cout<<"Euclidean Regularity between [0,"<<maxER<<"]"<<std::endl; break;
    case WILLMORE_ENERGY: std::cout<<"Willmore Energy between [0,"<<maxW<<"]"<<std::endl; break;
  }
  
  //mesh
  viewer.data_list[0].clear();
  viewer.data_list[0].set_face_based(true);
  viewer.data_list[0].show_lines=false;
  
  Eigen::MatrixXd C;  //color
  switch (viewingMode){
    case STANDARD: C = hedra::default_mesh_color(); break;
    case MOEBIUS_REGULARITY: hedra::scalar2RGB((meshMode==ORIGINAL_MESH ? MRData.origMR : MRData.deformMR), 0,maxMR,C); break;
    case EUCLIDEAN_REGULARITY: hedra::scalar2RGB((meshMode==ORIGINAL_MESH ? MRData.origER : MRData.deformER), 0,maxER, C); break;
    case WILLMORE_ENERGY: hedra::scalar2RGB((meshMode==ORIGINAL_MESH ? MRData.origW : MRData.deformW), 0,maxW, C); break;
  }
  
  viewer.data_list[0].set_mesh((meshMode==ORIGINAL_MESH ? VOrig : VRegular), T);
  if (viewingMode==EUCLIDEAN_REGULARITY){
    Eigen::MatrixXd TC(T.rows(),3);
    for (int i=0;i<T.rows();i++)
      TC.row(i)=C.row(TF(i));
    viewer.data_list[0].set_colors(TC);
  } else
    viewer.data_list[0].set_colors(C);
  
  viewer.data_list[0].set_edges((meshMode==ORIGINAL_MESH ? VOrig : VRegular),EV,hedra::default_edge_color());
  
  double radius = 0.25*igl::avg_edge_length(VOrig, T);
  hedra::point_spheres(constPoses, radius, (hedra::passive_handle_color()).replicate(F.rows(),1), 8, VHandles, FHandles, CHandles);
  viewer.data_list[1].clear();
  viewer.data_list[1].set_mesh(VHandles, FHandles);
  viewer.data_list[1].set_face_based(true);
  viewer.data_list[1].set_colors(CHandles);
  viewer.data_list[1].show_faces = true;
  viewer.data_list[1].show_lines = false;
}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  switch(key)
  {
    default:
      return false;
    case '1': {meshMode=(MeshModes)((meshMode+1)%2); break;}
    case '2': {viewingMode=(ViewingModes)((viewingMode+1)%4); break;}
    case '3': {
      ERCoeff+=0.1;
      hedra::compute_moebius_regular(MRData, MRCoeff, ERCoeff, constPoses, false, VRegular);
      std::cout<<"Euclidean Regularity Coefficient: "<<ERCoeff<<std::endl;
      break;
      
    }
    case '4': {
      ERCoeff=(ERCoeff-0.1>=0.0 ? ERCoeff-0.1 : 0.0);
      hedra::compute_moebius_regular(MRData, MRCoeff, ERCoeff, constPoses, false, VRegular);
      std::cout<<"Euclidean Regularity Coefficient: "<<ERCoeff<<std::endl;
      break;}
      
  }
  
  update_mesh(viewer);
  return true;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  
  cout<<"1  Show original/Moebius regular mesh"<<endl<<
  "2  Toggle measurements "<<endl<<
  "3  Increase Euclidean Regularity "<<endl<<
  "4  Decrease Euclidean Regularity "<<endl;
  hedra::polygonal_read_OFF(TUTORIAL_SHARED_PATH "/intersection.off", VOrig, D, F);
  hedra::polygonal_edge_topology(D, F, EV, FE, EF, EFi, FEs, innerEdges);
  hedra::triangulate_mesh(D,F,T, TF);
  
  //TODO: code to choose entire boundary
  vector<vector<int> > boundaryList;
  igl::boundary_loop(T, boundaryList);
  
  VectorXi boundaryMask=VectorXi::Zero(VOrig.rows());
  for (int i=0;i<boundaryList.size();i++)
    for (int j=0;j<boundaryList[i].size();j++)
      boundaryMask(boundaryList[i][j])=1;
  
  constIndices.resize(boundaryMask.sum());
  constPoses.resize(boundaryMask.sum(),3);
  int counter=0;
  for (int i=0;i<boundaryMask.size();i++){
    if (boundaryMask(i)){
      constIndices(counter)=i;
      constPoses.row(counter++)=VOrig.row(i);
    }
  }
  
  
  hedra::setup_moebius_regular(VOrig, D, F, T, EV, FE, EF, EFi, FEs, innerEdges, constIndices, MRData);
  hedra::compute_moebius_regular(MRData,  MRCoeff, ERCoeff, constPoses, false, VRegular);
  
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;

  viewer.core.background_color<<0.75,0.75,0.75,1.0;
  
  viewer.append_mesh();
  viewer.selected_data_index=0;
  update_mesh(viewer);
  viewer.launch();
  
}
