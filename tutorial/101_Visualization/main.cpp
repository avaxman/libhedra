#include <igl/opengl/glfw/Viewer.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/polyhedral_face_normals.h>
#include <hedra/polygonal_face_centers.h>
#include <hedra/point_spheres.h>
#include <hedra/polygonal_edge_lines.h>

Eigen::MatrixXd VMesh, VSpheres, VEdges, VNormals, FEs;
Eigen::MatrixXi F,EV, EF, FE,TMesh, TSpheres, TEdges, TNormals, EFi;
Eigen::VectorXi D,TF, innerEdges;
Eigen::MatrixXd faceCenters, faceNormals, normalColors, sphereColors;
Eigen::MatrixXd CNormals, CEdges, CMesh, CSpheres;

bool ShowPolygonalEdges=true;
bool ShowVertexSpheres=false;
bool ShowFaceNormals=false;

void update_meshes(igl::opengl::glfw::Viewer& viewer)
{
  
  if (ShowPolygonalEdges){
    std::cout<<"Showing polygonal lines"<<std::endl;
    viewer.data_list[1].show_faces=true;
  } else viewer.data_list[1].show_faces=false;
 
  if (ShowFaceNormals){
    viewer.data_list[2].show_faces=true;
    std::cout<<"Showing face normals"<<std::endl;
  }  else viewer.data_list[2].show_faces=false;
  
  if (ShowVertexSpheres){
    viewer.data_list[3].show_faces=true;
    std::cout<<"Showing vertex spheres"<<std::endl;
  } else viewer.data_list[3].show_faces=false;
  
}

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  switch(key){
    case '1': ShowPolygonalEdges=!ShowPolygonalEdges; break;
    case '2': ShowVertexSpheres=!ShowVertexSpheres; break;
    case '3': ShowFaceNormals=!ShowFaceNormals; break;
  }
  update_meshes(viewer);
  return false;
}


int main(int argc, char *argv[])
{
  using namespace std;
  
  hedra::polygonal_read_OFF(argv[1], VDOG, DDOG, FDOG);
  libigl::read_obj_(argv[2], VWedge, DWedge, F);
  

  sphereColors.resize(VMesh.rows(),3);
  for (int i=0;i<VMesh.rows();i++)
    sphereColors.row(i)<<(double)i/(double)VMesh.rows(), 1.0-(double)i/(double)VMesh.rows(), 0.0;
  
  hedra::triangulate_mesh(D, F,TMesh,TF);
  hedra::polygonal_edge_topology(D,F, EV,FE, EF, EFi, FEs, innerEdges);
  hedra::polygonal_face_centers(VMesh,D, F,faceCenters);
  hedra::polyhedral_face_normals(VMesh,D, F,faceNormals);
  
  normalColors.resize(F.rows(),3);
  normalColors.col(0)=Eigen::VectorXd::Constant(F.rows(),0.5);
  normalColors.col(1)=Eigen::VectorXd::Constant(F.rows(),0.5);
  normalColors.col(2)=Eigen::VectorXd::Constant(F.rows(),1.0);
  
  igl::opengl::glfw::Viewer viewer;
  
  //Polygonal triangulated mesh
  viewer.data_list[0].clear();
  viewer.data_list[0].set_mesh(VMesh,TMesh);
  viewer.data_list[0].set_colors(hedra::default_mesh_color());
  viewer.data_list[0].set_face_based(true);
  viewer.data_list[0].show_lines=false;
  
  //edges mesh
  viewer.append_mesh();
  hedra::polygonal_edge_lines(VMesh, F, TMesh, EV,  VEdges, TEdges, CEdges,0.25);
  viewer.data_list[1].set_mesh(VEdges, TEdges);
  viewer.data_list[1].set_colors(CEdges);
  viewer.data_list[1].show_lines=false;
  
  //spheres mesh
  viewer.append_mesh();
  hedra::point_spheres(VMesh, 0.02, sphereColors, 10, VSpheres, TSpheres, CSpheres);
  viewer.data_list[2].set_mesh(VSpheres, TSpheres);
  viewer.data_list[2].set_colors(CSpheres);
  viewer.data_list[2].show_lines=false;
  viewer.data_list[2].show_faces=false;
  
  //normals mesh
  viewer.append_mesh();
  hedra::line_cylinders(faceCenters, faceCenters+faceNormals, 0.02, normalColors, 10, VNormals, TNormals, CNormals);
  
  viewer.data_list[3].set_mesh(VNormals, TNormals);
  viewer.data_list[3].set_colors(CNormals);
  viewer.data_list[3].set_face_based(false);
  viewer.data_list[3].show_lines=false;
  viewer.data_list[3].show_faces=false;
  
  viewer.selected_data_index=0;

  viewer.callback_key_down = &key_down;
  update_meshes(viewer);
  viewer.launch();
  
  return 0;
}


