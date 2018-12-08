#include <igl/opengl/glfw/Viewer.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/edge_mesh.h>
#include <hedra/planarity.h>
#include <hedra/concyclity.h>
#include <hedra/scalar2RGB.h>
#include <hedra/polygonal_edge_lines.h>

enum ViewingMode{PLANARITY, CONCYCLITY, EDGE_LENGTHS} ViewingMode;

Eigen::MatrixXd VMesh, VEdgeMesh, CMesh;
Eigen::MatrixXd VEdges, CEdges;
Eigen::MatrixXi F, TMesh, TEdges, TEdgeMesh;
Eigen::VectorXi D, TF, TEEdgeMesh;
Eigen::MatrixXi EV, EF, FE, EFi;
Eigen::MatrixXd FEs;
Eigen::VectorXi innerEdges;
Eigen::VectorXd planarity;
Eigen::VectorXd concyclity;
Eigen::VectorXd edgeLengths;
double minEdgeLength, maxEdgeLength;


bool update_meshes(igl::opengl::glfw::Viewer & viewer)
{
  using namespace Eigen;
  using namespace std;
  
  if ((ViewingMode==PLANARITY)||(ViewingMode==CONCYCLITY)){
    Eigen::MatrixXd C;
    if (ViewingMode==PLANARITY)
      hedra::scalar2RGB(planarity, 0.0,1.0, C);
    else
      hedra::scalar2RGB(concyclity, 0.0,5.0, C);
    
    CMesh.conservativeResize(TMesh.rows(),3);
    for (int i=0;i<TMesh.rows();i++)
      CMesh.row(i)=C.row(TF(i));
    
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VMesh, TMesh);
    viewer.data_list[0].set_colors(CMesh);
    viewer.data_list[0].show_lines=false;
    
  }
  
  if (ViewingMode==EDGE_LENGTHS){
    Eigen::MatrixXd C;
    hedra::scalar2RGB(edgeLengths, minEdgeLength,maxEdgeLength, C);
    CMesh.conservativeResize(TEdgeMesh.rows(),3);
    for (int i=0;i<TEdgeMesh.rows();i++)
      CMesh.row(i)=C.row(TEEdgeMesh(i));
    
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VEdgeMesh, TEdgeMesh);
    viewer.data_list[0].set_colors(CMesh);
    viewer.data_list[0].show_lines=false;
    
  }
  
  return true;
}



bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  switch(key)
  {
    case '1': ViewingMode=PLANARITY; break;
    case '2': ViewingMode=CONCYCLITY; break;
    case '3': ViewingMode=EDGE_LENGTHS; break;
  }
  update_meshes(viewer);
  return false;
}


int main(int argc, char *argv[])
{
  
  // Load a mesh in OFF format
  using namespace std;
  using namespace Eigen;
  
  std::cout<<R"(
  1 Visualize planarity between [0,1]
  2 Visualize concyclity between [0,5]
  3 Visualize normalized edge lengths
  )";
  
  hedra::polygonal_read_OFF(TUTORIAL_SHARED_PATH "/intersection.off", VMesh, D, F);
  hedra::polygonal_edge_topology(D, F, EV, FE, EF,EFi,FEs,innerEdges);
  hedra::triangulate_mesh(D, F, TMesh, TF);
  hedra::edge_mesh(VMesh,D,F,EV, EF, VEdgeMesh, TEdgeMesh,TEEdgeMesh);
  
  hedra::planarity(VMesh,D,F,planarity);
  hedra::concyclity(VMesh,D,F,concyclity);
  
  edgeLengths.resize(EV.rows());
  for (int i=0;i<EV.rows();i++)
    edgeLengths(i)=(VMesh.row(EV(i,0))-VMesh.row(EV(i,1))).norm();
  
  minEdgeLength=edgeLengths.minCoeff();
  maxEdgeLength=edgeLengths.maxCoeff();
  
  igl::opengl::glfw::Viewer viewer;
  
  //Polygonal triangulated mesh
  viewer.data_list[0].clear();
  viewer.data_list[0].set_mesh(VMesh,TMesh);
  viewer.data_list[0].set_face_based(true);
  viewer.data_list[0].show_lines=false;
  viewer.data_list[0].invert_normals=true;
  
  //edges mesh
  viewer.append_mesh();
  hedra::polygonal_edge_lines(VMesh, F, TMesh, EV,  VEdges, TEdges, CEdges,0.25);
  viewer.data_list[1].set_mesh(VEdges, TEdges);
  viewer.data_list[1].set_colors(CEdges);
  viewer.data_list[1].show_lines=false;
  
  viewer.selected_data_index=0;
  
  ViewingMode=PLANARITY;
  viewer.callback_key_down=&key_down;
  viewer.core.background_color<<0.75,0.75,0.75,1.0;
  update_meshes(viewer);
  viewer.launch();
  
}
