#include <igl/unproject_onto_mesh.h>
#include <igl/viewer/Viewer.h>
#include <igl/readDMAT.h>
#include <igl/jet.h>
#include <igl/per_vertex_normals.h>
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
#include <hedra/complex_moebius_deform.h>
#include <algorithm>

using namespace Eigen;
using namespace std;


igl::viewer::Viewer viewer;
typedef enum {MC_ERROR, IA_ERROR, FACE_CONCYCLITY, QC_ERROR} ViewingModes;
ViewingModes viewingMode=FACE_CONCYCLITY;

bool showDeformed=true;
std::vector<int> handles;
std::vector<Eigen::RowVector3d> handlePoses;
int currHandle;

Eigen::MatrixXd origV, deformV, currV;
Eigen::MatrixXd edgeOrigV, edgeDeformV;
Eigen::MatrixXi F, edgeT, currT, polyT;
Eigen::VectorXi D, polyTF,edgeTE, currTF;
Eigen::MatrixXi EV, EF, FE, EFi, I4;
Eigen::MatrixXd FEs;
Eigen::VectorXi innerEdges;
Eigen::Vector3d spans;
bool isExactMC=false;
bool isExactIAP=false;

bool editing=false;
bool choosingHandleMode=false;
double currWinZ;
hedra::ComplexMoebiusData cmdata;

double rigidityFactor=0.1;      //inversion control ratio
double MCSensitivity=0.01;     //ratio of abs(cr)
double IASensitivity=1.0;      //circle IA degrees difference between faces
double QCSensitivity=0.5;      //abs(max_sing/min_sing)
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
            
        case MC_ERROR:
        case IA_ERROR:
            currV=(showDeformed ? edgeDeformV : edgeOrigV);
            currT=edgeT;
            currTF=edgeTE;
            break;
    }
}

bool UpdateCurrentView()
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
        case MC_ERROR:
            hedra::dc_error(origV, deformV, D, F, EV, EF, EFi, innerEdges, scalarMeasure);
            scalarMeasure/=MCSensitivity;
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
    
    double sphereRadius=spans.sum()/200.0;
    MatrixXd sphereGreens(handles.size(),3);
    sphereGreens.col(0).setZero();
    sphereGreens.col(1).setOnes();
    sphereGreens.col(2).setZero();
    
    MatrixXd bigV=currV;
    MatrixXi bigT=currT;
    MatrixXd bigTC=TC;
    
    hedra::point_spheres(bc, sphereRadius, sphereGreens, 10, false, true, bigV, bigT, bigTC);
    
    viewer.core.show_lines=false;
    Eigen::MatrixXd OrigEdgeColors(EV.rows(),3);
    OrigEdgeColors.col(0)=Eigen::VectorXd::Constant(EV.rows(),0.0);
    OrigEdgeColors.col(1)=Eigen::VectorXd::Constant(EV.rows(),0.0);
    OrigEdgeColors.col(2)=Eigen::VectorXd::Constant(EV.rows(),0.0);
    
    viewer.data.clear();
    viewer.data.set_mesh(bigV,bigT);
    viewer.data.set_colors(bigTC);
    viewer.data.compute_normals();
    viewer.data.set_edges(bigV,EV,OrigEdgeColors);
    return true;
}

bool mouse_move(igl::viewer::Viewer& viewer, int mouse_x, int mouse_y)
{
    if (!editing)
        return false;
    
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    Eigen::RowVector3f NewPos=igl::unproject<float>(Eigen::Vector3f(x,y,currWinZ),
                                                    viewer.core.view * viewer.core.model,
                                                    viewer.core.proj,
                                                    viewer.core.viewport);
    
    handlePoses[handlePoses.size()-1]=NewPos.cast<double>();
    
    Eigen::MatrixXd A(3*F.rows(),3);
    for (int i=0;i<F.rows();i++)
        A.block(3*i,0,3,3)=Eigen::Matrix3d::Identity();
    
    Eigen::MatrixXd bc(handles.size(),currV.cols());
    for (int i=0;i<handles.size();i++)
        bc.row(i)=handlePoses[i].transpose();
    
    hedra::complex_moebius_deform(cmdata, bc,150, deformV);
    hedra::edge_mesh(deformV,D,F,EV, EF, edgeDeformV, edgeT,edgeTE);
    UpdateCurrentView();
    return true;
    
}


bool mouse_up(igl::viewer::Viewer& viewer, int button, int modifier)
{
    if (((igl::viewer::Viewer::MouseButton)button==igl::viewer::Viewer::MouseButton::Left))
        return false;
    
    editing=false;
    
    return true;
}

bool mouse_down(igl::viewer::Viewer& viewer, int button, int modifier)
{
    if (((igl::viewer::Viewer::MouseButton)button==igl::viewer::Viewer::MouseButton::Left))
        return false;
    int vid, fid;
    Eigen::Vector3f bc;
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    if (!choosingHandleMode){
        editing=true;
        return false;
    }
    if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core.view * viewer.core.model,
                                viewer.core.proj, viewer.core.viewport, currV, F, fid, bc))
    {
        //add the closest vertex to the handles
        Eigen::MatrixXf::Index maxRow, maxCol;
        bc.maxCoeff(&maxRow);
        int CurrVertex=F(fid, maxRow);
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
                                                      viewer.core.view * viewer.core.model,
                                                      viewer.core.proj,
                                                      viewer.core.viewport);
        
        currWinZ=WinCoords(2);
        std::cout<<"Choosing Vertex :"<<CurrVertex<<std::endl;
        
        Eigen::VectorXi b(handles.size());
        for (int i=0;i<handles.size();i++)
            b(i)=handles[i];
        
        hedra::complex_moebius_precompute(b, isExactMC, isExactIAP, rigidityFactor, cmdata);
        //hedra::affine_maps_precompute(V,D,F,EV,EF,EFi, FE, b, 3, affine_data);
        
        UpdateCurrentView();
    }
    return true;
}

bool key_up(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
    switch(key)
    {
            
        case '1': choosingHandleMode=false;
            break;
    }
    return false;
}

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
    switch(key)
    {
        case '1': choosingHandleMode=true;
            break;
            
        case '2': viewingMode=(ViewingModes)(((int)viewingMode+1)%4);
            UpdateCurrentView();
            break;
      case '3': showDeformed=!showDeformed;
        UpdateCurrentView();
        break;
        
    }
    return false;
}


int main(int argc, char *argv[])
{
    
    cout<<"press 1+right button to select new handles"<<endl;
    cout<<"press 2 button to toggle between different quality measures"<<endl;
    cout<<"press 3 to toggle between original and deformed mesh"<<endl;
    cout<<"press the right button and drag the current handle for deformation"<<endl;
    
    // Load a mesh in OFF format
    hedra::polygonal_read_OFF(TUTORIAL_SHARED_PATH "/bar2d.off", origV, D, F);
    hedra::polygonal_edge_topology(D, F, EV, FE, EF,EFi,FEs,innerEdges);
    hedra::triangulate_mesh(D, F, polyT, polyTF);
    hedra::complex_moebius_setup(origV,D,F,polyTF,EV,EF,EFi,FE,FEs, innerEdges, cmdata);
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
    UpdateCurrentView();
    viewer.launch();
    
    
}
