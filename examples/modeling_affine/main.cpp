#include <hedra/hedra_read_off.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/viewer/Viewer.h>
#include <iostream>

Eigen::MatrixXd V, C;
Eigen::MatrixXi F;

bool LeftButtonDown=false;
std::vector<int> Handles;
int CurrentHandle;


void UpdateCurrentView(igl::viewer::Viewer& viewer){
    
}

bool mouse_up(igl::viewer::Viewer& viewer, int button, int modifier)
{
    if (((igl::viewer::Viewer::MouseButton)button==igl::viewer::Viewer::MouseButton::Left))
        LeftButtonDown=false;
    
    return true;
}

bool mouse_move(igl::viewer::Viewer& viewer, int mouse_x, int mouse_y)
{
    if (!LeftButtonDown)  //not actively deforming
        return false;
    
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    Vector3f NewPos=igl::unproject<float>(Vector3f(x,y,CurrWinZ),
                                          viewer.core.view * viewer.core.model,
                                          viewer.core.proj,
                                          viewer.core.viewport);
    
    
    if (ReorientationMode){
        RowVector3d CurrCRVec=(NewPos.transpose().cast<double>()-mm.DeformV.row(ConstCRIndices[CurrActiveCRHandle])).normalized();
        ConstCR[CurrActiveCRHandle]<<CurrCRVec;
    } else if (DeformationMode){
        ConstPoses[CurrActivePosHandle]=NewPos.cast<double>();
        
    }
    
    mm.UpdateMinimal(ConstPoses, ConstCR ,1000, MinimalMode, VRFactor, FRFactor);
    RecomputeInterpolation=true;
    UpdateCurrentView();
    return true;
}


bool mouse_down(igl::viewer::Viewer& viewer, int button, int modifier)
{
    if ((igl::viewer::Viewer::MouseButton)button==igl::viewer::Viewer::MouseButton::Left)
        return false;  //left button is only for the default viewer/camera stuff
    
    if (!(EditingMode==EDIT_MINIMAL))
        return false;  //nothing to do in other modes
    
    RightButtonDown=true;
    
    if ((!ChoosingHandleMode)&&(!ChoosingCRMode))
        return false;
    
    //finding out which vertex
    int ClosestVertex=0;
    
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = (double)viewer.core.viewport(3) - (double)viewer.current_mouse_y;
    double MinDistance=3276700.0;
    V=mm.DeformV;
    for (int i=0;i<V.rows();i++){
        Vector3f WinCoords=igl::project<float>(V.row(i).cast<float>(),
                                               viewer.core.view * viewer.core.model,
                                               viewer.core.proj,
                                               viewer.core.viewport);
        double Distance=(WinCoords(0)-x)*(WinCoords(0)-x)+(WinCoords(1)-y)*(WinCoords(1)-y)+WinCoords(2)*WinCoords(2);
        
        if (Distance<MinDistance){
            ClosestVertex=i;
            MinDistance=Distance;
            CurrWinZ=WinCoords(2);
        }
    }
    
    if (ChoosingHandleMode){
        //checking if handle already exists
        bool Found=false;
        for (int i=0;i<ConstPosIndices.size();i++){
            if (ConstPosIndices[i]==ClosestVertex){
                Found=true;
                CurrActivePosHandle=i;
                break;
            }
        }
        
        if (!Found){
            ConstPosIndices.push_back(ClosestVertex);
            ConstPoses.push_back(V.row(ClosestVertex));
            CurrActivePosHandle=(int)(ConstPosIndices.size()-1);
        }
        cout << "Picked positional (vertex): "<<CurrActivePosHandle<< ")" << endl;
    }
    
    if (ChoosingCRMode){
        //checking if handle already exists
        bool Found=false;
        for (int i=0;i<ConstCRIndices.size();i++){
            if (ConstCRIndices[i]==ClosestVertex){
                Found=true;
                CurrActiveCRHandle=i;
                break;
            }
        }
        
        if (!Found){
            ConstCRIndices.push_back(ClosestVertex);
            ConstCR.push_back(mm.VCR.row(ClosestVertex));
            CurrActiveCRHandle=(int)(ConstCRIndices.size()-1);
        }
        
        cout << "Picked CR (vertex): "<<CurrActiveCRHandle<< ")" << endl;
    }
    
    //cout<<"ConstPosIndicesMat: "<<ConstPosIndicesMat<<endl;
    //cout<<"ConstCRIndicesMat: "<<ConstCRIndicesMat<<endl;
    
    //mm.InitWillmore(ConstPosIndicesMat, ConstCRIndicesMat, isExactMC, isExactPlan);
    
    UpdateCurrentView(viewer);
    return true;
}


int main(int argc, char *argv[])
{

    // Load a mesh in OFF format
    hedra::hedra_read_OFF(DATA_PATH "/six.off", V, D, F);

    igl::viewer::Viewer viewer;
    viewer.callback_mouse_down = mouse_down;
    viewer.callback_mouse_up = mouse_up;
    viewer.callback_mouse_move = mouse_move;
  
    // Show mesh
    UpdateCurrentView(viewer);
    viewer.launch();
}
