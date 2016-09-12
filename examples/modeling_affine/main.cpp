#include <hedra/hedra_read_OFF.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/hedra_edge_topology.h>
#include <hedra/affine_maps_deform.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/viewer/Viewer.h>
#include <igl/readDMAT.h>


bool LeftButtonDown=false;
std::vector<int> Handles;
std::vector<Eigen::RowVector3d> HandlePoses;
int CurrentHandle;

const Eigen::RowVector3d sea_green(70./255.,252./255.,167./255.);
Eigen::MatrixXd V,U;
Eigen::MatrixXi F, T;
Eigen::VectorXi S,b, D;
Eigen::MatrixXi EV, EF, FE;
Eigen::RowVector3d mid;
double anim_t = 0.0;
double anim_t_dir = 0.03;
hedra::AffineData affine_data;


bool UpdateCurrentView(igl::viewer::Viewer & viewer)
{
    using namespace Eigen;
    using namespace std;
    MatrixXd bc(b.size(),V.cols());
    for(int i = 0;i<b.size();i++)
    {
        bc.row(i) = V.row(b(i));
        switch(S(b(i)))
        {
            case 0:
            {
                const double r = mid(0)*0.25;
                bc(i,0) += r*sin(0.5*anim_t*2.*igl::PI);
                bc(i,1) -= r+r*cos(igl::PI+0.5*anim_t*2.*igl::PI);
                break;
            }
            case 1:
            {
                const double r = mid(1)*0.15;
                bc(i,1) += r+r*cos(igl::PI+0.15*anim_t*2.*igl::PI);
                bc(i,2) -= r*sin(0.15*anim_t*2.*igl::PI);
                break;
            }
            case 2:
            {
                const double r = mid(1)*0.15;
                bc(i,2) += r+r*cos(igl::PI+0.35*anim_t*2.*igl::PI);
                bc(i,0) += r*sin(0.35*anim_t*2.*igl::PI);
                break;
            }
            default:
                break;
        }
    }
    Eigen::MatrixXd A;
    hedra::affine_maps_deform(affine_data,bc, U,U,A);
    viewer.data.set_vertices(U);
    viewer.data.compute_normals();
    if(viewer.core.is_animating)
    {
        anim_t += anim_t_dir;
    }
    return false;
}

bool key_down(igl::viewer::Viewer &viewer, unsigned char key, int mods)
{
    switch(key)
    {
        case ' ':
            viewer.core.is_animating = !viewer.core.is_animating;
            return true;
    }
    return false;
}


int main(int argc, char *argv[])
{
    
    // Load a mesh in OFF format
    using namespace std;
    using namespace Eigen;
    
    hedra::hedra_read_OFF(DATA_PATH "/six.off", V, D, F);
    igl::readDMAT(DATA_PATH "/six-selection.dmat",S);
    
    // vertices in selection
    igl::colon<int>(0,V.rows()-1,b);
    b.conservativeResize(stable_partition( b.data(), b.data()+b.size(),
                                          [](int i)->bool{return S(i)>=0;})-b.data());
    // Centroid
    mid = 0.5*(V.colwise().maxCoeff() + V.colwise().minCoeff());
    // Precomputation
    hedra::hedra_edge_topology(D, F, EV, FE, EF);
    hedra::affine_maps_precompute(V,D,F,EF,EV, b,affine_data);
    
    // Set color based on selection
    //TODO: replace with spheres and planarity plotting
    MatrixXd C(F.rows(),3);
    RowVector3d purple(80.0/255.0,64.0/255.0,255.0/255.0);
    RowVector3d gold(255.0/255.0,228.0/255.0,58.0/255.0);
    for(int f = 0;f<F.rows();f++)
    {
        if( S(F(f,0))>=0 && S(F(f,1))>=0 && S(F(f,2))>=0)
        {
            C.row(f) = purple;
        }else
        {
            C.row(f) = gold;
        }
    }
    
    // Plot the mesh with pseudocolors
    igl::viewer::Viewer viewer;
    viewer.data.set_mesh(U, F);
    viewer.data.set_colors(C);
    viewer.callback_pre_draw = &UpdateCurrentView;
    viewer.callback_key_down = &key_down;
    viewer.core.is_animating = false;
    viewer.core.animation_max_fps = 30.;
    cout<<
    "Press [space] to toggle animation"<<endl;
    viewer.launch();
}
