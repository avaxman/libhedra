#include <igl/viewer/Viewer.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/planarity.h>
#include <hedra/concyclity.h>
#include <hedra/scalar2RGB.h>

enum ViewingMode{PLANARITY, CONCYCLITY, DIHEDRAL} ViewingMode;

Eigen::MatrixXd V, midedgeV;
Eigen::MatrixXi F, T, midedgeT;
Eigen::VectorXi D, TF;
Eigen::MatrixXi EV, EF, FE, EFi;
Eigen::MatrixXd FEs;
Eigen::VectorXi innerEdges;
Eigen::Vector3d spans;
Eigen::VectorXd planarity;
Eigen::VectorXd concyclity;
Eigen::VectorXd dihedrals;


bool UpdateCurrentView(igl::viewer::Viewer & viewer)
{
    using namespace Eigen;
    using namespace std;
    
    Eigen::MatrixXd bigV;
    Eigen::MatrixXi bigT;
    Eigen::MatrixXd bigTC;
    
    if ((ViewingMode==PLANARITY)||(ViewingMode==CONCYCLITY)){
        bigV=V;
        bigT=T;
        Eigen::MatrixXd C;
        if (ViewingMode==PLANARITY)
            hedra::scalar2RGB(planarity, 0.0,1.0, C);
        else{
            cout<<"concyclity"<<concyclity<<endl;
            hedra::scalar2RGB(concyclity, 0.0,5.0, C);
            cout<<"C.rows() vs F.rows()"<<C.rows()<<","<<F.rows()<<endl;
        }
        bigTC.resize(T.rows(),3);
        for (int i=0;i<T.rows();i++)
            bigTC.row(i)=C.row(TF(i));
    }
    
    viewer.core.show_lines=false;
    Eigen::MatrixXd OrigEdgeColors(EV.rows(),3);
    OrigEdgeColors.col(0)=Eigen::VectorXd::Constant(EV.rows(),0.0);
    OrigEdgeColors.col(1)=Eigen::VectorXd::Constant(EV.rows(),0.0);
    OrigEdgeColors.col(2)=Eigen::VectorXd::Constant(EV.rows(),0.0);
    
    viewer.data.clear();
    viewer.data.set_mesh(bigV,bigT);
    viewer.data.set_colors(bigTC);
    viewer.data.compute_normals();
    viewer.data.set_edges(V,EV,OrigEdgeColors);
    return true;
}



bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
    switch(key)
    {
        case '1': ViewingMode=PLANARITY; break;
        case '2': ViewingMode=CONCYCLITY; break;
        case '3': ViewingMode=DIHEDRAL; break;
    }
    UpdateCurrentView(viewer);
    return false;
}


int main(int argc, char *argv[])
{
    
    // Load a mesh in OFF format
    using namespace std;
    using namespace Eigen;
    
    hedra::polygonal_read_OFF(DATA_PATH "/intersection.off", V, D, F);
    hedra::polygonal_edge_topology(D, F, EV, FE, EF,EFi,FEs,innerEdges);
    hedra::triangulate_mesh(D, F, T, TF);
    spans=V.colwise().maxCoeff()-V.colwise().minCoeff();
    
    hedra::planarity(V,D,F,planarity);
    hedra::concyclity(V,D,F,concyclity);
    
    ViewingMode=PLANARITY;
    igl::viewer::Viewer viewer;
    viewer.callback_key_down=&key_down;
    viewer.core.background_color<<0.75,0.75,0.75,1.0;
    UpdateCurrentView(viewer);
    viewer.launch();
    

}
