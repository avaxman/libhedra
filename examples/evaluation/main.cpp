#include <igl/viewer/Viewer.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/edge_mesh.h>
#include <hedra/planarity.h>
#include <hedra/concyclity.h>
#include <hedra/scalar2RGB.h>

enum ViewingMode{PLANARITY, CONCYCLITY, EDGE_LENGTHS} ViewingMode;

Eigen::MatrixXd V, edgeV;
Eigen::MatrixXi F, T, midedgeT;
Eigen::VectorXi D, TF, edgeTE;
Eigen::MatrixXi EV, EF, FE, EFi, edgeT;
Eigen::MatrixXd FEs;
Eigen::VectorXi innerEdges;
Eigen::Vector3d spans;
Eigen::VectorXd planarity;
Eigen::VectorXd concyclity;
Eigen::VectorXd edgeLengths;
double minEdgeLength, maxEdgeLength;


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
        else
            hedra::scalar2RGB(concyclity, 0.0,5.0, C);
        bigTC.resize(T.rows(),3);
        for (int i=0;i<T.rows();i++)
            bigTC.row(i)=C.row(TF(i));
    }
    
    if (ViewingMode==EDGE_LENGTHS){
        bigV=edgeV;
        bigT=edgeT;
        Eigen::MatrixXd C;
        hedra::scalar2RGB(edgeLengths, minEdgeLength,maxEdgeLength, C);
        bigTC.resize(edgeT.rows(),3);
        for (int i=0;i<edgeT.rows();i++)
            bigTC.row(i)=C.row(edgeTE(i));
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
        case '3': ViewingMode=EDGE_LENGTHS; break;
    }
    UpdateCurrentView(viewer);
    return false;
}


int main(int argc, char *argv[])
{
    
    // Load a mesh in OFF format
    using namespace std;
    using namespace Eigen;
    
    std::cout<<R"(
    1 Show planarity between [0,1]
    2 Show concyclity between [0,5]
    3 Show normalized edge lengths
    )";
    
    hedra::polygonal_read_OFF(DATA_PATH "/intersection.off", V, D, F);
    hedra::polygonal_edge_topology(D, F, EV, FE, EF,EFi,FEs,innerEdges);
    hedra::triangulate_mesh(D, F, T, TF);
    hedra::edge_mesh(V,D,F,EV, EF, edgeV, edgeT,edgeTE);
    spans=V.colwise().maxCoeff()-V.colwise().minCoeff();
    
    hedra::planarity(V,D,F,planarity);
    hedra::concyclity(V,D,F,concyclity);
    
    edgeLengths.resize(EV.rows());
    for (int i=0;i<EV.rows();i++)
        edgeLengths(i)=(V.row(EV(i,0))-V.row(EV(i,1))).norm();
    
    minEdgeLength=edgeLengths.minCoeff();
    maxEdgeLength=edgeLengths.maxCoeff();
    
    ViewingMode=PLANARITY;
    igl::viewer::Viewer viewer;
    viewer.callback_key_down=&key_down;
    viewer.core.background_color<<0.75,0.75,0.75,1.0;
    UpdateCurrentView(viewer);
    viewer.launch();
    

}
