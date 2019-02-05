// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_COPYLEFT_CGAL_EXTRACT_MESH_H
#define HEDRA_COPYLEFT_CGAL_EXTRACT_MESH_H
#include <igl/igl_inline.h>

#include <Eigen/Dense>
#include<CGAL/basic.h>
#include<CGAL/Cartesian.h>
#include<CGAL/Polygon_2.h>
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/gmpq.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <vector>


#define EPSILON 0.000001
#define MAXdouble       ((Number ) 3.40282346638528860e+38)
#define M_INFINITY        ((Number ) MAXdouble)
#define PI 3.14159265

//#define MAX_THREADS 2
//#define MAX_THREADS 6
#define MAX_THREADS 1
#define UNKNOWN_COLOR 100



typedef enum {MV_METHOD,LBD_METHOD,BV_METHOD} MethodType;

//typedef CGAL::Filtered_kernel< CGAL::Cartesian<CORE::Expr> > CKernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Simple_cartesian<CGAL::Gmpq>  EKernel;

typedef std::vector<int> SelectionCurve;

typedef Kernel::RT Number;
typedef EKernel::RT ENumber;
typedef Kernel::Point_2 Point2D;
typedef Kernel::Point_3 Point3D;
typedef Kernel::Ray_2 Ray2D;
typedef Kernel::Line_2 Line2D;
typedef Kernel::Line_3 Line3D;
typedef Kernel::Vector_2 Vector2D;
typedef Kernel::Vector_3 Vector3D;
typedef Kernel::Segment_2 Segment2D;
typedef Kernel::Segment_3 Segment3D;
typedef Kernel::Plane_3 Plane3D;
typedef Kernel::Triangle_3 Triangle3D;
typedef Kernel::Triangle_2 Triangle2D;
typedef Kernel::Circle_2 Circle2D;
typedef EKernel::Point_2 EPoint2D;
typedef EKernel::Point_3 EPoint3D;
typedef EKernel::Triangle_2 ETriangle2D;
typedef Kernel::Ray_3 Ray3D;
typedef Kernel::Direction_2 Direction2D;
typedef Kernel::Direction_3 Direction3D;
typedef Kernel::Sphere_3 Sphere3D;
typedef Kernel::Aff_transformation_3 Transform3D;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron3D;
typedef CGAL::Polygon_2<Kernel>    Polygon2D;

typedef std::pair<int, Point3D> PosConstraint;  //positional constraint


namespace hedra
{
  namespace copyleft
  {
    namespace cgal
    {
      
      using namespace boost;
      using namespace CGAL;
      
      //#define wx (-0.5)
      //#define wy (SQRT3/2.0)
      
      void Normalize(Vector3D& v){
        v=v/sqrt(squared_distance(Point3D(0,0,0)+v,Point3D(0,0,0)));
      }
      
      inline Number Norm(const Vector3D &v){return sqrt(v.x()*v.x()+v.y()*v.y()+v.z()*v.z());}
      
      
      struct MergeData {
        const bool operator()(const int& v1, const int& v2) const {return v1; }
      };
      
      struct EdgeData{
        int ID;
        bool isHex;
        int OrigHalfedge;
        bool isBoundary;
        
        EdgeData():ID(-1), isHex(false), OrigHalfedge(-1), isBoundary(false){}
        ~EdgeData(){};
      };
      
      class HexSegment{
      public:
        std::pair<int,int> Source;
        std::pair<int,int> Target;
        
        HexSegment(int su, int sv, int tu, int tv){
          std::pair<int,int> s(su,sv);
          std::pair<int,int> t(tu,tv);
          if (s<t){
            Source=s;
            Target=t;
          }else{
            Source=t;
            Target=s;
          }
        }
        
        ~HexSegment(){}
        
        const bool operator<(const HexSegment& h) const {
          if (Source<h.Source) return false;
          if (Source>h.Source) return true;
          
          if (Target<h.Target) return false;
          if (Target>h.Target) return true;
          
          return false; //both are equal
        }
        
      };
      
      class Vertex{
      public:
        int ID;
        static int RefCount;
        Point3D Coordinates;
        int AdjHalfedge;
        
        int AssocTriangle;  //for modifier vertices
        
        bool isHex;
        bool Valid;
        bool isHanging;  //for modifier vertices
        
        Vertex():ID(-1), AdjHalfedge(-1), isHex(false), Valid(true), isHanging(false){RefCount++;}
        ~Vertex(){RefCount--;}
      };
      
      class Halfedge{
      public:
        int ID;
        static int RefCount;
        int Origin;
        int Next;
        int Prev;
        int Twin;
        int AdjFace;
        Point2D TexCoords;
        bool isHex;
        bool Valid;
        
        int OrigHalfedge;
        
        Halfedge():ID(-1), Origin(-1), Next(-1), Prev(-1), Twin(-1), AdjFace(-1), isHex(false), Valid(true), OrigHalfedge(-1){RefCount++;}
        ~Halfedge(){RefCount--;}
      };
      
      
      class Face{
      public:
        int ID;
        static int RefCount;
        int AdjHalfedge;
        Vector3D Normal;
        Vector3D Centroid;
        
        int NumVertices;
        int Vertices[60];
        bool Valid;
        
        Face():ID(-1), AdjHalfedge(-1), NumVertices(-1), Valid(true){RefCount++;}
        ~Face(){RefCount--;}
      };
      
      class Mesh{
      public:
        
        std::vector<Vertex> Vertices;
        std::vector<Halfedge> Halfedges;
        std::vector<Face> Faces;
        
        std::vector<int> TransVertices;
        std::vector<int> InStrip;
        std::vector<std::set<int> > VertexChains;
        
        std::ofstream DebugLog;
        
        void JoinFace(int heindex);
        void UnifyEdges(int heindex);
        bool CheckMesh();
        void CleanMesh();
        void ComputeTwins();
        void WalkBoundary(int &CurrEdge);
        void RemoveEdge(int heindex);
        void RemoveFace(int findex, int heindex);
        void TestUnmatchedTwins();
        
        void GenerateMesh(Mesh& HexMesh, double EdgeLength);
        
        void SimplifyHexMesh();
        
        void CreateStrips();
        
        void Clear(){Vertices.clear(); Halfedges.clear(); Faces.clear();}
        
        void Allocate(int NumofVertices, int NumofFaces, int NumofHEdges)
        {
          Vertices.resize(NumofVertices);
          Faces.resize(NumofFaces);
          Halfedges.resize(NumofHEdges);
        }
        
        void ComputeNormals()
        {
          for (int i=0;i<Faces.size();i++){
            Vector3D Vector1=Vertices[Faces[i].Vertices[1]].Coordinates-Vertices[Faces[i].Vertices[0]].Coordinates;
            Vector3D Vector2=Vertices[Faces[i].Vertices[2]].Coordinates-Vertices[Faces[i].Vertices[0]].Coordinates;
            Faces[i].Normal=CGAL::cross_product(Vector1,Vector2);
            Normalize(Faces[i].Normal);
          }
        }
        
        Mesh(){}
        ~Mesh(){}
        
      };
      
      
      template <class ArrangementA, class ArrangementB, class ArrangementR>
      class Arr_mesh_generation_overlay_traits :
      public _Arr_default_overlay_traits_base<ArrangementA,ArrangementB,ArrangementR>
      {
      public:
        
        typedef typename ArrangementA::Face_const_handle    Face_handle_A;
        typedef typename ArrangementB::Face_const_handle    Face_handle_B;
        typedef typename ArrangementR::Face_handle          Face_handle_R;
        
        typedef typename ArrangementA::Vertex_const_handle    Vertex_handle_A;
        typedef typename ArrangementB::Vertex_const_handle    Vertex_handle_B;
        typedef typename ArrangementR::Vertex_handle          Vertex_handle_R;
        
        typedef typename ArrangementA::Halfedge_const_handle    Halfedge_handle_A;
        typedef typename ArrangementB::Halfedge_const_handle    Halfedge_handle_B;
        typedef typename ArrangementR::Halfedge_handle          Halfedge_handle_R;
        
      public:
        
        virtual void create_face (Face_handle_A f1,
                                  Face_handle_B f2,
                                  Face_handle_R f) const
        {
          // Overlay the data objects associated with f1 and f2 and store the result
          // with f.
          f->set_data (f1->data());
          return;
        }
        
        //-1 - triangle vertex (non-hex vertex)
        //-2 - hex vertex
        virtual void  create_vertex ( Vertex_handle_A v1, Vertex_handle_B v2, Vertex_handle_R v)
        {
          v->set_data(-2);
        }
        
        
        virtual void create_vertex ( Vertex_handle_A v1, Halfedge_handle_B e2, Vertex_handle_R v)
        {
          v->set_data(-1);
        }
        
        
        virtual void create_vertex ( Vertex_handle_A v1, Face_handle_B f2, Vertex_handle_R v)
        {
          v->set_data(-1);
        }
        
        virtual void create_vertex ( Halfedge_handle_A e1, Vertex_handle_B v2, Vertex_handle_R v)
        {
          v->set_data(-2);
        }
        
        virtual void create_vertex ( Face_handle_A f1, Vertex_handle_B v2, Vertex_handle_R v)
        {
          v->set_data(-2);
        }
        
        virtual void create_vertex ( Halfedge_handle_A e1, Halfedge_handle_B e2, Vertex_handle_R v)
        {
          v->set_data(-1);
        }
        
        
        virtual void create_edge ( Halfedge_handle_A e1, Halfedge_handle_B e2, Halfedge_handle_R e)
        {
          EdgeData data;
          data.ID=-2;
          data.OrigHalfedge=e1->data().OrigHalfedge;
          e->set_data(data);
          e->twin()->set_data(data);
        }
        
        
        virtual void create_edge ( Halfedge_handle_A e1, Face_handle_B f2, Halfedge_handle_R e)
        {
          EdgeData data;
          data.ID=-1;
          data.OrigHalfedge=e1->data().OrigHalfedge;
          e->set_data(data);
          e->twin()->set_data(data);
        }
        
        virtual void create_edge ( Face_handle_A f1, Halfedge_handle_B e2, Halfedge_handle_R e)
        {
          EdgeData data;
          data.ID=-2;
          e->set_data(data);
          e->twin()->set_data(data);
        }
      };
      
      typedef Arr_segment_traits_2<EKernel>                 Traits2;
      typedef Traits2::Point_2                              Point2;
      typedef Traits2::X_monotone_curve_2                   Segment2;
      typedef Arr_extended_dcel<Traits2, int,EdgeData,int>  Dcel;
      typedef Arrangement_2<Traits2, Dcel>                  Arr_2;
      typedef Arr_2::Face_iterator                          Face_iterator;
      typedef Arr_2::Face_handle                            Face_handle;
      typedef Arr_2::Edge_iterator                          Edge_iterator;
      typedef Arr_2::Halfedge_iterator                      Halfedge_iterator;
      typedef Arr_2::Vertex_iterator                        Vertex_iterator;
      typedef Arr_2::Vertex_handle                          Vertex_handle;
      typedef Arr_2::Halfedge_handle                        Halfedge_handle;
      typedef Arr_2::Ccb_halfedge_circulator                Ccb_halfedge_circulator;
      typedef Arr_mesh_generation_overlay_traits <Arr_2, Arr_2,Arr_2>  Overlay_traits;
      
      IGL_INLINE void generate_mesh(int N,
                                    const Eigen::MatrixXd& V,
                                    const Eigen::MatrixXi& F,
                                    const Eigen::MatrixXi& EV,
                                    const Eigen::MatrixXi& FE,
                                    const Eigen::MatrixXi& EF,
                                    const Eigen::MatrixXd& TC,
                                    const Eigen::MatrixXi& FTC,
                                    const double edgeLength,
                                    Eigen::MatrixXd& newV,
                                    Eigen::VectorXi& newD,
                                    Eigen::MatrixXi& newF){
        
        Mesh origMesh, genMesh;
        
        origMesh.Allocate(V.rows(), F.rows(), 3*F.rows());
        for (int i=0;i<V.rows();i++){
          origMesh.Vertices[i].Coordinates=Point3D(V(i,0),V(i,1),V(i,2));
          origMesh.Vertices[i].ID=i;
        }
        
        for (int i=0;i<F.rows();i++)
          for (int j=0;j<3;j++)
            origMesh.Faces[i].Vertices[j]=F(i,j);
        
        
        //origin, next, prev, twin, face, texu, texv
        for (int i=0;i<F.rows();i++){
          for (int j=0;j<3;j++){
            origMesh.Halfedges[3*i+j].ID=3*i+j;
            Eigen::RowVector2d uv=TC.row(FTC(i,j));
            origMesh.Halfedges[3*i+j].Origin=F(i,j);
            origMesh.Halfedges[3*i+j].Next=3*i+(j+1)%3;
            origMesh.Halfedges[3*i+j].Prev=3*i+(j+2)%3;
            origMesh.Halfedges[3*i+j].AdjFace=i;
            
            int otherface = (EF(FE(i,j),0)==i ? EF(FE(i,j),1) : EF(FE(i,j),0));
            if (otherface!=-1){
              for (int k=0;k<3;k++){
                if ((EV(otherface,(k+1)%3)==EV(i,j))&&(EV(otherface,k)==EV(i,(j+1)%3)))
                  origMesh.Halfedges[3*i+j].Twin=3*otherface+k;
                origMesh.Halfedges[3*otherface+k].Twin=3*i+j;
              }
            }
            origMesh.Halfedges[3*i+j].TexCoords=Point2D(uv(0), uv(1));
            origMesh.Faces[origMesh.Halfedges[3*i+j].AdjFace].AdjHalfedge=3*i+j;
            origMesh.Vertices[origMesh.Halfedges[3*i+j].Origin].AdjHalfedge=3*i+j;
          }
        }
        
        origMesh.GenerateMesh(genMesh, edgeLength);
        genMesh.SimplifyHexMesh();
        
        newV.conservativeResize(genMesh.Vertices.size(),3);
        for (int i=0;i<newV.rows();i++)
          newV.row(i)<<genMesh.Vertices[i].Coordinates.x(),genMesh.Vertices[i].Coordinates.y(),genMesh.Vertices[i].Coordinates.z();
        
        newD.conservativeResize(genMesh.Faces.size());
        for (int i=0;i<newD.size();i++){
          int ebegin=genMesh.Faces[i].AdjHalfedge;
          int ecurr=ebegin;
          newD(i)=0;
          do{
            newD(i)++;
            ecurr=genMesh.Halfedges[ecurr].Next;
          }while(ebegin!=ecurr);
        }
        
        newF.conservativeResize(newD.rows(),newD.maxCoeff());
        for (int i=0;i<newF.rows();i++){
          int ebegin=genMesh.Faces[i].AdjHalfedge;
          int ecurr=ebegin;
          int currIndex=0;
          do{
            newF(i,currIndex++)=genMesh.Halfedges[ecurr].Origin;
            ecurr=genMesh.Halfedges[ecurr].Next;
          }while(ebegin!=ecurr);
        }
        
      }
    }
  }
}
  
#include "generate_mesh.cpp"
  
#endif
  
  
