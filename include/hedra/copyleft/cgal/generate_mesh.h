// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2019 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_COPYLEFT_CGAL_EXTRACT_MESH_H
#define HEDRA_COPYLEFT_CGAL_EXTRACT_MESH_H
#include <igl/igl_inline.h>

#include <Eigen/Dense>
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <hedra/copyleft/cgal/basic_cgal_definitions.h>
#include <hedra/dcel.h>
#include <vector>



namespace hedra
{
  namespace copyleft
  {
    namespace cgal
    {
      
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
      
      /*class Vertex{
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
            normalize(Faces[i].Normal);
          }
        }
        
        Mesh(){}
        ~Mesh(){}
        
      };*/
      
      const int PARAM_LINE_VERTEX = -2;
      const int ORIGINAL_VERTEX = -1;
      
      struct ArrEdgeData{
        bool isParam;
        int origEdge;
        int newHalfedge;
      };
      
      //f1: original mesh
      //f2: parameter lines
      template <class ArrangementA, class ArrangementB, class ArrangementR>
      class Arr_mesh_generation_overlay_traits :
      public CGAL::_Arr_default_overlay_traits_base<ArrangementA,ArrangementB,ArrangementR>
      {
      public:
        
        typedef typename ArrangementA::Face_const_handle    Face_handle_A;
        typedef typename ArrangementB::Face_const_handle    Face_handle_B;
        typedef typename ArrangementR::Face_handle          Face_handle_C;
        
        typedef typename ArrangementA::Vertex_const_handle    Vertex_handle_A;
        typedef typename ArrangementB::Vertex_const_handle    Vertex_handle_B;
        typedef typename ArrangementR::Vertex_handle          Vertex_handle_C;
        
        typedef typename ArrangementA::Halfedge_const_handle    Halfedge_handle_A;
        typedef typename ArrangementB::Halfedge_const_handle    Halfedge_handle_B;
        typedef typename ArrangementR::Halfedge_handle          Halfedge_handle_C;
        
      public:
        
        virtual void create_face (Face_handle_A f1,
                                  Face_handle_B f2,
                                  Face_handle_C f) const
        {
          // Overlay the data objects associated with f1 and f2 and store the result
          // with f.
          f->set_data (f1->data());
        }
        
        virtual void  create_vertex ( Vertex_handle_A v1, Vertex_handle_B v2, Vertex_handle_C v)
        {
          v->set_data(PARAM_LINE_VERTEX);
        }
        
        virtual void create_vertex ( Vertex_handle_A v1, Halfedge_handle_B e2, Vertex_handle_C v)
        {
          v->set_data(ORIGINAL_VERTEX);
        }
        
        virtual void create_vertex ( Vertex_handle_A v1, Face_handle_B f2, Vertex_handle_C v)
        {
          v->set_data(ORIGINAL_VERTEX);
        }
        
        virtual void create_vertex ( Halfedge_handle_A e1, Vertex_handle_B v2, Vertex_handle_C v)
        {
          v->set_data(PARAM_LINE_VERTEX);
        }
        
        virtual void create_vertex ( Face_handle_A f1, Vertex_handle_B v2, Vertex_handle_C v)
        {
          v->set_data(PARAM_LINE_VERTEX);
        }
        
        virtual void create_vertex ( Halfedge_handle_A e1, Halfedge_handle_B e2, Vertex_handle_C v)
        {
          v->set_data(ORIGINAL_VERTEX);
        }
        
        virtual void create_edge ( Halfedge_handle_A e1, Halfedge_handle_B e2, Halfedge_handle_C e)
        {
          ArrEdgeData aed;
          aed.isParam=true;
          aed.origEdge = e1->data().origEdge;
          aed.newHalfedge = -1;
          e->set_data(aed);
          e->twin()->set_data(aed);
        }
        
        virtual void create_edge ( Halfedge_handle_A e1, Face_handle_B f2, Halfedge_handle_C e)
        {
          ArrEdgeData aed;
          aed.isParam=false;
          aed.origEdge = e1->data().origEdge;
          aed.newHalfedge = -1;
          e->set_data(aed);
          e->twin()->set_data(aed);
        }
        
        virtual void create_edge ( Face_handle_A f1, Halfedge_handle_B e2, Halfedge_handle_C e)
        {
          ArrEdgeData aed;
          aed.isParam=true;
          aed.origEdge = -1;
          aed.newHalfedge = -1;
          e->set_data(aed);
          e->twin()->set_data(aed);
        }
      };
      
      typedef CGAL::Arr_linear_traits_2<EKernel>                 Traits2;
      typedef Traits2::Point_2                              Point2;
      typedef Traits2::Segment_2                   Segment2;
      typedef Traits2::Line_2                   Line2;
      typedef Traits2::X_monotone_curve_2                  X_monotone_curve_2;
      
      typedef CGAL::Arr_extended_dcel<Traits2, int,ArrEdgeData,int>  Dcel;
      typedef CGAL::Arrangement_2<Traits2, Dcel>                  Arr_2;
      typedef Arr_2::Face_iterator                          Face_iterator;
      typedef Arr_2::Face_handle                            Face_handle;
      typedef Arr_2::Edge_iterator                          Edge_iterator;
      typedef Arr_2::Halfedge_iterator                      Halfedge_iterator;
      typedef Arr_2::Vertex_iterator                        Vertex_iterator;
      typedef Arr_2::Vertex_handle                          Vertex_handle;
      typedef Arr_2::Halfedge_handle                        Halfedge_handle;
      typedef Arr_2::Ccb_halfedge_circulator                Ccb_halfedge_circulator;
      typedef Arr_mesh_generation_overlay_traits <Arr_2, Arr_2,Arr_2>  Overlay_traits;
      
      
      //for now doing quad (u,v,-u -v) only!
      Point2 paramCoord2texCoord(Eigen::RowVectorXd paramCoord, int Resolution)
      {
        ENumber u=ENumber((int)(paramCoord(0)*(double)Resolution),Resolution);
        ENumber v=ENumber((int)(paramCoord(1)*(double)Resolution),Resolution);
        return Point2(u,v);
      }
      
      IGL_INLINE void stitch_boundaries(const Eigen::VectorXi triEF,
                                        const Eigen::VectorXi triInnerEdges,
                                        Eigen::MatrixXd& currV,
                                       Eigen::VectorXi& VH,
                                       Eigen::VectorXi& HV,
                                       Eigen::VectorXi& HF,
                                       Eigen::VectorXi& FH,
                                       Eigen::VectorXi& nextH,
                                       Eigen::VectorXi& prevH,
                                       Eigen::VectorXi& twinH,
                                       std::vector<bool>& isParamVertex,
                                       std::vector<int>& HE2origEdges,
                                       std::vector<bool>& isParamHE,
                                       std::vector<int>& overlayFace2Tri,
                                       const double closeTolerance)
      {

        using namespace Eigen;
        
        //TODO: tie all endpoint vertices to original triangles
        
        VectorXi old2NewV=VectorXi::Constant(currV.rows(),-1);
        
        std::vector<std::vector<int>> origEdges2HE(triEF.rows());
        for (int i=0;i<HE2origEdges.size();i++)
          origEdges2HE[HE2origEdges[i]].push_back(i);
        
        
        
        
        //for every original inner edge, stitching up boundary (original boundary edges don't have any action item)
        for (int i=0;i<triInnerEdges.size();i++){
          //first sorting to left and right edges according to faces
          int currEdge=triInnerEdges(i);
          
          int leftFace=triEF(currEdge,0);
          int rightFace=triEF(currEdge,1);
          
          std::vector<int> leftHE, rightHE;
          
          for (int k=0;k<origEdges2HE[currEdge].size();k++){
            if (overlayFace2Tri[HF(origEdges2HE[currEdge][k])]==leftFace)
              leftHE.push_back(origEdges2HE[currEdge][k]);
            else if (overlayFace2Tri[HF(origEdges2HE[currEdge][k])]==rightFace)
              rightHE.push_back(origEdges2HE[currEdge][k]);
            else
              int kaka=8;  //shouldn't happen
            
          }
          
          //if the parameterization is seamless, left and right halfedges should be perfectly matched, but it's not always the case
         
    
          
        }
        
      }
      
      IGL_INLINE void generate_mesh(int N,
                                    const Eigen::MatrixXd& V,
                                    const Eigen::MatrixXi& F,
                                    const Eigen::MatrixXi& EV,
                                    const Eigen::MatrixXi& FE,
                                    const Eigen::MatrixXi& EF,
                                    const Eigen::MatrixXi& EFi,
                                    const Eigen::VectorXi& innerEdges,
                                    const Eigen::MatrixXd& PC,
                                    const Eigen::MatrixXi& FPC,
                                    Eigen::MatrixXd& newV,
                                    Eigen::VectorXi& newD,
                                    Eigen::MatrixXi& newF){
        
        
        using namespace Eigen;
        using namespace std;
        VectorXi VH;
        MatrixXi EH;
        VectorXi HV,HE,HF,FH;
        VectorXi nextH,prevH,twinH;
        
        double minrange=(PC.colwise().maxCoeff()-PC.colwise().minCoeff()).minCoeff();
        int resolution=pow(10,ceil(log10(100000/minrange)));
        
       // hedra::DCEL(VectorXd::Constant(F.rows(),3),F,EV,EF,EFi,innerEdges,VH,EH,FH,HV,HE,HF,nextH,prevH,twinH);
        
        //creating an single-triangle arrangement
        
        //Intermediate growing DCEL
        std::vector<bool> isParamVertex;
        std::vector<int> DList;
        std::vector<int> HE2origEdges;
        std::vector<bool> isParamHE;
        std::vector<int> overlayFace2Triangle;
        
        MatrixXd currV(isParamVertex.size(),3);
        VH.resize(currV.rows());
        HV.resize(HE2origEdges.size());
        HF.resize(HE2origEdges.size());
        FH.resize(0);
        nextH.resize(HE2origEdges.size());
        prevH.resize(HE2origEdges.size());
        twinH.resize(HE2origEdges.size());
        
        for (int ti=0;ti<F.rows();ti++){
          
          Arr_2 paramArr, triangleArr, overlayArr;
          
          for (int j=0;j<3;j++){
            RowVectorXd PC1 = PC.row(FPC(ti,j));
            RowVectorXd PC2 = PC.row(FPC(ti,(j+1)%3));
            
            Halfedge_handle he=CGAL::insert_non_intersecting_curve(triangleArr, Segment2(paramCoord2texCoord(PC1,resolution),paramCoord2texCoord(PC2,resolution)));
            
            ArrEdgeData aed;
            aed.isParam=false;
            aed.origEdge =FE(ti,j);
            he->set_data(aed);
            he->twin()->set_data(aed);
          }
          
          for (Face_iterator fi= triangleArr.faces_begin(); fi != triangleArr.faces_end(); fi++){
            if (fi->is_unbounded())
              fi->data()=-1;
            else
              fi->data()=ti;
          }
          
          //creating an arrangement of parameter lines
          MatrixXd facePC(3, PC.cols());
          for (int i=0;i<3;i++)
              facePC.row(i)=PC.row(FPC(ti,i));
          
          for (int i=0;i<facePC.cols();i++){
            //inserting unbounded lines
            
            int coordMin = (int)std::floor(facePC.col(i).minCoeff()-1.0);
            int coordMax = (int)std::ceil(facePC.col(i).minCoeff()+1.0);
            vector<X_monotone_curve_2> lineCurves;
            for (int coordIndex=coordMin;coordIndex<=coordMax;coordIndex++){
              
              //The line coord = coordIndex
              RowVectorXd LineCoord1 = RowVectorXd::Zero(facePC.cols());
              RowVectorXd LineCoord2 = RowVectorXd::Ones(facePC.cols());
              LineCoord1(i)=coordIndex;
              LineCoord2(i)=coordIndex;
              lineCurves.push_back(Line2(paramCoord2texCoord(LineCoord1,resolution), paramCoord2texCoord(LineCoord2,resolution)));
            }
            insert(paramArr, lineCurves.begin(), lineCurves.end());
          }
          
          
          //Constructing the overlay arrangement
          Overlay_traits ot;
          overlay (triangleArr, paramArr, overlayArr, ot);
          
          //creating new halfedge structure from given mesh
          int formerNumVertices = currV.rows();
          int formerNumHalfedges =nextH.rows();
          int formerNumFaces =FH.rows();
          
          int currFace=0, currVertex=0, currHalfedge=0;
          for (Face_iterator fi=overlayArr.faces_begin();fi!=overlayArr.faces_end();fi++){
            if (fi->data()==-1)
              continue;  //one of the outer faces
            
            overlayFace2triangle.push_back(fi->data());
            fi->data()=formerNumFaces+currFace;
            currFace++;
            int DFace=0;
            Ccb_halfedge_circulator hebegin=fi->outer_ccb ();
            Ccb_halfedge_circulator heiterate=hebegin;
            do{
              DFace++;
              if (heiterate->source()->data()<0){  //new vertex
                
                isParamVertex.push_back(heiterate->source()->data()==PARAM_LINE_VERTEX);
                heiterate->source()->data()=formerNumVertices+currVertex;
                currVertex++;
              }
              
              if (heiterate->data().newHalfedge<0){  //new halfedge
                HE2origEdges.push_back(heiterate->data().origEdge);
                isParamHE.push_back(heiterate->data().isParam);
                heiterate->data().newHalfedge=formerNumHalfedges+currHalfedge;
                currHalfedge++;
              }
              heiterate++;
            }while(heiterate!=hebegin);
            
          }
          
          currV.conservativeResize(currV.rows()+currVertex,3);
          VH.conservativeResize(VH.size()+currVertex);
          HV.conservativeResize(HV.size()+currHalfedge);
          HF.conservativeResize(HF.size()+currHalfedge);
          FH.conservativeResize(FH.size()+currFace);
          nextH.conservativeResize(nextH.size()+currHalfedge);
          prevH.conservativeResize(prevH.size()+currHalfedge);
          twinH.conservativeResize(twinH.size()+currHalfedge);
          
          for (Face_iterator fi=overlayArr.faces_begin();fi!=overlayArr.faces_end();fi++){
            if (fi->data()==-1)
              continue;  //one of the outer faces
            
            Ccb_halfedge_circulator hebegin=fi->outer_ccb ();
            Ccb_halfedge_circulator heiterate=hebegin;
            //now assigning nexts and prevs
            do{
              nextH(heiterate->data().newHalfedge) = heiterate->next()->data().newHalfedge;
              prevH(heiterate->data().newHalfedge) = heiterate->prev()->data().newHalfedge;
              twinH(heiterate->data().newHalfedge) = heiterate->twin()->data().newHalfedge;
              if (heiterate->twin()->data().newHalfedge!=-1)
               twinH(heiterate->twin()->data().newHalfedge) =heiterate->data().newHalfedge;
              
              HV(heiterate->data().newHalfedge) = heiterate->source()->data();
              VH(heiterate->source()->data())=heiterate->data().newHalfedge;
              HF(heiterate->data().newHalfedge) = fi->data();
              FH(fi->data())=heiterate->data().newHalfedge;
              heiterate++;
            }while (heiterate!=hebegin);
          }
          
          //constructing the actual vertices
          for (Vertex_iterator vi=overlayArr.vertices_begin();vi!=overlayArr.vertices_end();vi++){
            
            if (vi->data()<0)
              continue;
            
            ENumber BaryValues[3];
            
            ENumber Sum=0;
            
            for (int i=0;i<3;i++){
              //finding out barycentric coordinates
              RowVectorXd PC2 = PC.row(FPC(ti,(i+1)%3));
              RowVectorXd PC3 = PC.row(FPC(ti,(i+2)%3));
              ETriangle2D t(vi->point(), paramCoord2texCoord(PC2,resolution),  paramCoord2texCoord(PC3,resolution));
              BaryValues[i]=t.area();
              Sum+=BaryValues[i];
            }
                            
            for (int i=0;i<3;i++)
              BaryValues[i]/=Sum;
            
            EPoint3D ENewPosition(0,0,0);
            for (int i=0;i<3;i++){
              
              EPoint3D vertexCoord(ENumber((int)(V(F(ti,i),0)*(double)resolution),resolution),
                                   ENumber((int)(V(F(ti,i),1)*(double)resolution),resolution),
                                   ENumber((int)(V(F(ti,i),2)*(double)resolution),resolution));
              ENewPosition=ENewPosition+(vertexCoord-CGAL::ORIGIN)*BaryValues[i];
            }
                          
            RowVector3d newPosition(to_double(ENewPosition.x()), to_double(ENewPosition.y()), to_double(ENewPosition.z()));
            currV.row(vi->data())=newPosition;
                            
            //DebugLog<<"Creating Vertex "<<vi->data()<<" with 2D coordinates ("<<vi->point().x()<<","<<vi->point().y()<<") "<<" and 3D Coordinates ("<<std::setprecision(10) <<NewPosition.x()<<","<<NewPosition.y()<<","<<NewPosition.z()<<")\n";
          }
        }
        
        //mesh unification
        
        stitch_boundaries(currV,VH,HV,HF,FH,nextH,prevH,twinH, isParamVertex, HE2origEdges, isParamHE, overlayFace2Triangle);
        
        //consolidation
        newV=currV;
      
        newD.conservativeResize(FH.rows());
        for (int i=0;i<newD.size();i++){
          int ebegin=FH(i);
          int ecurr=ebegin;
          newD(i)=0;
          do{
            newD(i)++;
            ecurr=nextH(ecurr);
          }while(ebegin!=ecurr);
        }
        
        newF.conservativeResize(newD.rows(),newD.maxCoeff());
        for (int i=0;i<newF.rows();i++){
          int ebegin=FH(i);
          int ecurr=ebegin;
          int currIndex=0;
          do{
            newF(i,currIndex++)=HV(ecurr);
            ecurr=nextH(ecurr);
          }while(ebegin!=ecurr);
        }
      }
    }
  }
}

#endif
  
  
