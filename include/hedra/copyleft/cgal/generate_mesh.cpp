#include <hedra/copyleft/cgal/generate_mesh.h>
#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <CGAL/Search_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <math.h>

//typedef CGAL::Search_traits_2<R> Traits;
//typedef CGAL::Euclidean_distance<Traits> Distance;
//typedef CGAL::Fair<Traits> Fair;
//typedef CGAL::Orthogonal_k_neighbor_search<Traits,Distance,Fair> Neighbor_search;
//typedef Neighbor_search::Tree Tree;

int hedra::copyleft::cgal::Vertex::RefCount=0;
int hedra::copyleft::cgal::Halfedge::RefCount=0;
int hedra::copyleft::cgal::Face::RefCount=0;

using namespace boost;
using namespace CGAL;




//#define wx (-0.5)
//#define wy (SQRT3/2.0)

#define SQRT3 sqrt(3.0)


struct MergeData {
    const bool operator()(const int& v1, const int& v2) const {return v1; }
};



typedef CGAL::Arr_segment_traits_2<EKernel>                     Traits2;
typedef Traits2::Point_2                                       Point2;
typedef Traits2::X_monotone_curve_2                            Segment2;
typedef CGAL::Arr_extended_dcel<Traits2, int,hedra::copyleft::cgal::EdgeData,int>          Dcel;
typedef CGAL::Arrangement_2<Traits2, Dcel>                     Arr_2;
typedef Arr_2::Face_iterator								   Face_iterator;
typedef Arr_2::Face_handle									   Face_handle;
typedef Arr_2::Edge_iterator								   Edge_iterator;
typedef Arr_2::Halfedge_iterator							   Halfedge_iterator;
typedef Arr_2::Vertex_iterator								   Vertex_iterator;
typedef Arr_2::Vertex_handle									Vertex_handle;
typedef Arr_2::Halfedge_handle									Halfedge_handle;
typedef Arr_2::Ccb_halfedge_circulator					       Ccb_halfedge_circulator;
typedef hedra::copyleft::cgal::Arr_mesh_generation_overlay_traits <Arr_2,
                                      Arr_2,
                                      Arr_2>  Overlay_traits;


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


Point2 Hex2Euc(double u, double v, int Resolution)
{
   ENumber cx=ENumber((int)(((double)u*3.0-(double)v*1.5)*Resolution),Resolution);
   ENumber cy=ENumber((int)(((double)v*SQRT3/2.0)*Resolution),Resolution);
   return Point2(cx,cy);
}

Point2D Hex2Euc(double u, double v)
{
   double cx=u*3.0-v*1.5;
   double cy=v*SQRT3/2.0;
   return Point2D(cx,cy);
}

Point2D Euc2Hex(ENumber x, ENumber y)
{
	double cx=x.to_double();
	double cy=y.to_double();
	double v=cy*2.0/SQRT3;
	double u=(cx+1.5*v)/3.0;
	return Point2D(u,v);
}

Point2D Euc2Hex(double cx, double cy)
{
	double v=cy*2.0/SQRT3;
	double u=(cx+1.5*v)/3.0;
	return Point2D(u,v);
}


bool CircleSegmentIntersect(double x1, double y1, double x2, double y2, double cx, double cy, double cr ) 
{
	Vector2D v1=Point2D(x1,y1)-Point2D(cx,cy);
	Vector2D v2=Point2D(x2,y2)-Point2D(cx,cy);

	double v1n=v1.squared_length();
	double v2n=v2.squared_length();
	double v1v2=v1*v2;

	double a=v1n+v2n-2*v1v2;
	double b=-2*v2n+2*v1v2;
	double c=v2n-cr;

	double delta= b * b - 4 * a * c;
	if (delta<0)
		return false;

	double t1=(-b-sqrt(delta))/(2*a);
	double t2=(-b+sqrt(delta))/(2*a);

	return (((t1>=-0.3)&&(t1<=1.3))||((t2>=-0.3)&&(t2<=1.3)));  //leaving some margins for conservative checks
}

bool DoIntersect(Triangle2D& t, Circle2D& c)
{
	//checking if actual intersection
	for (int i=0;i<3;i++)
		if (CircleSegmentIntersect(t.vertex(i).x(),t.vertex(i).y(),t.vertex((i+1)%3).x(),t.vertex((i+1)%3).y(), c.center().x(), c.center().y(), c.squared_radius()))
			return true;
	
	//if (CGAL::do_intersect(c,t))
		//return true;

	//checking if circle includes triangle
	bool insides[3];
	for (int i=0;i<3;i++)
		insides[i]=((t.vertex(i)-c.center()).squared_length()<=1.3*c.squared_radius());

	if (insides[0] && insides[1] && insides[2])
		return true;  //circle includes triangle

	//checking if triangle includes center
	double BaryCoords[3];
	double Area=t.area();
	for (int i=0;i<3;i++)
		BaryCoords[i]=Triangle2D(c.center(), t.vertex((i+1)%3), t.vertex((i+2)%3)).area()/Area;

	for (int i=0;i<3;i++)
		if ((BaryCoords[i]<-0.3)||(BaryCoords[i]>1.3))
			return false;   //center is outside

	return true;
	

}

double round(double number)
{
    return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}


void hedra::copyleft::cgal::Mesh::GenerateMesh(Mesh& HexMesh, double EdgeLength)
{
	HexMesh.Vertices.clear();
	HexMesh.Halfedges.clear();	
	HexMesh.Faces.clear();

	//DebugLog.open("Debugging.txt");
	
	//generating regular elliptic pattern 2D arrangement, and locating vertices of mesh in it
	
	int eID=0, vID=0;

	for (int findex=0;findex<Faces.size();findex++){
		//if (findex>100)
		//	continue;

		//building small face overlays of one triangle and a few roughly surrounding hexes to retrieve the structure in the face
	
	   int ebegin=Faces[findex].AdjHalfedge;
	   int eiterate=ebegin;
    std::vector<Point2D> TriPoints2D;
	   do{
		   Point2D Location=Halfedges[eiterate].TexCoords;
		   TriPoints2D.push_back(Location);
		   eiterate=Halfedges[eiterate].Next;
	   }while (eiterate!=ebegin);
	   Triangle2D CurrTri(TriPoints2D[0], TriPoints2D[1], TriPoints2D[2]);

	   //if triangle is degenerate, exit without doing anything
	   if (abs(CurrTri.area())<10e-12)
		   continue;

		//finding a cover of the triangle by flood filling hexagons
    std::set<std::pair<int, int> > Circles;
		std::queue<std::pair<int, int> > CircleQueue;
		//inserting adjacent circles to vertices
		ebegin=Faces[findex].AdjHalfedge;
		eiterate=ebegin;
		do{
			Point2D uv=Euc2Hex(Halfedges[eiterate].TexCoords.x(), Halfedges[eiterate].TexCoords.y());
			//double v=Halfedges[eiterate].TexCoords.y()/wy;
			//double u=(Halfedges[eiterate].TexCoords.x()+1.5*v)/3.0;
			
			//if it intersects, the neighbors should in line as well
			for (int i=-4;i<=4;i++)
				for (int j=-4;j<=4;j++)
					CircleQueue.push(std::pair<int,int>(floor(uv.x())+i, floor(uv.y())+j));

			eiterate=Halfedges[eiterate].Next;
		}while (eiterate!=ebegin);

		std::set<std::pair<int, int> > Rejects;

		while (!CircleQueue.empty())
		{
			std::pair<int, int> CurrPair=CircleQueue.front();
			CircleQueue.pop();
			if ((Circles.find(std::pair<int, int>(CurrPair.first, CurrPair.second))!=Circles.end())||(Rejects.find(std::pair<int, int>(CurrPair.first, CurrPair.second))!=Rejects.end()))
				continue;
			Circle2D CurrCircle(Hex2Euc(CurrPair.first,CurrPair.second),1);
			
			if (!DoIntersect(CurrTri, CurrCircle)){
				Rejects.insert(CurrPair);
				continue;
			} else {
				Circles.insert(CurrPair);
			}
		
			//if it intersects, the neighbors should in line as well
			for (int i=-4;i<=4;i++)
				for (int j=-4;j<=4;j++){
					if ((i==0)&&(j==0))
						continue;
					//if (((Circles.find(pair<int, int>(CurrPair.first+i, CurrPair.second+j))==Circles.end()))&&(Rejects.find(pair<int, int>(CurrPair.first+i, CurrPair.second+j))==Rejects.end()))
						CircleQueue.push(std::pair<int, int>(CurrPair.first+i, CurrPair.second+j));
				}
		}

		//DebugLog<<"Working on triangle "<<findex<<"\n";
		double fmaxu=-3276700.0, fmaxv=-327600.0;
		double fminu=3276700.0, fminv=327600.0;
		Arr_2 PrimalArr, DualArr, TriangleArr, FullArr;
		ebegin=Faces[findex].AdjHalfedge;
		eiterate=ebegin;
		do{
			if (Halfedges[eiterate].TexCoords.x()>fmaxu) fmaxu=Halfedges[eiterate].TexCoords.x();
			if (Halfedges[eiterate].TexCoords.y()>fmaxv) fmaxv=Halfedges[eiterate].TexCoords.y();
			if (Halfedges[eiterate].TexCoords.x()<fminu) fminu=Halfedges[eiterate].TexCoords.x();
			if (Halfedges[eiterate].TexCoords.y()<fminv) fminv=Halfedges[eiterate].TexCoords.y();
			eiterate=Halfedges[eiterate].Next;
		}while (eiterate!=ebegin);
	
	   double minrange=std::min(fmaxu-fminu,fmaxv-fminv);
       int Resolution=pow(10,ceil(log10(1000/minrange)));


		//building the one-triangle arrangement
	   ebegin=Faces[findex].AdjHalfedge;
	   eiterate=ebegin;
	   std::vector<Point2> TriPoints;
	   std::vector<EPoint3D> ETriPoints3D;
	   std::vector<EdgeData> EdgeDatas;
	   do{
		   Point2D Location=Halfedges[eiterate].TexCoords;
		   //rounding location by a bit to avoid numerical problems later
		   /*double newx=Location.x(), newy=Location.y();
		   if (abs(round((double)(Location.x()*1000.0*(double)Resolution))-Location.x()*1000.0*(double)Resolution)<(1.0/5000.0))
			   newx=round(Location.x()*1000.0*(double)Resolution)/(1000.0*(double)Resolution);

		   if (abs(round(Location.y()*1000.0*(double)Resolution)-Location.y()*1000.0*(double)Resolution)<(1.0/5000.0))
			   newy=round(Location.y()*1000.0*(double)Resolution)/(1000.0*(double)Resolution);*/


		   Point3D Position=Vertices[Halfedges[eiterate].Origin].Coordinates;
		   ENumber cx=ENumber((int)(Location.x()*Resolution),Resolution);
		   ENumber cy=ENumber((int)(Location.y()*Resolution),Resolution);
		   ENumber x=ENumber((int)(Position.x()*Resolution),Resolution);
		   ENumber y=ENumber((int)(Position.y()*Resolution),Resolution);
		   ENumber z=ENumber((int)(Position.z()*Resolution),Resolution);
		   TriPoints.push_back(Point2(cx,cy));
		   ETriPoints3D.push_back(EPoint3D(x,y,z));
		   int DomEdge;
		   if ((Halfedges[eiterate].Twin<0)||(Halfedges[eiterate].Twin>eiterate))
			   DomEdge=eiterate;
		   else
			   DomEdge=Halfedges[eiterate].Twin;
		   EdgeData ed; ed.OrigHalfedge=DomEdge;
		   ed.isBoundary=(Halfedges[eiterate].Twin<0);
		   EdgeDatas.push_back(ed);
		   eiterate=Halfedges[eiterate].Next;
	   }while(ebegin!=eiterate);

	   for (int i=0;i<3;i++){
		   Halfedge_handle he=CGAL::insert_non_intersecting_curve(TriangleArr, Segment2(TriPoints[i],TriPoints[(i+1)%3]));
		   he->set_data(EdgeDatas[i]);
		   if (EdgeDatas[i].isBoundary)
			   he->source()->data()=he->target()->data()=0;
		   else
			   he->source()->data()=he->target()->data()=1;

		   he->twin()->set_data(EdgeDatas[i]);
	   }
	   
	   for (Face_iterator fi= TriangleArr.faces_begin(); fi != TriangleArr.faces_end(); fi++){
		   if (fi->is_unbounded())
			   fi->data()=0;
		   else
			   fi->data()=1;
	   }

        //Face_handle  uf = PrimalArr.unbounded_face();

       //creating the primal hexagon grid (actually triangular grid)
       std::vector<Segment2> PrimalSegments;
	   std::set<HexSegment> HexSegments;
	   for (std::set<std::pair<int, int> >::iterator ci=Circles.begin();ci!=Circles.end();ci++)  {
		   HexSegments.insert(HexSegment(ci->first,ci->second,ci->first+1,ci->second+1));
		   HexSegments.insert(HexSegment(ci->first,ci->second,ci->first+1,ci->second+2));
		   HexSegments.insert(HexSegment(ci->first,ci->second,ci->first,ci->second+1));
		   HexSegments.insert(HexSegment(ci->first,ci->second,ci->first-1,ci->second-1));
		   HexSegments.insert(HexSegment(ci->first,ci->second,ci->first-1,ci->second-2));
		   HexSegments.insert(HexSegment(ci->first,ci->second,ci->first,ci->second-1));
	   }
	   for (std::set<HexSegment>::iterator hi=HexSegments.begin();hi!=HexSegments.end();hi++)  {
		   Point2 p1=Hex2Euc(hi->Source.first,hi->Source.second, Resolution);
		   Point2 p2=Hex2Euc(hi->Target.first,hi->Target.second, Resolution);

		   PrimalSegments.push_back(Segment2(p1,p2));
       }
       CGAL::insert_non_intersecting_curves(PrimalArr,PrimalSegments.begin(), PrimalSegments.end());

	   //creating dual arrangement
       //uf = DualArr.unbounded_face();

       int NumPrimVertices=PrimalArr.number_of_vertices();
       int NumPrimEdges=PrimalArr.number_of_edges();
       int NumPrimFaces=PrimalArr.number_of_faces();

       //faces - > dual vertices
       int ID=0;
       std::vector<Segment2> DualCurves;
       std::vector<Point2> FaceCenters(PrimalArr.number_of_faces());
       for (Face_iterator fi= PrimalArr.faces_begin(); fi != PrimalArr.faces_end(); fi++){
           if (fi->is_fictitious()||fi->is_unbounded())
               continue;
           Ccb_halfedge_circulator hc=fi->outer_ccb ();
           Point2 Location;
           for (int i=0;i<3;i++){
               Location=Location+(hc->source()->point()-CGAL::ORIGIN);
               hc++;
           }

           Location=CGAL::ORIGIN+(Location-CGAL::ORIGIN)/ENumber(3);

		   //removing potential small faces by projecting point with barycentric coordinates which are very close to 0 or 1 with respect to the triangle
		    //finding out barycentric coordinates
		   ENumber BaryValues[3];
		   ENumber Sum=0;
		   for (int i=0;i<3;i++){
			   ETriangle2D t(Location, TriPoints[(i+1)%3], TriPoints[(i+2)%3]);
			   BaryValues[i]=t.area();
			   Sum+=BaryValues[i];
		   }
		   for (int i=0;i<3;i++)
			   BaryValues[i]/=Sum;

		   double Bary1=to_double(BaryValues[0]);
		   double Bary2=to_double(BaryValues[1]);
		   double Bary3=to_double(BaryValues[2]);

		   for (int i=0;i<3;i++){
			   double dBary=to_double(BaryValues[i]);
			   if (abs(dBary)<0.00005){
				   BaryValues[i]=ENumber(0);
				   ENumber SumOthers=(BaryValues[(i+1)%3]+BaryValues[(i+2)%3]);
				   BaryValues[(i+1)%3]/=SumOthers;
				   BaryValues[(i+2)%3]/=SumOthers;
				   //DebugLog<<"Moving Location from ("<<to_double(Location.x())<<","<<to_double(Location.y())<<") to ";
				   Location=CGAL::ORIGIN+((TriPoints[(i+1)%3]-CGAL::ORIGIN)*BaryValues[(i+1)%3]+(TriPoints[(i+2)%3]-CGAL::ORIGIN)*BaryValues[(i+2)%3]);
				   //DebugLog<<"("<<to_double(Location.x())<<","<<to_double(Location.y())<<")\n";
			   }
		   }
		   
           FaceCenters[ID]=Location;
           fi->data()=ID++;

       }

       //edges -> dual edges
       for (Edge_iterator ei=PrimalArr.edges_begin();ei!=PrimalArr.edges_end();ei++){
           //counter++;
           if (ei->face()->is_fictitious()||ei->face()->is_unbounded()||ei->twin()->face()->is_fictitious()||ei->twin()->face()->is_unbounded())
               continue;

           Point2 p1=FaceCenters[ei->face()->data()];
           Point2 p2=FaceCenters[ei->twin()->face()->data()];
           Segment2 s2(p1,p2);
           DualCurves.push_back(s2);
           //Halfedge_handle eh=EllArr.insert_at_vertices(s2, v1,v2);
       }
       CGAL::insert_non_intersecting_curves(DualArr, DualCurves.begin(), DualCurves.end());

	   //DualArr=PrimalArr;

	   Overlay_traits ot;
	   overlay (TriangleArr, DualArr, FullArr, ot);


	   for (Face_iterator fi=FullArr.faces_begin();fi!=FullArr.faces_end();fi++){
		   if (!fi->data())
			   continue;  //not participating

		   Ccb_halfedge_circulator hebegin=fi->outer_ccb ();
		   Ccb_halfedge_circulator heiterate=hebegin;
		   do{
			   	   
			   if (heiterate->source()->data()<0){  //new vertex
				   Vertex NewVertex;
				   NewVertex.ID=HexMesh.Vertices.size();
				   NewVertex.isHex=(heiterate->source()->data()==-2);
				   HexMesh.Vertices.push_back(NewVertex);
				   heiterate->source()->data()=NewVertex.ID;
			   }

			   if (heiterate->data().ID<0){  //new halfedge
				   Halfedge NewHalfedge;
				   NewHalfedge.ID=HexMesh.Halfedges.size();
				   NewHalfedge.isHex=(heiterate->data().ID==-2);
				   NewHalfedge.Origin=heiterate->source()->data();
				   NewHalfedge.OrigHalfedge=heiterate->data().OrigHalfedge;
				   HexMesh.Vertices[heiterate->source()->data()].AdjHalfedge=NewHalfedge.ID;
				   HexMesh.Halfedges.push_back(NewHalfedge);
				   heiterate->data().ID=NewHalfedge.ID;
			   }
			   heiterate++;
		   }while(heiterate!=hebegin);

		   //now assigning nexts and prevs
		   do{
			   HexMesh.Halfedges[heiterate->data().ID].Next=heiterate->next()->data().ID;
			   HexMesh.Halfedges[heiterate->data().ID].Prev=heiterate->prev()->data().ID;
			   HexMesh.Halfedges[heiterate->data().ID].Twin=heiterate->twin()->data().ID;
			   if (heiterate->twin()->data().ID>=0)
				   HexMesh.Halfedges[heiterate->twin()->data().ID].Twin=heiterate->data().ID;

			   heiterate++;
		   }while (heiterate!=hebegin);
	   }
	   
	   //constructing the actual vertices
	   for (Vertex_iterator vi=FullArr.vertices_begin();vi!=FullArr.vertices_end();vi++){
		   if (vi->data()<0)
			   continue;

		   //finding out barycentric coordinates
		   ENumber BaryValues[3];
		   ENumber Sum=0;
		   for (int i=0;i<3;i++){
			   ETriangle2D t(vi->point(), TriPoints[(i+1)%3], TriPoints[(i+2)%3]);
			   BaryValues[i]=t.area();
			   Sum+=BaryValues[i];
		   }
		   for (int i=0;i<3;i++)
			   BaryValues[i]/=Sum;

		   EPoint3D ENewPosition(0,0,0);
		   for (int i=0;i<3;i++)
			   ENewPosition=ENewPosition+(ETriPoints3D[i]-CGAL::ORIGIN)*BaryValues[i];

		   Point3D NewPosition(to_double(ENewPosition.x()), to_double(ENewPosition.y()), to_double(ENewPosition.z()));
		   HexMesh.Vertices[vi->data()].Coordinates=NewPosition;

		   //DebugLog<<"Creating Vertex "<<vi->data()<<" with 2D coordinates ("<<vi->point().x()<<","<<vi->point().y()<<") "<<" and 3D Coordinates ("<<std::setprecision(10) <<NewPosition.x()<<","<<NewPosition.y()<<","<<NewPosition.z()<<")\n";
	   }

	   for (Face_iterator fi=FullArr.faces_begin();fi!=FullArr.faces_end();fi++){
		   if (!fi->data())
			   continue;
		      
		   int FaceSize=0;
		   Ccb_halfedge_circulator hebegin=fi->outer_ccb ();
		   Ccb_halfedge_circulator heiterate=hebegin;
		   do{ FaceSize++;  heiterate++; }while(heiterate!=hebegin);
		   int CurrPlace=0;
		   
		   Face NewFace;
		   NewFace.ID=HexMesh.Faces.size();
		   NewFace.NumVertices=FaceSize;
		   NewFace.AdjHalfedge=hebegin->data().ID;
		   
		   do{ 
			   NewFace.Vertices[CurrPlace++]=heiterate->source()->data();
			   HexMesh.Halfedges[heiterate->data().ID].AdjFace=NewFace.ID;
			   heiterate++; 
		   }while(heiterate!=hebegin);
		   HexMesh.Faces.push_back(NewFace);
	   }

	}

	//DebugLog.close();

}


struct PointPair{
	int Index1, Index2;
	double Distance;

	PointPair(int i1, int i2, double d):Index1(i1), Index2(i2), Distance(d){}
	~PointPair(){}

	const bool operator<(const PointPair& pp) const {
		if (Distance>pp.Distance) return false;
		if (Distance<pp.Distance) return true;

		if (Index1>pp.Index1) return false;
		if (Index1<pp.Index1) return true;

		if (Index2>pp.Index2) return false;
		if (Index2<pp.Index2) return true;

		return false;
		
	}
};


std::vector<std::pair<int,int>> FindVertexMatch(std::vector<Point3D>& Set1, std::vector<Point3D>& Set2)
{
  std::set<PointPair> PairSet;
	for (int i=0;i<Set1.size();i++)
    for (int j=0;j<Set2.size();j++){
      Vector3D vec =Set1[i]-Set2[j];
      PairSet.insert(PointPair(i,j,hedra::copyleft::cgal::Norm(vec)));
    }
  

	//adding greedily legal connections until graph is full
	std::vector<bool> Set1Connect(Set1.size());
	std::vector<bool> Set2Connect(Set2.size());

	std::vector<std::pair<int, int> > Result;

	for (int i=0;i<Set1.size();i++)
		Set1Connect[i]=false;

	for (int i=0;i<Set2.size();i++)
		Set2Connect[i]=false;

	if (Set1.size()!=Set2.size())
		int kaka=9;

	int NumConnected=0;
	for (std::set<PointPair>::iterator ppi=PairSet.begin();ppi!=PairSet.end();ppi++)
	{
		PointPair CurrPair=*ppi;
		//checking legality - if any of one's former are connected to ones latters or vice versa
		bool FoundConflict=false;
		for (int i=0;i<Result.size();i++){
			if (((Result[i].first>CurrPair.Index1)&&(Result[i].second<CurrPair.Index2))||
				((Result[i].first<CurrPair.Index1)&&(Result[i].second>CurrPair.Index2))){
					FoundConflict=true;
					break;
			}
		}
		if (FoundConflict)
			continue;  

		//otherwise this edge is legal, so add it
		Result.push_back(std::pair<int, int>(CurrPair.Index1, CurrPair.Index2));
		if (!Set1Connect[CurrPair.Index1]) NumConnected++;
		if (!Set2Connect[CurrPair.Index2]) NumConnected++;
		Set1Connect[CurrPair.Index1]=Set2Connect[CurrPair.Index2]=true;
		if (NumConnected==Set1.size()+Set2.size())
			break;  //all nodes are connected
	}

	if (NumConnected!=Set1.size()+Set2.size())
		int kaka=9;

	return Result;

}

struct sortercol
{
	int whichdim;

	int index;
	Point3D Coordinates;

	sortercol(int& w, int& i, Point3D& c):whichdim(w), index(i), Coordinates(c){}
    ~sortercol(){}

    const bool operator<(const sortercol& sc) const
	{
		switch(whichdim){
			case 0:
				if (Coordinates.x()<sc.Coordinates.x()) return true;
				if (Coordinates.x()>sc.Coordinates.x()) return false;

				if (Coordinates.y()<sc.Coordinates.y()) return true;
				if (Coordinates.y()>sc.Coordinates.y()) return false;

				if (Coordinates.z()<sc.Coordinates.z()) return true;
				if (Coordinates.z()>sc.Coordinates.z()) return false;

				return (index<sc.index);  //botyh equal
				break;

			case 1:
				if (Coordinates.y()<sc.Coordinates.y()) return true;
				if (Coordinates.y()>sc.Coordinates.y()) return false;

				if (Coordinates.z()<sc.Coordinates.z()) return true;
				if (Coordinates.z()>sc.Coordinates.z()) return false;

				if (Coordinates.x()<sc.Coordinates.x()) return true;
				if (Coordinates.x()>sc.Coordinates.x()) return false;

				return (index<sc.index);  //botyh equal
				break;

			case 2:
				if (Coordinates.z()<sc.Coordinates.z()) return true;
				if (Coordinates.z()>sc.Coordinates.z()) return false;

				if (Coordinates.x()<sc.Coordinates.x()) return true;
				if (Coordinates.x()>sc.Coordinates.x()) return false;

				if (Coordinates.y()<sc.Coordinates.y()) return true;
				if (Coordinates.y()>sc.Coordinates.y()) return false;

				return (index<sc.index);  //botyh equal
				break;

		}
      
    }
};

struct TwinFinder{
	int index;
	int v1,v2;

	TwinFinder(int i, int vv1, int vv2):index(i), v1(vv1), v2(vv2){}
	~TwinFinder(){}

	const bool operator<(const TwinFinder& tf) const
	{
		if (v1<tf.v1) return false;
		if (v1>tf.v1) return true;

		if (v2<tf.v2) return false;
		if (v2>tf.v2) return true;

		return false;
	}


};

void hedra::copyleft::cgal::Mesh::JoinFace(int heindex)
{
	if (Halfedges[heindex].Twin<0)
		return;  //there is no joining of faces are

	int Face1=Halfedges[heindex].AdjFace;
	int Face2=Halfedges[Halfedges[heindex].Twin].AdjFace;

	DebugLog<<"Merging Faces "<<Face1<<" and "<<Face2<<"\n";
	DebugLog<<"By Edges "<<Halfedges[heindex].Prev<<"->"<<heindex<<"->"<<Halfedges[heindex].Next<<"\n";
	DebugLog<<"And twins "<<Halfedges[Halfedges[heindex].Twin].Prev<<"->"<<Halfedges[heindex].Twin<<"->"<<Halfedges[Halfedges[heindex].Twin].Next<<"\n";
	DebugLog<<"Vertices "<<Halfedges[heindex].Origin<<" and "<<Halfedges[Halfedges[heindex].Next].Origin<<"\n";

	//check if a lonely edge
	if ((Halfedges[heindex].Prev==Halfedges[heindex].Twin)&&(Halfedges[heindex].Next==Halfedges[heindex].Twin)){
		Halfedges[heindex].Valid=Halfedges[Halfedges[heindex].Twin].Valid=false;
		Vertices[Halfedges[heindex].Origin].Valid=Vertices[Halfedges[Halfedges[heindex].Twin].Origin].Valid=false;
		if ((Faces[Halfedges[heindex].AdjFace].AdjHalfedge==heindex)||(Halfedges[Faces[Halfedges[heindex].AdjFace].AdjHalfedge].Twin==heindex)){  //should find another one somehow
			for (int i=0;i<Halfedges.size();i++){
				if ((!Halfedges[i].Valid))
					continue;
				if (Halfedges[i].AdjFace==Halfedges[heindex].AdjFace){
					Faces[Halfedges[heindex].AdjFace].AdjHalfedge=i;
					break;
				}
			}
		}
		return;
	}

	//check if spike edge
	if ((Halfedges[heindex].Prev==Halfedges[heindex].Twin)||(Halfedges[heindex].Next==Halfedges[heindex].Twin)){

		
		int CloseEdge=heindex;
		if (Halfedges[heindex].Prev==Halfedges[heindex].Twin)
			CloseEdge=Halfedges[heindex].Twin;

		Halfedges[CloseEdge].Valid=Halfedges[Halfedges[CloseEdge].Twin].Valid=false;
		Vertices[Halfedges[CloseEdge].Origin].AdjHalfedge=Halfedges[Halfedges[CloseEdge].Twin].Next;
		Faces[Face1].AdjHalfedge=Halfedges[CloseEdge].Prev;

		Halfedges[Halfedges[CloseEdge].Prev].Next=Halfedges[Halfedges[CloseEdge].Twin].Next;
		Halfedges[Halfedges[Halfedges[CloseEdge].Twin].Next].Prev=Halfedges[CloseEdge].Prev;

		Vertices[Halfedges[Halfedges[CloseEdge].Twin].Origin].Valid=false;

		return;
	}

	Faces[Face1].AdjHalfedge=Halfedges[heindex].Next;
	if (Face2!=Face1)
		Faces[Face2].Valid=false;

	//if ((Halfedges[heindex].Origin==123923)||(Halfedges[Halfedges[heindex].Next].Origin==123923))
	//	int kaka=9;

	//Faces[Face2].AdjHalfedge=Halfedges[Halfedges[heindex].Twin].Next;

	Halfedges[heindex].Valid=Halfedges[Halfedges[heindex].Twin].Valid=false;

	Halfedges[Halfedges[heindex].Next].Prev=Halfedges[Halfedges[heindex].Twin].Prev;
	Halfedges[Halfedges[Halfedges[heindex].Twin].Prev].Next=Halfedges[heindex].Next;

	Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Prev=Halfedges[heindex].Prev;
	Halfedges[Halfedges[heindex].Prev].Next=Halfedges[Halfedges[heindex].Twin].Next;

	Vertices[Halfedges[heindex].Origin].AdjHalfedge=Halfedges[Halfedges[heindex].Twin].Next;
	Vertices[Halfedges[Halfedges[heindex].Next].Origin].AdjHalfedge=Halfedges[heindex].Next;

	/*int hebegin=Halfedges[heindex].Next;
	int heiterate=hebegin;
	int IterationNum=0;
	do{
		Halfedges[heiterate].AdjFace=Face1;
		heiterate=Halfedges[heiterate].Next;
		IterationNum++;
		if (IterationNum>100)
			int akka=9;
	}while (heiterate!=hebegin);*/

	//all other floating halfedges should renounce this one
	for (int i=0;i<Halfedges.size();i++)
		if (Halfedges[i].AdjFace==Face2)
			Halfedges[i].AdjFace=Face1;

}


//edge has to be sourcing a 2-valence vertex!!
void hedra::copyleft::cgal::Mesh::UnifyEdges(int heindex)
{	
	//if (Halfedges[heindex].Twin<0)
	//	return;
	if (Halfedges[heindex].Twin<0){
		DebugLog<<"Unifying edge "<<heindex<<" with source "<<Halfedges[heindex].Origin<<"\n";
		DebugLog<<"new source "<<Halfedges[Halfedges[heindex].Prev].Origin<<"\n";

		DebugLog<<"old positions: ("<<Vertices[Halfedges[Halfedges[heindex].Prev].Origin].Coordinates<<")->("<<Vertices[Halfedges[heindex].Origin].Coordinates<<")->("<<Vertices[Halfedges[Halfedges[heindex].Next].Origin].Coordinates<<")\n";
	}
	//adjusting source
	Vertices[Halfedges[heindex].Origin].Valid=false;
	Halfedges[heindex].Origin=Halfedges[Halfedges[heindex].Prev].Origin;
	Vertices[Halfedges[heindex].Origin].AdjHalfedge=heindex;

	Faces[Halfedges[heindex].AdjFace].AdjHalfedge=Halfedges[heindex].Next;

	if (Halfedges[heindex].Twin<0)
		DebugLog<<"Removing edge "<<Halfedges[heindex].Prev<<" and connecting "<<Halfedges[Halfedges[heindex].Prev].Prev<<"->"<<heindex<<"\n";

	//adjusting halfedges
	Halfedges[Halfedges[heindex].Prev].Valid=false;
	Halfedges[heindex].Prev=Halfedges[Halfedges[heindex].Prev].Prev;
	Halfedges[Halfedges[heindex].Prev].Next=heindex;

	if (Halfedges[heindex].Twin<0)
		DebugLog<<"new positions: ("<<Vertices[Halfedges[heindex].Origin].Coordinates<<")->("<<Vertices[Halfedges[Halfedges[heindex].Next].Origin].Coordinates<<")\n";
	

	//adjusting twin, if exists
	if (Halfedges[heindex].Twin>=0){ 
		Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Valid=false;
		Halfedges[Halfedges[heindex].Twin].Next=Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Next;
		Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Prev=Halfedges[heindex].Twin;
		Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].AdjHalfedge=Halfedges[Halfedges[heindex].Twin].Next;
	} 
}

void hedra::copyleft::cgal::Mesh::CleanMesh()
{
	//removing nonvalid vertices
  std::vector<int> TransVertices(Vertices.size());
  std::vector<Vertex> NewVertices;
	for (int i=0;i<Vertices.size();i++){
		if (!Vertices[i].Valid)
			continue;

		Vertex NewVertex=Vertices[i];
		NewVertex.ID=NewVertices.size();
		NewVertices.push_back(NewVertex);
		TransVertices[i]=NewVertex.ID;
	}

	Vertices=NewVertices;
	for (int i=0;i<Halfedges.size();i++)
		Halfedges[i].Origin=TransVertices[Halfedges[i].Origin];

	for (int i=0;i<Faces.size();i++)
		for (int j=0;j<Faces[i].NumVertices;j++)
			Faces[i].Vertices[j]=TransVertices[Faces[i].Vertices[j]];

	//return;

	//removing nonvalid faces
  std::vector<Face> NewFaces;
  std::vector<int> TransFaces(Faces.size());
	for (int i=0;i<Faces.size();i++){
		if (!Faces[i].Valid)
			continue;

		Face NewFace=Faces[i];
		NewFace.ID=NewFaces.size();
		NewFaces.push_back(NewFace);
		TransFaces[i]=NewFace.ID;
	}
	Faces=NewFaces;
	for (int i=0;i<Halfedges.size();i++)
		Halfedges[i].AdjFace=TransFaces[Halfedges[i].AdjFace];

	//removing nonvalid halfedges
  std::vector<Halfedge> NewHalfedges;
  std::vector<int> TransHalfedges(Halfedges.size());
	for (int i=0;i<Halfedges.size();i++){
		if (!Halfedges[i].Valid)
			continue;

		Halfedge NewHalfedge=Halfedges[i];
		NewHalfedge.ID=NewHalfedges.size();
		NewHalfedges.push_back(NewHalfedge);
		TransHalfedges[i]=NewHalfedge.ID;
	}

	Halfedges=NewHalfedges;
	for (int i=0;i<Faces.size();i++)
		Faces[i].AdjHalfedge=TransHalfedges[Faces[i].AdjHalfedge];

	for (int i=0;i<Vertices.size();i++)
		Vertices[i].AdjHalfedge=TransHalfedges[Vertices[i].AdjHalfedge];

	for (int i=0;i<Halfedges.size();i++){
		if (Halfedges[i].Twin!=-1)
			Halfedges[i].Twin=TransHalfedges[Halfedges[i].Twin];
		Halfedges[i].Next=TransHalfedges[Halfedges[i].Next];
		Halfedges[i].Prev=TransHalfedges[Halfedges[i].Prev];
	}

}

bool hedra::copyleft::cgal::Mesh::CheckMesh()
{
	for (int i=0;i<Vertices.size();i++){
		if (!Vertices[i].Valid)
			continue;

		if (!Halfedges[Vertices[i].AdjHalfedge].Valid)
			return false;

		if (Halfedges[Vertices[i].AdjHalfedge].Origin!=i)
			return false;
	}

	for (int i=0;i<Halfedges.size();i++){
		if (!Halfedges[i].Valid)
			continue;

		if (Halfedges[Halfedges[i].Next].Prev!=i)
			return false;

		if (Halfedges[Halfedges[i].Prev].Next!=i)
			return false;

		if (!Vertices[Halfedges[i].Origin].Valid)
			return false;

		if (!Faces[Halfedges[i].AdjFace].Valid)
			return false;

		if (Halfedges[Halfedges[i].Next].Origin==Halfedges[i].Origin)  //a degenerate edge
			return false;

		if (Halfedges[i].Twin>=0)
			if (Halfedges[Halfedges[i].Twin].Twin!=i)
				return false;

		if (Halfedges[i].isHex){  //checking that it is not left alone
			if (Halfedges[i].Prev==Halfedges[i].Twin)
				return false;

			if (Halfedges[i].Next==Halfedges[i].Twin)
				return false;

		}
	}

	for (int i=0;i<Faces.size();i++){
		if (!Faces[i].Valid)
			continue;

		int hebegin=Faces[i].AdjHalfedge;
		int heiterate=hebegin;
		int NumEdges=0;
		do{
			if (!Halfedges[heiterate].Valid)
				return false;

			if (Halfedges[heiterate].AdjFace!=i)
				return false;

			heiterate=Halfedges[heiterate].Next;
			NumEdges++;
			if (NumEdges>1000)
				return false;

		}while (heiterate!=hebegin);
	}
}

void hedra::copyleft::cgal::Mesh::ComputeTwins()
{
	//twinning up edges 
  std::set<TwinFinder> Twinning;
	for (int i=0;i<Halfedges.size();i++){
		if (Halfedges[i].Twin>=0)
			continue;

		std::set<TwinFinder>::iterator Twinit=Twinning.find(TwinFinder(0,Halfedges[Halfedges[i].Next].Origin, Halfedges[i].Origin));
		if (Twinit!=Twinning.end()){
			Halfedges[Twinit->index].Twin=i;
			Halfedges[i].Twin=Twinit->index;
			Twinning.erase(*Twinit);
		} else {
			Twinning.insert(TwinFinder(i,Halfedges[i].Origin,Halfedges[Halfedges[i].Next].Origin));
		}
	}

}

void hedra::copyleft::cgal::Mesh::WalkBoundary(int &CurrEdge)
{
	do{
		CurrEdge=Halfedges[CurrEdge].Next;
		if (Halfedges[CurrEdge].Twin<0)
			break;  //next boundary over a 2-valence vertex
		CurrEdge=Halfedges[CurrEdge].Twin;
	}while (Halfedges[CurrEdge].Twin>=0);

	if (Halfedges[CurrEdge].Twin>=0)
		int kaka=9;
}

//typedef pair< vector<int>, vector<int> > VertexWalk;

typedef adjacency_list <vecS, vecS, undirectedS> Graph;

void hedra::copyleft::cgal::Mesh::TestUnmatchedTwins()
{
	std::vector<int> Untwinned;
	for (int i=0;i<Halfedges.size();i++)
		if ((Halfedges[i].Twin<0)&&(Halfedges[i].Valid))
			Untwinned.push_back(i);

	for (int i=0;i<Untwinned.size();i++){
		for (int j=0;j<Untwinned.size();j++){
			Vector3D diff1=Vertices[Halfedges[Untwinned[i]].Origin].Coordinates-Vertices[Halfedges[Halfedges[Untwinned[j]].Next].Origin].Coordinates;
			Vector3D diff2=Vertices[Halfedges[Untwinned[j]].Origin].Coordinates-Vertices[Halfedges[Halfedges[Untwinned[i]].Next].Origin].Coordinates;
			if ((sqrt(diff1.squared_length())<10e-4)&&(sqrt(diff2.squared_length())<10e-4)){
				DebugLog<<"Halfedge "<<Untwinned[i]<<":("<<Halfedges[Untwinned[i]].Origin<<","<<Halfedges[Halfedges[Untwinned[i]].Next].Origin<<") is untwinned to ";
				DebugLog<<"Halfedge "<<Untwinned[j]<<":("<<Halfedges[Untwinned[j]].Origin<<","<<Halfedges[Halfedges[Untwinned[j]].Next].Origin<<")\n";
				DebugLog<<Vertices[Halfedges[Untwinned[i]].Origin].Coordinates<<" and "<<Vertices[Halfedges[Untwinned[j]].Origin].Coordinates<<"\n";
			}
		}
	}
}

void hedra::copyleft::cgal::Mesh::RemoveEdge(int heindex)
{
	if ((Halfedges[heindex].Twin>=0)&&(Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].NumVertices<=3)){
		RemoveFace(Halfedges[Halfedges[heindex].Twin].AdjFace,Halfedges[heindex].Twin);
		return;
	}

	DebugLog<<"Removing Edge "<<heindex<<" with Twin "<<Halfedges[heindex].Twin<<"\n";
	Halfedges[heindex].Valid=false;
	Halfedges[Halfedges[heindex].Next].Prev=Halfedges[heindex].Prev;
	Halfedges[Halfedges[heindex].Prev].Next=Halfedges[heindex].Next;
	DebugLog<<"Connecting "<<Halfedges[heindex].Prev<<"->"<<Halfedges[heindex].Next<<"\n";
	Vertices[Halfedges[heindex].Origin].AdjHalfedge=Halfedges[heindex].Next;
	DebugLog<<"Vertex "<<Halfedges[heindex].Origin<<" points to "<<Halfedges[heindex].Next<<"\n";
	int LeftVertex=Halfedges[heindex].Origin;
	int RemoveVertex=Halfedges[Halfedges[heindex].Next].Origin;
	DebugLog<<"vertex "<<RemoveVertex<<" is removed"<<"\n";
	Vertices[RemoveVertex].Valid=false;
	Halfedges[Halfedges[heindex].Next].Origin=Halfedges[heindex].Origin;
	DebugLog<<"halfedge "<<Halfedges[heindex].Next<<" has vertex"<<Halfedges[heindex].Origin<<" as origin\n";
	Faces[Halfedges[heindex].AdjFace].AdjHalfedge=Halfedges[heindex].Next;
	Faces[Halfedges[heindex].AdjFace].NumVertices--;
	if (Faces[Halfedges[heindex].AdjFace].NumVertices==0)
		Faces[Halfedges[heindex].AdjFace].Valid=false;

	if (Halfedges[heindex].Twin>=0){
		Halfedges[Halfedges[heindex].Twin].Valid=false;
		Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Prev=Halfedges[Halfedges[heindex].Twin].Prev;
		Halfedges[Halfedges[Halfedges[heindex].Twin].Prev].Next=Halfedges[Halfedges[heindex].Twin].Next;
		DebugLog<<"Connecting "<<Halfedges[Halfedges[heindex].Twin].Prev<<"->"<<Halfedges[Halfedges[heindex].Twin].Next<<"\n";
		//Vertices[Halfedges[Halfedges[heindex].Twin].Origin].AdjHalfedge=Halfedges[Halfedges[heindex].Twin].Next;
		Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Origin=LeftVertex;
		Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].AdjHalfedge=Halfedges[Halfedges[heindex].Twin].Next;
		Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].NumVertices--;
		if (Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].NumVertices==0)
			Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].Valid=false;
	}

	for (int i=0;i<Halfedges.size();i++)
		if (Halfedges[i].Origin==RemoveVertex){
			Halfedges[i].Origin=LeftVertex;
			DebugLog<<"now halfedge "<<i<<" has vertex"<<LeftVertex<<" as origin\n";
		}

}

void hedra::copyleft::cgal::Mesh::RemoveFace(int findex, int heindex)
{
	//return;
	DebugLog<<"Removing Face "<<findex<<" with initial edge "<<heindex<<" and allged # vertices "<<Faces[findex].NumVertices<<"\n";
	
	int LeftVertex=Halfedges[heindex].Origin;
	int hebegin=heindex;
	int heiterate=hebegin;

	DebugLog<<"Leaving Vertex "<<LeftVertex<<"\n";
	Faces[findex].Valid=false;
	std::vector<int> ReplaceOrigins;
	do{
		DebugLog<<"Removing edge part "<<heiterate<<" with twin "<<Halfedges[heiterate].Twin<<" and origin "<<Halfedges[heiterate].Origin<<"\n";
		
		Halfedges[heiterate].Valid=false;
		if (Halfedges[heiterate].AdjFace!=findex)
			int kaka=9;
		if (Halfedges[heiterate].Origin!=LeftVertex){
			Vertices[Halfedges[heiterate].Origin].Valid=false; 
			ReplaceOrigins.push_back(Halfedges[heiterate].Origin);
		}
		if (Halfedges[heiterate].Twin<0){
			heiterate=Halfedges[heiterate].Next;
			continue;
		}
		int ReduceEdge=Halfedges[heiterate].Twin;
		Halfedges[ReduceEdge].Valid=false;
		Halfedges[Halfedges[ReduceEdge].Next].Prev=Halfedges[ReduceEdge].Prev;
		Halfedges[Halfedges[ReduceEdge].Prev].Next=Halfedges[ReduceEdge].Next;

		DebugLog<<"Connecting "<<Halfedges[ReduceEdge].Prev<<"->"<<Halfedges[ReduceEdge].Next<<"\n";
	
		Faces[Halfedges[ReduceEdge].AdjFace].AdjHalfedge=Halfedges[ReduceEdge].Next;
		Vertices[LeftVertex].AdjHalfedge=Halfedges[ReduceEdge].Next;
		DebugLog<<"Vertex "<<LeftVertex<<"("<<TransVertices[LeftVertex]<<") now points to "<<Halfedges[ReduceEdge].Next<<"\n";

		heiterate=Halfedges[heiterate].Next;
	}while (heiterate!=hebegin);

	for (int i=0;i<Halfedges.size();i++)
		for (int j=0;j<ReplaceOrigins.size();j++)
			if (Halfedges[i].Origin==ReplaceOrigins[j]){
				DebugLog<<"Now halfedge "<<i<<" has vertex "<<LeftVertex<<"("<<TransVertices[LeftVertex]<<") as origin instead of "<<ReplaceOrigins[j]<<"("<<TransVertices[ReplaceOrigins[j]]<<")\n";
				Halfedges[i].Origin=LeftVertex;
			}

	//in case there are leftovers
	if (!Halfedges[Vertices[LeftVertex].AdjHalfedge].Valid){
		DebugLog<<"Discarting Vertex "<<LeftVertex<<"("<<TransVertices[LeftVertex]<<")\n";
		Vertices[LeftVertex].Valid=false;
	}
}



void hedra::copyleft::cgal::Mesh::SimplifyHexMesh()
{
	//unifying vertices which are similar

	DebugLog.open("Debugging.txt");

	//CheckMesh();

	int MaxOrigHE=-3276700.0;
	for (int i=0;i<Halfedges.size();i++)
		MaxOrigHE=std::max(MaxOrigHE, Halfedges[i].OrigHalfedge);

	std::vector< std::vector<int> > BoundEdgeCollect1(MaxOrigHE+1);
	std::vector< std::vector<int> > BoundEdgeCollect2(MaxOrigHE+1);
	std::vector<bool> Marked(Halfedges.size());
	for (int i=0;i<Halfedges.size();i++) Marked[i]=false;
	//finding out vertex correspondence along twin edges of the original mesh by walking on boundaries
	for (int i=0;i<Halfedges.size();i++){
		if ((Halfedges[i].OrigHalfedge<0)||(Marked[i]))
			continue;

		//find the next beginning of a boundary
		int PrevOrig;
		int CurrEdge=i;
		do{
			PrevOrig=Halfedges[CurrEdge].OrigHalfedge;
			WalkBoundary(CurrEdge);
		}while(PrevOrig==Halfedges[CurrEdge].OrigHalfedge);

		//filling out strips of boundary with the respective attached original halfedges
		int BeginEdge=CurrEdge;
		std::vector<std::pair<int,int> > CurrEdgeCollect;
		do{
			CurrEdgeCollect.push_back(std::pair<int, int> (Halfedges[CurrEdge].OrigHalfedge, CurrEdge));
			Marked[CurrEdge]=true;
			WalkBoundary(CurrEdge);
		}while (CurrEdge!=BeginEdge);

		PrevOrig=-1000;
		bool In1;
		for (int j=0;j<CurrEdgeCollect.size();j++){
			if (CurrEdgeCollect[j].first!=PrevOrig)
				In1=BoundEdgeCollect1[CurrEdgeCollect[j].first].empty();

			if (In1) 
				BoundEdgeCollect1[CurrEdgeCollect[j].first].push_back(CurrEdgeCollect[j].second);
			else
				BoundEdgeCollect2[CurrEdgeCollect[j].first].push_back(CurrEdgeCollect[j].second);
			PrevOrig=CurrEdgeCollect[j].first;
		}
	}

	//editing the edges into two vector lists per associated original edge

	std::vector< std::vector<int> > VertexSets1(MaxOrigHE+1), VertexSets2(MaxOrigHE+1);
	for (int i=0;i<MaxOrigHE+1;i++){
		for (int j=0;j<BoundEdgeCollect1[i].size();j++)
			VertexSets1[i].push_back(Halfedges[BoundEdgeCollect1[i][j]].Origin);

		if (BoundEdgeCollect1[i].size()>0)
			VertexSets1[i].push_back(Halfedges[Halfedges[BoundEdgeCollect1[i][BoundEdgeCollect1[i].size()-1]].Next].Origin);

		for (int j=0;j<BoundEdgeCollect2[i].size();j++)
			VertexSets2[i].push_back(Halfedges[BoundEdgeCollect2[i][j]].Origin);

		if (BoundEdgeCollect2[i].size()>0)
			VertexSets2[i].push_back(Halfedges[Halfedges[BoundEdgeCollect2[i][BoundEdgeCollect2[i].size()-1]].Next].Origin);

		std::reverse(VertexSets2[i].begin(),VertexSets2[i].end());
	}

	//finding out vertex matches
	std::vector<std::pair<int, int> > VertexMatches;
	for (int i=0;i<MaxOrigHE+1;i++){
		std::vector<Point3D> PointSet1(VertexSets1[i].size());
		std::vector<Point3D> PointSet2(VertexSets2[i].size());
		for (int j=0;j<PointSet1.size();j++)
			PointSet1[j]=Vertices[VertexSets1[i][j]].Coordinates;

		for (int j=0;j<PointSet2.size();j++)
			PointSet2[j]=Vertices[VertexSets2[i][j]].Coordinates;


		std::vector<std::pair<int, int> > CurrMatches;
		if ((!PointSet1.empty())&&(!PointSet2.empty()))
			CurrMatches=FindVertexMatch(PointSet1, PointSet2);

		for (int j=0;j<CurrMatches.size();j++){
			CurrMatches[j].first =VertexSets1[i][CurrMatches[j].first];
			CurrMatches[j].second=VertexSets2[i][CurrMatches[j].second];
		}

		VertexMatches.insert( VertexMatches.end(), CurrMatches.begin(), CurrMatches.end() );
	}

	//finding connected components, and uniting every component into a random single vertex in it (it comes out the last mentioned)

    Graph MatchGraph;
	for (int i=0;i<Vertices.size();i++)
		add_vertex(MatchGraph);
	for (int i=0;i<VertexMatches.size();i++)
		add_edge(VertexMatches[i].first, VertexMatches[i].second, MatchGraph);

	double MaxDist=-327670000.0;
	for (int i=0;i<VertexMatches.size();i++)
		MaxDist=std::max(MaxDist, (Vertices[VertexMatches[i].first].Coordinates-Vertices[VertexMatches[i].second].Coordinates).squared_length());

	//vector<int> TransVertices(Vertices.size());
	TransVertices.resize(Vertices.size());
    int NumNewVertices = connected_components(MatchGraph, &TransVertices[0]);

	/*for (int i=0;i<Vertices.size();i++){
		if (TransVertices[i]==73)
			int kaka=8;
		if (TransVertices[i]==408)
			int kaka=9;
	}*/

	

	//seeing if there are faces that are going to degenerate
	/*vector<int> TransVertices2(NewVertices.size());
	for (int i=0;i<NewVertices.size();i++)
		TransVertices2[i]=i;*/

	//adding other vertex to the degeneration if needed
	bool ThereisChange;
	do{
		ThereisChange=false;
		for (int i=0;i<Halfedges.size();i++){	
			if ((!Halfedges[i].Valid)||(Halfedges[i].Twin>=0))
				continue;

			if (TransVertices[Halfedges[i].Origin]!=TransVertices[Halfedges[Halfedges[i].Next].Origin])
				continue;  //this edge is OK

			if (Faces[Halfedges[i].AdjFace].NumVertices<=3) {	
				for (int j=0;j<Faces[Halfedges[i].AdjFace].NumVertices;j++){
					if (TransVertices[Faces[Halfedges[i].AdjFace].Vertices[j]]!=TransVertices[Halfedges[i].Origin]){
						add_edge(Faces[Halfedges[i].AdjFace].Vertices[j], Halfedges[i].Origin, MatchGraph);
						ThereisChange=true;
					}
				}
			}
		}

		if (ThereisChange){
			TransVertices.clear();	
			TransVertices.resize(Vertices.size());
			NumNewVertices = connected_components(MatchGraph, &TransVertices[0]);
		}
	}while (ThereisChange);
	
	//removing edges (and consequent faces) which will degenerate
	for (int i=0;i<Halfedges.size();i++){
		//break;
		if ((!Halfedges[i].Valid)||(Halfedges[i].Twin>=0))
			continue;
		if (TransVertices[Halfedges[i].Origin]!=TransVertices[Halfedges[Halfedges[i].Next].Origin])
			continue;  //this edge is OK

		if (Faces[Halfedges[i].AdjFace].NumVertices<=3)
			RemoveFace(Halfedges[i].AdjFace, i);
		else
			RemoveEdge(i);

		//CheckMesh();
	}

	std::vector<Vertex> NewVertices(NumNewVertices);
	for (int i=0;i<Vertices.size();i++){  //redundant, but not terrible
		if (!Vertices[i].Valid)
			continue;
		Vertex NewVertex=Vertices[i];
		NewVertex.ID=TransVertices[i];
		NewVertices[TransVertices[i]]=NewVertex;
	}

	
	Vertices=NewVertices;

	for (int i=0;i<Faces.size();i++)
		for (int j=0;j<Faces[i].NumVertices;j++)
			Faces[i].Vertices[j]=TransVertices[Faces[i].Vertices[j]];

	for (int i=0;i<Halfedges.size();i++){
		Halfedges[i].Origin=TransVertices[Halfedges[i].Origin];
		if (Halfedges[i].Valid)
			Vertices[Halfedges[i].Origin].AdjHalfedge=i;
	}

		
	//CheckMesh();

	//twinning up edges 
	std::set<TwinFinder> Twinning;
	for (int i=0;i<Halfedges.size();i++){
		if ((Halfedges[i].Twin>=0)||(!Halfedges[i].Valid))
			continue;

		std::set<TwinFinder>::iterator Twinit=Twinning.find(TwinFinder(0,Halfedges[Halfedges[i].Next].Origin, Halfedges[i].Origin));
		if (Twinit!=Twinning.end()){
			Halfedges[Twinit->index].Twin=i;
			Halfedges[i].Twin=Twinit->index;
			DebugLog<<"Twinning "<<i<<" and "<<Twinit->index<<"\n";
			Twinning.erase(*Twinit);
		} else {
			Twinning.insert(TwinFinder(i,Halfedges[i].Origin,Halfedges[Halfedges[i].Next].Origin));
		}
	}

	//heckMesh();
	TestUnmatchedTwins();

	for (int i=0;i<Halfedges.size();i++){
		if (Halfedges[i].isHex)
			DebugLog<<"Hex edge "<<i<<"\n";
		else
			DebugLog<<"Triangle edge "<<i<<"\n";

		DebugLog<<"Origin: "<<Halfedges[i].Origin<<"\n";
		DebugLog<<"Prev: "<<Halfedges[i].Prev<<"\n";
		DebugLog<<"Next: "<<Halfedges[i].Next<<"\n";
		DebugLog<<"Twin: "<<Halfedges[i].Twin<<"\n";
		DebugLog<<"Face: "<<Halfedges[i].AdjFace<<"\n";
	}

	/*for (int i=0;i<Faces.size();i++)
	{
		if (!Faces[i].Valid){
			Faces[i].NumVertices=0;
			continue;
		}
		int hebegin=Faces[i].AdjHalfedge;
		int heiterate=hebegin;
		Faces[i].NumVertices=0;
		do{
			Faces[i].Vertices[Faces[i].NumVertices]=Halfedges[heiterate].Origin;
			Faces[i].NumVertices++;
			heiterate=Halfedges[heiterate].Next;
		}while (heiterate!=hebegin);
	}

	return;*/

	for (int i=0;i<Halfedges.size();i++){
		//if ((mod(i,1000)==0))
		//	CheckMesh();
		if ((!Halfedges[i].isHex)&&(Halfedges[i].Valid)){
			//if (i==12550)
			//	int kaka=8;
			JoinFace(i);	
		}
	}


	//CheckMesh();

	
	//unifying chains of edges

	//counting valences
	std::vector<int> Valences(Vertices.size());
	for (int i=0;i<Vertices.size();i++)
		Valences[i]=0;

	for (int i=0;i<Halfedges.size();i++){
		if (Halfedges[i].Valid){
			Valences[Halfedges[i].Origin]++;
			//Valences[Halfedges[Halfedges[i].Next].Origin]++;
			if (Halfedges[i].Twin<0)  //should account for the target as well
				Valences[Halfedges[Halfedges[i].Next].Origin]++;
		}
	}

	for (int i=0;i<Valences.size();i++)
		if ((Vertices[i].Valid)&&(Valences[i]<2))
			Vertices[i].Valid=false;

	for (int i=0;i<Vertices.size();i++)
		if ((Vertices[i].Valid)&&(Valences[i]<=2))//&&(Halfedges[Vertices[i].AdjHalfedge].Twin>=0))
			UnifyEdges(Vertices[i].AdjHalfedge);

	//CheckMesh();

	

	//re-assigning vertices on every face according to halfedges
	for (int i=0;i<Faces.size();i++)
	{
		if (!Faces[i].Valid){
			Faces[i].NumVertices=0;
			continue;
		}
		int hebegin=Faces[i].AdjHalfedge;
		int heiterate=hebegin;
		Faces[i].NumVertices=0;
		do{
			Faces[i].Vertices[Faces[i].NumVertices]=Halfedges[heiterate].Origin;
			Faces[i].NumVertices++;
			heiterate=Halfedges[heiterate].Next;
		}while (heiterate!=hebegin);
	}


	//remove non-valid components
	CleanMesh();

	//checking if mesh is valid
	CheckMesh();

	//computing centroids
	for (int i=0;i<Faces.size();i++){
		Faces[i].Centroid=Vector3D(0.0,0.0,0.0);
		for (int j=0;j<Faces[i].NumVertices;j++){
			Faces[i].Centroid=Faces[i].Centroid+(Vertices[Faces[i].Vertices[j]].Coordinates-CGAL::ORIGIN);
		}
		Faces[i].Centroid=Faces[i].Centroid/(double)Faces[i].NumVertices;
	}

}




void hedra::copyleft::cgal::Mesh::CreateStrips()
{
	InStrip.resize(Halfedges.size());
	for (int i=0;i<Halfedges.size();i++)
		InStrip[i]=-1;
	

	int CurrStripNum=0;
	for (int i=0;i<Halfedges.size();i++){
		if (InStrip[i]!=-1)
			continue;

		//building a strip starting from this face
		std::queue<int> CurrHalfedges;
		CurrHalfedges.push(i);
		bool Used=false;
		while (!CurrHalfedges.empty()){
			int CurrHalfedge=CurrHalfedges.front();
			CurrHalfedges.pop();

			if (InStrip[CurrHalfedge]!=-1)
				continue;

			InStrip[CurrHalfedge]=CurrStripNum;

			if (Halfedges[CurrHalfedge].Twin==-1)
				continue;

			InStrip[Halfedges[CurrHalfedge].Twin]=CurrStripNum;
			Used=true;

			//looking at the face of the twin
			int NextFace=Halfedges[Halfedges[CurrHalfedge].Twin].AdjFace;
			if (Faces[NextFace].NumVertices!=6)
				continue;  //a singularity (hopefully, or A BUG) which is the end of this strip

			//otherwise, continue on with the strip
			int NextEdge=Halfedges[Halfedges[Halfedges[Halfedges[CurrHalfedge].Twin].Next].Next].Next;
			CurrHalfedges.push(NextEdge);
			NextEdge=Halfedges[Halfedges[Halfedges[CurrHalfedge].Next].Next].Next;
			CurrHalfedges.push(NextEdge);		
		}

		if (Used)
			CurrStripNum++;
	}

	int NumDualStrips=CurrStripNum;
	VertexChains.resize(NumDualStrips);
	for (int i=0;i<NumDualStrips;i++){
		std::set<int> CurrVertexChain;
		for (int j=0;j<Halfedges.size();j++){
			if (InStrip[j]!=i)
				continue;
			CurrVertexChain.insert(Halfedges[Halfedges[j].Next].Origin);
			CurrVertexChain.insert(Halfedges[Halfedges[j].Prev].Origin);
			CurrVertexChain.insert(Halfedges[Halfedges[Halfedges[j].Next].Next].Origin);
			CurrVertexChain.insert(Halfedges[j].Origin);
		}
		VertexChains[i]=CurrVertexChain;
	}
	
}

