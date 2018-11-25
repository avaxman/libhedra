// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_WILLMORE_ENERGY_H
#define HEDRA_WILLMORE_ENERGY_H
#include <igl/igl_inline.h>
#include <hedra/dual_mesh.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
  // returns a mesh after vertex insertion, which creates a mesh by connected midedges inside each face and around each vertex
  // Inputs:
  //  V  eigen double matrix     #V by 3 - vertex coordinates
  //  D  eigen int vector        #F by 1 - face degrees
  //  F  eigen int matrix        #F by max(D) - vertex indices in face
  //  FE eign int matrix         #F by max(D) - edges by order in face
  
  // Outputs:
  //  newV  eigen double matrix  new vertices
  //  newD  eigen int vector    new valences
  //  newF eigen int matrix     new faces
  IGL_INLINE bool willmore_energy(const Eigen::MatrixXd& V,
                                  const Eigen::VectorXi& VH,
                                  const Eigen::VectorXi& HV,
                                  const Eigen::VectorXi& HE,
                                  const Eigen::VectorXi& HF,
                                  const Eigen::VectorXi& twinH,
                                  const Eigen::VectorXi& nextH,
                                  const Eigen::VectorXi& prevH,
                                  Eigen::VectorXd& W)
  {
    using namespace Eigen;
    using namespace std;
    W.conservativeResize(V.rows()); W.setZero();
    VectorXd MR(V.rows()); MR.setZero();
    for (int i=0;i<V.rows();i++){
      int beginH=VH(i);
      int currH=beginH;
      bool isBoundaryVertex=true;
      
      
      while ((twinH(currH)!=-1)){
        currH=nextH(twinH(currH));
        if (currH==beginH) {isBoundaryVertex=false; break;}
      }
      
      beginH=currH;
      
      vector<int> vertexStarList;
      MatrixXd vertexStar;
      do{
        vertexStarList.push_back(HV(nextH(currH)));
        if(twinH(prevH(currH))==-1){  //last edge on the boundary should be accounted for
          vertexStarList.push_back(HV(prevH(currH)));
        }
        currH=twinH(prevH(currH));
      }while ((beginH!=currH)&&(currH!=-1));
      
      vertexStar.conservativeResize(vertexStarList.size(),4);
      for (int j=0;j<vertexStarList.size();j++)
        vertexStar.row(j)<<0.0,V.row(vertexStarList[j]);
      
      //computing tangent polygon
      RowVector4d qu; qu<<0.0, V.row(i);
      MatrixXd tangentPolygon(vertexStar.rows(),4);
      for (int j=0;j<vertexStar.rows();j++)
        tangentPolygon.row(j)<<QInv(vertexStar.row(j)-qu);
      
      //finding minimum vertex which is by definition convex
      /*int minXIndex, maxXIndex;
      double minX = tangentPolygon.col(0).minCoeff(&minXIndex);
      double maxX = tangentPolygon.col(0).maxCoeff(&maxXIndex);
      int minYIndex, maxYIndex;
      double minY = tangentPolygon.col(1).minCoeff(&minYIndex);
      double maxY = tangentPolygon.col(1).maxCoeff(&maxYIndex);
      int minZIndex, maxZIndex;
      double minZ = tangentPolygon.col(2).minCoeff(&minZIndex);
      double maxZ = tangentPolygon.col(2).maxCoeff(&maxZIndex);
      RowVector3d span; span<<maxX-minX, maxY-minY, maxZ-minZ;
      RowVector3i indices; indices<<minXIndex, minYIndex, minZIndex;
      
      int domDim;
      span.maxCoeff(&domDim);
      int startIndex=indices(domDim);
      double angleSum=0.0;
      RowVector4d CR0=QMult(tangentPolygon.row((startIndex+2)%vertexStar.rows())-tangentPolygon.row((startIndex+1)%vertexStar.rows()), QInv(tangentPolygon.row(startIndex)-tangentPolygon.row((startIndex+1)%vertexStar.rows())));
      double CR0angle, CR0length;
      RowVector3d CR0vec;
      factorize_quaternion(CR0, CR0length, CR0angle, CR0vec);
      for (int j=0;j<vertexStar.rows();j++){
        int v0=(j+startIndex)%vertexStar.rows();
        int v1=(j+1+startIndex)%vertexStar.rows();
        int v2=(j+2+startIndex)%vertexStar.rows();
        RowVector4d CR=QMult(tangentPolygon.row(v2)-tangentPolygon.row(v1), QInv(tangentPolygon.row(v0)-tangentPolygon.row(v1)));
        double CRangle, CRlength;
        RowVector3d CRvec;
        factorize_quaternion(CR, CRlength, CRangle, CRvec);
        angleSum+=(CRvec.tail(3).dot(CR0vec)>=0 ? CRangle : M_PI-CRangle);
      }*/
      
      //W(i)=abs(M_PI-angleSum/((double)vertexStar.rows()-2.0));
      VectorXi Dtp(1); Dtp(0)=vertexStar.rows();
      RowVectorXi Ftp(vertexStar.rows());
      for (int j=0;j<vertexStar.rows();j++)
        Ftp(j)=j;
      VectorXd planarity, regularity;
      hedra::planarity(tangentPolygon.block(0,1,tangentPolygon.rows(),3),Dtp, Ftp, planarity);
      hedra::regularity(tangentPolygon.block(0,1,tangentPolygon.rows(),3),Dtp, Ftp, regularity);
      W(i)=planarity(0);
      MR(i)=planarity(0);
    }
    
    //cout<<"W: "<<W<<endl;
    
    return true;
    
  }
}


#endif


