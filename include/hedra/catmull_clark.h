// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_CATMULL_CLARK_H
#define HEDRA_CATMULL_CLARK_H
#include <igl/igl_inline.h>
#include <hedra/polygonal_face_centers.h>
#include <hedra/dual_mesh.h>
#include <hedra/linear_subdivision_basics.h>
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
  IGL_INLINE bool catmull_clark(const Eigen::MatrixXd& V,
                                const Eigen::VectorXi& D,
                                const Eigen::MatrixXi& F,
                                const Eigen::MatrixXi& EV,
                                const Eigen::MatrixXi& EF,
                                const Eigen::MatrixXi& EFi,
                                const Eigen::MatrixXi& FE,
                                const Eigen::VectorXi& innerEdges,
                                const hedra::SubdivisionType st
                                Eigen::MatrixXd& newV,
                                Eigen::VectorXi& newD,
                                Eigen::MatrixXi& newF)
  
  {
    using namespace Eigen;
    using namespace std;
    
    switch (st){
      case hedra::LINEAR_SUBDIVISION: vertex_insertion(V, D,F, EV, EF, EFi, FE, innerEdges, &blendCCEdgePoint, &computeDualFacePoint, &computeCCVertexEdgePoint,&blendCCEdgePoint, &computeThreePoints, &computeCCBoundaryVertex, newV, newD, newF);
      case hedra::CANONICAL_MOEBIUS_SUBDIVISION: vertex_insertion(level, deformV[level], &computeMoebiusFourPoints, &computeDualMoebiusFacePoint, &computeCCMoebiusVertexEdgePoint,&computeMoebiusFourPoints, &computeMoebiusThreePoints, &computeCCMoebiusBoundaryVertex);
        
      default: return false
    }
    
    return true;
  }
  
  
  IGL_INLINE bool vertex_insertion(const Eigen::MatrixXd& V,
                                     const Eigen::VectorXi& D,
                                     const Eigen::MatrixXi& F,
                                     const Eigen::MatrixXi& EV,
                                     const Eigen::MatrixXi& EF,
                                     const Eigen::MatrixXi& EFi,
                                     const Eigen::MatrixXi& FE,
                                     const Eigen::VectorXi& innerEdges,
                                     const std::function<RowVector3d(const RowVector3d&, const RowVector3d&, const RowVector3d&, const RowVector3d&)>& edgeBlendFunction,
                                     const std::function<RowVector3d(const MatrixXd&, const MatrixXd&, const VectorXd&, MatrixXd&)>& facePointFunction,
                                     const std::function<RowVector3d(const RowVector3d&, const MatrixXd&, const bool, const RowVector3d&, const double, const MatrixXd&, MatrixXd&)>& vertexEdgePointFunction,
                                     const std::function<MatrixXd(const RowVector3d&, const RowVector3d&, const RowVector3d&, const RowVector3d&)>& boundaryEdgePoint,
                                     const std::function<MatrixXd(const RowVector3d&, const RowVector3d&, const RowVector3d&)>& earExtrapolationPoint,
                                     const std::function<MatrixXd(const RowVector3d&, const RowVector3d&, const RowVector3d&, const RowVector3d&, const RowVector3d&)>& boundaryVertexPoint,
                                     Eigen::MatrixXd& newV,
                                     Eigen::VectorXi& newD,
                                     Eigen::MatrixXi& newF)
  {
    
    //creating 1-ring lists and completing them
    vector<MatrixXd> vertexStars(V.rows());
    vector<vector<int> > intVertexStars(V.rows());
    vector<MatrixXd> canonicalEmbeddings(V.rows());
    vector<MatrixXd> oppositeQuadPoints(V.rows());  //the points opposite the quads, where the point is always after the completeVertexStars of same quad along the halfedge
    
    //the vertex star ring edge and face indices
    vector<VectorXi> ringEdgeIndices(V.rows());
    vector<VectorXi> ringFaceIndices(V.rows());
    
    //the points computed in the respective canonical forms of the vertices, to be blended
    vector<MatrixXd> ringEdgePoints(V.rows());
    vector<MatrixXd> ringFacePoints(V.rows());
    
    vector<MatrixXd> candidateFacePoints(F[level].rows());
    
    MatrixXd tangSphereCenters(V.rows(),3);
    VectorXd tangRadii(V.rows());
    
    vector<bool> isBoundaryVertex(V.rows());
    
    MatrixXi faceVertex2ringFaceIndices;
    vector<VectorXi> ringFace2faceVertexIndices;
    
    ComputeCanonicalRings(level, V, D[level], F[level], vertexStars, ringEdgeIndices, ringFaceIndices,faceVertex2ringFaceIndices, ringFace2faceVertexIndices, tangSphereCenters, tangRadii, isBoundaryVertex);
    
    canonicalEmbeddings.resize(V.rows());
    for (int i=0;i<V.rows();i++)
      canonicalEmbeddings[i]=original2canonical(V.row(i), tangSphereCenters.row(i),tangRadii(i), vertexStars[i]);
    
    //face points
    MatrixXd facePoints(F[level].rows(),3);
    for (int i=0;i<F[level].rows();i++){
      MatrixXd faceVertices(D[level](i),3);
      MatrixXd faceCircles(D[level](i),3);
      VectorXd faceRadii(D[level](i));
      for (int j=0;j<D[level](i);j++){
        faceVertices.row(j)=V.row(F[level](i,j));
        faceCircles.row(j)=tangSphereCenters.row(F[level](i,j));
        faceRadii(j)=tangRadii(F[level](i,j));
      }
      
      facePoints.row(i)=facePointFunction(faceVertices, faceCircles, faceRadii, candidateFacePoints[i]);
      //cout<<"candidateFacePoints[i]: "<<candidateFacePoints[i]<<endl;
    }
    
    //assigning face points to rings
    for (int i=0;i<V.rows();i++){
      ringFacePoints[i].conservativeResize(ringFaceIndices[i].size(),3);
      for (int j=0;j<ringFaceIndices[i].size();j++)
        ringFacePoints[i].row(j)= facePoints.row(ringFaceIndices[i](j)); //candidateFacePoints[ringFaceIndices[i](j)].row(ringFace2faceVertexIndices[i](j));  // ;
    }
    
    
    //computing vertex points
    MatrixXd vertexPoints(V.rows(),3);
    for (int i=0;i<V.rows();i++)
      vertexPoints.row(i)=vertexEdgePointFunction(V.row(i),  vertexStars[i], isBoundaryVertex[i], tangSphereCenters.row(i),tangRadii(i),ringFacePoints[i], ringEdgePoints[i]);
    
    //computing edge points
    MatrixXd edgePoints(EV[level].rows(),3);
    MatrixXd edgeCandidatePoints(2*EV[level].rows(),3);
    
    //allocating edge candidates
    for (int i=0;i<V.rows();i++){
      for (int j=0;j<ringEdgeIndices[i].size();j++){
        int currEdge=ringEdgeIndices[i](j);
        for (int k=0;k<2;k++){
          if (EV[level](currEdge,k)==i)
            edgeCandidatePoints.row(2*currEdge+k)=ringEdgePoints[i].row(j);
        }
      }
    }
    
    for (int i=0;i<EV[level].rows();i++)
      edgePoints.row(i)=edgeBlendFunction(V.row(EV[level](i,0)),edgeCandidatePoints.row(2*i), edgeCandidatePoints.row(2*i+1), V.row(EV[level](i,1)));
    
    
    //special treatment for boundary vertices and edges
    
    //boundary edges
    for (int i=0;i<EH[level].rows();i++){
      int currH;
      if (EH[level](i,0)==-1)
        currH=EH[level](i,1);
      else currH=EH[level](i,0);
      int ix2=HV[level](currH);
      int ix3=HV[level](nextH[level](currH));
      
      RowVector3d x2; x2<<V.row(ix2);
      RowVector3d x3; x3<<V.row(ix3);
      
      if (twinH[level](currH)!=-1)
        continue;
      
      bool x1Ear=(twinH[level](prevH[level](currH))==-1);
      bool x4Ear=(twinH[level](nextH[level](currH))==-1);
      RowVector3d x1, x4;
      
      if (!x1Ear){
        int otherBoundH=currH;
        do{
          otherBoundH=twinH[level](prevH[level](otherBoundH));
        }while(twinH[level](prevH[level](otherBoundH))!=-1);
        int ix1=HV[level](prevH[level](otherBoundH));
        x1=V.row(ix1);
      }
      
      if (!x4Ear){
        int otherBoundH=currH;
        do{
          otherBoundH=twinH[level](nextH[level](otherBoundH));
        }while(twinH[level](nextH[level](otherBoundH))!=-1);
        int ix4=HV[level](nextH[level](nextH[level](otherBoundH)));
        x4=V.row(ix4);
      }
      
      if (x1Ear)
        x1=earExtrapolationPoint(x4, x3, x2);
      if (x4Ear)
        x4=earExtrapolationPoint(x1, x2, x3);
      
      
      edgePoints.row(i)=boundaryEdgePoint(x1,x2,x3,x4);
    }
    
    //boundary vertices
    for (int i=0;i<EH[level].rows();i++){
      int currH;
      if (EH[level](i,0)==-1)
        currH=EH[level](i,1);
      else currH=EH[level](i,0);
      
      if (twinH[level](currH)!=-1)
        continue;
      
      int ix2=HV[level](currH);
      int ix3=HV[level](nextH[level](currH));
      RowVector3d x2; x2<<V.row(ix2);
      RowVector3d x3; x3<<V.row(ix3);
      bool x1Ear=(twinH[level](prevH[level](currH))==-1);
      RowVector3d x1;
      if (!x1Ear){
        int otherBoundH=currH;
        do{
          otherBoundH=twinH[level](prevH[level](otherBoundH));
        }while(twinH[level](prevH[level](otherBoundH))!=-1);
        int ix1=HV[level](prevH[level](otherBoundH));
        x1=V.row(ix1);
        
        vertexPoints.row(ix2)=boundaryVertexPoint(x1,edgePoints.row(HE[level](prevH[level](otherBoundH))), x2,edgePoints.row(i), x3);
      } else {
        vertexPoints.row(ix2) = V.row(ix2);
      }
    }
    
    
    cout<<"(vertexPoints-V).lpNorm<Infinity>()" <<(vertexPoints-V).lpNorm<Infinity>()<<endl;
    
    
    MatrixXd VNew(vertexPoints.rows()+facePoints.rows()+edgePoints.rows(),3);
    VNew<<vertexPoints, edgePoints, facePoints;
    int numNewFaces=D[level].sum();
    D[level+1]=VectorXi::Constant(numNewFaces, 4);
    //VectorXi faceIndices; igl::colon(0,facePoints.rows()-1,faceIndices);
    F[level+1].conservativeResize(numNewFaces,4);
    int currNewFace=0;
    for (int i=0;i<D[level].rows();i++)
      for (int j=0;j<D[level](i);j++)
        F[level+1].row(currNewFace++)<<F[level](i,j),
        V.rows()+FE[level](i,j),
        V.rows()+edgePoints.rows()+i,
        V.rows()+FE[level](i,(j+D[level](i)-1)%D[level](i));
    
    /*F[level+1]<<F[level].col(0), V.rows()+FE[level].col(0).array(), V.rows()+edgePoints.rows()+faceIndices.array(), V.rows()+FE[level].col(3).array(),
     F[level].col(1), V.rows()+FE[level].col(1).array(), V.rows()+edgePoints.rows()+faceIndices.array(), V.rows()+FE[level].col(0).array(),
     F[level].col(2), V.rows()+FE[level].col(2).array(), V.rows()+edgePoints.rows()+faceIndices.array(), V.rows()+FE[level].col(1).array(),
     F[level].col(3), V.rows()+FE[level].col(3).array(), V.rows()+edgePoints.rows()+faceIndices.array(), V.rows()+FE[level].col(2).array();*/
    return VNew;
  }
}


#endif


