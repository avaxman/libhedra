// This file is part of libhedra, a library for polygonal mesh processing
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_VISUALIZATION_SCHEMES_H
#define HEDRA_VISUALIZATION_SCHEMES_H

#include <Eigen/Core>


//This file contains the default libdirectional visualization paradigms
namespace hedra
{
  
  //Mesh colors
  Eigen::RowVector3d IGL_INLINE default_mesh_color(){
    return Eigen::RowVector3d::Constant(1.0);
  }
  
  //Color for faces that are selected for editing and constraints
  Eigen::RowVector3d IGL_INLINE selected_face_color(){
    return Eigen::RowVector3d(0.7,0.2,0.2);
  }
  
  Eigen::RowVector3d IGL_INLINE passive_handle_color(){
    return Eigen::RowVector3d(1.0,0.5,0.0);
  }
  
  Eigen::RowVector3d IGL_INLINE active_handle_color(){
    return Eigen::RowVector3d(0.5,1.0,0.0);
  }
  
  Eigen::RowVector3d IGL_INLINE default_edge_color(){
    return Eigen::RowVector3d(0.0,0.0,0.0);
  }
}

#endif
