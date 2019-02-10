// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2019 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_COPYLEFT_CGAL_BASIC_DEFINTIONS_H
#define HEDRA_COPYLEFT_CGAL_BASIC_DEFINTIONS_H
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
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <vector>


namespace hedra
{
  namespace copyleft
  {
    namespace cgal
    {

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
      
      void normalize(Vector3D& v){
        v=v/sqrt(squared_distance(Point3D(0,0,0)+v,Point3D(0,0,0)));
      }
      
      inline Number norm(const Vector3D &v){return sqrt(squared_distance(Point3D(0,0,0)+v,Point3D(0,0,0)));}
      
    }
  }
}

#endif
  
  
