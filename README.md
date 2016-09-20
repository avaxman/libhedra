# libhedra

libhedra is a library that targets the processing of polygonal and polyhedral meshes, that are not necessarily triangular. It supports the generation, analysis, editing, and optimization of such meshes.

##Installation

libhedra is a header-only library, building on [libigl](http://libigl.github.io/libigl/) and consequently [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). They are bundled as submodules within libhedra. As such, minimal installation is needed. Some function would require extra dependencies, and they will be mentioned in context.

To get the library, use:

```bash
git clone --recursive https://github.com/avaxman/libhedra.git
```

to compile the examples, go into the respective library (e.g., "examples/basic") and enter:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

Using the library then amounts to including the relevant header files in the `include` directory.

##Mesh representation

[libigl](http://libigl.github.io/libigl/) represents triangle meshes by a pair of matrices: a matrix `V` that stores vertex coordinates, and a matrix ``F`` that stores oriented triangles as vertex indices into ``v``. Polygonal meshes are characterized by an arbitrary number of vertices per face, and so the mesh representation in libhedra provides an extension of this compact representation with two matrices and a vector:

```cpp
Eigen::MatrixXd V;
Eigen::VectorXi D;
Eigen::MatrixXi F;
```

Where `V` is a #V by 3 coordinates matrix, `D` is an #F by 1 vector that stores the degree (amount of vertices or edges) in a face, and `F` is a #F by `max(D)` array of oriented polygonal faces. In every row of `F`, i.e. `F.row(i)`, only the first `D(i)` elements are occupied. The rest are "don't care" and are usually set to `-1` to avoid errors.

This representation has a few advantages:

1. Triangle meshes are represented the same as in libigl, ignoring D (which is then a uniform array of 3's), allowing for full compatibility.

2. Uniform meshes are represented compactly.

3. Avoiding dynamic allocations inside array is memory friendly.

4. Simple and intuitive indexing.

The disadvantage is when a mesh contains relatively a few high-degree faces accompanied by many low-degree ones. Then, `max(D)` is high, and a major part of the memory occupied by `F` is not used. But this type of meshes are rare in practice.

###Mesh Combinatorics

While `V` and `F` represent the raw information about a polygonal mesh, combinatorial information is often required. It can be obtained with the following function:

```cpp
#include <hedra/hedra_edge_topology>

hedra::hedra_edge_topology(D, F, EV, FE, EF, EFi, FEs,InnerEdges);
```

where:

| Name                     | Description                                                                         |
| :----------------------- | :---------------------------------------------------------------------------------- |
| `EV`            | $\left|E\right| \times 2$ of edge endpoints (into `V`)                                             |
| `FE`               | $\left|F\right| \times max(D)$ consecutive edges per face, where `EV.row(FE(i,j))` contains `F(i,j)` as one of its endpoints.|
| `EF`              | $\left|E\right| \times 2$ of indices into `F` of the adjacent faces.|
| `EFi`            | The relative location of each edge $i$ in faces `EF.row(i)`. That is, `FE(EF(i,0),EFi(i,0))=i` and `FE(EF(i,1),EFi(i,1))=i`. In case of a boundary edge, only the former holds (also for EF).
| `FEs`             | `FEs(i,j)` holds the sign of edge `FE(i,j)` in face `i`. That is, if the edge is oriented positively or negatively relative to the face. |
| `innerEdges`             | a vector of indices into `EV` of inner (non-boundary) edges. |

Example: TBD.




##Loading and Visualization

Meshes can be loaded from OFF files, which have a similar data structure, with the following function:

```cpp
hedra_read_OFF(filename, V, D, F);
```

libhedra builds upon libigl viewer for visualization and interaction. libigl viewer currently only supports triangle meshes, and therefore we visualize our meshes by triangulating the faces using `hedra::triangulate_mesh`. Visualization can be worked out through the following code snippet from example "basic":

```cpp
hedra::hedra_read_OFF(DATA_PATH "/rhombitruncated_cubeoctahedron.off", V, D, F);
hedra::triangulate_mesh(D, F,T,TF);
hedra::hedra_edge_topology(D,F, EV,FE, EF, EFi, FEs, innerEdges);
igl::viewer::Viewer viewer;
viewer.data.clear();
viewer.data.set_mesh(V, T);
```

``T`` is the resulting triangle mesh, and `TF` is a vector relating every triangle to the face it tesselates. In this manner, one can easily adress the set of triangles belonging to a single face for purposes of visualization (say, a uniform color).

Nevertheless, the default edges in the overlay of libigl will show the triangulated edges. To see the polygonal edges and nothing else, the following code from "basic" can be used:

```cpp
Eigen::MatrixXd origEdgeColors(EV.rows(),3);
origEdgeColors.col(0)=Eigen::VectorXd::Constant(EV.rows(),0.0);
origEdgeColors.col(1)=Eigen::VectorXd::Constant(EV.rows(),0.0);
origEdgeColors.col(2)=Eigen::VectorXd::Constant(EV.rows(),0.0);
viewer.data.clear();
viewer.data.set_mesh(V, T);
viewer.data.set_edges(V,EV,OrigEdgeColors);
```
###Augmenting the viewer
libhedra adds some extra functionality to the libigl viewer with the following functions:

| Function   | Description       |
| :----------------------- | :---------------------------------------------------------------------------------- |
| `Scalar2RGB`            | Converts a value within $\left[0,1\right]$ into the cool/warm color map [#moreland_2009][], which can be fed to `igl::set_colors`. |               
| `point_spheres`               | creates spheres with configurable radius, resolution, and color that can be used, e.g., for visualizing deformation handles.|
| `line_cylinders`               | creates cylinders that can be used to visualize vectors and lines. |

`line_cylinders` and `point_spheres` create new meshes, and these need to be concatenated to a given mesh in order to be visualized, as the following code (from "basic"):

```cpp
hedra::point_spheres(bc, sphereRadius, sphereGreens, 10, false, sphereV, sphereT, sphereTC);
Eigen::MatrixXd bigV(V.rows()+sphereV.rows(),3);
Eigen::MatrixXi bigT(T.rows()+sphereT.rows(),3);
if (sphereV.rows()!=0){
    bigV<<V, sphereV;
    bigT<<T, sphereT+Eigen::MatrixXi::Constant(sphereT.rows(), sphereT.cols(), V.rows());
        bigTC<<TC, sphereTC;
} else{
    bigV<<V;
    bigT<<T;
    bigTC<<TC;
}
    
viewer.data.clear();
viewer.data.set_mesh(bigV,bigT);

```

**Note**: `sphereT` indices are relative to `sphereV`, and therefore need to be adjusted to indices in `bigV` before concatenation.

Example: TBD


    
##Evaluation

TBD

##Modelling with Affine Maps

TBD

[#moreland_2009]: Kenneth Moreland. [Diverging Color Maps for Scientific Visualization](http://www.kennethmoreland.com/color-maps).






