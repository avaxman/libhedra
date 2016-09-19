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

##Loading and Visualization

TBD

##Evaluation

TBD

##Modelling with Affine Maps

TBD






