# libhedra - Polygonal and Polyhedral Mesh Processing

<!---[![Build Status](https://travis-ci.org/libigl/libigl.svg?branch=master)](https://travis-ci.org/libigl/libigl)
[![Build status](https://ci.appveyor.com/api/projects/status/mf3t9rnhco0vhly8/branch/master?svg=true)](https://ci.appveyor.com/project/danielepanozzo/libigl-6hjk1/branch/master)
![](libigl-teaser.png)!--->

<https://github.com/avaxman/libhedra/>

libhedra is a C++ library for processing polygonal and polyhedral meshes, built as an extension to [libigl](https://www.github.com/libigl/libigl) on the foundations of [Eigen](http://eigen.tuxfamily.org/). Polygonal meshes are meshes whose faces are not necessarily triangular. Polyhedral meshes are polygonal meshes whose faces are additionally flat. These meshes are extensively in use in architectural and industrial geometry, and are important theoretical objects in discrete differential geometry.

## Installation
libhedra is a header-only library where each file generally includes one function. To use the library, simply add the _include_ directory to your include path and make sure Directional and its prerequisites are set up properly. After that you can include any files you need normally, using for example `#include <libhedra/affine_maps_deform.h>`.

To get the library, simply clone the repository using:
```git
git clone --recursive https://github.com/avaxman/libhedra.git
```

## Features
The current version is 1.0, comprising the following features:

1. Visualization of polygonal meshes with properties like normals and degrees.
2. Evaluation of measurements like face planarity, concyclity, and regularity.
3. Piecewise-affine handle-based deformation.
4. Complex piecewise-M\"obius deformations.
5. Polygonal subdivision meshes.

libhedra is **a header-only library**. You do not need to compile anything to use,
just include directional headers (e.g. `#include <libhedra/affine_maps_deform.h>`) and run.  Each
header file contains a single function (e.g. `hedra/dual_mesh.h` contains
`hedra::dual_mesh()`). 


## Tutorial
A [Tutorial](https://avaxman.github.io/libhedra/tutorial/) that walks through the entire functionality of libhedra is available. To compile it, go to the `tutorial` folder, open a shell and call:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```
This should properly set up the tutorial project, with the individual chapters as subprojects, and create project makefiles upon which you can build it using your favourite compiler. For windows, you should use `cmake-gui ..` and follow the instructions to create a compilable Visual Studio file.

## Coding Guidelines and Tips

libhedra inherits and follows the strict coding guidelines of libigl: please take a look [here](http://libigl.github.io/libigl/style-guidelines) before submitting your pull requests.


## How to Contribute

If you are interested in joining development, please fork the repository and submit a [pull request](https://help.github.com/articles/using-pull-requests/) with your changes.

## License
Directional is primarily [MPL2](http://www.mozilla.org/MPL/2.0/) licensed ([FAQ](http://www.mozilla.org/MPL/2.0/FAQ.html)). Some files contain third-party code under other licenses. 

## Attribution

If you use libhedra in your academic projects, please cite the implemented papers appropriately. To cite the library in general, you could use this BibTeX entry:

```bibtex
@misc{libhedra,
title = {{libhedra}: geometric processing and optimization of polygonal meshes,
author = {Amir Vaxman and others},
note = {https://github.com/avaxman/libhedra},
year = {2017},
}
```

If you would like to suggest further topics, would like to collaborate in implementation, complain about bugs or ask questions, please address [Amir Vaxman] (mailto:avaxman@gmail.com) (or open an issue in the repository)

## Contact

Directional is led by [Amir Vaxman](http://www.staff.science.uu.nl/~vaxma001/). Please [contact me](mailto:avaxman@gmail.com) if you have questions or comments. For troubleshooting, please post an [issue](https://github.com/avaxman/Directional/issues) on github.

If you're using libirectional in your projects, quickly [drop me a note](mailto:avaxman@gmail.com). Tell me who you are and what you're using it for. This helps justify spending time maintaining this library!

## Future Plans

The following functionality is planned for libhedra:

* Parallel and offset meshes, including evaluation of the Steiner formula for discrete curvature.
* Constrained optimization using augmented Lagrangians.
* Conformal Mesh Deformations with M&ouml;bius Transformations: integrating the working demo [MoebiusCode](https://github.com/avaxman/MoebiusCode) which already relies on libhedra.
* Polyhedral patterns parametrization and optimization.
* Polyhedral patterns
* Canonical M\"{o}bius subdivision
If you would like to suggest further topics, would like to collaborate in implementation, complain about bugs or ask questions, please address [Amir Vaxman](avaxman@gmail.com) (or open an issue in the repository).

## Copyright
2017 Amir Vaxman and others.

Please see individual files for appropriate copyright notices.
