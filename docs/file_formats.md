libhedra file formats
===================


### Data structures used from libigl:

- [.dmat](./dmat) uncompressed ASCII/binary files for dense matrices
- [.off](http://wias-berlin.de/software/tetgen/fformats.off.html) Geomview's polyhedral file format
- [.obj](http://en.wikipedia.org/wiki/Wavefront_.obj_file#File_format) Wavefront object file format. Usually unsafe to assume anything more than vertex positions and triangle indices are supported
- [.png](https://en.wikipedia.org/wiki/Portable_Network_Graphics) Portable Network Graphics image file. IGLLIB (in the libiglpng extra) supports png image files via the [yimg](https://github.com/yig/yimg) library. Alpha channels and compression are supported.


