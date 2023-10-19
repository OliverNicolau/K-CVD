# The order K Farthest Color Voronoi Diagram - 2023

Doctoral studies at [USI Lugano](http://inf.usi.ch): Implementation of Geometrical Algorithms -
The order K Farthest Color Voronoi Diagram (supervised by Prof.ssa Evanthia
Papadopoulou).

## Description

This are various ipelets to compute various kinds of Vornoi Diagrams.
Most of them are developed from existing code bases (CGAL libraries and demos), expanding them in various ways (see forked repository).

This repository includes the work in progress on the higher order Color Vornoi Diagram, or k-CVD.

## Prerequisites

- CMake: version >= 3.12 (can work with older if aproppiate tweaks are made)
- CGAL: version ~ 5.5.1-2 (sable release in debian 12, also available through macports. TODO test & update to new release)
- IPE: version >= 7 (stable release in debian 12)

--------------------------------------------------------------------------------

## Installation

### Download/clone this repository

Use your favourite git tools, or download the zip version from this github page.

### CMake + make

Generate the makefile using CMake, and run make:

```shell
cd /path_to_repository/kCVD_FCVD_HVD_Bisectors
cmake -DCMAKE_BUILD_TYPE="Release" .
make
```

### Install ipelets

Move the .so and .lua files into your IPE ipelets folder. An example for Mac computers:

```shell
cp libCGAL_bisectors.so /path_to_ipe-X.X.XX/Ipe.app/Contents/Resources/ipelets/libCGAL_bisectors.so
cp libCGAL_bisectors.lua /path_to_ipe-X.X.XX/Ipe.app/Contents/Resources/ipelets/libCGAL_bisectors.so
```

--------------------------------------------------------------------------------

## Notes

If in an Macbook with a "Apple silicon" (M1, M2 and other arm64 architecture CPU), you might have needed to compile Ipe from source.

In any case that you need to use a local version of IPE, not installed with a package manager (apt, brew, macports...), you will need to specify the relevant paths to the IPE library.

If in MacOS, you can use the following commands instead of the ones in the CMake section:

```shell
cd /path_to_repository/kCVD_FCVD_HVD_Bisectors
cmake -DIPE_INCLUDE_DIR="/path_to_ipe-X.X.XX/src/include" -DIPE_LIBRARIES="/path_to_ipe-X.X.XX/build/Ipe.app/Contents/Frameworks/libipe.X.X.XX.dylib -DCMAKE_BUILD_TYPE="Release" ."
make
```

If in a Linux system, set apropiately the IPE_INCLUDE_DIR and IPE_LIBRARIES variables.

- ```IPE_INCLUDE_DIR``` : should be the path to the directory containing ```ipebase.h```

- ```IPE_LIBRARIES``` : should be the path to the dynamic library file, usually ```libipe.X.X.XX.so``` (Linux) or ```libipe.X.X.XX.dylib``` (Macos)


--------------------------------------------------------------------------------

## Ipelet instructions

The k-CVD ipelet expects polygons as input (one per color), it then promts for the desired order k.

TODO: allow for detection of IPE colors as input instead of polygons.