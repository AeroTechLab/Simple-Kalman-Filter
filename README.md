# Simple-Kalman-Filter

Basic implementation of a Kalman filter for control state estimation and sensor fusion.

The algorithm is based on [the guide published at bzarg.com](http://www.bzarg.com/p/how-a-kalman-filter-works-in-pictures/#mjx-eqn-kalupdatefull)

### Building

This filter implementation depends on [Simple Matrix](https://github.com/LabDin/Simple-Matrix) project, which is added as a [git submodule](https://git-scm.com/docs/git-submodule).

For building this library e.g. with [GCC](https://gcc.gnu.org/) as a shared object, using reference **BLAS/LAPACK**, the following shell command (from root directory) would be required:

>$ gcc kalman_filters.c matrix/matrix.c -I. -I./matrix -shared -fPIC -o kalman_filters.so -lblas -llapack

### Documentation

Descriptions of how the functions and data structures work are available at the [Doxygen](http://www.stack.nl/~dimitri/doxygen/index.html)-generated [documentation pages](https://labdin.github.io/Simple-Kalman-Filter/kalman__filters_8h.html)

