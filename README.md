GeNuSys
-------------
C++ header-only library for computations with generalized number systems.

License: GPLv3+

Overview of functionality
-------------
#### Linear algebra
* Basic linear algebra (both sparse and dense representations) and polynomials
* Decompositions and standard matrix forms

#### Number systems
* Computing digit sets: j-canonical, j-symmetric, adjoint, dense
* Construction of generalized number systems from a base matrix, digit set and norm (operator norm, frobenius norm or p-norm)
* Computing cycles, orbits and expansions
* Speeding up computations by reducing the volume which needs to be checked for periodic elements (using a heuristic to find integer similarity transforms)
* Simultaneous number systems

Requirements
-------------
* C++11 compatible compiler
* [GMP](https://gmplib.org/)
* CMake (for building tests & examples)

Examples and documentation
-------------
For code examples, see the [examples](examples) subdirectory, and especially [examples/example0.cpp](examples/example0.cpp).

Automatically generated source code reference/overview: [doxygen/html](doxygen/html).

Mathematical background and description of algorithms:
* [Burcsi, P.: Algorithmic aspects of generalized number systems. Eötvös Loránd Tudományegyetem, Budapest, Ph.D. dissertation, 2008, 1–56.](http://www.tnkcs.inf.elte.hu/vedes/Burcsi_Peter_Tezisek_en.pdf)
* [Kovács, A.: Radix expansion in lattices. Eötvös Loránd Tudományegyetem, Budapest, Ph.D.
dissertation, 2001, 1–98.](http://compalg.inf.elte.hu/~attila/pub/phd.ps)
* Matula, D. W.: Basic digit sets for radix representation. Journal of the ACM (JACM). 29, 1982, 1131–1143.

[10_genusys.pdf](10_genusys.pdf) contains a detailed description of the library and the implemented algorithms (in hungarian).

Build options
-------------
To disable threading, define the GENUSYS_NO_THREADING macro before including GeNuSys headers.
