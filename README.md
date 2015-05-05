DOCI Exact
==========

This is an exact solver for double occupied configuration interaction (DOCI).
It builds the hamiltonian using compressed row storage (sparse) and uses
lanczos (or fully diagonalization) to find the lowest eigenvalue and
eigenvector. Using that eigenvector, it than builds the second order reduced
density matrix (rdm) and stores it in a HDF5 file. It only builds the
elements of the rdm that can be non-zero for DOCI.

Build
-----
To build this, you need a C++11 compiler (gcc 4.8 or newer, clang 3.3 or newer),
the HDF5 libraries, the blas and lapack libraries and the arpack library. The
Makefile is quite simple, adjust the compilers and header/libraries as needed
for your system.

Input
-----
The program needs molecular integrals from [PSI4](https://github.com/psi4/psi4public). 
You can extract them using a plugin: https://github.com/wpoely86/doci_sdp-atom/blob/master/extern/mointegrals.cc

Documentation
-------------
Everything is documented with doxygen. Run `make doc` to build the HTML docs.
Or read the docs [online](http://wpoely86.github.io/DOCI-Exact/). If something
is unclear, do not hesitate to contact me.

License
-------
The code is available under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.txt) license.
