#  IXP Code
Scientific code for surface parameterisation and optimisation of aerofoil using singular value decomposition of an aerofoil training library. Developed as part of the third year Individual Exploratory Project (IXP) for the title of MEng Aerospace Engineering from the University of Bristol.

The final submission of the IXP Report can be found in [GonzalezPerezAlejandro-IXPFinalReport.pdf](https://github.com/AlejandroGlezPerez/ixp-code/blob/master/GonzalezPerezAlejandro-IXPFinalReport.pdf).

## Required libraries
* IXP Code uses the [boost.uBLAS](https://github.com/boostorg/ublas) libraries, distributed as part of the [Boost C++ Libraries](https://www.boost.org/users/download/) for linear algebra.

* The [LAPACK](http://www.netlib.org/lapack/) library is required to perform SVD calculations.

* The (boost numeric bindings)[https://svn.boost.org/svn/boost/sandbox/numeric_bindings/] library provides the required linkage between the uBLAS containers.

* OpenMP is used for parallelisation of the calls to the flow solver. Most compilers provide native support for OpenMP by default. For compilation using LLVM (*clang*), the [OpenMP Library](https://openmp.llvm.org) must be installed separately.

## Additional requirements
### Flow solver
The implementation for the flow solver provided requires the **[rotorsim](#acknowledgments)** flow solver and **[UOBSMB](#acknowledgments)** mesh generator. Other suitable flow solvers and meshers may be used but their implementation is the responsibility of the user.

### Training library
A suitable training library must be provided. In order to run the provided example a folder must be created in the root directory named `library` containing the subdirectories `library/surfaces` and `library/meshes`. Surface files must be added to the `library/surfaces` folder in ascending numerical order using a `.dat` extension (*e.g.* `library/surfaces/0.dat`). The surface file should contain the aerofoil name, number of coordinate points and coordinates in the following format.

```
TITLE = "Eh 1.0-7.0 (DM_Smoothed)"
4
1.0000000000000000	0.0000000000000000
0.9998903416702609	-0.0000039999881933
0.9998903417581405	0.0000038760956449
1.0000000000000000	0.0000000000000000
```

Meshes are added to the `library/meshes` using the same name as their corresponding aerofoils but a `.blk` extension. The format used for storing the mesh information is dependent on the mesh type and compatibility with the chosen flow solver.

## Author
* **Alejandro González Pérez** - [AlejandroGlezPerez](https://github.com/AlejandroGlezPerez)

## Building
Compilation requires a compiler with support for the C++17 standard and new features, including the `<filesystem>` library. Using `clang`, to build the app:
```
$ clang++ -Xpreprocessor -fopenmp -lomp -Iinclude -std=c++17 src/*.cc -O3 -o IXP\ Code
```

## Acknowledgments
* The **rotorsim** flow solver and **UOBSMB** meshers were kindly provided by [Professor Christian B. Allen](mailto:c.b.allen@bristol.ac.uk).

* The training library used was retrieved from https://m-selig.ae.illinois.edu/ads/coord_database.html.
