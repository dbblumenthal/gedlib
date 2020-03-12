<!--------------------------------------------------------------------------
-                                                                          -
-   Copyright (C) 2018 by David B. Blumenthal                              -
-                                                                          -
-   This file is part of GEDLIB.                                           -
-                                                                          -
-   GEDLIB is free software: you can redistribute it and/or modify it      -
-   under the terms of the GNU Lesser General Public License as published  -
-   by the Free Software Foundation, either version 3 of the License, or   -
-   (at your option) any later version.                                    -
-                                                                          -
-   GEDLIB is distributed in the hope that it will be useful,              -
-   but WITHOUT ANY WARRANTY; without even the implied warranty of         -
-   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           -
-   GNU Lesser General Public License for more details.                    -
-                                                                          -
-   You should have received a copy of the GNU Lesser General Public       -
-   License along with GEDLIB. If not, see <http://www.gnu.org/licenses/>. -
-                                                                          -
--------------------------------------------------------------------------->

# GEDLIB (1.0)

## 1. About GEDLIB

GEDLIB is a C++ library for (suboptimally) computing edit distances between graphs using various state-of-the-art methods. GEDLIB allows you to build your graphs via the C++ API or load them from [GXL files](http://www.gupro.de/GXL/index.html). Several benchmark datasets are distributed with GEDLIB. For these datasets, GEDLIB provides predefined edit cost functions. You can easily extend GEDLIB by implementing new edit cost functions or new methods for computing the graph edit distance. An extensive Doxygen documentation is availabe [here](https://dbblumenthal.github.io/gedlib/).

## 2. License and Citing

The source code of GEDLIB is distributed under the [GNU Lesser General Public License](https://www.gnu.org/licenses/lgpl-3.0.en.html). If you want to use GEDLIB in a publication, please refer to the following papers:

- D. B. Blumenthal, S. Bougleux, J. Gamper, and L. Brun. &ldquo;GEDLIB: A C++ library for graph edit distance computation&rdquo;, GbRPR 2019, [https://doi.org/10.1007/978-3-030-20081-7_2](https://doi.org/10.1007/978-3-030-20081-7_2)
- D. B. Blumenthal, N. Boria, J. Gamper, S. Bougleux, and L. Brun. &ldquo;Comparing heuristics for graph edit distance computation&rdquo;, VLDB J. 29(1), pp. 419-458, 2020, [https://doi.org/10.1007/s00778-019-00544-1](https://doi.org/10.1007/s00778-019-00544-1)

## 3. Installation under Unix

GEDLIB uses the following external libraries:

- [CMake](https://cmake.org/), for compilation. For installation instructions, see [https://cmake.org/install/](https://cmake.org/install/). 
- [Doxygen](http://www.stack.nl/~dimitri/doxygen/), for creating the documentation. For installation instructions, see [https://www.stack.nl/~dimitri/doxygen/manual/install.html](https://www.stack.nl/~dimitri/doxygen/manual/install.html).
- [OpenMP](http://www.openmp.org/) compatible C++ compiler. Under Linux, OpenMP is supported by default. Under Max OS, please install [libomp](https://formulae.brew.sh/formula/libomp) using [Homebrew](https://brew.sh/).  After installing Homebrew, open a shell and execute `$ brew install libomp`.
- [Gurobi (version 8.01 or higher)](http://www.gurobi.com/), for solving mixed integer and linear programming problems. Gurobi is commercial software, but has a free academic licence. Download the binaries and header files into a directory `<GUROBI_ROOT>` and activate your Gurobi licence as described in the Gurobi documentation. If you cannot obtain a licence for Gurobi, you can install GEDLIB without it. In this case, the methods which use mixed integer or linear programming are not available.
- The following external libraries are distributed with GEDLIB:
    - [Eigen (version 3.3.4)](http://eigen.tuxfamily.org/index.php?title=Main_Page)
    - [NOMAD (version 3.8)](https://www.gerad.ca/nomad/Project/Home.html)
    - [LSAPE (version 5)](https://bougleux.users.greyc.fr/lsape/)
    - [LIBSVM (version 3.22)](https://www.csie.ntu.edu.tw/~cjlin/libsvm/)
    - [FANN (version 2.2.0)](http://leenissen.dk/fann/wp/)
    - [Boost (version 1.69.0)](https://www.boost.org/). Does not require compilation.


After having installed CMake, Doxygen, OpenMP (under Mac OS), and, possibly, Gurobi, execute the script `install.py` for installing GEDLIB and the external libraries distributed with GEDLIB:

```sh
python install.py [--help] [-h] [--doc] [--tests <arg>] [--gurobi <GUROBI_ROOT>] [--debug] [--clean] [--lib gxl|<indentifier>,<UserNodeID>,<UserNodeLabel>,<UserEdgeLabel>]
```

If you execute `install.py` without any arguments, only the external libraries distributed with GEDLIB are installed.

- `--help`, `-h`: Show help.  
- `--doc`: Build the [Doxygen documentation](https://dbblumenthal.github.io/gedlib/).
- `--tests <arg>`: Build test executables. Use the option  `--help` to display possible values of `<arg>`.
- `--gurobi <GUROBI_ROOT>`: Tell GEDLIB where to find your Gurobi installation.
- `--debug`: Build in debug mode. 
- `--clean`: Delete the build directoy and update the makefile before the build.
- `--lib gxl`: Build the shared library `lib/libgxlgedlib.so` for usage with graphs given in the [GXL file format](http://www.gupro.de/GXL/index.html).
- `--lib <indentifier>,<UserNodeID>,<UserNodeLabel>,<UserEdgeLabel>`: Build the shared library `lib/lib<indentifier>gedlib.so` for graphs with custom node ID, node label, and edge label types. For example, executing `$ python install.py --lib mytypes,int,double,double` builds the shared library `lib/libmytypesgedlib.so` for usage with graphs whose node IDs are of type `int` and whose node and edge labels are of type `double`.

## 4. Building an Application that Uses GEDLIB

### 4.1 Building an Application that Uses GEDLIB as Header-Only Library

For building an application that uses GEDLIB as a header-only library, it suffices to execute the installation script `install.py` without any options. Subsequently, carry out the following steps:

- Add the following directories to your include directories:
    - `<GEDLIB_ROOT>`
    - `<GEDLIB_ROOT>/ext/boost.1.69.0`
    - `<GEDLIB_ROOT>/ext/eigen.3.3.4/Eigen`
    - `<GEDLIB_ROOT>/ext/nomad.3.8.1/src`
    - `<GEDLIB_ROOT>/ext/nomad.3.8.1/ext/sgtelib/src`
    - `<GEDLIB_ROOT>/ext/lsape.5/include`
    - `<GEDLIB_ROOT>/ext/libsvm.3.22` 
    - `<GEDLIB_ROOT>/ext/fann.2.2.0/include`
    - If you want to install GEDLIB with Gurobi, additionally add the following include directory:
        - Under MacOS: `<GUROBI_ROOT>/mac64/include`
        - Under Linux: `<GUROBI_ROOT>/linux64/include`
- Add the following directories to your link directories:
    - `<GEDLIB_ROOT>/ext/nomad.3.8.1/lib`
    - `<GEDLIB_ROOT>/ext/libsvm.3.22`
    - `<GEDLIB_ROOT>/ext/fann.2.2.0/lib`
    - If you want to install GEDLIB with Gurobi, additionally add the following link directory:
        - Under MacOS: `<GUROBI_ROOT>/mac64/lib`
        - Under Linux: `<GUROBI_ROOT>/linux64/lib`
- Link your application against the following shared libraries:
    - `<GEDLIB_ROOT>/ext/fann.2.2.0/lib/libdoublefann.2.dylib`
    - `<GEDLIB_ROOT>/ext/libsvm.3.22/libsvm.so`
    - `<GEDLIB_ROOT>/ext/nomad.3.8.1/libnomad.so`
    - If you want to install GEDLIB with Gurobi, additionally link your application against the following libraries:
        - Under MacOS: `<GUROBI_ROOT>/mac64/lib/libgurobi80.so` and `<GUROBI_ROOT>/mac64/lib/libgurobi_c++.a`
        - Under Linux: `<GUROBI_ROOT>/linux64/lib/libgurobi80.so` and `<GUROBI_ROOT>/linux64/lib/libgurobi_c++.a`
- If you want to install GEDLIB with Gurobi, define `GUROBI` before including the header `src/env/ged_env.hpp`.

### 4.2 Building an Application that Uses GEDLIB as Shared Library

#### 4.2.1 Graphs Given as GXL Files

If you want to build an application that uses GEDLIB as a shared for graphs given as GXL files, make sure that you have installed GEDLIB with the option `--lib gxl`. Subsequently, carry out the folowing steps:

- Add the following directories to your include directories:
    - `<GEDLIB_ROOT>`
    - `<GEDLIB_ROOT>/ext/boost.1.69.0`
    - `<GEDLIB_ROOT>/ext/eigen.3.3.4/Eigen`
    - `<GEDLIB_ROOT>/ext/nomad.3.8.1/src`
    - `<GEDLIB_ROOT>/ext/nomad.3.8.1/ext/sgtelib/src`
    - `<GEDLIB_ROOT>/ext/lsape.5/include`
    - `<GEDLIB_ROOT>/ext/libsvm.3.22` 
    - `<GEDLIB_ROOT>/ext/fann.2.2.0/include`
    - If you want to install GEDLIB with Gurobi, additionally add the following include directory:
        - Under MacOS: `<GUROBI_ROOT>/mac64/include`
        - Under Linux: `<GUROBI_ROOT>/linux64/include`
- Add the directory `<GEDLIB_ROOT>/lib` to your link directories.
- Link your application against the shared library `libgxlgedlib.so`.

#### 4.2.2 Graphs with User-Defined Node ID, Node Label, and Edge Label Types

If you want to build an application that uses GEDLIB as a shared library for graphs with custom node ID, node label, and edge label types, make sure that you have installed GEDLIB with the option `--lib <indentifier>,<UserNodeID>,<UserNodeLabel>,<UserEdgeLabel>`. Subsequently, carry out the folowing steps:

- Add the directory `<GEDLIB_ROOT>/lib` to your link directories.
- Link your application against the shared library `lib<indentifier>gedlib.so`.

## 5. Using GEDLIB

### 5.1 General Usage

In your source file, include the header `src/env/ged_env.hpp` and create your environment object:

```cpp
#include "src/env/ged_env.hpp"
ged::GEDEnv<UserNodeID, UserNodeLabel, UserEdgeLabel> env;
```

All functionality of GEDLIB is accessible via the environment class `ged::GEDEnv`, whose template parameters `UserNodeID`, `UserNodeLabel`, and `UserEdgeLabel` must be set to the types of the node IDs, the node labels, and the edge labels of your graphs. Use the member functions `ged::GEDEnv::add_graph()`, `ged::GEDEnv::add_node()`, and `ged::GEDEnv::add_edge()`, and `ged::GEDEnv::set_edit_costs()` for adding graphs and edit costs to the environment. Once you're done with this, call `ged::GEDEnv::init()` to initialize the environment.

When you have initialized the environment, you are ready to run the GED methods. To select the method, call `ged::GEDEnv::set_method()`. You can then initialize the selected method via `ged::GEDEnv::init_method()`, and run it via `ged::GEDEnv::run_method()`. Note that most method can be run also without initialization, but doing so typically has a negative impact on both runtime and accuracy. Use the member functions `ged::GEDEnv::get_lower_bound()`, `ged::GEDEnv::get_upper_bound()`, `ged::GEDEnv::get_node_map()`, `ged::GEDEnv::get_runtime()`, `ged::GEDEnv::get_init_time()`, `ged::GEDEnv::get_graph_class()`, and `ged::GEDEnv::get_graph_name()` for accessing the results of your GED computations. Use `ged::GEDEnv::quasimetric_costs()` to check if your edit costs are quasimetric on the graphs contained in your environment.
    
### 5.2 Additional Functionality for Graphs Given as GXL Files
    
If you use GEDLIB with the template parameter `UserNodeID` set to `ged::GXLNodeID` a.k.a. `std::string` and the template parameters `UserNodeLabel` and `UserEdgeLabel` set to `ged::GXLLabel` a.k.a. `std::map<std::string,std::string>`, GEDLIB offers additional functionality for loading graphs given in the GXL file format. For those graphs, you do not have to use the member functions `ged::GEDEnv::add_graph()`, `ged::GEDEnv::add_node()`, and `ged::GEDEnv::add_edge()` for adding graphs to your environment. Instead, you can simply load all of them at once by one call to `ged::GEDEnv::load_gxl_graphs()`.

### 5.3 Using GEDLIB as a Shared Library

#### 5.3.1 Graphs Given as GXL Files

If you want to use GEDLIB as a shared library for graphs given as GXL files, make sure that you have installed GEDLIB with the option `--lib gxl`. In your source file, you have to define `GXL_GEDLIB_SHARED` before including `src/env/ged_env.hpp`:

```cpp
  #define GXL_GEDLIB_SHARED  
  #include "src/env/ged_env.hpp"  
  ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env  
  // ...
```

#### 5.3.2 Graphs with User-Defined Node ID, Node Label, and Edge Label Types
  
If you want to use GEDLIB as a shared library for graphs with custom node ID, node label, and edge label types, make sure that you have installed GEDLIB with the option `--lib <indentifier>,<UserNodeID>,<UserNodeLabel>,<UserEdgeLabel>`. In your source file, you have to define `<IDENTIFIER>_GEDLIB_SHARED` before including `src/env/ged_env.hpp`, where `<IDENTIFIER>` is the upper case transformation of `<identifier>`. For example, if you have installed GEDLIB by executing `$ python install.py --lib mytypes,int,double,double`, you have to do the following:

```cpp
  #define MYTYPES_GEDLIB_SHARED  
  #include "src/env/ged_env.hpp"  
  ged::GEDEnv<int, double, double> env  
  // ...
```

### 5.4 Examples

For an extensively commented example of how to use GEDLIB, have a look at `median/tests/median_letter_demo.cpp`. For more examples, see the `.cpp` files contained in the sub-directories of the `tests/` directory. 

## 6. Reproducability Packages

GEDLIB has been used for several research papers. For reproducing the experiments reported in these papers, follow the instructions below.

##### D. B. Blumenthal, S. Bougleux, J. Gamper, and L. Brun. &ldquo;Ring based approximation of graph edit distance&rdquo;, S+SSPR 2018, vol. 11004 of LNCS, pp. 293-303, [https://doi.org/10.1007/978-3-319-97785-0_28](https://doi.org/10.1007/978-3-319-97785-0_28)

In order to reproduce the experiments reported in this paper, install GEDLIB with the option `--tests sspr2018`. After installation, open a shell and execute the following commands:

```sh
$ cd <GEDLIB_ROOT>/tests/sspr2018/bin
$ ./learn_ring_params
$ ./learn_subgraph_depths
$ ./learn_walks_depth
$ ./test_lsape_based_methods
```

After having executed these commands, the results of the experiments are contained in the folder `tests/sspr2018/output/`.

##### D. B. Blumenthal, N. Boria, J. Gamper, S. Bougleux, and L. Brun. &ldquo;Comparing heuristics for graph edit distance computation&rdquo;, VLDB J. 29(1), pp. 419-458, 2020, [https://doi.org/10.1007/s00778-019-00544-1](https://doi.org/10.1007/s00778-019-00544-1)

In order to reproduce the experiments reported in this paper, install GEDLIB with the options `--tests vldbj2020` and `--gurobi <GUROBI_ROOT>`. After installation, open a shell and execute the following commands:

```sh
$ cd <GEDLIB_ROOT>/tests/vldbj2020/bin
$ ./vldbj_train_subgraph
$ ./vldbj_train_walks
$ ./vldbj_train_ring
$ ./vldbj_train_ml
$ ./vldbj_test_lsape_based_methods
$ ./vldbj_test_lp_based_methods
$ ./vldbj_test_ls_based_methods
$ ./vldbj_test_misc_methods
$ ./vldbj_test_best_methods
```

After having executed these commands, the results of the experiments are contained in the folder `tests/vldbj2020/results/`. For creating TikZ figures and tables that visualize the results, run the script `process_results.py`.

# 7. Datasets

GEDLIB comes with several datassets which contain graphs given in the [GXL file format](http://www.gupro.de/GXL/index.html). They are contained in the following subdirectories of the directory `data/datasets/`:

- <b>`AIDS`, `Fingerprint`, `GREC`, `Letter`, `Mutagenicity`, `Protein`:</b> These datasets are taken from the [IAM Graph Database](http://www.fki.inf.unibe.ch/databases/iam-graph-database). 
  You can use them for scientific work, but are requested to include the following reference to your paper:
    - K. Riesen, H. Bunke:
      &ldquo;IAM graph database repository for graph based pattern recognition and machine learning&rdquo;,
      [https://doi.org/10.1007/978-3-540-89689-0\_33](https://doi.org/10.1007/978-3-540-89689-0_33)
- <b>`CMU-GED`:</b> This datasets is taken from the [Graph Data Repository for Graph Edit Distance](http://www.rfai.li.univ-tours.fr/PublicData/GDR4GED/home.html). 
  You can use it for scientific work, but are requested to include the following reference to your paper:
    - Z. Abu-Aisheh, R. Raveaux, J.-Y. Ramel:
      &ldquo;A graph database repository and performance evaluation metrics for graph edit distance&rdquo;,
      [https://doi.org/10.1007/978-3-319-18224-7\_14](https://doi.org/10.1007/978-3-319-18224-7_14)
- <b>`acyclic`, `alkane`, `mao`, `pah`:</b> These datasets are taken from [GREYC's Chemistry Dataset](https://brunl01.users.greyc.fr/CHEMISTRY/).
- <b> `S-MOL`:</b> Synthetically generated graphs with varying number of node labels whose structure is similar to the structure of `pah` graphs.
- <b> `S-MOL-5`:</b> Synthetically generated graphs with 5 node labels whose structure is similar to the structure of `pah` graphs. 

For each dataset, the directory `data/collections/` contains an XML file which lists the contained graphs' GXL files along with their classes. These files match the document type definition `data/collections/GraphCollection.dtd` and can hence be used as input for `ged::GEDEnv::load_gxl_graphs()`. The Python script `data/collections/sample.py` can be used to generate samples of the datasets.

## 8. Directory Structure

After executing `install.py`, the directoy `<GEDLIB_ROOT>` has the following internal structure:

```
.
├── CMakeLists.txt --------------- used for building GEDLIB
├── doxyfile.in ------------------ used for creating the documentation
├── LICENSE.md ------------------- a copy of GNU LGPL
├── README.md -------------------- this file
├── install.py ------------------- installs GEDLIB
├── _data  
|   ├── _datasets ---------------- contains several datasets
|   |   └── ... 
|   └── _collections ------------- contains XML files that list the graphs in the datasets
|       ├── GraphCollection.dtd -- document type definition for graph collections 
|       ├── sample.py ------------ generates a sample from a collection file
|       └── ... 
├── _ext ------------------------- contains external libraries
|   └── ...
├── _docs ------------------------ contains documentation if GEDLIB has been installed with option --doc
|   └── ...
├── _lib ------------------------- contains shared libraries if GEDLIB has been installed with option --lib
|   └── ...
├── _src ------------------------- contains the sources of GEDLIB  
|   ├── CMakeLists.txt ----------- used for building GEDLIB
|   ├── _edit_costs -------------- contains edit costs for several datasets
|   |   └── ... 
|   ├── _env --------------------- contains the architecture of GEDLIB
|   |   ├── ged_env.hpp ---------- include this header in your application
|   |   └── ... 
|   ├── _methods ----------------- contains various methods for computing GED 
|   |   └── ... 
|   └── _util -------------------- contains utility classes and functions used by the GED methods
|       └── ... 
├── _median ---------------------- contains median graph computation, clustering, and indexing
|   ├── CMakeLists.txt ----------- used for building the executables
|   ├── _bin --------------------- contains the executables if GEDLIB has been built with option --median
|   |   └── ...
|   ├── _collections ------------- contains graph collections used by the median graph computation
|   |   └── ...
|   ├── _output ------------------ contains the median graphs once the executables have been run
|   |   └── ...
|   ├── _src --------------------- contains the median graphs once the executables have been run
|   |   └── ...
|   └── _tests ------------------- contains tests
|       └── ...
└── _tests ----------------------- contains various tests
    └── ...
```
