# ParBLiSS kmerind 
[![Apache 2.0 License](https://img.shields.io/badge/license-Apache%20v2.0-blue.svg)](LICENSE)
[![Build Status](https://img.shields.io/travis/ParBLiSS/kmerind.svg)](https://travis-ci.org/ParBLiSS/kmerind)
[![Build Status](https://travis-ci.org/ParBLiSS/kmerind.svg?branch=master)](https://travis-ci.org/ParBLiSS/kmerind)
[![Test Coverage](https://img.shields.io/codecov/c/github/ParBLiSS/kmerind.svg)](http://codecov.io/github/ParBLiSS/kmerind?branch=master)

`Kmerind` is a library in the **Par**allel **B**ioinformatics **Li**brary for **S**hort **S**equences project (ParBLiSS).

Kmerind provides **k**-**mer** **ind**exing capability for biological sequence data.

Please take a look at our [Wiki](https://github.com/ParBLiSS/kmerind/wiki).

## Overview
ParBLiSS is a C++ library for distributed and multi-core bioinformatics algorithms.  It requires C++ 11 features and MPI (OpenMP not required).  The library is implemented as a set of templated classes.  As such, most of the code is in header form, and are incorporated into the user code via `#include`.  

K-merind provides basic parallel sequence file access and k-mer index construction and query.  Currently, it supports indices for frequency, position, and quality of kmers from short reads and whole genomes, using FASTQ and FASTA formats.

## Build

### Dependencies

Required:

- c++11 supporting compiler
-- `g++` (version 4.8.1+ due to "decltype" and other c++11 features) or
-- `icpc` (version 16+ due to constexpr functions and initializers) or
-- `clang` (version 3.5+ - cmake generated make file has problems with prior versions. or 3.7+ if openmp is used)
- `cmake` (version 2.8+)
- an MPI implementation, one of the following
-- `openmpi` (version 1.7+ due to use of MPI_IN_PLACE)
-- `mpich2` (version 1.5 +)
-- `mvapich` (tested with version 2.1.7)
-- `intel mpi library` (poorly tested)

See 
http://en.cppreference.com/w/cpp/compiler_support


Optional libraries are:

`boost_log`, `boost_system`, `boost_thread`, `boost_program-options`

These are only needed if you intend to turn on boost log engine.


Optional tools include:

- `ccmake` (for graphical cmake configuration)
- `perl`, and perl packages `Term::ANSIColor`, `Getopt::ArgvFile`, `Getopt::Long`, `Regexp::Common` (for g++ error message formatting)


### Getting the source

```sh
git clone https://github.com/ParBLiSS/kmerind.git
cd kmerind
git submodule init
git submodule update
```

### Configuring for build

```sh
mkdir kmerind-build
cd kmerind-build
cmake ../kmerind
```

alternatively, instead of `cmake ../kmerind`, you can use

```sh
ccmake ../kmerind
```


### CMake Parameters

The following are important parameters:

- `CMAKE_BUILD_TYPE`:  defaults to `Release`.
- `ENABLE_TESTING`:  `On` allows `BUILD_TEST_APPLICATIONS` to show, which enables building the test applications
- `BUILD_EXAMPLE_APPLICATION`: `On` allows applications in the `examples` directory to be built
- `LOG_ENGINE`: chooses which log engine to use.
- `LOGGER_VERBOSITY`: chooses the type of messages to prin.

- `ENABLE_SANITIZER`: turns on g++'s address or thread sanitizer.  Use `SANITIZER_STYLE` to configure.  This is for debugging
- `ENABLE_STLFILT`:  turns on g++ error message post processing to make them human readable.  Control verbosity via `STLFIL_VERBOSITY`

It is highly recommended that `ccmake` be used until you've become familiar with the available CMake options and how to specify them on the commandlinie.


### Compiling

```sh
make
```

Important for developers using Intel Compilers, please see the "Intel Compiler Specific Issues" section at the end of the document. 


### Running the tests

```sh
ctest -T Test
```

or 

```sh
make test
```


### Building the documentation

```sh
make doc
```

### Using Bliss
Please see  [Wiki](https://github.com/ParBLiSS/kmerind/wiki).



## Configuring Eclipse:

Cmake typically uses a out-of-source build.  to generate eclipse compatible `.project` and `.cproject` files, supply
```
-G"Eclipse CDT4 - Unix Makefiles"
```
to cmake.

Recommend that `ptp`, `egit`, and `cmake ed` also be installed.


### Intel Compiler Specific Issues
With Intel C Compiler (icc) version 15, the following compilation error is observed:
```sh
Internal error: assertion failed at: "shared/cfe/edgcpfe/il.c", line 18295
```

While there is very little information to be found on the internet related to this error, we have theorized that this is a compiler bug related to auto type deduction in templated function instantiation.  It appears that ICC is unable to auto deduce the data type and size of a statically sized array of the form
```sh
datatype x[len]
```
which as a function parameter is specified as
```sh
datatype (&x)[len]
```
with datatype and len being template parameters for the function.

This error appears only for bitgroup_ops.hpp.  Attempts to replicate the error in a separate test code was not successful.  The workaround is to fully specify the template parameters for the function so to avoid automatic type deduction in this case.

It is not clear if other function parameter forms also cause this error.
 
