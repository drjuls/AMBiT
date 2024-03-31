# Compiling AMBiT

We’ve tested AMBiT and know it works on the Linux and macOS operating
systems. It will probably work on MS Windows or other Unix operating
systems (e.g. FreeBSD, Illumos, etc), but we haven't tried it so we
can't say for sure.

In order to compile AMBiT you’ll need the following software libraries
and tools:

- A C++ compiler with support for C++11, such as GCC, Clang, or the
  Intel C++ compiler.

- [GSL](https://www.gnu.org/software/gsl/) - The GNU Scientific
  Library.

- The [Boost](https://www.boost.org/) C++ libraries.

- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) (v3) -
  C++ linear algebra package.

- [LAPACK](http://www.netlib.org/lapack/) and
  [BLAS](http://www.netlib.org/blas/) - linear algebra subroutines.
  Can be substituted for internal libraries in the (proprietary)
  [Intel Math Kernel Library](https://software.intel.com/en-us/mkl)
  (MKL). Compiling against MKL allows expensive linear-algebra
  operations (including generating angular data) to be automatically
  parallelised at run-time.

- [Google Abseil](https://github.com/abseil/abseil-cpp)

- [CMake](https://cmake.org/) software build tool.

- (Optional) OpenMP and MPI - parallel programming libraries to run
    AMBiT on high-performance computing clusters. These will generally be
    cluster-specific, although OpenMP is supported by default on both GCC
    and the Intel C++ compiler.

Most of these dependencies can be installed from your operating system's
pacakge manager (e.g. apt-get in Ubuntu, dnf in Fedora or Homebrew in
macOS). Some package managers have separate "development" packages, which
you must also install to be able to compile AMBiT against them. These
development packages tend to have names ending with some variation on "dev",
such as `liblapack-dev` on Ubuntu.

It is important that Boost dependencies are built with the same compiler as
the rest of AMBiT. Some package managers (such as Homebrew or the module
systems on HPC clusters) do not guarantee this, which can cause crashes and
unpredictable behaviour at run-time. Consult your system's documentation
for instructions on how to ensure AMBiT and Boost are built with compatible
compilers.

## Compilers
AMBiT makes extensive use of modern C++ features, and requires a compiler which supports at least
C++17. Notably, AMBiT requires a toolchain (C++ compiler and standard library) which supports the
`filesystem` library introduced in the C++17 standard; a table showing the minimum versions of
major compilers which support this feature can be found at this link:
<https://en.cppreference.com/w/cpp/compiler_support/17#C.2B.2B17_library_features>.

The following compilers are known to work with AMBiT:
- GCC with `libstdc++` v8.x or later
- Intel "Classic" C++ compiler (`icpc`, not the newer LLVM-based compiler `icpx`) v2023 or later
- Clang with `libc++` v7.x or later

There is a bug in the newer LLVM-based Intel C++ compiler (`icpx`, released in 2024) which causes the
compiler to crash when compiling AMBiT; use the "Classic" Intel C++ compiler or GCC instead.

## Building with CMake
### Quick start
Do this to build AMBiT in `./build` and install it to `./installed/bin/ambit`:

```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=./installed
cmake --build build
cmake --install build
```

The AMBiT build process is based around a build tool called CMake, which handles both build
configuration and the actual compilation. CMake is a powerful, highly-customisable tool; we provide
a brief overview here, but see the official 
[Kitware CMake FAQ](https://gitlab.kitware.com/cmake/community/-/wikis/FAQ)
and [documentation](https://cmake.org/documentation/) for more detailed instructions.

A CMake build has three conceptual stages:

1. Configuration, including setting compiler options and finding external libraries,
2. Compilation, which can optionally be done in parallel on multi-core machines, 
3. Installation, where CMake copies the final `ambit` binary executable to another directory.

We recommend compiling AMBiT in a dedicated folder, as this keeps everything neat and tidy, and
makes it easy to quickly clean up things like shared objects by just deleting the build folder.
For concreteness, we'll use `build` as the name of the build folder for the examples throughout
this guide.

The first stage, configuration, is the most complex and is CMake's "special sauce". It lets us
specify a choice of compiler (and compiler flags), toggle different modes of parallelism, and set
up (or automatically infer) the paths and options for linking to external libraries. To configure
the build, make sure you're in the base-level AMBiT directory and run:

```bash
cmake -B build
```

This creates a folder called `build` in the current directory and populates it with basic CMake
scaffolding. We can provide additional options to CMake via command line options of the form:

```bash
cmake -B build -D<OPTION>=<VALUE>
```

For example, to build AMBiT with MPI and OpenMP, we need to set `USE_MPI=true` and 
`USE_OPENMP=true` during configuration, we'd do:

```bash
cmake -B build -DUSE_MPI=true -DUSE_OPENMP=true
```

We can leave most build options unspecified and CMake will attempt to automatically infer
sensible defaults, but these inferences can be explicitly overridden if required. It may be
necessary to explicitly specify the paths to some of the libraries required by AMBiT, since CMake
may not be able to automatically find the libraries on all systems (especially if there are
multiple versions available).

CMake can often automatically find and configure the external libraries needed by AMBiT, but may
sometimes struggle on systems with non-standard installation directories (especially on HPC
systems using non-standard environment modules). The CMake configuration step will fail if it can't
find any of the required libraries, in which case you'll need to manually set the paths using one
of the CMake options listed below.

Here is a short list of CMake options for AMBiT:

- `CMAKE_INSTALL_PREFIX`: location to install AMBiT (you must have read/write permissions for this
  directory)
- `CMAKE_CXX_COMPILER`: C++ compiler to use when compiling AMBiT. For MPI-enabled builds, this must
  be manually set to the correct MPI compiler wrappers (e.g. `mpiCC` for OpenMPI)
- `CMAKE_CXX_FLAGS`: Flags to pass to the C++ compiler
- `CMAKE_Fortran_COMPILER` (case sensitive): Fortran compiler to use when compiling Davidson
  eigensolver
- `CMAKE_Fortran_FLAGS` (case sensitive): Flags to pass to the Fortran compiler
- `CMAKE_PREFIX_PATH`: additional paths for CMake configuration files for external libraries. Only
  necessary when CMake's automatic package configuration fails
- `EIGEN_INCLUDE_DIR`: directory containing Eigen header files
- `USE_OPENMP`: toggle whether to use OpenMP multithreading parallelism
- `USE_MPI`: toggle whether to use MPI parallelism
- `USE_MKL`: toggle whether to use MKL for linear algebra operations
- `MKL_INTERFACE_FULL`: MKL library interface to link against (see below for more details)
- `MKL_THREADING`: backend library to use for multithreaded linear algebra operations (see below
  for more details)
- `MKL_MPI`: choice of MPI library to use for MKL (must be the same as used by AMBiT)
- `CMAKE_BUILD_TYPE`: toggle whether to build an optimised (`CMAKE_BUILD_TYPE=Release`, the
  default) or debug (`CMAKE_BUILD_TYPE=Debug`) binary
-`BUILD_TESTING`: (for developers) toggle whether to build unit testing suite (see below for more
  details)

Once we've configured the build, we can compile AMBiT by doing

```bash
cmake --build build
```

Additionally, you can pass the `-j` option to tell `scons` to use
multiple cores when compiling AMBiT. For example, we can use all available CPU cores on a system by
doing:

```bash
cmake --build build -j $(nproc)
```

which will considerably speed up the build process by using compiling multiple files at a time in
parallel.

Finally, we can install AMBiT by running:

```bash
cmake --install build
```

By default, this will try to install AMBiT somewhere like `/usr/bin/ambit`. This is often not
writable on shared systems like HPC clusters, so we usually want to set `CMAKE_INSTALL_PREFIX`
during configuration to some directory we have write-access to. In this case, CMake will install
AMBiT to `${CMAKE_INSTALL_PREFIX}/bin/abit`, so for example:

```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=./installed
cmake --build build
cmake --install build
```

will install AMBiT to `./installed/bin/ambit`.

## Building with MKL
Intel's Math Kernel Library (MKL) can provide substantial speedups for linear algebra operations,
but it has many configuration options and can be difficult to link against. Fortunately, MKL
releases greater than `2021.3` have good support for CMake, so most of the necessary configuration
can by done by setting CMake variables. CMake can usually find MKL automatically, but if it fails
then you may need to manually specify the path to the MKL CMake configuration (i.e. the folder
containing `MKLConfig.cmake`) by setting the `CMAKE_PREFIX_PATH` configuration variable.

It's important to note that MKL supports two different linking interfaces for its dynamic
libraries: 64-bit integers (the default) and 32-bit integers.
These can be set via the `MKL_INTERFACE_FULL` CMake option. Eigen (which AMBiT uses for linear
algebra operations) can use MKL as a backend, but it only works with the 32-bit interface and will
crash at run-time (usually with a segmentation fault) if compiled with the default 64-bit
interface.

As a result, you must explicitly request the 32-bit interface when compiling AMBiT by setting
`MKL_INTERFACE_FULL` in the configuration stage. The exact value this variable must be set to can
be different depending on your configuration, but it's usually one of either
`-DMKL_INTERFACE_FULL=lp64` or `-DMKL_INTERFACE_FULL=intel_lp64` (for Intel compilers), or
`-DMKL_INTERFACE_FULL=gf_lp64` (for GNU). If `MKL_INTERFACE_FULL` is not explicitly set, then
CMake will print the possible options in the configuration output. Consult the documentation for
your version of MKL, as well as for Eigen, if you are unsure:

- [CMake Config for
oneMKL](https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-linux/2024-0/cmake-config-for-onemkl.html),
- [Using the ILP64 Interface vs. LP64 Interface](https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-linux/2024-0/using-the-ilp64-interface-vs-lp64-interface.html),
- [Using Intel MKL from Eigen](https://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html)

Finally, MKL supports using multithreading to do linear algebra operations in parallel, which can
provide substantial speedups when enabled. There are multiple supported backends, but it's
important to select one that's compatible with the compiler you are using to build AMBiT. You can
specify which backend MKL should use by setting the `MKL_THREADING` CMake variable to one of the
following values:

- `sequential`: no multithreading
- `intel_thread`: Intel's *OpenMP* implementation (default)
- `tbb_thread`: Intel's *Threading Building Blocks* library
- `gnu_thread`: The GNU/GCC *OpenMP* implementation

While it is technically possible to "mix and match" libraries, the safest choice is to set
`MKL_THREADING=gnu_thread` if you're using GCC, or `MKL_THREADING=intel_thread` if you're using the
Intel compilers. If no value is supplied then it will default to using
`MKL_THREADING=intel_thread`

### Building the test suite
AMBiT has a set of automated unit and regression tests via 
[GoogleTest](https://github.com/google/googletest) and CMake's inbuilt `ctest` runner program.
Running tests can often catch bugs well before they show up in production, so we we recommend
running the test suite early and often if you're developing AMBiT or adding your own functionality.

To enable tests, set `ENABLE_TESTING=true` when configuring with CMake, then compile as normal -
CMake will automatically download GoogleTest and build the test suite. Once the tests have been
built, change into the build directory (e.g. `./build`) and run `ctest`, e.g.:

```bash
cmake -B build -DENABLE_TESTING=true
cmake --build build
cd build
ctest
```

Note that the unit tests do not currently support MPI parallelism, so you must also compile with
`USE_MPI=no` to build the test suite.

Use `ctest -V` to print "verbose" diagnostic output to the terminal for the whole test suite, 
and `ctest --output-on-failure` to print diagnostics for only those tests which fail. 
Additionally, it can sometimes be useful to only run a subset of the test suite while debugging, 
which can by providing a test or group of tests via `ctest -R` e.g.:

```bash
ctest -R AngularDataTester
```

Will run just the AngularData tests and exit (`ctest -R` also accepts test IDs in the form of
regular expressions).

You can also run `ctest --help` for a full list of possible options when running unit tests.

## Troubleshooting common build errors

### CMake fails to automatically find the path to external dependencies
**Example error output**:
```
CMake Error at src/CMakeLists.txt:16 (find_package):
  Could not find a package configuration file provided by "absl" with any of
  the following names:

    abslConfig.cmake
    absl-config.cmake

  Add the installation prefix of "absl" to CMAKE_PREFIX_PATH or set
  "absl_DIR" to a directory containing one of the above files.  If "absl"
  provides a separate development package or SDK, be sure it has been
  installed.
```
This error occurs when you have installed on of AMBiT's external dependencies in a non-standard
location (e.g. not in `/usr/lib` or similar global locations). 

**Fix**: All of AMBiT's external dependencies have native support for CMake, so ensure the
directories where the problem package is installed has one of the `.cmake` files mentioned in the
error message (they'll usually be installed in something like `/path/to/installation/lib/cmake`)
and modify the `CMAKE_PREFIX_PATH` build option to include this path.

### AMBiT crashes on a segmentation fault when calling into MKL functions

Eigen only works with MKL as a backend via MKL's 32-bit integer interface
(which Intel calls `LP64`) (see: <https://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html>), which
we can control by setting the `-DMKL_INTERFACE_FULL` CMake variable. The exact value required
varies, but it's usually one of either `-DMKL_INTERFACE_FULL=lp64` or
`-DMKL_INTERFACE_FULL=intel_lp64` (for Intel compilers), or `-DMKL_INTERFACE_FULL=gf_lp64` (for GNU
compilers).

### Building with MKL and GCC fails due to missing `libiomp*.so` library
Intel MKL has multiple ways of enabling multithreading for libear algebra operations, which are
controlled by setting the `MKL_THREADING` CMake variable to one of the following choices:
- `sequential`: no multithreading
- `intel_thread`: Intel's *OpenMP* implementation (default)
- `tbb_thread`: Intel's *Threading Building Blocks* library
- `gnu_thread`: The GNU/GCC *OpenMP* implementation

If no value is supplied to `MKL_THREADING`, it will default to using `intel_thread`, which requires
the `libiomp` library to be available. This will usually be fine if you're using the Intel compilers,
but may not be present when compiling AMBiT with GCC so you'll need to set
`MKL_THREADING=gnu_thread`.

### Compilation fails when using Intel LLVM C++ compiler `icpx`
**Example error output:**
```
icpx: error: clang frontend command failed with exit code 139 (use -v to see invocation)
Intel(R) oneAPI DPC++/C++ Compiler 2024.1.0 (2024.1.0.20240308)
Target: x86_64-unknown-linux-gnu
Thread model: posix
InstalledDir: /opt/intel/oneapi/compiler/2024.1/bin/compiler
```

**Fix:** 
There is a bug in the newer LLVM-based Intel C++ compiler (`icpx`, released in 2024) which causes the
compiler to crash when compiling AMBiT; use the "Classic" Intel C++ compiler (`icpc`) or GCC instead.

## Make

We also supply some rudimentary Makefiles for systems where CMake is
unavailable. Plain `make` can't automatically figure out linker paths like
CMake, so these files need to be hand-tailored to your specific system
configuration. The system-specific parameters like compiler, compile flags
and linker paths are set in `make.machine`; once this is set run `make`,
`make debug` or `make test` in the root directory.
