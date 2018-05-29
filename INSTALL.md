## Compiling AMBiT
We’ve tested AMBiT and know it works on the Linux and macOS operating
systems. It will probably work on MS Windows or other Unix operating
systems (e.g. FreeBSD, Illumos, etc), but we haven't tried it so we 
can't say for sure.

In order to compile AMBiT you’ll need the following software libraries
and tools:

-   A C++ compiler with support for C++11, such as GCC, Clang, or the
    Intel C++ compiler.

-   [GSL](https://www.gnu.org/software/gsl/) - The GNU Scientific
    Library.

-   The [Boost](https://www.boost.org/) filesystem and system C++
    libraries (boost\_filesystem and boost\_system).

-   [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) (v3) -
    C++ linear algebra package.

-   [LAPACK](http://www.netlib.org/lapack/) and
    [BLAS](http://www.netlib.org/blas/) - linear algebra subroutines.
    Can be substituted for internal libraries in the (proprietary) 
    [Intel Math Kernel Library](https://software.intel.com/en-us/mkl) 
    (MKL). Compiling against MKL allows expensive linear-algebra 
    operations (including generating angular data) to be automatically
    parallelised at run-time.

-   [Google Sparsehash](https://github.com/sparsehash/sparsehash)

-   Python v2.7.x

-   [SCons](http://scons.org/) - a python based software build tool.

-   The Unix `pkg-config` tool if you want the build system to
    automatically find and link against libraries. This should be
    installed by default on macOS and most Linux distributions.

-   (Optional) OpenMP and MPI - parallel programming libraries to run 
    AMBiT on high-performance computing clusters. These will generally be 
    cluster-specific, although OpenMP is supported by default on both GCC 
    and the Intel C++ compiler.

Most of these dependencies can be installed from your operating system's 
pacakge manager (e.g. apt-get in Ubuntu, dnf in Fedora or Homebrew in 
macOS). Some package managers have separate "development" packages, which 
you must also install to be able to compile AMBiT against them. These 
development packages tend to have names ending with some variation on "dev",
such as `libsparsehash-dev` on Ubuntu or `sparsehash-devel` on Fedora.

The AMBiT build process is based around a build tool called SCons, whose
build-control scripts (similar to Makefiles) are full Python scripts. To
simplify the process, we specify build options such as compiler and
flags, locations of required libraries, and which (if any) methods of
parallelism to employ in a plain-text file `config.ini`. This file has a
standard key-value structure and its filename is hard-coded into the
build system – you cannot specify a different file to use when
controlling the build. A minimal “skeleton” file is included with the
source code as `config\_template.ini`, which we recommend you copy (not
rename) this file to `config.ini` and fill in the desired compilation
options. If no file named `config.ini` exists, then the build system
will make a copy of the minimal `config\_template.ini` and use that. See
the AMBiT user guide for an example of this configuration file.

Most build options can be left unspecified and the build system will
attempt to automatically infer sensible defaults, but these inferences
can be explicitly overridden if required. It may be necessary to explicitly
specify the paths to some of the libraries required by AMBiT, since SCons
may not be able to automatically find the libraries on all systems. If 
you need to do this, the Unix `locate` command can be useful to find the
full path to the libraries (it may be necessary to refres `locate`'s 
database with the `updatedb` if you've just installed a package). 
For example, to find `sparsehash`, do:

    $ locate sparsehash

Each library should have an "include" directory, similar to 
`/usr/local/include/sparsehash` - this is the directory which should go
in your `config.ini`.

Once the `config.ini` file has
been filled as required, the software executable can be built by
navigating to the top-level AMBiT directory and issuing the `scons`
command. The minimal version of this (allowing the build system to
automatically find the required libraries) would look like:

    $ cd /path/to/ambit/
    $ cp config_template.ini config.ini
    $ scons

Additionally, you can pass the `-j` option to tell `scons` to use
multiple cores when compiling AMBiT. For example, on a quad-core system
we can do:
    
    $ cd /path/to/ambit
    $ scons -j 4

which will considerably speed up the build process by using compiling 4
files at a time in parallel.

Once you have a `config.ini` file working, it should continue to work
unless dependencies and libraries get moved (very rare, even considering
system upgrades). We strive to maintain backwards compatibility in the
build system and try our best not to make changes which break working
configurations.

Finally, `scons` takes an optional `target` to build, which can be
either `release` (the default), `debug` (to enable debugging options and
disable compiler optimisations), or `test` to make the test-suite. For
example:

    $ scons release
    $ scons debug
    $ scons test

### Make
We also supply some rudimentary Makefiles for systems where SCons is 
unavailable. Plain `make` can't automatically figure out linker paths like
SCons, so these files need to be hand-tailored to your specific system 
configuration. The system-specific parameters like compiler, compile flags 
and linker paths are set in `make.machine`; once this is set run `make`,
`make debug` or `make test` in the root directory.
