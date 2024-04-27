# AMBiT
AMBiT is a high-precision atomic structure software package based on the Configuration Interaction + Many-Body Perturbation Theory (CI+MBPT) method and is developed and maintained by the Julian Berengut Group at the University of New South Wales in Sydney, Australia. While we are contactable by email, we would prefer for bug reports and feature requests to be made on this GitHub page.

### Features
Briefly, AMBiT is:

* Highly-accurate - AMBiT typically gives better than 1% accuracy for simple atoms and highly-charged ions, while more complicated open-shell systems (e.g. Lanthanides and Actinides) AMBiT agrees with experiment to approximately 10% (note that such systems would be either extremely difficult or outright impossible to simulate using other codes)
* Flexible - The particle-hole extension to CI+MBPT allows us to treat a wide range of atoms with many valence electrons which would be otherwise intractable
* Ab-initio - no fitting or semi-empirical fudge-factors are necessary to get good results
* Fast and scalable - AMBiT is highly parallelised and designed to take full advantage of modern high-performance computing architecture, while still being possible to run on a personal workstation (or even a notebook computer!)
* Modern - Written in C++11, AMBiT makes heavy use of modern computer science paradigms which allows for rapid development and integration of new functionality
* User friendly - we consider usability to be a first-class constraint on our development. 

### Installation
AMBiT is currently in active development: the ''dev'' branch is regularly updated with new features and performance improvements, while the ''master'' branch is more stable and only updated once we are sure it is correct and bug-free. If you are installing AMBiT globally (i.e. not in your personal directory) on a cluster then we recommend you use the master branch - otherwise use either dev or master depending on your tolerance for software changes. 

AMBiT uses the [CMake](https://cmake.org/) build-system. It also requires the following external libraries:
- [GSL](https://www.gnu.org/software/gsl/) - The GNU Scientific Library.
- [Boost](https://www.boost.org/)
- [Eigen v3](http://eigen.tuxfamily.org/index.php?title=Main_Page) - C++ linear algebra package.
- [LAPACK](http://www.netlib.org/lapack/) and [BLAS](http://www.netlib.org/blas/) - linear algebra subroutines. Can be substituted for the [Intel Math Kernel Library (MKL)](https://software.intel.com/en-us/mkl).
- [Google Abseil](https://github.com/abseil/abseil-cpp)

See [INSTALL.md](https://github.com/drjuls/AMBiT/blob/master/INSTALL.md) for more detailed
instructions on how to compile AMBiT.

### Contributing
We welcome contributions, bug reports and feature requests. See [CONTRIBUTING.md](https://github.com/drjuls/AMBiT/blob/master/CONTRIBUTING.md) for contribution and bug-report guidlines.

### Licensing
This project is licensed under the GNU GPL v3.0. See [LICENSE](https://github.com/drjuls/AMBiT/blob/master/LICENSE) for more details
