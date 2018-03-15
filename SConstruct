import subprocess # Required to call ./getGitInfo.sh
import os # Required to get the current working dir
import shutil # Required for copying files
import re
import ConfigParser

def get_build_target(src, build):
    # Constructs the path to the build object for some C++ source file, e.g.: 
    # Atom/ambit.cpp -> Atom/$BUILD/ambit.o

    # Strip the directory/module name and ".cpp" extension from the file
    basename = src[:src.find("cpp")].split('/')[1]
    obj_filename = basename + 'o'
    target = "./build/{}/{}".format(build, obj_filename)
    return(target)
    
def find_omp_flags(env):
    """ Searches for OpenMP compiler flags and, if not found, attempts to automatically add them.

    The OpenMP flags are compiler specific, so only certain compilers are supported here:

    gcc               ->      -fopenmp
    clang (Linux)     ->      -fopenmp=libiomp5
    Intel             ->      -qopenmp or -openmp
    Portland          ->      -mp
    clang (macOS)     ->      Not supported by default under XCode
    
    If the OpenMP flags were not directly supplied then the compiler choice is used to make a guess and
    append the correct set of flags. This is not capable of dealing with MPI compiler wrappers like
    mpiCC, since they provide no information about the underlying compiler.
    
    """

    # This pattern should match all OpenMP flags for gcc, clang, Intel and Portland compilers
    print("Checking for OpenMP flags...")
    omp_pattern = re.compile(r"(\-mp|\-openmp|\-qopenmp|\-fopenmp)\b")
    
    # Do the checks twice, in case the compiler and linker are different (this is usually unnecessary,
    # but sometimes NEEDS to happen)
    for case in ["CXX", "LINK"]:
        # Search both the compiler and flags, in case CXX or LINK is something like "g++ -fopenmp"
        match = omp_pattern.search(str(env[case + "FLAGS"]) + " " + str(env[case])) 
        if match:
            pass

        else:
            # If there were no openmp flags supplied, then attempt to deduce the correct one from the 
            # choice of compiler (note that sometimes the compiler is specified like "g++ -std=c++11", 
            # use a regular expression here)
            print("No OpenMP {}FLAGS supplied. Attempting to automatically generate...".format(case))
            compiler_pattern = re.compile(r"\b(g\+\+|clang\+\+|icpc|pgcc)")
            omp_flag_choices = {"g++": "-fopenmp", "clang++": "-fopenmp=libiomp5", "icpc": "-qopenmp", 
                             "pgcc": "-mp"}
            match = compiler_pattern.search(str(env[case]))
            if match:
                omp_flag = omp_flag_choices[match.group(1)]
                if case == "CXX":
                    env.AppendUnique(CXXFLAGS = omp_flag)
                elif case == "LINK":
                    env.AppendUnique(LINKFLAGS = omp_flag)
                print("Adding {} to {}FLAGS".format(omp_flag, case))
            else:
                print("Warning: could not automatically determine correct OpenMP compiler flags.")

    return(env)

def expand_flags(env, flags):
    """ Expands environment variables in supplied linker flags and passes them into the compilation env.
    
    """
    
    # Split the string of flags so we can parse each individually
    flags = flags.split()
    
    # Regex pattern: matches environment variables like ${FOO} or $FOO
    var_pattern = re.compile(r"(\$\{\w+\}|\$\w+\b)")
    expanded_flags = []
    for flag in flags:
        match = var_pattern.search(flag)
        
        if match:

            # Get the substring containing the variable and remove the $ and {} surrounding it
            span = match.span()
            stripped_var = match.group(0).strip("${}")
            try:
                flag = flag[:span[0]] + env["ENV"][stripped_var] + flag[span[1]:]
                expanded_flags.append(flag)
            except KeyError:
                flag = flag[:span[0]] + flag[span[1]:] # Remove the variable if it's not defined
                expanded_flags.append(flag)
        else:
            expanded_flags.append(flag)

    return(expanded_flags)


### Custom configuration/autoconf-like functions for AMBiT ###
def check_pkgconfig(context):
    # Checks that pkg-config exists and is configured
    context.Message("Checking pkg-config version... ")
    ret = context.TryAction("pkg-config --version")[0]
    context.Result(ret)
    return(ret)

def check_pkg(context, lib):
    # Checks that the specified library exists (and is findable by pkg-config)
    context.Message("Attempting to automatically find {}... ".format(lib))
    ret = context.TryAction("pkg-config --exists {}".format(lib))[0]
    context.Result(ret)
    return(ret)

def check_pkg_version(context, lib, version):
    # Checks the specified package is at the specified version or above
    context.Message("Checking {} is at least version {}...".format(lib, version))
    ret = context.TryAction("pkg-config --atleast-version={} {}".format(version, lib))[0]
    context.Result(ret)
    return(ret)

### Configuration functions ###
def read_config(env, conf, ambit_conf):
    """ Configures the supplied compilation environment according to the configuration options in conf.
    
        env - Pre-existing compilation environment, which should already have been initialised
        conf - ConfigParser object generated by parsing the build configuration file (config.ini)
        ambit_conf - ConfigParser object generated by parsing the dependencies configuration 
                     file (ambit_dependencies.ini) 
    """

    # The AMBiT directory will default to the current working directory if not specified in the config
    # file
    ambit_dir = conf.get("AMBiT options", "AMBiT path")
    if not ambit_dir:
        ambit_dir = dir_path = os.getcwd()
        print("AMBiT path not specified. Defaulting to {}".format(ambit_dir))

    # NOTE: These have their own, specially defined paths because they're usually in lots of places 
    # and don't play nicely with pkg-config and/or SCons builtin configuration checks
    try:
        gtest_libs_dir = conf.get("Dependency paths","gtest include path")
    except ConfigParser.NoOptionError:
        gtest_libs_dir = '' # gtest functionality can be safely left out in release builds

    try:
        eigen_dir = conf.get("Dependency paths", "Eigen path")
    except ConfigParser.NoOptionError:
        eigen_dir = ''

    try:
        sparsehash_dir = conf.get("Dependency paths", "Sparsehash path")
    except ConfigParser.NoOptionError:
        sparsehash_dir = ''

    libs = [s.strip() for s in conf.get("Dependencies", "Libs").split(',')]
    
    # Make sure to expand any environment variables in the paths
    lib_path = expand_flags(env, conf.get("Dependency paths", "Lib path"))
    custom_include = expand_flags(env, conf.get("Dependency paths", "Include path"))
    
    # Angular data should go in the AMBiT directory if not specified in the config file
    angular_data_dir = conf.get("AMBiT options", "Angular data")
    if not angular_data_dir:
        angular_data_dir = ambit_dir + "/AngularData"
        print("Angular data directory not specified. Defaulting to {}".format(ambit_dir))
    env.Append(CXXFLAGS = '-DANGULAR_DATA_DIRECTORY={}'.format(angular_data_dir))

    header_path =  [ambit_dir, eigen_dir, sparsehash_dir, custom_include]

    env.Append(CPPPATH=header_path)
    env.Append(LIBPATH=lib_path)
    env.Append(LIBPATH=gtest_libs_dir)
    env.Append(LIBS=libs)

    # Also grab the user specified compiler and compiler flags. This section is optional, so the keys may not 
    # exist
    try:
        custom_cxx_flags = conf.get("Compiler options", "CXXFLAGS")
        env.Append(CXXFLAGS = custom_cxx_flags)
    except ConfigParser.NoOptionError:
        pass
    # Link flags
    try:
        custom_linkflags = conf.get("Compiler options", "LINKFLAGS")
        env.Append(LINKFLAGS=custom_linkflags)
    except ConfigParser.NoOptionError:
        pass


    # Only replace the default compilers if the user has specified an alternative
    try:
        custom_cxx = conf.get("Compiler options", "CXX")
        if custom_cxx:
            env.Replace(CXX = custom_cxx)
    except ConfigParser.NoOptionError:
        pass

    try: 
        link = conf.get("Compiler options", "LINK")
        if link:
            env.Replace(LINK = link)
        else:
            env.Replace(LINK = env["CXX"]) # Try to fallback to the requested C++ compiler, if possible
    except ConfigParser.NoOptionError:
        env.Replace(LINK = env["CXX"]) 

    # Do the same for the Fortran compiler
    try:
        custom_F77 = conf.get("Compiler options", "F77")
        if custom_F77:
            env.Replace(FORTRAN = custom_F77)
    except ConfigParser.NoOptionError:
        pass

    # Final step before configuring: check if either OpenMP or MPI have been requested
    # NOTE: we only setup the -DAMBIT_USE_* flags here. The compiler and flags required by OpenMP and MPI
    # are system dependent, so must be specified in CXX and CXXFLAGS
    try:
        use_openmp = conf.getboolean("HPC options", "Use OpenMP") 
    except ValueError:
        use_openmp = False

    try:
        use_mpi = conf.getboolean("HPC options", "Use MPI") 
    except ValueError:
        use_mpi = False

    try:
        use_mkl = conf.getboolean("HPC options", "Use MKL")
    except ValueError:
        use_mkl = False

    if use_openmp:
        env.Append(CXXFLAGS = "-DAMBIT_USE_OPENMP")
        find_omp_flags(env) # Look for OpenMP flags/try to automagically infer them
    if use_mpi:
        env.Append(CXXFLAGS = "-DAMBIT_USE_MPI")
    if use_mkl:
        env.Append(CXXFLAGS = "-DEIGEN_USE_MKL_ALL")
        # Also read and parse the MKL link flags. These can just be of the form: -l<whatever> -I<path>
        try:
            mkl_flags = conf.get("HPC options", "MKL flags")
            if mkl_flags:
                mkl_flags = expand_flags(env, mkl_flags)
                env.MergeFlags(mkl_flags)
                
            else:
                print("Error: MKL flags must be explicitly defined in config.ini to use MKL")
                exit(-1)

        except ConfigParser.NoOptionError:
            print("Error: MKL flags must be explicitly defined in config.ini to use MKL")
            exit(-1)

    return(env)

def configure_environment(env, conf, ambit_conf):
    """ Runs autoconf-like checks to ensure that the compilation environment is properly configured.

        env - Pre-existing compilation environment, which should already have been initialised
        conf - ConfigParser object generated by parsing the build configuration file (config.ini)
        ambit_conf - ConfigParser object generated by parsing the dependencies configuration 
                     file (ambit_modules.ini) 
    """


    env_conf = Configure(env, custom_tests = {"check_pkgconfig": check_pkgconfig, 
                                              "check_pkg": check_pkg,
                                              "check_pkg_version": check_pkg_version})

    pkgconfig_exists = env_conf.check_pkgconfig()

    # Check the C++ and Fortran compilers exist and are properly configured
    if not env_conf.CheckCXX():
        print("Error: C++ compiler improperly installed/configured. Aborting")
        exit(-1)
    try:
      if not env_conf.CheckProg(env_conf.env['FORTRAN']):
          print("Error: Fortran compiler improperly installed/configured. Aborting.")
          exit(-1)
    except AttributeError: # Skip this check if we're using an old version of SCons
      pass

    # Run through the required libs, two approaches to find libraries:
    #   1) Look for it in the path specified in the configuration file
    #   2) Fallback: Try to find it automatically with pkg-config
    # If neither of these work then warn the user and bail out. Don't check for Boost here, since it doesn't support pkg-config
    for lib in [l for l in env["LIBS"] if l.find("boost") == -1]:
        if not env_conf.CheckLib(lib): # Can't find in current path
            if(not pkgconfig_exists or not env_conf.check_pkg(lib)):
                print("Warning: could not find library {}. Check [Dependency paths] in config.ini".format(lib))
                exit(-1)
            else:
               env_conf.env.ParseConfig("pkg-config --libs --cflags {}".format(lib))

    # Check the requested HPC libraries exist
    if "-DAMBIT_USE_OPENMP" in env_conf.env["CXXFLAGS"]:
        if not env_conf.CheckCXXHeader("omp.h"):
            print("Warning: failed to locate OpenMP headers. Aborting.")
            exit(-1)

    if "-DAMBIT_USE_MPI" in env_conf.env["CXXFLAGS"]:
        if not env_conf.CheckCXXHeader("mpi.h"):
            print("Warning: failed to locate MPI headers. Aborting.")
            exit(-1)

    # Check for Boost, Eigen and Sparsehash (these don't follow the nice pkg-config conventions)
    if not env_conf.CheckCXXHeader("boost/version.hpp"):
        print("Warning: could not find Boost headers. Check [Dependency paths] in config.ini")
        exit(-1)

    try:
        if not conf.get("Dependency paths", "Sparsehash path"):
            print("Sparsehash directory not specified...")
            if env_conf.check_pkg("libsparsehash"):
                env_conf.env.ParseConfig("pkg-config --libs --cflags libsparsehash")
            else:
                print("Failed to automatically locate Sparsehash headers. Specify Sparsehash path in config.ini")
                exit(-1)
    except ConfigParser.NoOptionError:
        print("Sparsehash directory not specified...")
        if env_conf.check_pkg("libsparsehash"):
            env_conf.env.ParseConfig("pkg-config --libs --cflags libsparsehash")
        else:
            print("Failed to automatically locate Sparsehash headers. Specify Sparsehash path in config.ini")
            exit(-1)
 

    try:
        if not conf.get("Dependency paths", "Eigen path"):
            print("Eigen directory not specified...")
            if env_conf.check_pkg("eigen3"):
                env_conf.env.ParseConfig("pkg-config --libs --cflags eigen3")
            else:
                print("Failed to automatically locate Eigen headers. Specify Eigen path in config.ini")
                exit(-1)
    except ConfigParser.NoOptionError:
        print("Eigen directory not specified...")
        if env_conf.check_pkg("eigen3"):
            env_conf.env.ParseConfig("pkg-config --libs --cflags eigen3")
        else:
            print("Failed to automatically locate Eigen headers. Specify Eigen path in config.ini")
            exit(-1) 

    # Also check the version of Eigen is at least 3.2 and gsl is version 2.x (or above) 
    if pkgconfig_exists:
        env_conf.check_pkg_version("gsl", "2.0")
        env_conf.check_pkg_version("eigen3", "3.2.0")
            
    env = env_conf.Finish()
    print("Finished configuring build environment.\n")
    return(env)



### Begin SCons commands ###

# Before we construct anything, make sure to update the gitinfo header file
with open("gitInfo.h", 'w') as fp:
    subprocess.call("./getGitInfo.sh", stdout=fp)

# First, grab the type of build from the command line and set up the compiler environment (default is gcc)
# NOTE: We need to explicitly import the shell (bash) from the environment so compiler checks work
env = Environment(CXX = 'g++', CC = 'gcc', LINK = 'g++', \
    FORTRAN = 'gfortran', \
    F77FLAGS = '-O2',
    SHELL = "/bin/bash",
    ENV = os.environ)

build = "Release" # Default build target (i.e. don't do debug unless specifically asked for it)

if 'debug' in COMMAND_LINE_TARGETS:
    build = "Debug"
    env.Append(CXXFLAGS = '-std=c++11 -g -Wno-deprecated-register -Wno-unused-result -O0')
elif 'test' in COMMAND_LINE_TARGETS:
    build = "Test"
    env.Append(CXXFLAGS = '-std=c++11')
    env.Append(LIBS = ['gtest', 'pthread'])
else:
    env.Append(CXXFLAGS = '-std=c++11 -O3')

# Now open the ini files containing the user config options, and the AMBiT requirements
conf = ConfigParser.SafeConfigParser(allow_no_value = True)
conf.optionxform = str # Necessary to make keys case-sensitive
if not conf.read("config.ini"):
    print("Warning: configuration file config.ini not found. Copying configuration template ...")
    shutil.copy2("config_template.ini", "config.ini")
    conf.read("config.ini")

ambit_conf = ConfigParser.SafeConfigParser(allow_no_value = True)
ambit_conf.optionxform = str # Necessary to make keys case-sensitive
# Bail out if we can't read config.ini
if not ambit_conf.read("ambit_modules.ini"):
    print("Error: could not read dependencies file ambit_modules.ini. Aborting build.")
    exit(-1)

# Configure the build environment and run autoconf-like checks, but only if not cleaning targets
env = read_config(env, conf, ambit_conf)
if not env.GetOption('clean'):
    env = configure_environment(env, conf, ambit_conf)

# And build all of our common libraries/modules
modules = ambit_conf.items("Modules")
common_libs = []
test_objs = []
for module, files in modules:

    # First, split the comma separated list of files and strip out whitespace
    files = [s.strip() for s in files.split(',')]

    # Next, compile all the object files for our module

    srcs = [module + '/' + item for item in files]
    module_lib_target = './build/' + build + '/' + module
    
    # Compile the object files
    objs = [env.Object(target = get_build_target(src, build), source = src) 
            for src in srcs]

    # Compile the unit tests into objects if requested
    if 'test' in COMMAND_LINE_TARGETS:
        test_objs += [env.Object(target = get_build_target(str(src), build), source = src) 
            for src in Glob("{}/*.test.cpp".format(module))]
    
    # Now link them all into one library
    common_libs.append(Library(target = module_lib_target, source = objs))

# Finally, put it all together in one executable (note that debug and release both alias to the ambit 
# executable. Eventually I'll add a test target in as well)
if 'test' in COMMAND_LINE_TARGETS:
    env.Program(target = 'ambit_test', source = test_objs + common_libs)
env.Program(target = 'ambit', source = common_libs)
env.Alias('release', 'ambit')
env.Alias('debug', 'ambit')
env.Alias('test', 'ambit_test')
