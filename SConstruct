import yaml
from collections import OrderedDict

# Function to force yaml to preserve the order we specify in the config file, since that's actually important at the linking stage
def ordered_load(stream):

    class OrderedLoader(yaml.SafeLoader):
        # Temporary loader class to pass to yaml.load(), derived from the base yaml loader
        pass

    def make_ordered_mapping(loader, node):
        # Converts a yaml node to an ordered dict
        loader.flatten_mapping(node)
        return(OrderedDict(loader.construct_pairs(node)))

    # Now feed the ordered mapping converter to the loader object and pass it along to yaml
    OrderedLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, make_ordered_mapping)
    return(yaml.load(stream, OrderedLoader))
        

# First, set up the compiler environment (cc and cflags and such) for the
# moment we'll just do a copy of the defult config in make.machine
Build = "Release"
defaultEnv = Environment(CXX = 'g++', CC = 'gcc', LINK = 'g++', \
    FORTRAN = 'gfortran', \
    CXXFLAGS = '-std=c++11')

# Now list the external libraries we need to link to, as well as their path
eigen_dir = ['/usr/local/include/Eigen3']
gtest_libs_dir = ['./gtest']
sparsehash_dir = ['/usr/local/include/sparsehash']
ambit_dir = ['/home/emily/Work/Atoms/Code/ambit']

libs = ['gsl', 'boost_filesystem', 'boost_system', 'lapack', 'blas']
lib_path = ['/usr/local/include'] + ambit_dir
header_path = ['/usr/local/include'] +  ambit_dir + eigen_dir + gtest_libs_dir + sparsehash_dir

# And append the paths to the current compilation environment
defaultEnv.Append(CPPPATH=header_path)
defaultEnv.Append(LIBPATH=lib_path)
defaultEnv.Append(LIBS=libs)

# Now open our YAML configuration file and set up our individual modules/objects. Note that this method is explicitly unsafe if the config file is dodgy
with open("config.yml", 'r') as fp:
    modules = ordered_load(fp)

# And build all of our common libraries
srcs = []
common_objs = []
common_libs = []
for module, files in modules.items():

    # First, compile all the object files for our module
    srcs = [module + '/' + item for item in files]
    module_lib_dir = module + '/' + Build
    obj = defaultEnv.Object(source = srcs)
    common_objs += obj
    
    # Now link them all into one library
    common_libs.append(Library(target = module_lib_dir + '/' + module, source = obj))
    #libs.append(module)
    #lib_path.append(module_lib_dir)

#defaultEnv.Append(LIBPATH=lib_path)
#defaultEnv.Append(LIBS=libs)

#ambit = defaultEnv.Object(["Atom/ambit.cpp", "Atom/ambit_recursive.cpp"])

# Finally, put it all together in one executable
defaultEnv.Program(target = 'ambit', source = common_libs)
