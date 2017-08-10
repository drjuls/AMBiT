import subprocess
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

def get_build_target(src, module, build):
    # Constructs the path to the build object for some C++ source file, e.g.: 
    # Atom/ambit.cpp -> Atom/$BUILD/ambit.o

    # Strip the directory/module name and ".cpp" extension from the file
    basename = src[:src.find("cpp")].split('/')[1]
    obj_filename = basename + 'o'
    target = "{}/{}/{}".format(module, build, obj_filename)
    return(target)

### Begin SCons commands ###

# Before we construct anything, make sure to update the gitinfo header file
with open("gitInfo.h", 'w') as fp:
    subprocess.call("./getGitInfo.sh", stdout=fp)

# First, grab the type of build from the command line and set up the compiler environment (cc and cflags 
# and such). Default build mode is release
env = Environment(CXX = 'g++', CC = 'gcc', LINK = 'g++', \
    FORTRAN = 'gfortran', \
    F77FLAGS = '-O2')
if 'debug' in COMMAND_LINE_TARGETS:
    build = "Debug"
    env.Append(CXXFLAGS = '-std=c++11 -g -Wno-deprecated-register -Wno-unused-result -O0')
else:
    build = "Release"
    env.Append(CXXFLAGS = '-std=c++11 -O3')

# Now open our YAML configuration file and set up our individual modules/objects. Note that this method 
# can only read basic Python data structures from the config file so we (hopefully) don't get owned if 
# someone supplies a dodgy yaml file
with open("config.yml", 'r') as fp:
    conf = ordered_load(fp)

modules = conf['Modules']
build_info = conf['Build Info']

# Now grab the external libraries we need to link to, as well as their path
ambit_dir = build_info["AMBiT"]
gtest_libs_dir = build_info["gtest"]
eigen_dir = build_info["Eigen"]
sparsehash_dir = build_info["Sparsehash"]

libs = build_info["Libs"]
lib_path = build_info["Lib path"]
header_path =  [ambit_dir, eigen_dir , gtest_libs_dir , sparsehash_dir]
angular_data_dir= build_info["Angular data"]
env.Append(CXXFLAGS = '-DANGULAR_DATA_DIRECTORY={}'.format(angular_data_dir))

# And append the paths to the current compilation environment
env.Append(CPPPATH=header_path)
env.Append(LIBPATH=lib_path)
env.Append(LIBS=libs)

# Now grab the user specified compiler and compiler flags. This section is optional, so the keys may not 
# exist
if "CXXFLAGS" in build_info.keys():
    custom_cxx_flags = build_info["CXXFLAGS"]
    env.Append(CXXFLAGS = custom_cxx_flags)

    # Only replace the default compiler if the user has specified an alternative
if "CXX" in build_info.keys():
    custom_cxx = build_info["CXX"]
    if custom_cxx:
        env.Replace(CXX = custom_cxx)

# And build all of our common libraries
common_libs = []
for module, files in modules.items():

    # First, compile all the object files for our module
    srcs = [module + '/' + item for item in files]
    module_lib_target = module + '/' + build + '/' + module
    
    # Compile the object files
    objs = [env.Object(target = get_build_target(src, module, build), source = src) 
            for src in srcs]
    
    # Now link them all into one library
    common_libs.append(Library(target = module_lib_target, source = objs))

# Finally, put it all together in one executable (note that debug and release both alias to the ambit 
# executable. Eventually I'll add a test target in as well)
env.Program(target = 'ambit', source = common_libs)
env.Alias('release', 'ambit')
env.Alias('debug', 'ambit')
