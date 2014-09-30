# Run with gmake.
# To compile a debug version, define 'debug' in the command line, ie:
#     gmake debug=1
# To compile a google test, use
#     gmake test=1
# 'gmake test' and 'gmake debug' will also work, but not with 'clean'.

coremodules = Basis Configuration ExternalField HartreeFock MBPT Universal
modules = Atom $(coremodules)

ifdef test
  exes = ambit_test

  # The *.test.o objects need to be linked in separately because google tests are
  # left out of the linking of the static library (because there are no unresolved
  # symbols in the main program that resolve to that object - see googletest wiki).

  gtestcxx = $(foreach module, $(modules), $(wildcard $(module)/*.test.cpp))
  gtesttmp = $(join $(dir $(gtestcxx)),$(addprefix $(BUILD)/,$(notdir $(gtestcxx))))
  gtestobj = $(gtesttmp:%.cpp=%.o)
else
  exes = ambit
endif

packname = AtomPack
# other files that should be included in package
packfiles = Include.h CustomBasis.txt Atom/GetPot getGitInfo.sh

include make.machine

corelibnames = $(foreach module, $(coremodules), $(module)/$(BUILD)/$(module)$(LIBSUFFIX))
libnames = $(foreach module, $(modules), $(module)/$(BUILD)/$(module)$(LIBSUFFIX))

$(exes): $(libnames)
	$(LINK) $(LINKFLAGS) -o $@ $^ $(gtestobj) $(addprefix -L, $(LIBDIR)) \
                $(LIBOBJ) $(addprefix -l, $(LIBS))
	-@mkdir -p $(ANGULAR_DATA_DIRECTORY)

.EXPORT_ALL_VARIABLES:

$(libnames): %$(LIBSUFFIX): gitInfo.h FORCE
	-@mkdir -p $(notdir $*)/$(BUILD)
	$(MAKE) -C $(notdir $*) -f make.$(notdir $*) module=$(notdir $*)

# Always rebuild gitInfo.h for the latest compile info
gitInfo.h: FORCE
	./getGitInfo.sh > gitInfo.h

#if a command fails, delete the target
.DELETE_ON_ERROR:

.PHONY: FORCE clean veryclean pack unpack unpackbackend

# FORCE is used to always rebuild subdirectories.
FORCE:

test:
	$(MAKE) test=1

debug:
	$(MAKE) debug=1

clean:
	-rm $(exes); \
	$(foreach module, $(modules), $(MAKE) -C $(module) -f make.$(module) clean;)

veryclean:
	-rm ambit ambit_test
	-rm -r $(foreach module, $(modules), $(module)/$(BUILD))

pack:
	-tar cf $(packname).tar $(packfiles) \
	     Makefile make.instructions make.dependencies make.machine\
	     $(foreach module, $(modules), $(module)/*.cpp $(module)/*.c $(module)/*.h $(module)/*.f $(module)/make.$(module))
	gzip $(packname).tar

unpack:
	gunzip $(packname).tar.gz
	-tar xvf $(packname).tar
