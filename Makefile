# Run with gmake.
# To compile a debug version, define 'debug' in the command line, ie:
#     make debug=1
#
coremodules = Configuration MBPT Basis HartreeFock Universal
modules = Atom $(coremodules)
exes = ambit

packname = AtomPack
# other files that should be included in package
packfiles = Include.h CustomBasis.txt Atom/GetPot getGitInfo.sh

include make.machine

corelibnames = $(foreach module, $(coremodules), $(module)/$(BUILD)/$(module)$(LIBSUFFIX))
libnames = $(foreach module, $(modules), $(module)/$(BUILD)/$(module)$(LIBSUFFIX))

ambit: $(libnames)
	$(LINK) $(LINKFLAGS) -o $@ $^ $(addprefix -L, $(LIBDIR)) \
	     $(LIBOBJ) $(addprefix -l, $(LIBS))

.EXPORT_ALL_VARIABLES:

$(libnames): %$(LIBSUFFIX): gitInfo.h FORCE
	-@mkdir $(notdir $*)/$(BUILD)
	$(MAKE) -C $(notdir $*) -f make.$(notdir $*) module=$(notdir $*)

# Always rebuild gitInfo.h for the latest compile info
gitInfo.h: FORCE
	sh getGitInfo.sh > gitInfo.h

#if a command fails, delete the target
.DELETE_ON_ERROR:

.PHONY: FORCE clean veryclean pack unpack unpackbackend

# FORCE is used to always rebuild subdirectories.
FORCE:

clean:
	-rm $(exes); \
	$(foreach module, $(modules), $(MAKE) -C $(module) -f make.$(module) clean;)

veryclean:
	-rm $(exes)
	-rm -r $(foreach module, $(modules), $(module)/$(BUILD))

pack:
	-tar cf $(packname).tar $(packfiles) \
	     Makefile make.instructions make.dependencies make.machine\
	     $(foreach module, $(modules), $(module)/*.cpp $(module)/*.c $(module)/*.h $(module)/*.f $(module)/make.$(module))
	gzip $(packname).tar

unpack:
	gunzip $(packname).tar.gz
	-tar xvf $(packname).tar
