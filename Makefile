# Run with gmake.
# To compile a debug version, define 'debug' in the command line, ie:
#     make debug=1
#
coremodules = Configuration MBPT Basis HartreeFock Universal
modules = Atom $(coremodules)
exes = atom rap.exe

packname = AtomPack
# other files that should be included in package
packfiles = Include.h config.txt mendel.tab CustomBasis.txt Atom/make.RMatrixPrimer

include make.machine

corelibnames = $(foreach module, $(coremodules), $(module)/$(BUILD)/$(module)$(LIBSUFFIX))
libnames = $(foreach module, $(modules), $(module)/$(BUILD)/$(module)$(LIBSUFFIX))

atom: $(libnames)
	$(LINK) $(LINKFLAGS) -o $@ $^ $(addprefix -L, $(LIBDIR)) \
	     $(LIBOBJ) $(addprefix -l, $(LIBS))

rap.exe: Atom/$(BUILD)/RMatrixPrimer$(LIBSUFFIX) $(corelibnames)
	$(LINK) $(LINKFLAGS) -o $@ $^ $(addprefix -L, $(LIBDIR)) \
	     $(addprefix -l, $(LIBS))

.EXPORT_ALL_VARIABLES:

$(libnames): %$(LIBSUFFIX): FORCE
	-@mkdir $(notdir $*)/$(BUILD)
	$(MAKE) -C $(notdir $*) -f make.$(notdir $*) module=$(notdir $*)

# RMatrixPrimer is special because it resides in directory "Atom".
Atom/$(BUILD)/RMatrixPrimer$(LIBSUFFIX): %$(LIBSUFFIX): FORCE
	-@mkdir Atom/$(BUILD)
	$(MAKE) -C Atom -f make.$(notdir $*) module=$(notdir $*)

#if a command fails, delete the target
.DELETE_ON_ERROR:

.PHONY: FORCE clean veryclean pack unpack unpackbackend

# FORCE is used to always rebuild subdirectories.
# Could also use .PHONY, but this is easier for, e.g., RMatrixPrimer.
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

unpackbackend:
	gunzip $(packname).tar.gz
	cp Atom/Atom.cpp Atom/Atom.cpp.save
	cp Atom/Atom_Open.cpp Atom/Atom_Open.cpp.save
	cp Atom/RMatrixPrimer.cpp Atom/RMatrixPrimer.cpp.save
	-tar xvf $(packname).tar
	mv Atom/Atom.cpp.save Atom/Atom.cpp
	mv Atom/Atom_Open.cpp.save Atom/Atom_Open.cpp
	mv Atom/RMatrixPrimer.cpp.save Atom/RMatrixPrimer.cpp
