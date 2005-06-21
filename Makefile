# To compile a debug version, define 'debug' in the command line, ie:
#     make debug=1
#
modules = Atom Configuration MBPT Basis HartreeFock Universal
exe = atom

packname = AtomPack
#other files that should be included in package
packfiles = Include.h config.txt mendel.tab CustomBasis.txt

include make.machine

libnames = $(foreach module, $(modules), $(module)/$(BUILD)/$(module)$(LIBSUFFIX))

$(exe): $(libnames)
	$(LINK) -o $(exe) $(libnames) $(addprefix -L, $(LIBDIR)) \
             $(addprefix -l, $(EXTRALIBS) $(LAPACKLIB) $(BLASLIB) m c)

.EXPORT_ALL_VARIABLES:

$(libnames): %$(LIBSUFFIX):
	-@mkdir $(notdir $*)/$(BUILD)
	$(MAKE) -C $(notdir $*) -f make.$(notdir $*) module=$(notdir $*)

#if a command fails, delete the target
.DELETE_ON_ERROR:

.PHONY: clean veryclean $(libnames)
clean:
	-rm $(exe); \
	$(foreach module, $(modules), $(MAKE) -C $(module) -f make.$(module) clean;)

veryclean:
	-rm $(exe)
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
	-tar xvf $(packname).tar
	mv Atom/Atom.cpp.save Atom/Atom.cpp
	mv Atom/Atom_Open.cpp.save Atom/Atom_Open.cpp
