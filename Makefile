# To compile a debug version, define 'debug' in the command line, ie:
#     make debug=1
#
modules = Atom Configuration MBPT Basis HartreeFock Universal
exe = atom

packname = AtomPack
#other files that should be included in package
packfiles = Include.h config.txt mendel.tab CustomBasis.txt

# change this stuff only if you know what you are doing
ifdef debug
 BUILD = Debug
 CXXFLAGS += -g
 CFLAGS += -g
else
 BUILD = Release
 CXXFLAGS += -O3
 CFLAGS += -O3
endif

SHELL = /bin/sh
CC = gcc
CXX = g++
SRCDIR = $(shell pwd)

INCLUDE = $(SRCDIR) /usr/local/include
CXXFLAGS += $(addprefix -I, $(INCLUDE))
CFLAGS += $(addprefix -I, $(INCLUDE))

RANLIB = ranlib
AR = ar

LIBSUFFIX = lib.a
libnames = $(foreach module, $(modules), $(module)/$(BUILD)/$(module)$(LIBSUFFIX))

PLAT	     = _FREEBSD
PREF         = /usr/local/lib/lib
BLASLIB      = $(PREF)blas$(PLAT).a
LAPACKLIB    = $(PREF)clapack$(PLAT).a
LIBF2C       = $(PREF)f2c.a

$(exe): $(libnames)
	$(CXX) -o $(exe) $(libnames) $(LAPACKLIB) $(BLASLIB) $(LIBF2C) -lm -lc

.EXPORT_ALL_VARIABLES:

$(libnames): %$(LIBSUFFIX):
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
             Makefile make.instructions make.dependencies \
             $(foreach module, $(modules), $(module)/*.cpp $(module)/*.c $(module)/*.h $(module)/make.$(module))
	gzip $(packname).tar

unpack:
	gunzip $(packname).tar.gz
	-tar xvf $(packname).tar
