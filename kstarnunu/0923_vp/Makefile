# Generated automatically from Makefile.in by configure.
#
# File: Makefile.in
# Description : Makefile Template for fsim_analysis
#
# 1. User Specifications
OBJS = ./UserUtility.o ./shape.o ./k_sfw.o ./ExKFitter.o ana_fullreconk.o #change here
#OBJS = ana_d0pi.o
MODULE = ./ana_fullreconk.so  #change here
LIBS = -L$(BELLE_RUN_DIR)/lib/so -ltuple -lmdst -lkfitter -lip -ltagv -lkid -lparticle -lhamlet -leid -lbenergy -lfullrecon -ltables 

# 2. System Specifications
#    --- Do not change without knowledge
# Compiler Setup with machine dependence

FC = g77
CC = gcc
CXX = g++

DEFS =  -DWORDS_BIGENDIAN=1 -DHAVE_LIBCURSES=1 -DHAVE_LIBREADLINE=1 -DHAVE_LIBNSL=1 -DHAVE_LIBSOCKET=1 -DHAVE_LIBDL=1 -DFORTRAN_PPU=1  -DCERNLIB_TYPE
CPPFLAGS = 
DEPCPPFLAGS = -MM

FFLAGS = -O -g -PIC
CFLAGS = -DBELLE_TARGET_H=\"belle-x86_64-unknown-linux-gnu-g++.h\" -g -O2
CXXFLAGS = -DHEP_SHORT_NAMES -DBELLE_SHORT_NAMES -DDSTXX_NOINLINE -DBELLE_TARGET_H=\"belle-x86_64-unknown-linux-gnu-g++.h\" -g -fPIC
SOFLAGS = -shared -Wl,-export-dynamic
LDFLAGS = 
SYSLIB = -lcrypt -lm -lgcc -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -lg2c
CLHEPINC = /belle/local/include
CLHEPLIB = -L$(BELLE_LIB_DIR)  -lbelleCLHEP
CLHEPINC = 
CLHEPLIB = 

# Include directories
INCLUDES_C = $(MY_TOP_DIR)/include $(BELLE_TOP_DIR)/include $(CLHEPINC) ../
INCLUDES_FORTRAN = $(BELLE_TOP_DIR)/inc

# Dependence description
include $(BELLE_RUN_DIR)/src/config/Makefile.panther

COMPILE_FCPP := $(FC) -c $(PANTHER_FMACROS) $(INCLUDES_FORTRAN:%=-I%) $(CPPFLAGS) $(FFLAGS)
COMPILE_FC := $(FC) -c  $(INCLUDES_FORTRAN:%=-I%) $(FFLAGS)
COMPILE_CC := $(CC) -c  $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(CFLAGS)
COMPILE_CXX := $(CXX) -c  $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(CXXFLAGS)

LINK_FCPP := $(FC)
LINK_FC := $(FC)
LINK_CC := $(CC)
LINK_CXX := $(CXX)

DEPEND_FCPP := /lib/cpp -M $(DEFS) $(PANTHER_FMACROS) $(INCLUDES_FORTRAN:%=-I%) $(CPPFLAGS) $(FFLAGS)
DEPEND_CC := $(CC) $(DEPCPPFLAGS) $(DEFS) $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(CFLAGS)
DEPEND_CXX := $(CXX) $(DEPCPPFLAGS) $(DEFS) $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(CXXFLAGS)

LINUX_G77_BUG = -e 's/\.F\.o:/.o:/' 

%.o:%.c
	$(COMPILE_CC) $<

%.d:%.c
	$(SHELL) -ec '$(DEPEND_CC) $< | sed -e "s/$*.o[ :]*/$@ &/g" -e 's/\.[12][0-9][0-9][0-9][0-9][0-9][0-9][0-9][a-z]\.h/.tdf/g' > $@'

%.o:%.cc
	$(COMPILE_CXX) $<

%.d:%.cc
	$(SHELL) -ec '$(DEPEND_CXX) $< | sed -e "s/$*.o[ :]*/$@ &/g" -e 's/\.[12][0-9][0-9][0-9][0-9][0-9][0-9][0-9][a-z]\.h/.tdf/g'> $@'


%.o:%.F
	$(COMPILE_FCPP) $<

%.d:%.F
	$(SHELL) -ec '$(DEPEND_FCPP) $< | sed $(LINUX_G77_BUG) -e "s/$*.o[ :]*/$@ &/g" -e 's/\.[12][0-9][0-9][0-9][0-9][0-9][0-9][0-9][a-z]\.inc/.tdf/g'> $@'


# CERNLIB
ifeq "$(CERN)/$(CERN_LEVEL)" "/"
  CERNLIB_LIB_DIR = /belle/cern/96a/lib
else
  CERNLIB_LIB_DIR = $(CERN)/$(CERN_LEVEL)/lib
endif


# Dependencies
all::	$(OBJS)
	$(LINK_CXX) -o $(MODULE) $(SOFLAGS) $(OBJS) $(LIBS) \
	$(CLHEPLIB) $(CERNLIB) $(SYSLIB)
	cp $(MODULE) /gpfs/group/belle/users_old/boyuan/special_topic/skim/scriptkstar
#	cp $(MODULE) /gpfs/fs02/belle/users/boyuan/special_topic/skim/scriptkstar/
#	cp $(MODULE) /gpfs/fs02/belle/users/boyuan/special_topic/skim/scriptkstar/
check:
	./run.csh $(BELLE_LEVEL) > ./out 2>&1

clean::
	rm -f $(OBJS) $(MODULE) *~ *.o fort.* fpda_pid.* test.hbk out

