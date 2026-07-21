###############################################################################
#                                                                             #
# Makefile for glm                                                            #
#                                                                             #
#  Part of GLM (General Lake Model)                                           #
#                                                                             #
#  Developed by :                                                             #
#      AquaticEcoDynamics (AED) Group                                         #
#      School of Agriculture and Environment                                  #
#      The University of Western Australia                                    #
#                                                                             #
#      http://aquatic.science.uwa.edu.au/                                     #
#                                                                             #
#  Copyright 2013-2026 : The University of Western Australia                  #
#                                                                             #
#   GLM is free software: you can redistribute it and/or modify               #
#   it under the terms of the GNU General Public License as published by      #
#   the Free Software Foundation, either version 3 of the License, or         #
#   (at your option) any later version.                                       #
#                                                                             #
#   GLM is distributed in the hope that it will be useful,                    #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#   GNU General Public License for more details.                              #
#                                                                             #
#   You should have received a copy of the GNU General Public License         #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

ifeq ($(shell uname),Linux)
  OSTYPE=$(shell uname -s)
else ifeq ($(shell uname),Darwin)
  OSTYPE=$(shell uname -s)
else ifeq ($(shell uname),FreeBSD)
  OSTYPE=$(shell uname -s)
else
  OSTYPE=$(shell uname -o)
endif

BUILDDATE=$(shell date -u +%Y%m%d-%H%MUTC)

ifeq ($(MDEBUG),true)
  DEBUG=true
endif
ifeq ($(WITH_PLOTS),)
  WITH_PLOTS=true
  ifeq ($(WITH_XPLOTS),)
    WITH_XPLOTS=true
  endif
  ifeq ($(PLOTDIR),)
    PLOTDIR=../libplot
  endif
endif

ifeq ($(UTILDIR),)
  UTILDIR=../libutil
endif

ifeq ($(WITH_CHECKS),)
  ifeq ($(DEBUG),true)
    WITH_CHECKS=true
  else
    WITH_CHECKS=false
  endif
endif

srcdir=src
incdir=src
objdir=obj
moddir=mod
libdir=lib

TARGETS=glm libglm.${so_ext}
DEFINES=
FINCLUDES+=-I$(UTILDIR)/include
CINCLUDES+=-I$(UTILDIR)/include -I../ancillary/include
LIBS=-L$(UTILDIR)/lib -lutil
GLM_DEPS=$(UTILDIR)/lib/libutil.a
ifeq ($(WITH_PLOTS),true)
  DEFINES+=-DPLOTS
  ifeq ($(WITH_XPLOTS),true)
    DEFINES+=-DXPLOTS
  endif
  GLM_DEPS+=$(PLOTDIR)/lib/libplot.a
  CINCLUDES+=-I$(PLOTDIR)/include
endif

ifeq ($(OSTYPE),Darwin)
  ifeq ($(HOMEBREW),true)
     FINCLUDES+=-I/usr/local/include -I/opt/homebrew/include
     CINCLUDES+=-I/usr/local/include -I/opt/homebrew/include
     LIBS+=-L/usr/local/lib -L/opt/homebrew/lib
  else
     FINCLUDES+=-I/opt/local/include
     CINCLUDES+=-I/opt/local/include
     LIBS+=-L/opt/local/lib
  endif
  #EXTRALINKFLAGS=-Wl,-no_compact_unwind
  EXTRALINKFLAGS=-Wl,-no_compact_unwind,-headerpad_max_install_names
  SHARED=-dynamiclib -undefined dynamic_lookup
  so_ext=dylib
else ifeq ($(OSTYPE),Msys)
  ifeq ("$(wildcard '../ancillary/windows')", "")
    LIBS+=-L../ancillary/lib
  else
    CINCLUDES+=-I../ancillary/windows/include
    LIBS+=-L../ancillary/windows/lib
  endif
  SHARED=-shared
  so_ext=dll
  # do_model_non_avg (and siblings) declare several MaxInf/MaxOut/MaxVars-sized
  # arrays (WQNew[MaxInf*MaxVars] alone is 10000*60 doubles, ~4.8MB) as plain
  # locals rather than heap-allocating them. That fits fine under Linux/macOS's
  # ~8MB default thread stack, but silently overflows Windows' 2MB default
  # link-time stack reserve - crashing with no error message right on entry to
  # the function, before any of its own code runs. Reserve more stack at link
  # time until these arrays are heap-allocated properly upstream.
  EXTRALINKFLAGS=-Wl,--stack,67108864
else ifeq ($(OSTYPE),FreeBSD)
  FINCLUDES+=-I/usr/local/flang/include -I/usr/local/include
  CINCLUDES+=-I/usr/local/include
  LIBS+=-L/usr/local/flang/lib -L/usr/local/lib -lm
  ifeq ($(MDEBUG),true)
    DBG_LIBS=-fsanitize=address -static-libsan
  endif
  SHARED=-shared
  so_ext=so
else
  CINCLUDES+=-I/usr/local/include
  EXTRALINKFLAGS=-Wl,-z,relro,--export-dynamic
  ifeq ($(MDEBUG),true)
    DBG_LIBS=-fsanitize=address -static-libasan
  endif
  SHARED=-shared
  so_ext=so
endif

ifeq ($(FABM),true)
  ifeq ($(FABMDIR),)
    FABMDIR=../fabm-git
  endif
  DEFINES+=-DFABM

  FABMLIB=fabm
  ifeq ($(DEBUG),true)
    WITH_CHECKS=true
  endif

  FINCLUDES+=-I$(FABMDIR)/include -I$(FABMDIR)/src/drivers/glm -I$(FABMDIR)/build/modules
  FABMLIBS=-L$(FABMDIR)/build -l$(FABMLIB)

  ifeq ($(USE_DL),true)
    FABMTARGETS=libglm_wq_fabm.${so_ext}
  endif
endif

EXTFFLAGS=
ifeq ($(WITH_AED),true)
  DEFINES+=-DAED
endif
ifeq ($(WITH_API),true)
  DEFINES+=-DAPI
  WITH_AED=true
endif

ifeq ($(WITH_AED),true)
  AEDWATDIR=../libaed-water
  FINCLUDES+=-I$(AEDWATDIR)/include -I$(AEDWATDIR)/mod
  CINCLUDES+=-I$(AEDWATDIR)/include
  AEDLIBS=-L$(AEDWATDIR)/lib -laed-water
  ifeq ($(WITH_API),true)
    AEDAPIDIR=../libaed-api
    FINCLUDES+=-I$(AEDAPIDIR)/include -I$(AEDAPIDIR)/mod
    CINCLUDES+=-I$(AEDAPIDIR)/include
    AEDLIBS+=-L$(AEDAPIDIR)/lib -laed-api
  endif
  ifdef AEDBENDIR
    AEDLIBS+=-L$(AEDBENDIR)/lib -laed-benthic
  else
    EXTFFLAGS+=-DNO_BENTHIC
  endif
  ifdef AEDDMODIR
    AEDLIBS+=-L$(AEDDMODIR)/lib -laed-demo
  else
    EXTFFLAGS+=-DNO_DEMO
  endif

  # the plus version libs
  ifdef AEDRIPDIR
    AEDLIBS+=-L$(AEDRIPDIR)/lib -laed-riparian
  else
    EXTFFLAGS+=-DNO_RIPARIAN
  endif
  ifdef AEDLGTDIR
    AEDLIBS+=-L$(AEDLGTDIR)/lib -laed-lighting
  else
    EXTFFLAGS+=-DNO_LGT
  endif
  ifdef AEDDEVDIR
    EXTRA_FLAGS+=-DWITH_AED_PLUS
    AEDLIBS+=-L$(AEDDEVDIR)/lib -laed-dev
    ifdef PHREEQDIR
       AEDLIBS+=-L$(PHREEQDIR)/build -lPhreeqcRM
    endif
    ifdef MODFLOWDIR
       AEDLIBS+=-L$(MODFLOWDIR)/build/srcbmi -lmf6
    endif
  else
    EXTFFLAGS+=-DNO_DEV
  endif

  ifeq ($(USE_DL),true)
    AEDTARGETS=libglm_wq_aed.${so_ext}
  endif

  GLM_DEPS+=$(AEDWATDIR)/lib/libaed-water.a
endif

ifeq ($(AED2),true)
  DEFINES+=-DAED2

  ifeq ($(AED2DIR),)
    AED2DIR=../libaed2
  endif

  FINCLUDES+=-I$(AED2DIR)/include -I$(AED2DIR)/mod
  AED2LIBS=-L$(AED2DIR)/lib -laed2
  ifneq ("$(wildcard ${AED2PLS}/Makefile)","")
    AED2PLBS=-L${AED2PLS}/lib -laed2+
  endif

  ifeq ($(USE_DL),true)
    AED2TARGETS=libglm_wq_aed2.${so_ext}
  endif

  GLM_DEPS+=$(AED2DIR)/lib/libaed2.a
endif

FLIBS=
# Select specific compiler bits
ifeq ($(F90),ifort)
  LINK=$(CC)
  FINCLUDES+=-I/opt/intel/include
  DEBUG_FFLAGS=-g -traceback -DDEBUG=1
  OMPFLAG=-qopenmp
  OPT_FFLAGS=-O3
  FFLAGS=-warn all -module ${moddir} -static-intel -mp1 -stand f08 -warn nounused $(DEFINES) $(FINCLUDES)
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-check bounds -check noarg_temp_created
  endif
  FFLAGS+=-real-size 64
  FLIBS+=-L/opt/intel/lib
  FLIBS+=-lifcore -lsvml -lifport
  FLIBS+=-limf -lintlc -liomp5 -lifport
  #EXTFFLAGS=-warn-no-unused-dummy-argument
else ifeq ($(F90),ifx)
  DEBUG_FFLAGS=-g -traceback -DDEBUG=1 -O0
  OMPFLAG=-qopenmp
  OPT_FFLAGS=-O3
  ifeq ($(OSTYPE),Msys)
    LINK=$(FC) /nofor_main
    FFLAGS=-fpp -warn=all -module=${moddir} -static-intel -mp1 -stand=f23 -warn=nounused $(DEFINES) $(FINCLUDES)
    ifeq ($(WITH_CHECKS),true)
      FFLAGS+=-check=all -check=noarg_temp_created
    endif
    FFLAGS+=-real-size=64 -fpscomp
  else
    FFLAGS=-fpp -warn all -module ${moddir} -static-intel -mp1 -stand f23 -warn nounused $(DEFINES) $(FINCLUDES)
    LINK=$(FC) -nofor-main
    ifeq ($(WITH_CHECKS),true)
      FFLAGS+=-check all -check=noarg_temp_created
    endif
    FFLAGS+=-real-size 64 -fpscomp
    FLIBS+=-L/opt/intel/lib
    FLIBS+=-lifcore -lsvml -lifport
    FLIBS+=-limf -lintlc -liomp5 -lifport
  endif
else ifeq ($(F90),flang)
  LINK=$(FC)
  DEBUG_FFLAGS=-g -DDEBUG=1
  OPT_FFLAGS=-O3
  FFLAGS=-module-dir ${moddir} $(DEFINES) $(FINCLUDES)
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-Mbounds
  endif
  FFLAGS+=-fdefault-real-8 -fdefault-double-8
  FLIBS+=-L../ancillary/lib
else
  LINK=$(FC)
  DEBUG_FFLAGS=-g -fbacktrace -DDEBUG=1
  OPT_FFLAGS=-O3
  FFLAGS=-Wall -J ${moddir} -Wno-c-binding-type -ffree-line-length-none -std=f2018 -fall-intrinsics $(DEFINES) $(FINCLUDES)
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-fcheck=all,no-array-temps
  endif
  FFLAGS+=-fdefault-real-8 -fdefault-double-8
  OMPFLAG=-fopenmp
  FLIBS+=-lgfortran -lgomp
 #EXTFFLAGS+=-Wno-unused-dummy-argument -Wno-unused-value
endif

ifneq ($(USE_DL),true)
  WQLIBS=$(AEDLIBS) $(FABMLIBS) $(AED2LIBS)
endif

ifeq ($(DEBUG),true)
  ifeq ($(OSTYPE),FreeBSD)
    DEBUG_CFLAGS=-g -DDEBUG=1
  else
    DEBUG_CFLAGS=-g -fbounds-check -DDEBUG=1
  endif
  OPT_CFLAGS=
  OPT_FFLAGS=
else
  DEBUG_FFLAGS=
  DEBUG_CFLAGS=
  # OPT_CFLAGS=-O4 -Ofast -frounding-math
  OPT_CFLAGS=-O3
  # OPT_CFLAGS=
  # OPT_FFLAGS=
endif

LIBS+=-lnetcdf
# If variable NETCDFLIB is not empty, use it to
# set path to the library
ifneq ($(NETCDFLIB),)
  LIBS+=-L$(NETCDFLIB)
endif

ifeq ($(PLOTDIR),)
  PLOTDIR=../../libplot
endif

PLOTLIBS=-L$(PLOTDIR)/lib -lplot -lgd
ifeq ($(OSTYPE),Darwin)
  XLIBS+=-framework Cocoa
else ifeq ($(OSTYPE),Msys)
  XLIBS+=-lgdi32
else
  XLIBS+=-lX11
endif
ifdef AEDDEVDIR
  LIBS+=$(PLOTLIBS) $(XLIBS)
else ifeq ($(WITH_PLOTS),true)
  LIBS+=$(PLOTLIBS)
  ifeq ($(WITH_XPLOTS),true)
    LIBS+=$(XLIBS)
  endif
endif

ifeq ($(OSTYPE),Msys)
  RES=${objdir}/glm_rc.o
  RESP=${objdir}/glm+_rc.o
else
  RES=
  RESP=
endif

ifeq ($(FENCE),true)
  LIBS+=-lefence
endif
ifeq ($(MDEBUG),true)
  DEBUG_CFLAGS+=-fsanitize=address
  LIBS+=$(DBG_LIBS)
endif

CFLAGS=-Wall -I$(UTILDIR) -I$(PLOTDIR) $(CINCLUDES) $(DEFINES) $(DEBUG_CFLAGS) $(OPT_CFLAGS)
CFLAGS+=-fPIC
FFLAGS+=$(DEBUG_FFLAGS) $(OPT_FFLAGS)
FFLAGS+=-fPIC

OBJS=${objdir}/glm_globals.o \
     ${objdir}/glm_util.o \
     ${objdir}/glm_csv.o \
     ${objdir}/glm_mobl.o \
     ${objdir}/glm_mixu.o \
     ${objdir}/glm_wqual.o \
     ${objdir}/glm_ptm.o \
     ${objdir}/glm_layers.o \
     ${objdir}/glm_surface.o \
     ${objdir}/glm_input.o \
     ${objdir}/glm_plot.o \
     ${objdir}/glm_output.o \
     ${objdir}/glm_ncdf.o \
     ${objdir}/glm_restart.o \
     ${objdir}/glm_lnum.o \
     ${objdir}/glm_init.o \
     ${objdir}/glm_flow.o \
     ${objdir}/glm_mixer.o \
     ${objdir}/glm_deep.o \
     ${objdir}/glm_stress.o \
     ${objdir}/glm_bird.o \
     ${objdir}/glm_model.o \
     ${objdir}/glm_types.o \
     ${objdir}/glm_const.o \
     ${objdir}/glm_debug.o \
     ${objdir}/glm_balance.o \
     ${objdir}/glm_heatexchange.o \
     ${objdir}/glm_oxygenation.o

GLMOBJS=${objdir}/glm_main.o

LIBOBJS=${objdir}/glm_lib.o

ifeq ($(USE_DL),true)
  LIBS+=-ldl
  CFLAGS+=-DUSE_DL_LOADER=1
  FFLAGS+=-DUSE_DL_LOADER=1
  TARGETS+=$(AEDTARGETS) $(FABMTARGETS)
else
  ifeq ($(WITH_AED),true)
    OBJS+=${objdir}/glm_zones.o
  endif
  ifeq ($(WITH_API),true)
    OBJS+=${objdir}/glm_api_zones.o \
          ${objdir}/glm_api_aed.o
  endif
  ifeq ($(WITH_AED),true)
    OBJS+=${objdir}/glm_aed.o \
          ${objdir}/aed_external.o
  endif
endif

all: $(TARGETS)

${libdir}:
	@mkdir ${libdir}

${objdir}:
	@mkdir ${objdir}

${moddir}:
	@mkdir ${moddir}

glm: ${objdir} ${moddir} $(OBJS) $(GLMOBJS) $(GLM_DEPS) $(RES)
	$(LINK) -o $@ $(EXTRALINKFLAGS) $(OBJS) $(GLMOBJS) $(RES) $(WQLIBS) $(LIBS) $(FLIBS)

glm+: ${objdir} ${moddir} $(OBJS) $(GLMOBJS) $(GLM_DEPS) $(RESP)
	$(LINK) -o $@ $(EXTRALINKFLAGS) $(OBJS) $(GLMOBJS) $(RESP) $(WQLIBS) $(LIBS) $(FLIBS)

# Shared library for Python/ctypes (libglm.so or libglm.dylib)
# Build after: make glm (or glm+). Requires AED libs. Use: make libglm.so WITH_PLOTS=false
libglm.${so_ext}: ${objdir} ${moddir} $(OBJS) $(LIBOBJS) $(GLM_DEPS)
	$(FC) ${SHARED} -o $@ $(OBJS) $(LIBOBJS) $(RES) $(WQLIBS) $(LIBS) $(FLIBS)

libglm+.${so_ext}: ${objdir} ${moddir} $(OBJS) $(LIBOBJS) $(GLM_DEPS)
	$(FC) ${SHARED} -o $@ $(OBJS) $(LIBOBJS) $(RESP) $(WQLIBS) $(LIBS) $(FLIBS)

clean: ${objdir} ${moddir}
	@touch ${objdir}/1.o ${moddir}/1.mod 1.t 1__genmod.f90 glm 1.${so_ext} glm_test_bird macos/glm.app macos/glm+.app
	@touch debian/.debhelper debian/files debian/control
	@touch debian/glm debian/glm.debhelper.log debian/glm.substvars
	@touch debian/glm+ debian/glm+.debhelper.log debian/glm+.substvars
	@/bin/rm -f glm glm+ libglm.${so_ext}
	@/bin/rm ${moddir}/*.mod ${objdir}/*.o *.t *__genmod.f90 *.${so_ext} glm_test_bird
	@/bin/rm -rf debian/.debhelper debian/files debian/control
	@/bin/rm -rf debian/glm debian/glm.debhelper.log debian/glm.substvars
	@/bin/rm -rf debian/glm+ debian/glm+.debhelper.log debian/glm+.substvars
	@echo Made clean

distclean: clean
	@/bin/rm -rf ${objdir} ${moddir} glm glm+ libglm.${so_ext} macos/glm.app macos/glm+.app

${objdir}/%.o: ${srcdir}/%.c ${incdir}/glm.h
	$(CC) $(CFLAGS) $(EXTRA_FLAGS) -c $< -o $@

${objdir}/%.o: ${srcdir}/%.F90 ${incdir}/glm.h
	$(FC) $(FFLAGS) $(EXTFFLAGS) -c $< -o $@

${objdir}/glm_main.o: ${srcdir}/glm_main.c ${incdir}/glm.h
	$(CC) -DBUILDDATE=\"${BUILDDATE}\" $(CFLAGS) $(EXTRA_FLAGS) -c $< -o $@

${objdir}/glm_rc.o: win/glm.rc
	windres $< -o $@

${objdir}/glm+_rc.o: win/glm+.rc
	windres $< -o $@

%.${so_ext}:
	$(LD) ${SHARED} $(LDFLAGS) \
                        -o $@ $^ -L/opt/intel/lib/intel64/ $(LIBS)

#                        -E -Bdynamic -undefined suppress -o $@ $^ -L/opt/intel/lib/intel64/ $(LIBS)

# Build rules

libglm_wq_aed.${so_ext}: ${objdir}/glm_zones.o ${objdir}/glm_aed.o ${objdir}/glm_plugin.o
	$(CC) ${SHARED} $(LDFLAGS) -o $@ $^ $(AEDLIBS)

libglm_wq_aed+.${so_ext}: ${objdir}/glm_zones.o ${objdir}/glm_aed.o ${objdir}/glm_plugin.o
	$(CC) ${SHARED} $(LDFLAGS) -o $@ $^ $(AEDPLBS)

libglm_wq_fabm.${so_ext}: ${objdir}/glm_zones.o ${objdir}/glm_fabm.o ${objdir}/ode_solvers.o ${objdir}/glm_plugin.o
	$(CC) ${SHARED} $(LDFLAGS) -o $@ $^ $(FABMLIBS)

# special needs dependancies

${objdir}/aed_external.o: ../libaed-water/src/aed_external.F90
	$(FC) $(FFLAGS) $(EXTFFLAGS) $(OMPFLAG) -c $< -o $@

${objdir}/glm_globals.o: ${srcdir}/glm_globals.c ${incdir}/glm_globals.h ${incdir}/glm.h
${objdir}/glm_plugin.o: ${srcdir}/glm_plugin.c ${incdir}/glm_plugin.h ${incdir}/glm.h
${objdir}/glm_mixer.o: ${srcdir}/glm_mixer.c ${incdir}/glm_mixer.h ${incdir}/glm.h ${incdir}/glm_debug.h
