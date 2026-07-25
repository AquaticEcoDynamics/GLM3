#!/bin/sh

# CURDIR should be the directory of the project we are building
export CURDIR=`pwd`
# CWD should be the tools directory in which CURDIR lives
export CWD=`dirname ${CURDIR}`

#
# These are defaults for glm
#
export WITH_AED=true
export AED=true
export WITH_AED_PLUS=false
if [ -d ${CWD}/libaed-dev ] ; then
  export WITH_AED_PLUS=true
fi
export WITH_API=true
export API=true
export USE_DL=false
export WITH_PLOTS=true
export WITH_XPLOTS=true
export WITH_MPI=false


export PLOTDIR=${CWD}/libplot
export UTILDIR=${CWD}/libutil

case `uname` in
  "Darwin"|"Linux"|"FreeBSD")
    export OSTYPE=`uname -s`
    ;;
  MINGW*)
    export OSTYPE="Msys"
    ;;
esac

case $OSTYPE in
  Darwin)
    export LIB_EXT=dylib
    ;;
  Msys)
    export LIB_EXT=dll
    ;;
  *)
    export LIB_EXT=so
    ;;
esac

if [ "$OSTYPE" = "FreeBSD" ] ; then
  export FC=flang
  export CC=clang
  export MAKE=gmake
else
  export FC=gfortran
  export CC=gcc
  export MAKE=make
fi

ARGS=""
while [ $# -gt 0 ] ; do
  ARGS="$ARGS $1"
  case $1 in
    --debug)
      export DEBUG=true
      ;;
    --checks)
      export WITH_CHECKS=true
      ;;
    --mdebug)
      export MDEBUG=true
      ;;
    --fence)
      export FENCE=true
      ;;
    --with-aed)
      export WITH_AED=true
      ;;
    --without-aed)
      export WITH_AED=false
      ;;
    --with-aed-plus)
      export WITH_AED_PLUS=true
      ;;
    --without-aed-plus)
      export WITH_AED_PLUS=false
      ;;
    --with-lib)
      export WITH_LIB=true
      ;;
    --without-lib)
      export WITH_LIB=false
      ;;
    --gfort)
      export FC=gfortran
      ;;
    --ifx)
      export FC=ifx
      ;;
    --ifort)
      export FC=ifort
      ;;
    --clang)
      export CC=clang
      ;;
    --flang)
      export FC=flang
      ;;
    --no-gui)
      export WITH_PLOTS=false
      export WITH_XPLOTS=false
      ;;
    --auto-prereq)
      export AUTO_PREQ=true
      ;;
    --help)
      echo "build_glm accepts the following flags:"
      echo "  --debug            : build with debugging symbols"
      echo "  --gfort            : use the gfortran compiler"
      echo "  --ifort            : use the older intel fortran compiler"
      echo "  --ifx              : use the newer intel fortran compiler"
      echo "  --clang            : use the clang C/C++ compiler"
      echo "  --flang            : use the flang fortran compiler"
      echo
      echo "  --with-aed-plus    : build with aed and aed-plus enabled"
      echo "  --without-aed-plus : build without aed and aed-plus enabled"
      echo "  --with-lib         : build library (libglm) as well"
      echo "  --without-lib      : dont build libglm (default)"
      echo
      echo "  --auto-prereq      : if needed, also build ancillary pre-requirments"
      echo

      exit 0
      ;;
    *)
      ;;
  esac
  shift
done

export F77=$FC
export F90=$FC
export F95=$FC

export MPI=OPENMPI

. ${CWD}/build_env.inc

if [ "$PLOTDIR" = "" ] ; then
  export PLOTDIR=../libplot
fi
if [ "$UTILDIR" = "" ] ; then
  export UTILDIR=../libutil
fi

if [ "$FABM" = "true" ] ; then
  if [ ! -d $FABMDIR ] ; then
    echo "FABM directory not found"
    export FABM=false
  else
    which cmake > /dev/null 2>&1
    if [ $? -ne 0 ] ; then
      echo "cmake not found - FABM cannot be built"
      export FABM=false
    fi
  fi
  if [ "$FABM" = "false" ] ; then
    echo build with FABM requested but FABM cannot be built
    exit 1
  fi

  export FABMHOST=glm
  cd ${FABMDIR}
  if [ ! -d build ] ; then
    mkdir build
  fi
  cd build
# export FFLAGS="$FFLAGS -fPIC"
  if [ "${USE_DL}" = "true" ] ; then
    cmake ${FABMDIR} -DBUILD_SHARED_LIBS=1 || exit 1
  else
    cmake ${FABMDIR} || exit 1
  fi
  ${MAKE} || exit 1
fi

if [ "${WITH_AED}" = "true" ] || [ "${WITH_API}" = "true" ] ; then
  . ${CWD}/build_aedlibs.inc
fi

if [ -d "${UTILDIR}" ] ; then
  echo "making libutil"
  cd "${UTILDIR}"
  ${MAKE} || exit 1
  cd "${CWD}"
fi

if [ "$OSTYPE" = "Msys" ] ; then
  if [ ! -d ancillary/lib ] ; then
    echo making windows ancillary extras
    cd ancillary
    ./build.sh || exit 1
  fi
fi

if [ "$WITH_PLOTS" = "true" ] ; then
  echo "making libplot"
  cd "${PLOTDIR}"
  ${MAKE} || exit 1
fi

cd "${CURDIR}"
if [ -f obj/aed_external.o ] ; then
  /bin/rm obj/aed_external.o
fi

# Update versions in resource files
if [ -f src/glm.h ] ; then
  VERSION=`grep GLM_VERSION src/glm.h | cut -f2 -d\"`
else
  VERSION=`grep GLM_VERSION include/glm.h | cut -f2 -d\"`
fi
cd "${CURDIR}/win"
${CURDIR}/vers.sh $VERSION
#cd ${CURDIR}/win-dll
#${CURDIR}/vers.sh $VERSION
cd "${CURDIR}"
get_commit_id >> ${CWD}/cur_state.log

export LIBRARY_PATH=$LIB
if [ "$WITH_LIB" = "true" ] ; then
  LIBTARG="libglm.${LIB_EXT}"
else
  LIBTARG=""
fi
${MAKE} glm $LIBTARG AEDBENDIR=$DAEDBENDIR AEDDMODIR=$DAEDDMODIR || exit 1
if [ "${DAEDDEVDIR}" != "" ] ; then
  if [ -d "${DAEDDEVDIR}" ] ; then
    echo now build plus version
    /bin/rm obj/aed_external.o
    /bin/rm obj/glm_main.o
    if [ "$WITH_LIB" = "true" ] ; then
      LIBTARG="libglm+.${LIB_EXT}"
    else
      LIBTARG=""
    fi
    ${MAKE} glm+ $LIBTARG WITH_AED_PLUS=1 \
                     AEDBENDIR=$DAEDBENDIR AEDDMODIR=$DAEDDMODIR \
                     AEDRIPDIR=$DAEDRIPDIR AEDLGTDIR=$DAEDLGTDIR \
                     AEDDEVDIR=$DAEDDEVDIR PHREEQDIR=$PHREEQDIR || exit 1
  fi
fi

cd "${CWD}"

# =====================================================================

exit 0
