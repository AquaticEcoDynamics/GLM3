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
    --gfort)
      export FC=gfortran
      ;;
    --ifx)
      export FC=ifx
      ;;
    --ifort)
      export FC=ifort
      ;;
    --flang-new)
      export FC=flang-new
      ;;
    --flang)
      export FC=flang
      ;;
    --no-gui)
      export WITH_PLOTS=false
      export WITH_XPLOTS=false
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
    if [ $? != 0 ] ; then
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

if [ "$OSTYPE" = "FreeBSD" ] ; then
  echo not making flang extras
  # cd ancillary/freebsd
  # ./fetch.sh
  # ${MAKE} || exit 1
elif [ "$OSTYPE" = "Msys" ] ; then
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
${MAKE} AEDBENDIR=$DAEDBENDIR AEDDMODIR=$DAEDDMODIR || exit 1
if [ "${DAEDDEVDIR}" != "" ] ; then
  if [ -d "${DAEDDEVDIR}" ] ; then
    echo now build plus version
    /bin/rm obj/aed_external.o
    /bin/rm obj/glm_main.o
    ${MAKE} glm+ WITH_AED_PLUS=1 AEDBENDIR=$DAEDBENDIR AEDDMODIR=$DAEDDMODIR \
                                 AEDRIPDIR=$DAEDRIPDIR AEDLGTDIR=$DAEDLGTDIR \
                                 AEDDEVDIR=$DAEDDEVDIR PHREEQDIR=$PHREEQDIR || exit 1
  fi
fi

cd "${CWD}"

# =====================================================================

exit 0
