#!/bin/sh

if [ "$GLM_CONFIGURED" != "true" ] ; then
  . ./GLM_CONFIG
fi

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
    --fabm)
      export FABM=true
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

. ${CWD}/build_env.inc

export F77=$FC
export F90=$FC
export F95=$FC

export MPI=OPENMPI

if [ "$AED2DIR" = "" ] ; then
  export AED2DIR=../libaed2
fi
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

if [ "${AED2}" = "true" ] ; then
  cd "${AED2DIR}"
  ${MAKE} || exit 1
  cd ..
  if [ "${AED2PLS}" != "" ] ; then
    if [ -d "${AED2PLS}" ] ; then
      cd "${AED2PLS}"
      ${MAKE} || exit 1
      cd ..
    fi
  fi
fi

if [ "${AED}" = "true" ] || [ "${API}" = "true" ] ; then
  export WITH_AED_PLUS='true'
  . ${CWD}/build_aedlibs.inc
fi

if [ -d "${UTILDIR}" ] ; then
  echo "making libutil"
  cd "${UTILDIR}"
  ${MAKE} || exit 1
  cd "${CURDIR}/.."
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

${MAKE} AEDBENDIR=$DAEDBENDIR AEDDMODIR=$DAEDDMODIR || exit 1
if [ "${DAEDDEVDIR}" != "" ] ; then
  if [ -d "${DAEDDEVDIR}" ] ; then
    echo now build plus version
    /bin/rm obj/aed_external.o
    ${MAKE} glm+ AEDBENDIR=$DAEDBENDIR AEDDMODIR=$DAEDDMODIR AEDRIPDIR=$DAEDRIPDIR AEDLGTDIR=$DAEDLGTDIR AEDDEVDIR=$DAEDDEVDIR PHREEQDIR=$PHREEQDIR || exit 1
  fi
fi

cd "${CURDIR}/.."

# =====================================================================
# Package building bit
#
# Disabled: building a .deb here requires dh_testdir/fakeroot (debhelper)
# to be installed, which a plain build environment doesn't have by
# default. Uncomment if you actually want to produce a .deb from this
# script; otherwise `make` alone already leaves glm/glm+ built above.

# ***************************** Linux *********************************
# if [ "$OSTYPE" = "Linux" ] ; then
#   if [ $(lsb_release -is) = Ubuntu ] ; then
#     BINPATH=binaries/ubuntu/$(lsb_release -rs)
#     if [ ! -d "${BINPATH}" ] ; then
#       mkdir -p "${BINPATH}"/
#     fi
#     cd ${CURDIR}
#     if [ -x glm+ ] ; then
#        /bin/cp debian/control-with+ debian/control
#     else
#        /bin/cp debian/control-no+ debian/control
#     fi
#     VERSDEB=`head -1 debian/changelog | cut -f2 -d\( | cut -f1 -d-`
#     echo debian version $VERSDEB
#     if [ "$VERSION" != "$VERSDEB" ] ; then
#       echo updating debian version
#       dch --newversion ${VERSION}-0 "new version ${VERSION}"
#     fi
#     VERSRUL=`grep 'version=' debian/rules | cut -f2 -d=`
#     if [ "$VERSION" != "$VERSRUL" ] ; then
#       sed -i "s/version=$VERSRUL/version=$VERSION/" debian/rules
#     fi
#
#     fakeroot ${MAKE} -f debian/rules binary || exit 1
#
#     cd ..
#
#     mv glm*.deb ${BINPATH}/
#   else
#     BINPATH="binaries/$(lsb_release -is)/$(lsb_release -rs)"
#     echo "No package build for $(lsb_release -is)"
#   fi
# fi

exit 0
