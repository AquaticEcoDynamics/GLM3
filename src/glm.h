/******************************************************************************
 *                                                                            *
 * glm.h                                                                      *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 2013-2026 : The University of Western Australia                  *
 *                                                                            *
 *  This file is part of GLM (General Lake Model)                             *
 *                                                                            *
 *  GLM is free software: you can redistribute it and/or modify               *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  GLM is distributed in the hope that it will be useful,                    *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                            *
 ******************************************************************************/
#ifndef _GLM_H_
#define _GLM_H_

#ifndef DEBUG
#define DEBUG         0
#endif
#define USE_FILLVALUE 1

#define GLM_VERSION  "4.0.0"

#define POINT         0
#define Z_SHAPE       1
#define T_SHAPE       2
#define XYT_SHAPE     4
#define XYZT_SHAPE    5
#define XYNT_SHAPE    6
#define R_SHAPE       1

#ifndef PATH_MAX
#define PATH_MAX  1024
#endif

#ifndef __STDC__
!-------------------------------------------------------------------------------
!! Fortran version

!#  define surfLayer NumLayers
#  define botmLayer 1

#  ifndef AED_REAL
#    define AED_REAL REAL(kind=C_DOUBLE)
#  endif
#  define CAED_REAL REAL(kind=C_DOUBLE)
#  define NF90_REALTYPE NF90_DOUBLE
#  define NC_FILLER NC_FILL_DOUBLE

#  ifndef FLOAT
#    define FLOAT(x) (x)
#  endif
#  define DOUBLETYPE double precision
#  define CINTEGER INTEGER(kind=C_INT32_T)
#  define CSIZET   INTEGER(kind=C_SIZE_T)
#  define FLOGICAL LOGICAL(kind=C_BOOL)
#  define CLOGICAL LOGICAL(kind=C_BOOL)
#  define CCHARACTER CHARACTER(C_CHAR)

#  define stdin  5
#  define stdout 6
#  define stderr 0

!-------------------------------------------------------------------------------
#else
//------------------------------------------------------------------------------
  // C Version of header

# include <stdint.h>

  #define surfLayer (NumLayers-1)
  #define botmLayer 0

  #define onshoreLayer (surfLayer+2)
  #define offshoreLayer (surfLayer+1)

// This is how we sort out the differences in multidimensional array indexing
// the order of indeces has been reversed in this version of GLM, but generally
//  the only lace you will know is in macros like this and when allocating arrays
//#define WQ_INF_(a,i,j) a[_IDX_2d(MaxInf,MaxVars,i,j)]
  #define WQ_INF_(a,i,j) a[_IDX_2d(MaxVars,MaxInf,j,i)]

  #define AMOD fmod
  typedef double AED_REAL;
  #define NC_REALTYPE NC_DOUBLE
  #define NC_FILLER NC_FILL_DOUBLE
  typedef double DOUBLETYPE;
  typedef int CLOGICAL;
  typedef _Bool FLOGICAL;
  typedef int32_t CINTEGER;
  typedef int8_t  CCHARACTER;

  #define TRUE  1
  #define FALSE 0

  #ifdef _WIN32
    #if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
       #define _VISUAL_C_
    #endif
    #define snprintf _snprintf
    #define strcasecmp _stricmp
    #define strncasecmp _strnicmp
    double fmod(double x, double y);
    char *ctime_r(const time_t *t, char *b);
  #endif

  #if DEBUG
    #define CRASH(s) ( { int *x = (int*)1; fputs(s, stderr); *x = 1; } )
  #endif

//------------------------------------------------------------------------------
#endif

#if USE_FILLVALUE
#  define PARAM_FILLVALUE   ,NC_FILLER
#else
#  define PARAM_FILLVALUE
#endif

#define sqr(x)  ((x)*(x))
#define gprime(d1,d2) (((d2)-(d1))*g/(((d1)+(d2))/2.0))
#define combine_vol(c1,v1,c2,v2) (((c1)*(v1)+(c2)*(v2))/((v1)+(v2)))


#define MISVAL -9999.
#ifndef _ZERO_
#define _ZERO_ 0.
#endif

#define LW_CC    1
#define LW_IN    2
#define LW_NET   3

#endif
