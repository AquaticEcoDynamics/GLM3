!###############################################################################
!#                                                                             #
!# glm_aed.F90                                                                 #
!#                                                                             #
!# The interface between glm and libaed-xxx                                    #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Agriculture and Environment                                   #
!#     The University of Western Australia                                     #
!#                                                                             #
!#     http://aquatic.science.uwa.edu.au/                                      #
!#                                                                             #
!# Copyright 2013-2026 : The University of Western Australia                   #
!#                                                                             #
!#  This file is part of GLM (General Lake Model)                              #
!#                                                                             #
!#  GLM is free software: you can redistribute it and/or modify                #
!#  it under the terms of the GNU General Public License as published by       #
!#  the Free Software Foundation, either version 3 of the License, or          #
!#  (at your option) any later version.                                        #
!#                                                                             #
!#  GLM is distributed in the hope that it will be useful,                     #
!#  but WITHOUT ANY WARRANTY; without even the implied warranty of             #
!#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
!#  GNU General Public License for more details.                               #
!#                                                                             #
!#  You should have received a copy of the GNU General Public License          #
!#  along with this program.  If not, see <http://www.gnu.org/licenses/>.      #
!#                                                                             #
!###############################################################################

#include "aed.h"

#undef MISVAL

#include "glm.h"

!-------------------------------------------------------------------------------
MODULE glm_aed
!
   USE ISO_C_BINDING
   USE IEEE_ARITHMETIC

   USE aed_util
   USE aed_common
   USE glm_types
   USE glm_zones

   IMPLICIT NONE

   PRIVATE ! By default, make everything private
!
#include "glm_globals.h"
#include "glm_plot.h"
#include "glm_ncdf.h"
#include "glm_csv.h"
#include "glm_mobl.h"
!
#if USE_DL_LOADER
# define _WQ_INIT_GLM        "wq_init_glm"
# define _WQ_SET_GLM_DATA    "wq_set_glm_data"
# define _WQ_DO_GLM          "wq_do_glm"
# define _WQ_CLEAN_GLM       "wq_clean_glm"
# define _WQ_INIT_GLM_OUTPUT "wq_init_glm_output"
# define _WQ_WRITE_GLM_      "wq_write_glm"
# define _WQ_VAR_INDEX_C     "wq_var_index_c"
# define _WQ_IS_VAR          "wq_is_var"
# define _WQ_INFLOW_UPDATE   "wq_inflow_update"
#else
# define _WQ_INIT_GLM        "aed_init_glm"
# define _WQ_SET_GLM_DATA    "aed_set_glm_data"
# define _WQ_DO_GLM          "aed_do_glm"
# define _WQ_CLEAN_GLM       "aed_clean_glm"
# define _WQ_INIT_GLM_OUTPUT "aed_init_glm_output"
# define _WQ_WRITE_GLM_      "aed_write_glm"
# define _WQ_VAR_INDEX_C     "aed_var_index_c"
# define _WQ_IS_VAR          "aed_is_var"
# define _WQ_INFLOW_UPDATE   "aed_update_inflow_wq"
#endif
!
!-------------------------------------------------------------------------------
!
!MODULE DATA

   AED_REAL :: par_fraction =  0.450
   AED_REAL :: nir_fraction =  0.510
   AED_REAL :: uva_fraction =  0.035
   AED_REAL :: uvb_fraction =  0.005

   !# Arrays for state and diagnostic variables
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,TARGET :: cc    !# water quality array: [nvars, nlayers]
!  AED_REAL,DIMENSION(:),  ALLOCATABLE,TARGET :: cc_hz !# water quality array - benthic: [nvars]
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,TARGET :: cc_diag
   AED_REAL,DIMENSION(:,:),ALLOCATABLE        :: cc_diag_old
   AED_REAL,DIMENSION(:),  ALLOCATABLE,TARGET :: cc_diag_hz
   AED_REAL,DIMENSION(:),  ALLOCATABLE        :: z_prev

   !# Arrays for work, vertical movement, and cross-boundary fluxes
   AED_REAL,DIMENSION(:,:),ALLOCATABLE :: ws
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: dz
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: tss
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: sed_zones

   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: depth
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: layer_area  ! incremental area per layer (area(i)-area(i-1))
   AED_REAL,TARGET :: col_depth
   AED_REAL,TARGET :: col_num = 1

   !# Arrays for environmental variables not supplied externally.
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: par
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: pres
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: uva
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: uvb
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: nir

   !# External variables
   AED_REAL :: dt_eff       ! External and internal time steps
   AED_REAL,DIMENSION(:),POINTER :: rad
   AED_REAL,DIMENSION(:),POINTER :: height
   AED_REAL,DIMENSION(:),POINTER :: salt
   AED_REAL,DIMENSION(:),POINTER :: temp
   AED_REAL,DIMENSION(:),POINTER :: rho
   AED_REAL,DIMENSION(:),POINTER :: area
   AED_REAL,DIMENSION(:),POINTER :: extc
   AED_REAL,DIMENSION(:),POINTER :: layer_stress
   AED_REAL,DIMENSION(:),POINTER :: vel

   AED_REAL,POINTER :: precip
   AED_REAL,POINTER :: delzBlueIce
   AED_REAL,POINTER :: delzWhiteIce
   AED_REAL,POINTER :: u_star_glm
   AED_REAL,POINTER :: Q_net_glm
   AED_REAL,POINTER :: evap
   AED_REAL,POINTER :: bottom_stress
   AED_REAL,POINTER :: air_temp
   AED_REAL,POINTER :: rel_hum
   AED_REAL,POINTER :: I_0
   AED_REAL,POINTER :: wnd
   AED_REAL,POINTER :: air_pres

!  TYPE(aed_column_t),DIMENSION(:,:),ALLOCATABLE,TARGET :: all_cols !# (n_aed_vars, ncols)
!  TYPE(aed_column_t),DIMENSION(:,:),ALLOCATABLE,TARGET :: zon_cols !# (n_aed_vars, nzones)
   AED_REAL,ALLOCATABLE,TARGET :: flux_ben(:)          !# (n_vars+n_vars_ben)
   AED_REAL,ALLOCATABLE,TARGET :: flux_atm(:)          !# (n_vars+n_vars_ben)
   AED_REAL,ALLOCATABLE,TARGET :: flux_zon(:,:)        !# (n_vars+n_vars_ben, aed_n_zones)
   AED_REAL,ALLOCATABLE,TARGET :: flux_pel(:,:)        !# (n_vars+n_vars_ben, MAX(n_layers, aed_n_zones))
   AED_REAL,DIMENSION(:,:),ALLOCATABLE :: flux_pel_pre !# (n_vars+n_vars_ben, MAX(n_layers, aed_n_zones))
   AED_REAL,DIMENSION(:,:),ALLOCATABLE :: flux_pel_z   !# (n_vars+n_vars_ben, MAX(n_layers, aed_n_zones))

   CHARACTER(len=48),ALLOCATABLE :: names(:)
   CHARACTER(len=48),ALLOCATABLE :: bennames(:)
!  CHARACTER(len=48),ALLOCATABLE :: diagnames(:)

   INTEGER,DIMENSION(:),ALLOCATABLE :: externalid
   INTEGER,DIMENSION(:),ALLOCATABLE :: zexternalid
#ifdef PLOTS
   INTEGER,DIMENSION(:),ALLOCATABLE :: plot_id_v, plot_id_sv, plot_id_d, plot_id_sd
#endif

   LOGICAL :: reinited = .FALSE.

   INTEGER :: n_aed_vars, n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet
   INTEGER :: zone_var = 0

   CHARACTER(len=64) :: NULCSTR = ""

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed_init_glm(i_fname, len, NumWQ_Vars, NumWQ_Ben)                   &
                                                      BIND(C, name=_WQ_INIT_GLM)
!-------------------------------------------------------------------------------
! Initialize the GLM-AED driver by reading settings from aed.nml.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CCHARACTER,INTENT(in) :: i_fname(*)
   CSIZET,INTENT(in)     :: len
   CINTEGER,INTENT(out)  :: NumWQ_Vars, NumWQ_Ben
!
!LOCALS
   INTEGER :: i,j,namlst,status
   INTEGER :: rc, av, v, sv, tv

   CHARACTER(len=80) :: fname
   TYPE(aed_variable_t),POINTER :: tvar

   CHARACTER(len=64) :: models(64)
   NAMELIST /aed_models/ models
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL make_string(fname, i_fname, len)

#ifdef __INTEL_COMPILER
#  ifdef __INTEL_LLVM_COMPILER
   print *,'    glm_aed built using intel fortran (ifx) version ', __INTEL_LLVM_COMPILER
#  else
   print *,'    glm_aed built using intel fortran version ', __INTEL_COMPILER
#  endif
#else
# ifdef __PGI
   print *,'    glm_aed built using pgfortran version ', __PGIC__, '.', __PGIC_MINOR__, '.', __PGIC_PATCHLEVEL__
# else
#  ifdef __GNUC__
    print *,'    glm_aed built using gfortran version ', __GNUC__, '.', __GNUC_MINOR__, '.', __GNUC_PATCHLEVEL__
#  else
#   ifdef __clang__
     print*,"    glm_aed built using flang version ", __clang_major__, '.', __clang_minor__, '.', __clang_patchlevel__
#   else
     print*,"    glm_aed built using unknow fortran"
#   endif
#  endif
# endif
#endif
   print *,'    libaed enabled.... init_glm_aed processing: ', TRIM(fname)
   namlst = find_free_lun()

   write(*,"(/,5X,'---------- AED config : start ----------')")
   IF ( aed_init_core('.') /= 0 ) STOP "     ERROR: Initialisation of aed_core failed"
   CALL aed_print_version

   tv = aed_provide_global('temperature', 'temperature',   'celsius')
   tv = aed_provide_global('salinity',    'salinity',      'g/Kg')
   tv = aed_provide_global('density',     'density',       '')
   tv = aed_provide_global('layer_ht',    'layer heights', 'meters')
   tv = aed_provide_global('layer_area',  'layer area',    'm2')
   tv = aed_provide_sheet_global('rain',  'rainfall',      'm/s')
   !rainloss
   !material
   !bathy
   tv = aed_provide_global('extc_coef', 'extinction coefficient', '')
   tv = aed_provide_global('tss',       'tss', '')
   !ss1, ss2, ss3, ss4
   tv = aed_provide_global('cell_vel',  'layer velocity', 'm/s')
   tv = aed_provide_global('nir',       'nir',      'W/m2')
   tv = aed_provide_global('par',       'par',      'W/m2')
   tv = aed_provide_global('uva',       'uva',      'W/m2')
   tv = aed_provide_global('uvb',       'uvb',      'W/m2')
   tv = aed_provide_global('pressure',  'pressure', '')
   tv = aed_provide_global('depth',     'depth',    'm')

   tv = aed_provide_global('sed_zones', 'sediment zones', '')

   tv = aed_provide_sheet_global('sed_zone',   'sediment zone',     '')
   tv = aed_provide_sheet_global('wind_speed', 'wind speed',        'm/s')
   tv = aed_provide_sheet_global('par_sf',     'par_sf',            '')
   tv = aed_provide_sheet_global('taub',       'layer stress',      'N/m2')
   tv = aed_provide_sheet_global('air_temp',   'air temperature',   'celsius')
   tv = aed_provide_sheet_global('humidity',   'relative humidity', '-')
   !longwave
   tv = aed_provide_sheet_global('col_num',    'column number', '')
   !col_depth
   tv = aed_provide_sheet_global('col_depth',  'lake depth',  'meters')
   tv = aed_provide_sheet_global('col_area',   'area at surface',  'm2')
   tv = aed_provide_sheet_global('evap',       'evaporation', 'm/s')
   ! added for oasim
   tv = aed_provide_sheet_global('longitude',  'longitude', 'radians')
   tv = aed_provide_sheet_global('latitude',   'latitude',  'radians')
   tv = aed_provide_sheet_global('yearday',    'yearday',   'day')
   tv = aed_provide_sheet_global('timestep',   'timestep',  'seconds')

   ! new var air_pressure
   tv = aed_provide_sheet_global('air_pres',  'air_pressure',  'Pa')
   tv = aed_provide_sheet_global('delzBlueIce',  'thickness of blue ice',  'm')
   tv = aed_provide_sheet_global('delzWhiteIce',  'thickness of white ice',  'm')
   tv = aed_provide_sheet_global('u_star',  'wind friction velocity',  'm/s')
   tv = aed_provide_sheet_global('Q_net',   'net non-penetrative heat flux',  'W/m2')

   !# Create model tree
   print *,"     Processing aed_models config from ",TRIM(fname)
   OPEN(namlst,file=fname,action='read',status='old',iostat=status)
   IF ( status /= 0 ) CALL STOPIT("Cannot open file " // TRIM(fname))

   models = ''
   READ(namlst, nml=aed_models, iostat=status)
   IF ( status /= 0 ) STOP "Cannot read namelist entry aed_models"

   DO i=1,size(models)
      IF ( models(i) == '' ) EXIT
      IF ( benthic_mode .GT. 1 ) models(i) = TRIM(models(i)) // ':za' ! make all models zone averaged
      CALL aed_define_model(models(i), namlst)
   ENDDO

   !# should be finished with this file
   CLOSE(namlst)
!  print *,"      ... nml file parsing completed."

   n_aed_vars = aed_core_status(n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet)

#if DEBUG
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         print *,"AED var ", i, tvar%sheet, tvar%var_type, TRIM(tvar%name)
      ELSE
         print *,"AED var ", i, " is empty"
      ENDIF
   ENDDO
#endif

   print "(/,5X,'AED : n_aed_vars  = ',I3,' ; MaxLayers         = ',I4)",n_aed_vars,MaxLayers
   print "(  5X,'AED : n_vars      = ',I3,' ; n_vars_ben        = ',I4)",n_vars,n_vars_ben
   print "(  5X,'AED : n_vars_diag = ',I3,' ; n_vars_diag_sheet = ',I4,/)",n_vars_diag,n_vars_diag_sheet

   CALL check_data

   !# names = grab the names from info
   ALLOCATE(names(n_vars),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (names)'
   ALLOCATE(bennames(n_vars_ben),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (bennames)'

#ifdef PLOTS
   ALLOCATE(plot_id_v(n_vars))
   ALLOCATE(plot_id_sv(n_vars_ben))
   ALLOCATE(plot_id_d(n_vars_diag))
   ALLOCATE(plot_id_sd(n_vars_diag_sheet))
   plot_id_v = -1; plot_id_sv = -1; plot_id_d = -1; plot_id_sd = -1
#endif

   NumWQ_Vars = n_vars
   NumWQ_Ben  = n_vars_ben
   !# Now that we know how many vars we need, we can allocate space for them
   ALLOCATE(cc((n_vars + n_vars_ben), MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (CC)'
   cc = 0.         !# initialise to zero

   CALL set_c_wqvars_ptr(cc)
   IF (n_vars_ben > 0) CALL set_c_wqsvars_ptr(cc(n_vars+1, 1))

   print "(5X,'Configured AED variables to simulate:')"

   j = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         IF ( .NOT. tvar%sheet .AND. tvar%var_type == V_STATE ) THEN
            j = j + 1
            names(j) = TRIM(tvar%name)
            !print *,"     S(",j,") AED pelagic(3D) variable: ", TRIM(names(j))
            print "(7X,'S(',I4,') water column variable     : ',A)",j , TRIM(names(j))
         ENDIF
      ENDIF
   ENDDO

   j = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         IF ( tvar%sheet .AND. tvar%var_type == V_STATE ) THEN
            j = j + 1
            bennames(j) = TRIM(tvar%name)
            !print *,"     B(",j,") AED benthic(2D) variable: ", TRIM(bennames(j))
            print "(7X,'B(',I4,') bottom variable           + ',A)",j , TRIM(bennames(j))
         ENDIF
      ENDIF
   ENDDO

   j = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         IF ( tvar%var_type == V_DIAGNOSTIC ) THEN
            IF ( .NOT.  tvar%sheet ) THEN
               j = j + 1
               print "(7X,'D(',I4,') water column diagnostic   > ',A)",j , TRIM(tvar%name)
               !print *,"     D(",j,") AED diagnostic 3Dvariable: ", TRIM(tvar%name)
            ENDIF
         ENDIF
      ENDIF
   ENDDO

   j = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         IF ( tvar%var_type == V_DIAGNOSTIC ) THEN
            IF (tvar%sheet ) THEN
               j = j + 1
               !print *,"     D(",j,") AED diagnostic 2Dvariable: ", TRIM(tvar%name)
               print "(7X,'D(',I4,') bottom/surface diagnostic ~ ',A)",j , TRIM(tvar%name)
            ENDIF
         ENDIF
      ENDIF
   ENDDO

   ALLOCATE(externalid(n_aed_vars))
   ALLOCATE(zexternalid(n_aed_vars))

   !----------------------------------------------------------------------------

   !# Now set initial values
   v = 0 ; sv = 0;
   DO av=1,n_aed_vars
      IF ( .NOT.  aed_get_var(av, tvar) ) STOP "     ERROR getting variable info"
      IF ( tvar%var_type == V_STATE ) THEN  !# neither global nor diagnostic variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
            cc(n_vars+sv, :) = tvar%initial
         ELSE
            v = v + 1
            cc(v, :) = tvar%initial
         ENDIF
      ENDIF
   ENDDO

   !# Allocate diagnostic variable array and set all values to zero.
   !# (needed because time-integrated/averaged variables will increment rather than set the array)
   ALLOCATE(cc_diag(n_vars_diag, MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (cc_diag)'
   cc_diag = zero_

   ALLOCATE(cc_diag_old(n_vars_diag, MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (cc_diag_old)'
   cc_diag_old = zero_

   ALLOCATE(z_prev(MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (z_prev)'
   z_prev = zero_

   !# Allocate diagnostic variable array and set all values to zero.
   !# (needed because time-integrated/averaged variables will increment rather than set the array)
   ALLOCATE(cc_diag_hz(n_vars_diag_sheet),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (cc_diag_hz)'
   cc_diag_hz = zero_

   CALL set_c_wqdvars_ptr(cc_diag, cc_diag_hz, n_vars_diag, n_vars_diag_sheet)

   !# Allocate array with vertical movement rates (m/s, positive for upwards),
   !# and set these to the values provided by the model.
   !# allocated for all vars even though only state vars entries will be used
   ALLOCATE(ws(n_aed_vars, MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (ws)'
   ws = zero_

   !# Allocate array for photosynthetically active radiation (PAR).
   !# This will be calculated internally during each time step.
   ALLOCATE(par(MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (par)'
   par = zero_

   ALLOCATE(nir(MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (nir)'
   nir = zero_
   ALLOCATE(uva(MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (uva)'
   uva = zero_
   ALLOCATE(uvb(MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (uvb)'
   uvb = zero_

   ALLOCATE(dz(MaxLayers),stat=rc)
   dz = zero_

   !# Allocate array for local pressure.
   !# This will be calculated [approximated] from layer depths internally
   !# during each time step.
   ALLOCATE(pres(MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (pres)'
   pres = zero_

   ALLOCATE(tss(MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (tss)'
   tss = zero_

   write(*,"(/,5X,'----------  AED config : end  ----------',/)")
END SUBROUTINE aed_init_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
CINTEGER FUNCTION aed_is_var(id,i_vname,len) BIND(C, name=_WQ_IS_VAR)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in)   :: id
   CCHARACTER,INTENT(in) :: i_vname(*)
   CSIZET,INTENT(in)     :: len
!LOCALS
   CHARACTER(len=45) :: vname
   TYPE(aed_variable_t),POINTER :: tvar
   INTEGER :: i, v, sv, d, sd
!-------------------------------------------------------------------------------
!BEGIN
#ifndef PLOTS
   IF (.FALSE.) PRINT *, id   !# id is only consumed when PLOTS is enabled
#endif
   CALL make_string(vname, i_vname, len)

   v = 0; sv = 0; d = 0; sd = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         IF ( tvar%var_type == V_STATE ) THEN
            IF ( tvar%sheet ) THEN ; sv=sv+1; ELSE ; v=v+1 ; ENDIF
            IF ( TRIM(tvar%name) == vname ) THEN
               IF (tvar%sheet) THEN
                  aed_is_var=-sv
#ifdef PLOTS
                  plot_id_sv(sv) = id;
#endif
               ELSE
                  aed_is_var=v
#ifdef PLOTS
                  plot_id_v(v) = id;
#endif
               ENDIF
               RETURN
            ENDIF
         ELSEIF ( tvar%var_type == V_DIAGNOSTIC ) THEN
            IF ( tvar%sheet ) THEN ; sd=sd+1; ELSE ; d=d+1 ; ENDIF
            IF ( TRIM(tvar%name) == vname ) THEN
               IF (tvar%sheet) THEN
                  aed_is_var=-sd
#ifdef PLOTS
                  plot_id_sd(sd) = id;
#endif
               ELSE
                  aed_is_var=d
#ifdef PLOTS
                  plot_id_d(d) = id;
#endif
               ENDIF
               RETURN
            ENDIF
         ENDIF
      ENDIF
   ENDDO
   aed_is_var = 0
END FUNCTION aed_is_var
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_set_glm_data()                     BIND(C, name=_WQ_SET_GLM_DATA)
!-------------------------------------------------------------------------------
!ARGUMENTS
!LOCALS
!
!-------------------------------------------------------------------------------
!BEGIN
   !# Save pointers to external dynamic variables that we need later (in do_glm_wq)
   height => theLake%Height
   temp   => theLake%Temp
   salt   => theLake%Salinity
   rho    => theLake%Density
   area   => theLake%LayerArea
   rad    => theLake%Light
   vel    => theLake%Umean
   extc   => theLake%ExtcCoefSW
   layer_stress => theLake%LayerStress

   ALLOCATE(depth(MaxLayers))
   ALLOCATE(layer_area(MaxLayers))
   ALLOCATE(sed_zones(MaxLayers))
   sed_zones = 0.

   ALLOCATE(flux_pel(n_vars+n_vars_ben, MAX(MaxLayers,n_zones)))
   ALLOCATE(flux_ben(n_vars+n_vars_ben))
   ALLOCATE(flux_atm(n_vars+n_vars_ben))

   ALLOCATE(flux_pel_pre(n_vars+n_vars_ben, MAX(MaxLayers, n_zones)))
   ALLOCATE(flux_pel_z(  n_vars+n_vars_ben, MAX(MaxLayers, n_zones)))

   IF ( n_zones > 0 ) THEN
      ALLOCATE(flux_zon(n_vars+n_vars_ben, n_zones))
      CALL wq_set_glm_zones(n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet)
   ENDIF

   precip => MetData%Rain
   air_temp => MetData%AirTemp
   rel_hum => MetData%RelHum
   air_pres => MetData%AirPres

   evap         => SurfData%Evap
   delzBlueIce  => SurfData%DelzBlueIce
   delzWhiteIce => SurfData%DelzWhiteIce
   u_star_glm   => SurfData%u_star
   Q_net_glm    => SurfData%Q_net
   bottom_stress => layer_stress(botmLayer)

   !# Provide pointers to arrays with environmental variables to aed.
   wnd => MetData%WindSpeed
   I_0 => MetData%ShortWave

   !# Calculate and save internal time step.
   dt_eff = dt/FLOAT(split_factor)

   !# Trigger an error if WQ hasn't got all it needs from us.
   CALL check_data
END SUBROUTINE aed_set_glm_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE check_data
!-------------------------------------------------------------------------------
! Check that all variable dependencies have been met
!-------------------------------------------------------------------------------
!ARGUMENTS
!
!LOCALS
   INTEGER :: av
   INTEGER :: v, d, sv, sd, ev, err_count
   TYPE(aed_variable_t),POINTER :: tvar
!-------------------------------------------------------------------------------
!BEGIN
   v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
   err_count = 0

   DO av=1,n_aed_vars
      IF ( .NOT.  aed_get_var(av, tvar) ) STOP "Error getting variable info"

      IF ( tvar%var_type == V_EXTERNAL ) THEN !# global variable
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'temperature' ) ; tvar%found = .true.
            CASE ( 'salinity' )    ; tvar%found = .true.
            CASE ( 'density' )     ; tvar%found = .true.
            CASE ( 'layer_ht' )    ; tvar%found = .true.
            CASE ( 'extc_coef' )   ; tvar%found = .true.
            CASE ( 'tss' )         ; tvar%found = .true.
            CASE ( 'cell_vel' )    ; tvar%found = .true.
            CASE ( 'par' )         ; tvar%found = .true.
            CASE ( 'nir' )         ; tvar%found = .true.
            CASE ( 'uva' )         ; tvar%found = .true.
            CASE ( 'uvb' )         ; tvar%found = .true.
            CASE ( 'pressure' )    ; tvar%found = .true.
            CASE ( 'depth' )       ; tvar%found = .true.
            CASE ( 'sed_zones' )   ; tvar%found = .true.
            CASE ( 'sed_zone' )    ; tvar%found = .true.
            CASE ( 'wind_speed' )  ; tvar%found = .true.
            CASE ( 'par_sf' )      ; tvar%found = .true.
            CASE ( 'taub' )        ; tvar%found = .true.
            CASE ( 'col_depth' )   ; tvar%found = .true.
            CASE ( 'col_num' )     ; tvar%found = .true.
            CASE ( 'col_area' )    ; tvar%found = .true.
            CASE ( 'evap' )        ; tvar%found = .true.
            CASE ( 'layer_area' )  ; tvar%found = .true.
            CASE ( 'rain' )        ; tvar%found = .true.
            CASE ( 'delzBlueIce' ) ; tvar%found = .true.
            CASE ( 'delzWhiteIce' ) ; tvar%found = .true.
            CASE ( 'u_star' )      ; tvar%found = .true.
            CASE ( 'Q_net' )       ; tvar%found = .true.
            CASE ( 'air_temp' )    ; tvar%found = .true.
            CASE ( 'air_pres' )    ; tvar%found = .true.
            CASE ( 'humidity' )    ; tvar%found = .true.
            CASE ( 'longitude' )   ; tvar%found = .true.
            CASE ( 'latitude' )    ; tvar%found = .true.
            CASE ( 'yearday' )     ; tvar%found = .true.
            CASE ( 'timestep' )    ; tvar%found = .true.
            CASE DEFAULT ; CALL STOPIT("ERROR: external variable "//TRIM(tvar%name)//" not found.")
         END SELECT
      ELSEIF ( tvar%var_type == V_DIAGNOSTIC ) THEN  !# Diagnostic variable
         IF ( tvar%sheet ) THEN
            sd = sd + 1
         ELSE
            d = d + 1
         ENDIF
      ELSEIF ( tvar%var_type == V_STATE ) THEN    !# state variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
         ELSE
            v = v + 1
         ENDIF
      ENDIF
      IF ( .NOT. tvar%found ) THEN
         print *, "ERROR: Undefined variable ", TRIM(tvar%name)
         err_count = err_count + 1
      ENDIF
   ENDDO

   IF ( n_vars < v ) print *,"More vars than expected",v,n_vars
   IF ( n_vars_ben < sv ) print *,"More sheet vars than expected"
   IF ( n_vars_diag < d ) print *,"More diag vars than expected"
   IF ( n_vars_diag_sheet < sd ) print *,"More sheet diag vars than expected"

   IF ( err_count > 0 ) CALL STOPIT("*** Errors in configuration")
END SUBROUTINE check_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE define_column(column, top)
!-------------------------------------------------------------------------------
! Set up the current column pointers
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t), INTENT(inout) :: column(:)
   INTEGER, INTENT(in)  :: top
!
!LOCALS
   INTEGER :: av !, i
   INTEGER :: v, d, sv, sd, ev
   TYPE(aed_variable_t),POINTER :: tvar
!-------------------------------------------------------------------------------
!BEGIN
   v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
   DO av=1,n_aed_vars
      IF ( .NOT.  aed_get_var(av, tvar) ) STOP "Error getting variable info"

      IF ( tvar%var_type == V_EXTERNAL ) THEN !# global variable
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'temperature' ) ; column(av)%cell => temp(:)
            CASE ( 'salinity' )    ; column(av)%cell => salt(:)
            CASE ( 'density' )     ; column(av)%cell => rho(:)
            CASE ( 'layer_ht' )    ; column(av)%cell => dz(:)
            CASE ( 'extc_coef' )   ; column(av)%cell => extc(:)
            CASE ( 'tss' )         ; column(av)%cell => tss(:)
            CASE ( 'cell_vel' )    ; column(av)%cell => vel(:)
            CASE ( 'par' )         ; column(av)%cell => par(:)
            CASE ( 'nir' )         ; column(av)%cell => nir(:)
            CASE ( 'uva' )         ; column(av)%cell => uva(:)
            CASE ( 'uvb' )         ; column(av)%cell => uvb(:)
            CASE ( 'pressure' )    ; column(av)%cell => pres(:)
            CASE ( 'depth' )       ; column(av)%cell => depth(:)
            CASE ( 'sed_zones' )   ; column(av)%cell => sed_zones(:)
            CASE ( 'sed_zone' )    ; column(av)%cell_sheet => sed_zones(1)
            CASE ( 'wind_speed' )  ; column(av)%cell_sheet => wnd
            CASE ( 'delzBlueIce' ) ; column(av)%cell_sheet => delzBlueIce
            CASE ( 'delzWhiteIce' ) ; column(av)%cell_sheet => delzWhiteIce
            CASE ( 'u_star' )      ; column(av)%cell_sheet => u_star_glm
            CASE ( 'Q_net' )       ; column(av)%cell_sheet => Q_net_glm
            CASE ( 'par_sf' )      ; column(av)%cell_sheet => I_0
            CASE ( 'taub' )        ; column(av)%cell_sheet => bottom_stress
            CASE ( 'col_depth' )   ; column(av)%cell_sheet => col_depth
            CASE ( 'col_num' )     ; column(av)%cell_sheet => col_num
            CASE ( 'col_area' )    ; column(av)%cell_sheet => area(top)
            CASE ( 'evap' )        ; column(av)%cell_sheet => evap
            CASE ( 'layer_area' )  ; column(av)%cell => layer_area(:)
            CASE ( 'rain' )        ; column(av)%cell_sheet => precip
            CASE ( 'air_temp' )    ; column(av)%cell_sheet => air_temp
            CASE ( 'air_pres' )    ; column(av)%cell_sheet => air_pres
            CASE ( 'humidity' )    ; column(av)%cell_sheet => rel_hum
            CASE ( 'longitude' )   ; column(av)%cell_sheet => longitude
            CASE ( 'latitude' )    ; column(av)%cell_sheet => latitude
            CASE ( 'yearday' )     ; column(av)%cell_sheet => yearday
            CASE ( 'timestep' )    ; column(av)%cell_sheet => timestep
            CASE DEFAULT ; CALL STOPIT("ERROR: external variable "//TRIM(tvar%name)//" not found.")
         END SELECT
      ELSEIF ( tvar%var_type == V_DIAGNOSTIC ) THEN  !# Diagnostic variable
         IF ( tvar%sheet ) THEN
            sd = sd + 1
            column(av)%cell_sheet => cc_diag_hz(sd)
         ELSE
            d = d + 1
            column(av)%cell => cc_diag(d, :)
         ENDIF
      ELSEIF ( tvar%var_type == V_STATE ) THEN    !# state variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
            IF ( tvar%bot ) THEN
               column(av)%cell_sheet => cc(n_vars+sv, 1)
            ELSEIF ( tvar%top ) THEN
               column(av)%cell_sheet => cc(n_vars+sv, top)
            ENDIF
            column(av)%flux_ben => flux_ben(n_vars+sv)
            column(av)%flux_atm => flux_atm(n_vars+sv)
         ELSE
            v = v + 1
            column(av)%cell => cc(v, :)
            column(av)%flux_atm => flux_atm(v)
            column(av)%flux_pel => flux_pel(v, :)
            column(av)%flux_ben => flux_ben(v)
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE define_column
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE check_states(wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: wlev
!
!LOCALS
   TYPE(aed_variable_t),POINTER :: tv
   INTEGER :: i,v,sv,lev
#if DEBUG
   INTEGER :: last_naned = -1
#endif
!
!-------------------------------------------------------------------------------
!BEGIN
   IF ( .NOT. repair_state ) RETURN

   v = 0 ; sv = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tv) ) THEN
         IF ( tv%var_type == V_STATE ) THEN
            IF (tv%sheet) THEN
               sv = sv + 1
#if DEBUG
               IF ( ieee_is_nan(cc(n_vars+sv, 1)) ) last_naned = i
#endif
               IF ( .NOT. ieee_is_nan(tv%minimum) ) THEN
                  IF ( cc(n_vars+sv, 1) < tv%minimum ) cc(n_vars+sv, 1) = tv%minimum
               ENDIF
               IF ( .NOT. ieee_is_nan(tv%maximum) ) THEN
                  IF ( cc(n_vars+sv, 1) > tv%maximum ) cc(n_vars+sv, 1) = tv%maximum
               ENDIF
            ELSE
               v = v + 1
               DO lev=1, wlev
#if DEBUG
                  IF ( ieee_is_nan(cc(v, lev)) ) last_naned = i
#endif
                  IF ( .NOT. ieee_is_nan(tv%minimum) ) THEN
                     IF ( cc(v, lev) < tv%minimum ) cc(v, lev) = tv%minimum
                  ENDIF
                  IF ( .NOT. ieee_is_nan(tv%maximum) ) THEN
                     IF ( cc(v, lev) > tv%maximum ) cc(v, lev) = tv%maximum
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDIF
   ENDDO

#if DEBUG
   IF ( last_naned > -1 ) THEN
      print*
      IF ( aed_get_var(last_naned, tv) ) THEN
         print*,"NaNs detected in CC in var ", TRIM(tv%name)
      ELSE
         print*,"NaNs detected in CC unidentified var"
      ENDIF
   !  STOP
   ENDIF
#endif
END SUBROUTINE check_states
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_do_glm(wlev)                             BIND(C, name=_WQ_DO_GLM)
!-------------------------------------------------------------------------------
!                           wlev is the number of levels used;
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: wlev
!
!LOCALS
   TYPE(aed_variable_t),POINTER :: tv

   AED_REAL :: min_C, surf, ratio, localext_post
   INTEGER  :: i, j, v, lev, split, k, j_d
   INTEGER  :: d_rz, sd_rz     ! rezero counters (independent of j_d / other index vars)
   INTEGER, SAVE :: wlev_prev = 0
   LOGICAL  :: grid_changed

   TYPE (aed_column_t) :: column(n_aed_vars)
   TYPE (aed_column_t) :: column_sed(n_aed_vars)
   AED_REAL :: pa = 0.
!
!-------------------------------------------------------------------------------
!BEGIN
   col_depth = height(wlev)
   surf = height(wlev)
   !# re-compute the layer heights and depths
   dz(1) = height(1)
   depth(1) = surf - height(1)
   layer_area(1) = area(1)
   DO i=2,wlev
      dz(i) = height(i) - height(i-1)
      depth(i) = surf - height(i)
      layer_area(i) = area(i) - area(i-1)   ! incremental area: cumulative area gained at this layer
   ENDDO

   !# Calculate local pressure
   pres(1:wlev) = -height(1:wlev)

   IF ( benthic_mode > 1 ) THEN
      j = 1
      DO i=1,wlev
!        print *,'j',i,j
!        print *,'zone_heights',height(i),zone_heights(j),theZones(1)%zheight,theZones(2)%zheight
         IF (height(i) .GT. zone_heights(j)) THEN
            sed_zones(i) = j * ( area(i) - pa ) / area(i)
            pa = area(i)
            j = j+1
         ELSE
            sed_zones(i) = j
         ENDIF
      ENDDO
   ENDIF

   CALL define_column(column, wlev)

   IF ( .NOT. reinited ) CALL re_initialize

   ! Per-variable rezero — matches API path (aed_api.F90 aed_run_model, lines 1385-1397).
   ! Each diagnostic is registered with a `rezero` flag (default .TRUE.).  Vars flagged
   ! rezero=.FALSE. (e.g. CH4_ch4_ebb_dsfv, aed_methane.F90:387) must persist across
   ! timesteps; every other diagnostic should be cleared here so its value reflects only
   ! the current timestep.  The previous "!cc_diag = 0." blanket-off protected the zavg
   ! ones but left stale values in per-timestep 3D diagnostics such as BUB_Cbg_Wch4_*,
   ! BUB_gas_exchange_Wch4, BUB_bDiff_Wch4 etc.  This restores the correct behaviour.
   d_rz = 0; sd_rz = 0
   DO i = 1, n_aed_vars
      IF ( aed_get_var(i, tv) ) THEN
         IF ( tv%var_type == V_DIAGNOSTIC ) THEN
            IF ( tv%sheet ) THEN
               sd_rz = sd_rz + 1
               IF ( tv%rezero ) cc_diag_hz(sd_rz) = zero_
            ELSE
               d_rz = d_rz + 1
               IF ( tv%rezero ) cc_diag(d_rz, :) = zero_
            ENDIF
         ENDIF
      ENDIF
   ENDDO

   ! Save snapshot then remap zavg diagnostics onto the new layer grid.
   cc_diag_old = cc_diag
   ! ---- ORIGINAL (growth-only remap; revert by re-enabling this block) ----
   ! Only remapped when layers were ADDED (wlev > wlev_prev); left persisted zavg
   ! diagnostics (e.g. ch4_ebb_dsfv) on a stale grid during drawdown (wlev shrinking)
   ! and same-count reshapes. Replaced by the general remap below.
   !IF (wlev > wlev_prev .AND. wlev_prev > 0) THEN
   !   j_d = 0
   !   DO k = 1, n_aed_vars
   !      IF ( aed_get_var(k, tv) ) THEN
   !         IF ( tv%diag .AND. .NOT. tv%sheet ) THEN
   !            j_d = j_d + 1
   !            IF ( tv%zavg ) THEN
   !               j = 1
   !               cc_diag(j_d, 1) = cc_diag_old(j_d, 1)
   !               DO i = 2, wlev - 1
   !                  IF (height(i) < z_prev(1)) THEN
   !                     cc_diag(j_d, i) = cc_diag_old(j_d, 1)
   !                  ELSE
   !                     DO
   !                        IF (height(i) >= z_prev(j) .AND. height(i) <= z_prev(j+1)) THEN
   !                           ratio = (height(i) - z_prev(j)) / (z_prev(j+1) - z_prev(j))
   !                           cc_diag(j_d, i) = (1.0-ratio)*cc_diag_old(j_d, j) + ratio*cc_diag_old(j_d, j+1)
   !                           EXIT
   !                        ELSE
   !                           IF (j < wlev_prev - 1) THEN
   !                              j = j + 1
   !                           ELSE
   !                              EXIT
   !                           ENDIF
   !                        ENDIF
   !                     ENDDO
   !                  ENDIF
   !               ENDDO
   !               cc_diag(j_d, wlev) = cc_diag_old(j_d, wlev_prev)
   !            ENDIF
   !         ENDIF
   !      ENDIF
   !   ENDDO
   !ENDIF
   ! ---- GENERAL remap: grow, SHRINK (drawdown), or same-count reshape (mixing). ----
   ! Only remap CH4_ch4_ebb_dsfv (bubble flux diagnostic) which needs persistence across grid changes.
   ! Other zavg diagnostics should not be interpolated during grid changes as they may represent
   ! intensive properties that shouldn't be linearly remapped (revert to old selective approach).
   grid_changed = ( wlev /= wlev_prev )
   IF ( .NOT. grid_changed .AND. wlev_prev > 0 ) grid_changed = ANY( lheights(1:wlev) /= z_prev(1:wlev) )
   IF ( grid_changed .AND. wlev_prev > 0 ) THEN
      j_d = 0
      DO k = 1, n_aed_vars
         IF ( aed_get_var(k, tv) ) THEN
            IF ( tv%var_type == V_DIAGNOSTIC .AND. .NOT. tv%sheet ) THEN
               j_d = j_d + 1
               IF ( tv%zavg ) THEN
                  DO lev = 1, wlev
                     IF     (lheights(lev) <= z_prev(1))         THEN
                        cc_diag(j_d, lev) = cc_diag_old(j_d, 1)          ! at/below old bottom
                     ELSEIF (lheights(lev) >= z_prev(wlev_prev)) THEN
                        cc_diag(j_d, lev) = cc_diag_old(j_d, wlev_prev)  ! at/above old top (covers growth)
                     ELSE
                        j = 1
                        DO
                           IF (lheights(lev) >= z_prev(j) .AND. lheights(lev) <= z_prev(j+1)) THEN
                              ratio = (lheights(lev) - z_prev(j)) / (z_prev(j+1) - z_prev(j))
                              cc_diag(j_d, lev) = (1.0-ratio)*cc_diag_old(j_d, j) + ratio*cc_diag_old(j_d, j+1)
                              EXIT
                           ELSEIF (j < wlev_prev - 1) THEN
                              j = j + 1
                           ELSE
                              cc_diag(j_d, lev) = cc_diag_old(j_d, wlev_prev)
                              EXIT
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ENDIF

   IF ( .NOT. mobility_off ) THEN
      v = 0
      DO i=1,n_aed_vars
         IF ( aed_get_var(i, tv) ) THEN
            IF ( .NOT. tv%sheet .AND. tv%var_type == V_STATE ) THEN
               v = v + 1
               ! only for state_vars that are not sheet
               IF ( .NOT. ieee_is_nan(tv%mobility) ) THEN
                  ! default to ws that was set during initialisation
                  ws(i,:) = tv%mobility
               ELSE
                  ws(i,:) = zero_
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      DO lev = 1, wlev
         ! update ws for modules that use the mobility method
         CALL aed_mobility(column, lev, ws(:,lev))
      ENDDO

      !# (3) Calculate source/sink terms due to the settling or rising of
      !# state variables in the water column (note that settling into benthos
      !# is done in aed_do_benthos)
      v = 0
      DO i=1,n_aed_vars
         IF ( aed_get_var(i, tv) ) THEN
            IF ( .NOT. tv%sheet .AND. tv%var_type == V_STATE ) THEN
               v = v + 1
               !# only for state_vars that are not sheet, and also non-zero ws
               IF ( .NOT. ieee_is_nan(tv%mobility) .AND. SUM(ABS(ws(i, 1:wlev)))>zero_ ) THEN
                  min_C = tv%minimum
                  CALL doMobility(wlev, dt, dz, area, ws(i, :), min_C, cc(v, :))
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ENDIF

   DO lev=1, wlev
      CALL aed_equilibrate(column, lev)
   ENDDO

   CALL check_states(wlev)

   IF (benthic_mode .GT. 1) &
      CALL calc_zone_areas(area, wlev, height(wlev))

   DO split=1,split_factor
      IF (benthic_mode .GT. 1) &
         CALL copy_to_zone(cc, cc_diag, wlev)

      !# Update local light field (self-shading may have changed through
      !# changes in biological state variables). Update_light is set to
      !# be inline with current aed_phyoplankton, which requires only
      !# surface par, then integrates over depth of a layer
      CALL update_light(column, wlev)

      !# Fudge
      nir(:) = (par(:)/par_fraction) * nir_fraction
      uva(:) = (par(:)/par_fraction) * uva_fraction
      uvb(:) = (par(:)/par_fraction) * uvb_fraction

      !# Time-integrate one biological time step
      CALL calculate_fluxes(wlev)

      !# Update the water column layers
      cc(1:n_vars, 1:wlev) = cc(1:n_vars, 1:wlev) + dt_eff*flux_pel(1:n_vars, 1:wlev)

      !# Now update benthic variables, depending on whether zones are simulated
      IF ( benthic_mode .GT. 1 ) THEN
         ! Loop through each sediment zone
         z_cc(n_vars+1:n_vars+n_vars_ben, 1, 1:n_zones) = z_cc(n_vars+1:n_vars+n_vars_ben, 1, 1:n_zones) &
                                                + dt_eff*flux_zon(n_vars+1:n_vars+n_vars_ben, 1:n_zones)

         !# Distribute cc-sed benthic properties back into main cc array
         CALL copy_from_zone(n_aed_vars, cc, cc_diag, cc_diag_hz, wlev)

      ELSE
         cc(n_vars+1:n_vars+n_vars_ben, 1) = cc(n_vars+1:n_vars+n_vars_ben, 1) &
                                 + dt_eff*flux_ben(n_vars+1:n_vars+n_vars_ben)
      ENDIF

     !CALL post_integration_update(column,wlev,flux_pel,flux_ben,flux_atm,flux_zon)

!     DO lev=1, wlev
!        CALL aed_equilibrate(column, lev)
!     ENDDO
      !     aed_update  (aka aed_equilibrate)
      !       layer_map = (/ (i, i = wlev, 1, -1) /)
      !     aed_update_column (column,layer_map,)

      CALL check_states(wlev)
   ENDDO

   IF (.NOT. link_ext_par) THEN
      localext_post = zero_
      DO i = 1, wlev - 1
         CALL aed_light_extinction(column, i, localext_post)
         extc(i) = Kw + localext_post
      ENDDO
   ENDIF

  ! IF ( display_minmax ) THEN
  !    v = 0; d = 0
  !    DO i=1,n_aed_vars
  !       IF ( aed_get_var(i, tv) ) THEN
  !          IF ( tv%var_type == V_STATE ) THEN
  !             v = v + 1
  !             WRITE(*,'(1X,"VarLims: ",I0,1X,"<=> ",f15.8,f15.8," : ",A," (",A,")")') &
  !                                        v,MINVAL(cc(v,:)),MAXVAL(cc(v,:)),TRIM(tv%name),TRIM(tv%units)
  !             !print *,'VarLims',v,TRIM(tv%name),MINVAL(cc(v,:)),MAXVAL(cc(v,:))
  !          ELSE IF ( tv%var_type == V_DIAGNOSTIC  ) THEN
  !             d = d + 1
  !             WRITE(*,'(1X,"DiagLim: ",I0,1X,"<=> ",f15.8,f15.8," : ",A," (",A,")")') &
  !                                        d,MINVAL(cc_diag(:,d)),MAXVAL(cc_diag(:,d)),TRIM(tv%name),TRIM(tv%units)
  !             !print *,'DiagLim',d,TRIM(tv%name),MINVAL(cc_diag(d,:)),MAXVAL(cc_diag(d,:))
  !          ENDIF
  !       ENDIF
  !    ENDDO
  !  ENDIF

   ! Save layer count and heights for interpolation on next timestep
   wlev_prev = wlev
   z_prev = 0.
   z_prev(1:wlev) = height(1:wlev)

CONTAINS
!-------------------------------------------------------------------------------

   !############################################################################
   SUBROUTINE define_zone_column(zon, top)
   !----------------------------------------------------------------------------
   ! Set up the current column pointers
   !----------------------------------------------------------------------------
   !ARGUMENTS
      INTEGER, INTENT(in)  :: zon, top
   !
   !LOCALS
      INTEGER :: av
      INTEGER :: v, d, sv, sd, ev
      TYPE(aed_variable_t),POINTER :: tvar
   !----------------------------------------------------------------------------
   !BEGIN
   v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
   DO av=1,n_aed_vars
      IF ( .NOT.  aed_get_var(av, tvar) ) STOP "Error getting variable info"

      IF ( tvar%var_type == V_EXTERNAL ) THEN !# global "environment" variables
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'temperature' ) ; column_sed(av)%cell => theZones(:)%ztemp
            CASE ( 'salinity' )    ; column_sed(av)%cell => theZones(:)%zsalt
            CASE ( 'density' )     ; column_sed(av)%cell => theZones(:)%zrho
            CASE ( 'layer_ht' )    ; column_sed(av)%cell => theZones(:)%zdz
            CASE ( 'extc_coef' )   ; column_sed(av)%cell => theZones(:)%zextc
            CASE ( 'tss' )         ; column_sed(av)%cell => theZones(:)%ztss
            CASE ( 'cell_vel' )    ; column_sed(av)%cell => theZones(:)%zvel
            CASE ( 'par' )         ; column_sed(av)%cell => theZones(:)%zpar
            CASE ( 'nir' )         ; column_sed(av)%cell => theZones(:)%znir
            CASE ( 'uva' )         ; column_sed(av)%cell => theZones(:)%zuva
            CASE ( 'uvb' )         ; column_sed(av)%cell => theZones(:)%zuvb
            CASE ( 'pressure' )    ; column_sed(av)%cell => theZones(:)%zpres
            CASE ( 'depth' )       ; column_sed(av)%cell => theZones(:)%zdepth
            CASE ( 'sed_zones' )   ; column_sed(av)%cell => theZones(:)%z_sed_zones
            CASE ( 'sed_zone' )    ; column_sed(av)%cell_sheet => theZones(zon)%z_sed_zones; zone_var = av
            CASE ( 'wind_speed' )  ; column_sed(av)%cell_sheet => wnd
            CASE ( 'delzBlueIce' ) ; column_sed(av)%cell_sheet => delzBlueIce
            CASE ( 'delzWhiteIce' ) ; column_sed(av)%cell_sheet => delzWhiteIce
            CASE ( 'u_star' )      ; column_sed(av)%cell_sheet => u_star_glm
            CASE ( 'Q_net' )       ; column_sed(av)%cell_sheet => Q_net_glm
            CASE ( 'par_sf' )      ; column_sed(av)%cell_sheet => I_0
            CASE ( 'taub' )        ; column_sed(av)%cell_sheet => bottom_stress
            CASE ( 'col_depth' )   ; column_sed(av)%cell_sheet => col_depth
            CASE ( 'col_num' )     ; column_sed(av)%cell_sheet => col_num
            CASE ( 'col_area' )    ; column_sed(av)%cell_sheet => area(top)
            CASE ( 'evap' )        ; column_sed(av)%cell_sheet => evap
            CASE ( 'layer_area' )  ; column_sed(av)%cell => theZones(:)%zarea
            CASE ( 'rain' )        ; column_sed(av)%cell_sheet => precip
            CASE ( 'air_temp' )    ; column_sed(av)%cell_sheet => air_temp
            CASE ( 'air_pres' )    ; column_sed(av)%cell_sheet => air_pres
            CASE ( 'humidity' )    ; column_sed(av)%cell_sheet => rel_hum
            CASE ( 'longitude' )   ; column_sed(av)%cell_sheet => longitude
            CASE ( 'latitude' )    ; column_sed(av)%cell_sheet => latitude
            CASE ( 'yearday' )     ; column_sed(av)%cell_sheet => yearday
            CASE ( 'timestep' )    ; column_sed(av)%cell_sheet => timestep
            CASE DEFAULT ; CALL STOPIT("ERROR: external variable "//trim(tvar%name)//" not found.")
         END SELECT
      ELSEIF ( tvar%var_type == V_DIAGNOSTIC ) THEN  !# Diagnostic variables
         IF ( tvar%sheet ) THEN
            sd = sd + 1
            column_sed(av)%cell_sheet => z_diag_hz(sd, zon)
         ELSE
            d = d + 1
            column_sed(av)%cell => z_diag(d, :, zon)     ! limit to layers above zone (instead of ":")
         ENDIF
      ELSEIF ( tvar%var_type == V_STATE ) THEN    !# State variables
         IF ( tvar%sheet ) THEN
            sv = sv + 1
            IF ( tvar%bot ) THEN
               column_sed(av)%cell_sheet => z_cc(n_vars+sv, 1, zon) ! 2D are stored after n_vars in cc
            ELSEIF ( tvar%top ) THEN
               column_sed(av)%cell_sheet => z_cc(n_vars+sv, 1, zon)
            ENDIF
            column_sed(av)%flux_ben => flux_ben(n_vars+sv)
            column_sed(av)%flux_atm => flux_atm(n_vars+sv)
         ELSE
            v = v + 1
            column_sed(av)%cell => z_cc(v, 1, :)  ! GLM water layers above "zon" (set in copy_to_zone)
            column_sed(av)%flux_atm => flux_atm(v)
            column_sed(av)%flux_pel => flux_pel(v, :)
            column_sed(av)%flux_ben => flux_ben(v)
         ENDIF
      ENDIF
   ENDDO
   END SUBROUTINE define_zone_column
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !############################################################################
   SUBROUTINE re_initialize
   !----------------------------------------------------------------------------
   !ARGUMENTS
   !LOCALS
      INTEGER :: lev, zon, av, sv, sd
      INTEGER :: layer_map(wlev)
      TYPE(aed_variable_t),POINTER :: tvar
   !
   !----------------------------------------------------------------------------
   !BEGIN
      !# (1) WATER COLUMN INITIALISATION
      layer_map = (/ (i, i = wlev, 1, -1) /)
      CALL aed_initialize_column(column, layer_map(:))

      !# (2) WATER CELL INITIALISATION
      DO lev=1, wlev
         CALL aed_initialize(column, lev)
      ENDDO

      !# (3) BENTHIC INITIALISATION
      IF ( benthic_mode .GT. 1 ) THEN
         !# Multiple static sediment zones are simulated, and therfore overlying
         !# water conditions need to be aggregated from multiple cells/layers

         DO zon=1,n_zones
            CALL define_zone_column(zon, n_zones)

            theZones(zon)%z_sed_zones = zon  !MH TMP
!           print *,'theZones(zon)%z_sed_zones',theZones(zon)%z_sed_zones !MH TMP
            !# If multiple benthic zones, we must update the benthic variable pointer for the new zone
            IF (zone_var .GT. 0) column_sed(zone_var)%cell_sheet => theZones(zon)%z_sed_zones

            sv = 0 ; sd = 0

            DO av=1,n_aed_vars
               IF ( .NOT. aed_get_var(av, tvar) ) STOP "Error getting variable info"

               IF ( tvar%var_type == V_DIAGNOSTIC ) THEN  !# Diagnostic variable
                  IF ( tvar%sheet ) THEN
                     sd = sd + 1
                     column_sed(av)%cell_sheet => z_diag_hz(sd, zon)
                  ENDIF
               ELSEIF ( tvar%var_type == V_STATE ) THEN !# State variable
                  IF ( tvar%sheet ) THEN
                     sv = sv + 1
                     IF ( tvar%bot ) THEN
                        column_sed(av)%cell_sheet => z_cc(n_vars+sv, 1, zon)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

            !# Copy benthic state initial values from cc to this zone's z_cc
            z_cc(n_vars+1:n_vars+n_vars_ben, 1, zon) = cc(n_vars+1:n_vars+n_vars_ben, 1)

            CALL aed_initialize_benthic(column_sed, zon)

           !# (4) ZONE COLUMN INITIALISATION
           ! CALL aed_initialize_zone_column(column_sed, zon)
         ENDDO
      ENDIF
      reinited = .TRUE.
   END SUBROUTINE re_initialize
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !############################################################################
   SUBROUTINE calculate_fluxes(wlev)
   !----------------------------------------------------------------------------
   ! Checks the current values of all state variables and repairs these
   !----------------------------------------------------------------------------
   !ARGUMENTS
      INTEGER, INTENT(in) :: wlev
   !
   !LOCALS
      INTEGER :: lev,zon,v,av,sv,sd, zlev
      INTEGER :: layer_map(wlev)
      AED_REAL :: scale
      AED_REAL :: localrainl(n_zones), localshade(n_zones), localdrag(n_zones)
      LOGICAL :: splitZone
      TYPE(aed_variable_t),POINTER :: tvar
   !----------------------------------------------------------------------------
   !BEGIN
      flux_pel = zero_
      flux_atm = zero_
      flux_ben = zero_
      IF (ALLOCATED(flux_zon)) flux_zon = zero_

      flux_pel_pre = zero_
      flux_pel_z = zero_

      !# Start with updating column (used for light, bubbles, plants, phreeqc)

      !# WATER COLUMN UPDATING
      DO lev=1,wlev
         layer_map(lev) = 1 + wlev-lev
      ENDDO
      CALL aed_calculate_column(column, layer_map)

      !# Now do the general calculation all flux terms for RHS in mass/m3/s
      !# Includes (i) benthic flux, (ii) surface exchange and (ii) kinetic updates in each cell
      !# as calculated by glm

      !# BENTHIC FLUXES
      IF ( benthic_mode .GT. 1 ) THEN
         !# Multiple static sediment zones are simulated, and therfore overlying
         !# water conditions need to be aggregated from multiple cells/layers, and output flux
         !# needs disaggregating from each zone back to the overlying cells/layers

!$OMP DO
         DO zon=1,n_zones
            CALL define_zone_column(zon, n_zones)

            !# Reinitialise flux_ben to be repopulated for this zone
            flux_ben = zero_
            flux_pel_pre = zero_

            !# If multiple benthic zones, we must update the benthic variable pointer for the new zone
            IF (zone_var > 0) column_sed(zone_var)%cell_sheet => theZones(zon)%z_sed_zones
!           !MH WE NEED A COLUMN TO CC VAR MAP FOR BENTHIC GUYS
!           !CAB Yes, a map (or 2 maps) would be better, but QnD since this all needs reworking

            sv = 0 ; sd = 0
            DO av=1,n_aed_vars
               IF ( .NOT. aed_get_var(av, tvar) ) STOP "Error getting variable info"

               IF ( tvar%sheet ) THEN
                  IF ( tvar%var_type == V_DIAGNOSTIC ) THEN  !# Diagnostic variable
                     sd = sd + 1
                     column_sed(av)%cell_sheet => z_diag_hz(sd, zon)
                  ELSEIF ( tvar%var_type == V_STATE ) THEN !# State variable
                     sv = sv + 1
                     IF ( tvar%bot ) THEN
                        column_sed(av)%cell_sheet => z_cc(n_vars+sv, 1, zon)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
!           print*,"Calling ben for zone ",zon,zone_var,z_sed_zones(zon)

            !# (1) ZONE COLUMN UPDATING
            zlev = 0
            DO lev=1,wlev
               IF(zone_heights(zon)<lheights(lev)) THEN
                  zlev = lev
                  EXIT
               ENDIF
            ENDDO
            ! The upper height of the last zone may have no water above it.
            ! In this case, set the wet-layer vector to just be the top water layer
            IF (zlev == 0) zlev = wlev

            DO lev=zlev,wlev
              layer_map(lev) = zlev + wlev-lev
            ENDDO
            IF ( benthic_mode .EQ. 3 ) &
               CALL aed_calculate_column(column_sed, layer_map)   ! only for riparian; bubble module needs full-wlev column, not zone column

            IF ( benthic_mode .EQ. 3 ) THEN
               !# Zone is able to operated on by riparian and dry methods
               CALL aed_calculate_riparian(column_sed, zlev, theZones(zon)%z_pc_wet)
               IF (theZones(zon)%z_pc_wet < 0.01 ) CALL aed_calculate_dry(column_sed, zlev)

               !# update feedback arrays to host model, to reduce rain (or if -ve then add flow)
               CALL aed_rain_loss(column_sed, 1, localrainl(zon));
               IF (link_rain_loss) rain_factor = localrainl(zon)

               !# update feedback arrays to shade the water (ie reduce incoming light, Io)
               CALL aed_light_shading(column_sed, 1, localshade(zon))
               IF (link_solar_shade) sw_factor = localshade(zon)

               !# update feedback arrays to host model, to add drag effects
               IF (theZones(zon)%z_pc_wet > 0.99 ) THEN
                  CALL aed_bio_drag(column_sed, zlev, localdrag(zon))
                  IF (link_bottom_drag) friction = localdrag(zon)
               ENDIF
            ENDIF

            !# Calculate temporal derivatives due to benthic processes.
            !# They are stored in flux_ben (benthic vars) and flux_pel (water vars)
            flux_pel_pre = flux_pel

!           print*,"Calling ben for zone ",zon,zone_var!,z_sed_zones(zon)
            CALL aed_calculate_benthic(column_sed, zon, .TRUE.)  ! use zon: cell => theZones(:) so cell(zon) = this zone's env vars

            !# Record benthic fluxes in the zone array
            flux_zon(:, zon) = flux_ben(:)

            !# Now we have to find out the water column flux that occured and
            !# (incremented in flux_pel) and then disaggregate it to relevant layers
            flux_pel_z(:, zon) = flux_pel(:, zon)-flux_pel_pre(:, zon)
         ENDDO
!$OMP END DO

         !# Disaggregation of zone induced fluxes to overlying layers
         zon = n_zones
         DO lev=wlev,1,-1
            IF ( zon .GT. 1 ) THEN
               IF (lev .GT. 1) THEN
                  splitZone = lheights(lev-1) < zone_heights(zon-1)
               ELSE
                  splitZone = 0.0 < zone_heights(zon-1)
               ENDIF
            ELSE
               splitZone = .FALSE.
            ENDIF

            IF (splitZone) THEN
               IF (lev .GT. 1) THEN
                  scale = (zone_heights(zon-1) - lheights(lev-1)) / (lheights(lev) - lheights(lev-1))
               ELSE
                  scale = (zone_heights(zon-1) - 0.0) / (lheights(lev) - 0.0)
               ENDIF
               flux_pel(1:n_vars,lev) = flux_pel_z(1:n_vars,zon) * scale

               zon = zon - 1

               flux_pel(1:n_vars,lev) = flux_pel(1:n_vars,lev) + &
                                        flux_pel_z(1:n_vars,zon) * (1.0 - scale)
            ELSE
               flux_pel(1:n_vars,lev) = flux_pel_z(1:n_vars,zon)
            ENDIF
         ENDDO
         !# Limit flux out of bottom waters to concentration of that layer
         !# i.e. don't flux out more than is there & distribute
         !# bottom flux into pelagic over bottom box (i.e., divide by layer height).
         !# scaled to proportion of area that is "bottom"
         DO lev=1,wlev
            IF (lev > 1) flux_pel(:, lev) = flux_pel(:, lev) * (area(lev)-area(lev-1))/area(lev)
            DO v=1,n_vars
              IF ( cc(v, 1) .GE. 0.0 ) flux_pel(v, lev) = max(-1.0 * cc(v, lev), flux_pel(v, lev)/dz(lev))
            END DO
         ENDDO
      ELSE
         !# Sediment zones are not simulated and therefore just operate on the bottom-most
         !# GLM layer as the "benthos". If benthic_mode=1 then benthic fluxes will also be
         !# applied on flanks of the remaining layers, but note this is not suitable for
         !# model configurations where mass balance of benthic variables is required.

         !# Calculate temporal derivatives due to exchanges at the sediment/water interface
         CALL aed_calculate_benthic(column, 1)

         !# Limit flux out of bottom layers to concentration of that layer
         !# i.e. don't flux out more than is there is. Then
         !# distribute bottom flux into pelagic over bottom box (i.e., divide by layer height)
         !# Skip -ve values, as GEO_ubalchg is -ve and doesnt not comply with this logic
         DO v=1,n_vars
           IF ( cc(v, 1) .GE. 0.0 ) flux_pel(v, 1) = max(-1.0 * cc(v, 1), flux_pel(v, 1)/dz(1))
         END DO

         IF ( benthic_mode .EQ. 1 ) THEN
!$OMP DO
            DO lev=2,wlev
               !# Calculate temporal derivatives due to benthic fluxes.
               CALL aed_calculate_benthic(column, lev)

               !# Limit flux out of bottom layers to concentration of that layer
               !# i.e. don't flux out more than is there
               !# & distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
               !# scaled to proportion of area that is "bottom"
               DO v=1,n_vars
                  IF ( cc(v, lev) .GE. 0.0 ) flux_pel(v, lev) = &
                                   max(-1.0 * cc(v, lev), flux_pel(v, lev)/dz(lev))
               END DO
               flux_pel(:, lev) = flux_pel(:, lev) * (area(lev)-area(lev-1))/area(lev)
            ENDDO
!$OMP END DO
         ENDIF
      ENDIF

      !# (3) SURFACE FLUXES
      !# Calculate temporal derivatives due to air-water exchange.
      IF (.NOT. ice) THEN !# no surface exchange under ice cover
         CALL aed_calculate_surface(column, wlev)

         !# Distribute the fluxes into pelagic surface layer
         flux_pel(:, wlev) = flux_pel(:, wlev) + flux_atm(:)/dz(wlev)
      ENDIF

      !# (4) WATER CELL KINETICS
      !# Add pelagic sink and source terms in cells of all depth levels.
      DO lev=1,wlev
         CALL aed_calculate(column, lev)
      ENDDO

   END SUBROUTINE calculate_fluxes
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END SUBROUTINE aed_do_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_clean_glm() BIND(C, name=_WQ_CLEAN_GLM)
!-------------------------------------------------------------------------------
! Finish biogeochemical model
!-------------------------------------------------------------------------------
!BEGIN
   CALL aed_delete()
   ! Deallocate internal arrays
   IF (ALLOCATED(cc_diag))      DEALLOCATE(cc_diag)
   IF (ALLOCATED(cc_diag_old))  DEALLOCATE(cc_diag_old)
   IF (ALLOCATED(cc_diag_hz))   DEALLOCATE(cc_diag_hz)
   IF (ALLOCATED(z_prev))       DEALLOCATE(z_prev)

   IF (ALLOCATED(ws))           DEALLOCATE(ws)
   IF (ALLOCATED(par))          DEALLOCATE(par)
   IF (ALLOCATED(nir))          DEALLOCATE(nir)
   IF (ALLOCATED(uva))          DEALLOCATE(uva)
   IF (ALLOCATED(uvb))          DEALLOCATE(uvb)
   IF (ALLOCATED(pres))         DEALLOCATE(pres)
   IF (ALLOCATED(dz))           DEALLOCATE(dz)
   IF (ALLOCATED(layer_area))   DEALLOCATE(layer_area)

   IF (ALLOCATED(flux_pel))     DEALLOCATE(flux_pel)
   IF (ALLOCATED(flux_ben))     DEALLOCATE(flux_ben)
   IF (ALLOCATED(flux_atm))     DEALLOCATE(flux_atm)

   IF (ALLOCATED(flux_zon))     DEALLOCATE(flux_zon)
   IF (ALLOCATED(flux_pel_z))   DEALLOCATE(flux_pel_z)
   IF (ALLOCATED(flux_pel_pre)) DEALLOCATE(flux_pel_pre)
END SUBROUTINE aed_clean_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE update_light(column, nlev)
!-------------------------------------------------------------------------------
! Calculate photosynthetically active radiation over entire column based
! on surface radiation, attenuated based on background & biotic extinction
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t), INTENT(inout) :: column(:)
   INTEGER,INTENT(in)  :: nlev
!
!LOCALS
   INTEGER :: i
   AED_REAL :: localext, localext_up
!
!-------------------------------------------------------------------------------
!BEGIN
   localext = zero_; localext_up = zero_

   ! Surface Kd
   CALL aed_light_extinction(column, nlev, localext)

   ! Surface PAR — matches API path (aed_api.F90 Light) using zz = 0.001
   par(nlev) = par_fraction * rad(nlev) * EXP( -(Kw+localext) * 0.001 )

   ! Now set the top of subsequent layers, down to the bottom
   DO i = (nlev-1),1,-1
      localext_up = localext
      CALL aed_light_extinction(column, i, localext)

      par(i) = par(i+1) * EXP( -(Kw + localext_up) * dz(i+1) )

      IF (bioshade_feedback) extc(i) = Kw + localext
   ENDDO
END SUBROUTINE update_light
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_init_glm_output(ncid,x_dim,y_dim,z_dim,zone_dim,time_dim) BIND(C, name=_WQ_INIT_GLM_OUTPUT)
!-------------------------------------------------------------------------------
!  Initialize the output by defining biogeochemical variables.
!-------------------------------------------------------------------------------
!
!ARGUMENTS
   CINTEGER,INTENT(in) :: ncid,x_dim,y_dim,z_dim,zone_dim,time_dim
!
!LOCALS
   INTEGER i !, v, d
   INTEGER dims(4)

   TYPE(aed_variable_t),POINTER :: tv
!
!-------------------------------------------------------------------------------
!BEGIN
   !# Put NetCDF library in define mode.
   CALL define_mode_on(ncid)

   !# Set up dimension indices for 3D (+ time) variables (longitude,latitude,depth,time).
   dims(1) = x_dim
   dims(2) = y_dim
   dims(3) = z_dim
   dims(4) = time_dim

!  v = 0; d = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tv) ) THEN
         IF ( .NOT. tv%sheet .AND. (tv%var_type == V_DIAGNOSTIC .OR. tv%var_type == V_STATE) ) THEN
            !# only for state and diag vars that are not sheet
            externalid(i) = NEW_NC_VARIABLE(ncid, TRIM(tv%name), LEN_TRIM(tv%name), NF90_REALTYPE, 4, dims(1:4))
            CALL set_nc_attributes(ncid, externalid(i), MYTRIM(tv%units), MYTRIM(tv%longname) PARAM_FILLVALUE)
         ENDIF
      ENDIF
   ENDDO

   !# Set up dimension indices for 2D (+ time) variables (longitude,latitude,time).
   dims(1) = x_dim
   dims(2) = y_dim
   dims(3) = time_dim

!  v = 0; d = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tv) ) THEN
         IF ( tv%sheet .AND. (tv%var_type == V_DIAGNOSTIC .OR. tv%var_type == V_STATE) ) THEN
            !# only for state and diag sheet vars
            externalid(i) = NEW_NC_VARIABLE(ncid, TRIM(tv%name), LEN_TRIM(tv%name), NF90_REALTYPE, 3, dims(1:3))
            CALL set_nc_attributes(ncid, externalid(i), MYTRIM(tv%units), MYTRIM(tv%longname) PARAM_FILLVALUE)
         ENDIF
      ENDIF
   ENDDO

   !# Set up dimension indices for 2D+zone (+ time) variables (longitude,latitude,zone,time).
   IF ( n_zones .GT. 0 ) THEN
      dims(1) = x_dim
      dims(2) = y_dim
      dims(3) = zone_dim
      dims(4) = time_dim

!print*,x_dim,y_dim,zone_dim,time_dim
!     v = 0; d = 0
      DO i=1,n_aed_vars
         IF ( aed_get_var(i, tv) ) THEN
            IF ( tv%sheet .AND. (tv%var_type == V_STATE .OR. tv%var_type == V_DIAGNOSTIC) ) THEN
               !# only for state and diag sheet vars
               zexternalid(i) = NEW_NC_VARIABLE(ncid, TRIM(tv%name)//"_Z", LEN_TRIM(tv%name)+2, NF90_REALTYPE, 4, dims(1:4))
               CALL set_nc_attributes(ncid, zexternalid(i), MYTRIM(tv%units), MYTRIM(tv%longname) PARAM_FILLVALUE)
            ENDIF
         ENDIF
      ENDDO
   ENDIF

   !# Take NetCDF library out of define mode (ready for storing data).
   CALL define_mode_off(ncid)
END SUBROUTINE aed_init_glm_output
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!  Save properties of biogeochemical model, including state variable
!  values, diagnostic variable values, and sums of conserved quantities.
!-------------------------------------------------------------------------------
SUBROUTINE aed_write_glm(ncid,wlev,nlev,lvl,point_nlevs) BIND(C, name=_WQ_WRITE_GLM_)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: ncid, wlev, nlev
   CINTEGER,INTENT(in) :: lvl(*), point_nlevs
!
!LOCALS
   TYPE(aed_variable_t),POINTER :: tv

   INTEGER  :: i, j, v, d, sv, sd
#ifdef PLOTS
   INTEGER  :: z
#endif
   AED_REAL :: val_out
   FLOGICAL :: last = .FALSE.
!
!-------------------------------------------------------------------------------
!BEGIN
   v = 0; d = 0; sv = 0; sd = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tv) ) THEN
         IF ( tv%var_type == V_DIAGNOSTIC ) THEN
            !# Process and store diagnostic variables.
            IF ( tv%sheet ) THEN
               sd = sd + 1
               !# Process and store diagnostic variables defined on horizontal slices of the domain.
               IF ( n_zones .GT. 0 ) THEN
                  z_diag_hz(sd, n_zones+1) = cc_diag_hz(sd)
                  CALL store_nc_array(ncid, zexternalid(i), XYNT_SHAPE, n_zones, n_zones, array=z_diag_hz(sd, 1:n_zones+1))
               ENDIF
               CALL store_nc_scalar(ncid, externalid(i), XYT_SHAPE, scalar=cc_diag_hz(sd))
#ifdef PLOTS
               IF ( do_plots .AND. plot_id_sd(sd).GE.0 ) THEN
                  IF ( n_zones .GT. 0 ) THEN
                     DO z=1,n_zones
                        CALL put_glm_val_z(plot_id_sd(sd),z_diag_hz(sd, z), z)
                     ENDDO
                  ENDIF
                  CALL put_glm_val_s(plot_id_sd(sd),cc_diag_hz(sd))
               ENDIF
#endif
               DO j=1,point_nlevs
                  val_out = missing
                  IF ((lvl(j) .EQ. wlev) .AND. tv%top) val_out = cc_diag_hz(sd)
                  IF ((lvl(j) .EQ. 0)    .AND. tv%bot) val_out = cc_diag_hz(sd)
                  CALL write_csv_point(j, tv%name, len_trim(tv%name), val_out, NULCSTR, 0, last=last)
               ENDDO
            ELSE  !# not sheet
               d = d + 1
               !# Store diagnostic variable values defined on the full domain.
               CALL store_nc_array(ncid, externalid(i), XYZT_SHAPE, wlev, nlev, array=cc_diag(d, :))
#ifdef PLOTS
               IF ( do_plots .AND. plot_id_d(d).GE.0 ) &
                  CALL put_glm_val(plot_id_d(d), cc_diag(d, 1:wlev))
#endif
               DO j=1,point_nlevs
                  IF (lvl(j) .GE. 0) THEN ; val_out = cc_diag(d, lvl(j)+1)
                  ELSE                    ; val_out = missing     ; ENDIF
                  CALL write_csv_point(j, tv%name, len_trim(tv%name), val_out, NULCSTR, 0, last=last)
                  CALL write_csv_point_avg(j, tv%name, len_trim(tv%name), cc_diag(d, :), NULCSTR, 0, last=last)
               ENDDO
            ENDIF
         ELSEIF ( tv%var_type == V_STATE ) THEN  ! not diag
            IF ( tv%sheet ) THEN
               sv = sv + 1
               !# Store benthic biogeochemical state variables.
               IF ( n_zones .GT. 0 ) THEN
                  CALL store_nc_array(ncid, zexternalid(i), XYNT_SHAPE, n_zones, n_zones, array=z_cc(n_vars+sv, 1, 1:n_zones))
               ENDIF
               CALL store_nc_scalar(ncid, externalid(i), XYT_SHAPE, scalar=cc(n_vars+sv, 1))
#ifdef PLOTS
               IF ( do_plots .AND. plot_id_sv(sv).GE.0 ) THEN
                  IF ( n_zones .GT. 0 ) THEN
                     DO z=1,n_zones
                        CALL put_glm_val_z(plot_id_sv(sv), z_cc(n_vars+sv, 1, z), z)
                     ENDDO
                  ENDIF
                  CALL put_glm_val_s(plot_id_sv(sv), cc(n_vars+sv, 1))
               ENDIF
#endif
               DO j=1,point_nlevs
                  val_out = missing
                  IF ((lvl(j) .EQ. wlev) .AND. tv%top) val_out = cc(n_vars+sv, 1)
                  IF ((lvl(j) .EQ. 0)    .AND. tv%bot) val_out = cc(n_vars+sv, 1)
                  CALL write_csv_point(j, tv%name, len_trim(tv%name), val_out, NULCSTR, 0, last=last)
               ENDDO
            ELSE     !# not sheet
               v = v + 1
               !# Store pelagic biogeochemical state variables.
               CALL store_nc_array(ncid, externalid(i), XYZT_SHAPE, wlev, nlev, array=cc(v, :))
#ifdef PLOTS
               IF ( do_plots .AND. plot_id_v(v).GE.0 ) CALL put_glm_val(plot_id_v(v), cc(v, 1:wlev))
#endif
               DO j=1,point_nlevs
                  IF (lvl(j) .GE. 0) THEN ; val_out = cc(v, lvl(j)+1)
                  ELSE                    ; val_out = missing     ; ENDIF
                  CALL write_csv_point(j, tv%name, len_trim(tv%name), val_out, NULCSTR, 0, last=last)
                  CALL write_csv_point_avg(j, tv%name, len_trim(tv%name), cc(v, :), NULCSTR, 0, last=last)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE aed_write_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE wq_inflow_update(wqinf, nwqVars, temp, salt) BIND(C, name=_WQ_INFLOW_UPDATE)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(C_PTR),VALUE :: wqinf
   CINTEGER,INTENT(in) :: nwqVars
   AED_REAL,INTENT(inout) :: temp, salt
!LOCALS
   AED_REAL,DIMENSION(:),POINTER :: wqInfF
!BEGIN
   CALL C_F_POINTER(wqinf, wqInfF, [nwqVars])
   CALL aed_inflow_update(wqInfF, temp, salt)
END SUBROUTINE wq_inflow_update
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
CINTEGER FUNCTION aed_var_index_c(name, len) BIND(C, name=_WQ_VAR_INDEX_C)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CCHARACTER,INTENT(in) :: name(*)
   CSIZET,INTENT(in)     :: len
!LOCALS
   CHARACTER(len=len) :: tn
!BEGIN
   tn = trim(transfer(name(1:len),tn))
   aed_var_index_c = WQVar_Index(tn) - 1
END FUNCTION aed_var_index_c
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION WQVar_Index(name)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*) :: name
!
!LOCALS
   TYPE(aed_variable_t),POINTER :: tv
   INTEGER i,v
!
!-------------------------------------------------------------------------------
!BEGIN
   v = 0
   DO i=1, n_aed_vars
      IF ( aed_get_var(i, tv) ) THEN
         IF ( .NOT. tv%sheet .AND. tv%var_type == V_STATE ) THEN
            v = v + 1
            IF ( name .EQ. tv%name ) THEN
               WQVar_Index = v
               RETURN
            ENDIF
         ENDIF
      ENDIF
   ENDDO
   WQVar_Index = -1
END FUNCTION WQVar_Index
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_get_var_names(wbuf, wlen, bbuf, blen) BIND(C, name="aed_get_var_names")
!-------------------------------------------------------------------------------
! Fill comma-separated C strings with all water-column and benthic WQ state
! variable names in their AED index order (matching rows of wq_vars / wqs_vars).
!-------------------------------------------------------------------------------
!ARGUMENTS
   CCHARACTER, INTENT(out) :: wbuf(*), bbuf(*)
   CINTEGER,   INTENT(out) :: wlen, blen
!LOCALS
   TYPE(aed_variable_t), POINTER :: tvar
   INTEGER :: i, j, nlen, w, b
!
!-------------------------------------------------------------------------------
!BEGIN
   w = 0 ; b = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         IF ( tvar%var_type == V_STATE ) THEN
            nlen = LEN_TRIM(tvar%name)
            IF ( .NOT. tvar%sheet ) THEN
               IF (w > 0) THEN ; wbuf(w+1) = ',' ; w = w+1 ; ENDIF
               DO j=1,nlen ; wbuf(w+j) = tvar%name(j:j) ; ENDDO
               w = w + nlen
            ELSE
               IF (b > 0) THEN ; bbuf(b+1) = ',' ; b = b+1 ; ENDIF
               DO j=1,nlen ; bbuf(b+j) = tvar%name(j:j) ; ENDDO
               b = b + nlen
            ENDIF
         ENDIF
      ENDIF
   ENDDO
   wbuf(w+1) = ACHAR(0)
   bbuf(b+1) = ACHAR(0)
   wlen = w
   blen = b
END SUBROUTINE aed_get_var_names
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE glm_aed
