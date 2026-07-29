!###############################################################################
!#                                                                             #
!# glm_api_aed.F90                                                             #
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
!# Copyright 2024-2026 : The University of Western Australia                   #
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

#include "aed_api.h"

#undef MISVAL

#include "glm.h"

#define _UNUSED(x) if (.FALSE.) print*,shape(x)

!-------------------------------------------------------------------------------
MODULE glm_api_aed
!
   USE ISO_C_BINDING
   USE IEEE_ARITHMETIC

   USE aed_util
   USE aed_common
   USE glm_types
   USE aed_api
   USE aed_ptm
   USE aed_zones
   USE glm_api_zones

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
# define _WQ_SET_GLM_PTM     "api_set_glm_ptm"
# define _WQ_DO_GLM          "wq_do_glm"
# define _WQ_CLEAN_GLM       "wq_clean_glm"
# define _WQ_INIT_GLM_OUTPUT "wq_init_glm_output"
# define _WQ_WRITE_GLM_      "wq_write_glm"
# define _WQ_VAR_INDEX_C     "wq_var_index_c"
# define _WQ_IS_VAR          "wq_is_var"
# define _WQ_INFLOW_UPDATE   "wq_inflow_update"
#else
# define _WQ_INIT_GLM        "api_init_glm"
# define _WQ_SET_GLM_DATA    "api_set_glm_data"
# define _WQ_SET_GLM_PTM     "api_set_glm_ptm"
# define _WQ_DO_GLM          "api_do_glm"
# define _WQ_CLEAN_GLM       "api_clean_glm"
# define _WQ_INIT_GLM_OUTPUT "api_init_glm_output"
# define _WQ_WRITE_GLM_      "api_write_glm"
# define _WQ_VAR_INDEX_C     "api_var_index_c"
# define _WQ_IS_VAR          "api_is_var"
# define _WQ_INFLOW_UPDATE   "api_update_inflow_wq"
#endif




!
!-------------------------------------------------------------------------------
!MODULE DATA

  !AED_REAL :: par_fraction = 0.450
  !AED_REAL :: nir_fraction = 0.510
  !AED_REAL :: uva_fraction = 0.035
  !AED_REAL :: uvb_fraction = 0.005

   !# Arrays for state and diagnostic variables
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,TARGET :: cc    !# water quality array - water  : nvars, nlayers
!  AED_REAL,DIMENSION(:),  ALLOCATABLE,TARGET :: cc_hz !# water quality array - benthic: nvars
   AED_REAL,DIMENSION(:),  POINTER :: cc_hz !# water quality array - benthic: nvars
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,TARGET :: cc_diag
   AED_REAL,DIMENSION(:,:),ALLOCATABLE        :: cc_diag_old
   AED_REAL,DIMENSION(:),  ALLOCATABLE        :: z_prev
   AED_REAL,DIMENSION(:),  ALLOCATABLE        :: area_prev
   AED_REAL,DIMENSION(:),  ALLOCATABLE,TARGET :: cc_diag_hz

   LOGICAL,TARGET :: actv = .TRUE.
   INTEGER,TARGET :: botidx = 1, topidx = 0

   !# Arrays for work, vertical movement, and cross-boundary fluxes
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: dz
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: tss

   !# Arrays for environmental variables not supplied externally.
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: depth
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: layer_area  ! incremental area per layer (area(i)-area(i-1))

   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: pres
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: par
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: uva
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: uvb
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: nir

   !# External variables
   AED_REAL,DIMENSION(:),POINTER :: temp
   AED_REAL,DIMENSION(:),POINTER :: salt
   AED_REAL,DIMENSION(:),POINTER :: rad
   AED_REAL,DIMENSION(:),POINTER :: rho
   AED_REAL,DIMENSION(:),POINTER :: area
   AED_REAL,DIMENSION(:),POINTER :: extc
   AED_REAL,DIMENSION(:),POINTER :: layer_stress
   AED_REAL,DIMENSION(:),POINTER :: cvel

   AED_REAL,TARGET :: col_area

   AED_REAL,TARGET :: col_depth
   INTEGER,TARGET  :: col_num = 1

   AED_REAL,DIMENSION(:),POINTER :: sed_zones
   AED_REAL,TARGET :: matz

   AED_REAL,POINTER :: rain
   AED_REAL,POINTER :: evap
   AED_REAL,POINTER :: air_temp
   AED_REAL,POINTER :: humidity
   AED_REAL,POINTER :: I_0
   AED_REAL,POINTER :: longwave
   AED_REAL,POINTER :: wind
   AED_REAL,POINTER :: air_pres
   AED_REAL,POINTER :: u_star_api
   AED_REAL,POINTER :: Q_net_api
   AED_REAL,POINTER :: delzBlueIce_api
   AED_REAL,POINTER :: delzWhiteIce_api

   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: feedback

!  CHARACTER(len=48),ALLOCATABLE :: names(:)
!  CHARACTER(len=48),ALLOCATABLE :: bennames(:)
!  CHARACTER(len=48),ALLOCATABLE :: diagnames(:)

   INTEGER,DIMENSION(:),ALLOCATABLE :: externalid
   INTEGER,DIMENSION(:),ALLOCATABLE :: zexternalid
   INTEGER,DIMENSION(:),ALLOCATABLE :: plot_id_v, plot_id_sv, plot_id_d, plot_id_sd

   INTEGER :: n_aed_vars, n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet, n_ptm_vars

   CHARACTER(len=64) :: NULCSTR = ""

   TYPE(aed_part_group_t),DIMENSION(:),ALLOCATABLE :: ptm_bla


!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE api_init_glm(i_fname, len, NumWQ_Vars, NumWQ_Ben)                   &
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
   CHARACTER(len=80) :: fname

   TYPE(aed_coupling_t) :: conf

   !LOGICAL :: do_particle_bgc = .true. ! REMVE THIS LATER AND GET FROM GLOBALS
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL make_string(fname, i_fname, len)

   conf%glm_style_zones = .TRUE.

   conf%par_fraction =  0.450
   conf%nir_fraction =  0.510
   conf%uva_fraction =  0.035
   conf%uvb_fraction =  0.005

   conf%mobility_off = mobility_off
   conf%bioshade_feedback = bioshade_feedback
   conf%repair_state = repair_state
   conf%link_rain_loss = link_rain_loss
   conf%link_solar_shade = link_solar_shade
   conf%link_bottom_drag = link_bottom_drag

   conf%split_factor = split_factor
   conf%benthic_mode = benthic_mode

   conf%rain_factor => rain_factor
   conf%sw_factor => sw_factor
   conf%friction => friction

   conf%Kw => Kw

   conf%do_particle_bgc = do_particle_bgc

   conf%link_ext_par = link_ext_par

   CALL aed_set_coupling(conf)

   CALL api_set_glm_env()

   n_aed_vars = aed_configure_models(fname, n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet, n_ptm_vars)

   NumWQ_Vars  = n_vars
   NumWQ_Ben   = n_vars_ben
   CALL set_c_num_ptm_vars(n_ptm_vars)

   ALLOCATE(plot_id_v(n_vars))
   ALLOCATE(plot_id_sv(n_vars_ben))
   ALLOCATE(plot_id_d(n_vars_diag))
   ALLOCATE(plot_id_sd(n_vars_diag_sheet))
   plot_id_v = -1; plot_id_sv = -1; plot_id_d = -1; plot_id_sd = -1

   ALLOCATE(externalid(n_aed_vars))
   ALLOCATE(zexternalid(n_aed_vars))
END SUBROUTINE api_init_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_set_glm_env()
!-------------------------------------------------------------------------------
!ARGUMENTS
!
!LOCALS
   INTEGER :: status

   TYPE(aed_env_t) :: env(1)

   PROCEDURE(aed_mobility_fn_t),POINTER :: doMobilityP
!
!-------------------------------------------------------------------------------
!BEGIN
   !# Provide pointers to arrays with meteorological variables
   air_temp => MetData%AirTemp
   air_pres => MetData%AirPres
   humidity => MetData%RelHum
   wind     => MetData%WindSpeed
   rain     => MetData%Rain
   evap             => SurfData%Evap
   I_0              => MetData%ShortWave
   longwave         => MetData%LongWave
   u_star_api       => SurfData%u_star
   Q_net_api        => SurfData%Q_net
   delzBlueIce_api  => SurfData%delzBlueIce
   delzWhiteIce_api => SurfData%delzWhiteIce

   !# Set pointers to GLMs dynamic variables that will be updated later (in do_glm_wq)
   lheights => theLake%Height
   area     => theLake%LayerArea
   temp     => theLake%Temp
   salt     => theLake%Salinity
   rho      => theLake%Density
   cvel     => theLake%Umean
   rad      => theLake%Light
   extc     => theLake%ExtcCoefSW
   layer_stress => theLake%LayerStress

   !# Allocate arrays for various data constructs
   !# These will be calculated [approximated] from layers internally
   !# during each time step.

   ALLOCATE(depth(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (depth)'
   depth = one_

   ALLOCATE(layer_area(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (layer_area)'
   layer_area = zero_

   ALLOCATE(dz(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (dz)'
   dz = zero_

   ALLOCATE(sed_zones(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (sed_zones)'
   sed_zones = zero_

   ALLOCATE(pres(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (pres)'
   pres = zero_

   ALLOCATE(tss(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (tss)'
   tss = zero_

   ALLOCATE(par(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (par)'
   par = zero_

   ALLOCATE(nir(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (nir)'
   nir = zero_

   ALLOCATE(uva(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (uva)'
   uva = zero_

   ALLOCATE(uvb(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (uvb)'
   uvb = zero_

   ALLOCATE(feedback(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (feedback)'
   feedback = zero_

   matz = zero_

   timestep = dt

   env(1)%n_layers  =  MaxLayers
   env(1)%bot_idx   => botidx
   env(1)%top_idx   => topidx

   env(1)%timestep  => timestep
   env(1)%yearday   => yearday

   env(1)%longitude => longitude
   env(1)%latitude  => latitude
   env(1)%col_num   => col_num

   env(1)%active    => actv !# .true.

   env(1)%longwave      => longwave
   env(1)%air_temp      => air_temp
   env(1)%air_pres      => air_pres
   env(1)%humidity      => humidity
   env(1)%wind          => wind
   env(1)%rain          => rain
   env(1)%evap          => evap
   env(1)%I_0           => I_0

   env(1)%col_depth     => col_depth
   env(1)%col_area      => col_area
   env(1)%height        => lheights
   env(1)%depth         => depth
   env(1)%area          => layer_area   ! incremental area; cumulative theLake%LayerArea still available via 'area'
   env(1)%dz            => dz

   env(1)%temp          => temp
   env(1)%salt          => salt
   env(1)%cvel          => cvel
   env(1)%pres          => pres
   env(1)%rho           => rho
   env(1)%rad           => rad

   env(1)%extc          => extc
   env(1)%par           => par
   env(1)%nir           => nir
   env(1)%uva           => uva
   env(1)%uvb           => uvb

   env(1)%tss           => tss
   env(1)%ss1           =>  tss !ss1
   env(1)%ss2           =>  tss !ss2
   env(1)%ss3           =>  tss !ss3
   env(1)%ss4           =>  tss !ss4

   env(1)%ustar_bed     => feedback !ustar_bed
   env(1)%u_star        => u_star_api
   env(1)%Q_net         => Q_net_api
   env(1)%delzBlueIce   => delzBlueIce_api
   env(1)%delzWhiteIce  => delzWhiteIce_api
   env(1)%wv_uorb       => feedback !wv_uorb
   env(1)%wv_t          => feedback !wv_t
   env(1)%layer_stress  => layer_stress(1)
!  env(1)%dz_benthic    =>  ...

   env(1)%sed_zones     => sed_zones
   env(1)%sed_zone      => sed_zones(1)
   env(1)%mat_id        => matz

   env(1)%bathy         =>  feedback(1) !bathy(1)
!  env(1)%datum         =>  ...

   env(1)%biodrag       =>  feedback !bio_drag
   env(1)%bioextc       =>  feedback
   env(1)%solarshade    =>  feedback(1)
   env(1)%windshade     =>  feedback(1)
   env(1)%rainloss      =>  feedback(1) !rainloss(1)

   CALL aed_set_model_env(env, 1, MaxLayers)

   doMobilityP => doMobilityF
   CALL aed_set_mobility_fn(doMobilityP)
END SUBROUTINE api_set_glm_env
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_set_glm_data()                     BIND(C, name=_WQ_SET_GLM_DATA)
!-------------------------------------------------------------------------------
!ARGUMENTS
!
!LOCALS
   INTEGER :: status

   TYPE(aed_data_t) :: dat(1)
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (n_zones .GT. 0) &
      CALL api_set_glm_zones(n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet, n_aed_vars)

   !# Now that we know how many vars we need, we can allocate space for them
   ALLOCATE(cc((n_vars+n_vars_ben), MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (cc)'
   cc = zero_         !# initialise to zero

   CALL set_c_wqvars_ptr(cc)

   cc_hz => cc(n_vars+1:n_vars+n_vars_ben, 1)
   IF (n_vars_ben > 0) CALL set_c_wqsvars_ptr(cc_hz)

   !# Allocate diagnostic variable array and set all values to zero.
   !# (needed because time-integrated/averaged variables will increment rather than set the array)
   ALLOCATE(cc_diag(n_vars_diag, MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (cc_diag)'
   cc_diag = zero_
   ALLOCATE(cc_diag_old(n_vars_diag, MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (cc_diag_old)'
   cc_diag_old = zero_
   ALLOCATE(z_prev(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (z_prev)'
   z_prev = zero_
   ALLOCATE(area_prev(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (area_prev)'
   area_prev = zero_

   !# Allocate diagnostic variable array and set all values to zero.
   !# (needed because time-integrated/averaged variables will increment rather than set the array)
   ALLOCATE(cc_diag_hz(n_vars_diag_sheet),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (cc_diag_hz)'
   cc_diag_hz = zero_

   CALL set_c_wqdvars_ptr(cc_diag, cc_diag_hz, n_vars_diag, n_vars_diag_sheet)

   dat(1)%cc => cc(1:n_vars, :)
   dat(1)%cc_hz => cc(n_vars+1:n_vars+n_vars_ben, 1)
   dat(1)%cc_diag => cc_diag
   dat(1)%cc_diag_hz => cc_diag_hz

   CALL aed_set_model_data(dat, 1, MaxLayers)

   CALL aed_check_model_setup
END SUBROUTINE api_set_glm_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_set_glm_ptm(num_particle_groups,num_particles)    BIND(C, name=_WQ_SET_GLM_PTM)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: num_particle_groups,num_particles
!LOCALS
 ! INTEGER :: status
!
!-------------------------------------------------------------------------------
!BEGIN

   ALLOCATE(ptm_bla(num_particle_groups))

   CALL aed_ptm_init(num_particle_groups,num_particles,ptm_bla,n_ptm_vars,MaxLayers)

  ! CALL set_c_ptmvars_ptr(cc)
   CALL set_c_ptmstat_ptr(ptm_istat)
   CALL set_c_ptmenv_ptr(ptm_env)

END SUBROUTINE api_set_glm_ptm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
CINTEGER FUNCTION api_is_var(id,i_vname,len)            BIND(C, name=_WQ_IS_VAR)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in)   :: id
   CCHARACTER,INTENT(in) :: i_vname(*)
   CSIZET,INTENT(in)     :: len
!
!LOCALS
   CHARACTER(len=45) :: vname
   TYPE(aed_variable_t),POINTER :: tvar
   INTEGER :: i, v, sv, d, sd
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL make_string(vname, i_vname, len)

   v = 0; sv = 0; d = 0; sd = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         IF ( tvar%var_type == V_STATE ) THEN
            IF ( tvar%sheet ) THEN ; sv=sv+1; ELSE ; v=v+1 ; ENDIF
            IF ( TRIM(tvar%name) == vname ) THEN
               IF (tvar%sheet) THEN
                  api_is_var=-sv
                  plot_id_sv(sv) = id;
               ELSE
                  api_is_var=v
                  plot_id_v(v) = id;
               ENDIF
               RETURN
            ENDIF
         ELSEIF ( tvar%var_type == V_DIAGNOSTIC ) THEN
            IF ( tvar%sheet ) THEN ; sd=sd+1; ELSE ; d=d+1 ; ENDIF
            IF ( TRIM(tvar%name) == vname ) THEN
               IF (tvar%sheet) THEN
                  api_is_var=-sd
                  plot_id_sd(sd) = id;
               ELSE
                  api_is_var=d
                  plot_id_d(d) = id;
               ENDIF
               RETURN
            ENDIF
         ENDIF
      ENDIF
   ENDDO
   api_is_var = 0
END FUNCTION api_is_var
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE doMobilityF(N,dt,h,A,ww,min_C,mcc)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in)    :: N       !# number of vertical layers
   AED_REAL,INTENT(in)    :: dt      !# time step (s)
   AED_REAL,INTENT(in)    :: h(:)    !# layer thicknesses (m)
   AED_REAL,INTENT(in)    :: A(:)    !# layer areas (m2)
   AED_REAL,INTENT(in)    :: ww(:)   !# vertical speed (m/s)
   AED_REAL,INTENT(in)    :: min_C   !# minimum allowed cell concentration
   AED_REAL,INTENT(inout) :: mcc(:)  !# cell concentration
!
!LOCALS
!
!-------------------------------------------------------------------------------
!BEGIN
   _UNUSED(A)

!## The "Woods" example shows the problem.  Either version, "A" or "area" work with ifx
!## but, with gfortran - using local copy "area" works, using passed array "A" doesnt
!## Woods lake example Oxy and Oxy_sat NaN out with the A passed in.
!  CALL doMobility(N,dt,h,A,ww,min_C,mcc)
   CALL doMobility(N,dt,h,area,ww,min_C,mcc)
END SUBROUTINE doMobilityF
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_do_glm(wlev)                             BIND(C, name=_WQ_DO_GLM)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: wlev
!
!LOCALS
   LOGICAL :: doSurface
   INTEGER :: lev, zon, j, k, j_d
   INTEGER, SAVE :: wlev_prev = 0
   AED_REAL :: surf, pa, ratio
   AED_REAL :: I_old, I_new
   TYPE(aed_variable_t),POINTER :: tv
   LOGICAL :: grid_changed
!
!-------------------------------------------------------------------------------
!BEGIN
   pa = zero_
   topidx = wlev ; botidx = 1
   col_area = area(wlev)

   col_depth = lheights(wlev)
   surf = lheights(wlev)
   !# re-compute the layer heights and depths
   dz(1) = lheights(1)
   depth(1) = surf - lheights(1)
   layer_area(1) = area(1)
   DO lev=2,wlev
      dz(lev) = lheights(lev) - lheights(lev-1)
      depth(lev) = surf - lheights(lev)
      layer_area(lev) = area(lev) - area(lev-1)   ! incremental area: cumulative gained at this layer
   ENDDO

   !# Calculate local pressure
   pres(1:wlev) = -lheights(1:wlev)

   IF ( benthic_mode .GT. 1 ) THEN
      zon = 1
      DO lev=1,wlev
         IF (lheights(lev) .GT. aedZones(zon)%z_env%z_height) THEN
            sed_zones(lev) = zon * ( area(lev) - pa ) / area(lev)
            pa = area(lev)
            zon = zon+1
         ELSE
            sed_zones(lev) = zon
         ENDIF
      ENDDO
   ENDIF

   !# Reset benthic sheet diagnostics each timestep — api_copy_from_zone accumulates into cc_diag_hz
   cc_diag_hz = 0.

   !# Remap the PERSISTED (rezero=.FALSE.) zavg diagnostics onto the current
   !# layer grid. The host zeroes the layer-indexed cc_diag during between-step
   !# layer restructuring (whenever the grid changes), so the value visible at the
   !# top of this routine is unreliable. cc_diag_old, snapshotted at the END of the
   !# previous step (see bottom of routine) on the previous grid, holds the true
   !# post-benthic value and is the remap source. This lets the bubble module read
   !# the previous step's ebullition, correctly remapped to the current layers.
   grid_changed = ( wlev /= wlev_prev )
   IF ( .NOT. grid_changed .AND. wlev_prev > 0 ) grid_changed = ANY( lheights(1:wlev) /= z_prev(1:wlev) )
   IF ( grid_changed .AND. wlev_prev > 0 ) THEN
      !# Grid changed: interpolate cc_diag_old (previous grid) onto the new grid.
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

                  !# Conserve the area-integral for persisted zavg flux diags
                  !# (e.g. ch4_ebb_dsfv, mmol/m2/s). Pointwise height interpolation
                  !# preserves the per-area shape but not Sum(diag*dArea) -- the
                  !# total the bubble module integrates -- so rescale to the
                  !# previous step's total. (rezero=.TRUE. diags are wiped in
                  !# aed_run_model, so only persisted ones reach a reader.)
                  IF ( .NOT. tv%rezero ) THEN
                     I_old = zero_
                     DO lev = 1, wlev_prev
                        IF (lev == 1) THEN
                           I_old = I_old + cc_diag_old(j_d,lev)*area_prev(lev)
                        ELSE
                           I_old = I_old + cc_diag_old(j_d,lev)*(area_prev(lev)-area_prev(lev-1))
                        ENDIF
                     ENDDO
                     I_new = zero_
                     DO lev = 1, wlev
                        IF (lev == 1) THEN
                           I_new = I_new + cc_diag(j_d,lev)*area(lev)
                        ELSE
                           I_new = I_new + cc_diag(j_d,lev)*(area(lev)-area(lev-1))
                        ENDIF
                     ENDDO
                     IF (I_new > 1.0e-30) cc_diag(j_d,:) = cc_diag(j_d,:) * (I_old/I_new)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ELSEIF ( wlev_prev > 0 ) THEN
      !# Grid unchanged: restore persisted zavg diags from the snapshot (exact).
      j_d = 0
      DO k = 1, n_aed_vars
         IF ( aed_get_var(k, tv) ) THEN
            IF ( tv%var_type == V_DIAGNOSTIC .AND. .NOT. tv%sheet ) THEN
               j_d = j_d + 1
               IF ( tv%zavg .AND. .NOT. tv%rezero ) cc_diag(j_d,1:wlev) = cc_diag_old(j_d,1:wlev)
            ENDIF
         ENDIF
      ENDDO
   ENDIF

   doSurface = .NOT. ice
   CALL aed_run_model(1, wlev, doSurface)

   !# Save layer count, heights and cumulative areas for interpolation on next timestep
   wlev_prev = wlev
   z_prev(1:wlev) = lheights(1:wlev)
   area_prev(1:wlev) = area(1:wlev)

   !# Snapshot the persisted zavg diag AFTER the benthic update, before the host
   !# can zero it during layer restructuring. This is the remap source next step.
   cc_diag_old = cc_diag

END SUBROUTINE api_do_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_clean_glm()                           BIND(C, name=_WQ_CLEAN_GLM)
!-------------------------------------------------------------------------------
! Finish biogeochemical model
!-------------------------------------------------------------------------------
!BEGIN
   CALL aed_clean_model()
   ! Deallocate internal arrays
!  IF (ALLOCATED(cc))         DEALLOCATE(cc)   !# left allocated: shared with C side
!   IF (ASSOCIATED(cc_hz))     NULLIFY(cc_hz)
!   IF (ALLOCATED(cc_diag))    DEALLOCATE(cc_diag)
!   IF (ALLOCATED(cc_diag_hz)) DEALLOCATE(cc_diag_hz)
   IF (ALLOCATED(cc_diag_old)) DEALLOCATE(cc_diag_old)
   IF (ALLOCATED(z_prev))      DEALLOCATE(z_prev)
   IF (ALLOCATED(par))        DEALLOCATE(par)
   IF (ALLOCATED(nir))        DEALLOCATE(nir)
   IF (ALLOCATED(uva))        DEALLOCATE(uva)
   IF (ALLOCATED(uvb))        DEALLOCATE(uvb)
   IF (ALLOCATED(pres))       DEALLOCATE(pres)
   IF (ALLOCATED(dz))         DEALLOCATE(dz)
   IF (ALLOCATED(tss))        DEALLOCATE(tss)
   IF (ALLOCATED(depth))      DEALLOCATE(depth)
   IF (ALLOCATED(layer_area)) DEALLOCATE(layer_area)
   IF (ASSOCIATED(sed_zones)) DEALLOCATE(sed_zones)
   IF (ALLOCATED(feedback))   DEALLOCATE(feedback)
   IF (ALLOCATED(externalid)) DEALLOCATE(externalid)
   IF (ALLOCATED(zexternalid))DEALLOCATE(zexternalid)
   IF (ALLOCATED(plot_id_v))  DEALLOCATE(plot_id_v)
   IF (ALLOCATED(plot_id_sv)) DEALLOCATE(plot_id_sv)
   IF (ALLOCATED(plot_id_d))  DEALLOCATE(plot_id_d)
   IF (ALLOCATED(plot_id_sd)) DEALLOCATE(plot_id_sd)
   IF (ALLOCATED(ptm_bla))    DEALLOCATE(ptm_bla)
END SUBROUTINE api_clean_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_init_glm_output(ncid,x_dim,y_dim,z_dim,zone_dim,time_dim)       &
                                               BIND(C, name=_WQ_INIT_GLM_OUTPUT)
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
            IF ( tv%sheet .AND. (tv%var_type == V_DIAGNOSTIC .OR. tv%var_type == V_STATE) ) THEN
               !# only for state and diag sheet vars
               zexternalid(i) = NEW_NC_VARIABLE(ncid, TRIM(tv%name)//"_Z", LEN_TRIM(tv%name)+2, NF90_REALTYPE, 4, dims(1:4))
               CALL set_nc_attributes(ncid, zexternalid(i), MYTRIM(tv%units), MYTRIM(tv%longname) PARAM_FILLVALUE)
            ENDIF
         ENDIF
      ENDDO
   ENDIF

   !# Take NetCDF library out of define mode (ready for storing data).
   CALL define_mode_off(ncid)
END SUBROUTINE api_init_glm_output
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!  Save properties of biogeochemical model, including state variable
!  values, diagnostic variable values, and sums of conserved quantities.
!-------------------------------------------------------------------------------
SUBROUTINE api_write_glm(ncid,wlev,nlev,lvl,point_nlevs) BIND(C, name=_WQ_WRITE_GLM_)
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
   FLOGICAL :: last
!
!-------------------------------------------------------------------------------
!BEGIN
   last = .FALSE.
   v = 0; d = 0; sv = 0; sd = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tv) ) THEN
         IF ( tv%var_type == V_DIAGNOSTIC ) THEN
            !# Process and store diagnostic variables.
            IF ( tv%sheet ) THEN
               sd = sd + 1
               !# Process and store diagnostic variables defined on horizontal slices of the domain.
               IF ( n_zones .GT. 0 ) THEN
                  z_diag_hz(sd,n_zones+1) = cc_diag_hz(sd)
                  CALL store_nc_array(ncid, zexternalid(i), XYNT_SHAPE, n_zones, n_zones, array=z_diag_hz(sd,1:n_zones+1))
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
         ELSE IF ( tv%var_type == V_STATE ) THEN  ! not diag
            IF ( tv%sheet ) THEN
               sv = sv + 1
               !# Store benthic biogeochemical state variables.
               IF ( n_zones .GT. 0 ) THEN
                  CALL store_nc_array(ncid, zexternalid(i), XYNT_SHAPE, n_zones, n_zones, array=z_cc_hz(sv, 1:n_zones))
               ENDIF
               CALL store_nc_scalar(ncid, externalid(i), XYT_SHAPE, scalar=cc_hz(sv))
#ifdef PLOTS
               IF ( do_plots .AND. plot_id_sv(sv).GE.0 ) THEN
                  IF ( n_zones .GT. 0 ) THEN
                     DO z=1,n_zones
                        CALL put_glm_val_z(plot_id_sv(sv), z_cc_hz(sv, z), z)
                     ENDDO
                  ENDIF
                  CALL put_glm_val_s(plot_id_sv(sv), cc_hz(sv))
               ENDIF
#endif
               DO j=1,point_nlevs
                  val_out = missing
                  IF ((lvl(j) .EQ. wlev) .AND. tv%top) val_out = cc_hz(sv)
                  IF ((lvl(j) .EQ. 0)    .AND. tv%bot) val_out = cc_hz(sv)
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
END SUBROUTINE api_write_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE wq_inflow_update(wqinf, nwqVars, temp, salt)                        &
                                                 BIND(C, name=_WQ_INFLOW_UPDATE)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(C_PTR),VALUE :: wqinf
   CINTEGER,INTENT(in) :: nwqVars
   AED_REAL,INTENT(inout) :: temp, salt
!
!LOCALS
   AED_REAL,DIMENSION(:),POINTER :: wqInfF
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL C_F_POINTER(wqinf, wqInfF, [nwqVars])
   CALL aed_inflow_update(wqInfF, temp, salt)
END SUBROUTINE wq_inflow_update
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
CINTEGER FUNCTION aed_var_index_c(name, len)       BIND(C, name=_WQ_VAR_INDEX_C)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CCHARACTER,INTENT(in) :: name(*)
   CSIZET,INTENT(in)     :: len
!
!LOCALS
   CHARACTER(len=len) :: tn
!
!-------------------------------------------------------------------------------
!BEGIN
   tn = trim(transfer(name(1:len),tn))
   aed_var_index_c = aed_var_index(tn) - 1
END FUNCTION aed_var_index_c
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_get_var_names(wbuf, wlen, bbuf, blen) BIND(C, name="api_get_var_names")
!-------------------------------------------------------------------------------
! Fill comma-separated C strings with all water-column and benthic WQ state
! variable names in their AED index order (matching rows of wq_vars / wqs_vars).
! Caller convention (strlen-style): wlen/blen are char counts excluding the
! NUL terminator; wbuf/bbuf must be allocated to at least wlen+1 / blen+1 bytes
! (the routine writes ACHAR(0) at wbuf(wlen+1) and bbuf(blen+1)).
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
END SUBROUTINE api_get_var_names
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE glm_api_aed
