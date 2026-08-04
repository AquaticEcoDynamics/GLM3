!###############################################################################
!#                                                                             #
!# glm_zones.F90                                                               #
!#                                                                             #
!# The sediment zone processing bit for WQ                                     #
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

#undef MISVAL

#include "glm.h"

MODULE glm_zones

   USE ISO_C_BINDING

   USE glm_types
   USE aed_common

   IMPLICIT NONE

   PRIVATE ! By default, make everything private

   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: z_cc      ! (n_vars, n_levs, nsed_zones)
!  AED_REAL,DIMENSION(:,:),  ALLOCATABLE,TARGET :: z_cc_hz   ! (n_vars,         nsed_zones)
   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: z_diag    ! (n_vars, n_levs, nsed_zones)
   AED_REAL,DIMENSION(:,:),  ALLOCATABLE,TARGET :: z_diag_hz ! (n_vars,         nsed_zones)

   AED_REAL,DIMENSION(:),POINTER :: lheights

   AED_REAL,DIMENSION(:),POINTER :: zone_heights
   INTEGER :: nvars, nbenv, nvdiag, nvdiag_hz

   PUBLIC n_zones, zone_heights, lheights
   PUBLIC wq_set_glm_zones, copy_from_zone, copy_to_zone, calc_zone_areas

!  PUBLIC z_cc, z_cc_hz, z_diag, z_diag_hz, theZones, theLake
   PUBLIC z_cc, z_diag, z_diag_hz, theZones, theLake

CONTAINS

!###############################################################################
SUBROUTINE wq_set_glm_zones(numVars, numBenV, numDiagV, numDiagHzV)            &
                                                BIND(C, name="wq_set_glm_zones")
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: numVars, numBenV, numDiagV, numDiagHzV
!
!LOCALS
!
!-------------------------------------------------------------------------------
!BEGIN
   nvars = numVars
   nbenv = numBenV
   nvdiag = numDiagV
   nvdiag_hz = numDiagHzV

   lheights => theLake%Height
   zone_heights => theZones%zheight

   ALLOCATE(z_cc(numVars+numBenV, MaxLayers, n_zones)) ; z_cc = 0.
!  ALLOCATE(z_cc_hz(numVars+numBenV, n_zones))         ; z_cc_hz = 0.
   ALLOCATE(z_diag(numDiagV, MaxLayers, n_zones))      ; z_diag = 0.
   ALLOCATE(z_diag_hz(numDiagHzV, n_zones+1))          ; z_diag_hz = 0.
   theZones%zarea = 0.
END SUBROUTINE wq_set_glm_zones
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE calc_zone_areas(areas, wlev, surf)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,DIMENSION(:),INTENT(in) :: areas
   INTEGER,INTENT(in) :: wlev
   AED_REAL,INTENT(in) :: surf
!
!LOCALS
   INTEGER  :: lev, zon
   LOGICAL  :: w_zones
   AED_REAL :: scale
!
!-------------------------------------------------------------------------------
!BEGIN
   theZones%zarea = 0.  ; theZones%z_pc_wet = 0.
   w_zones = .FALSE.

   zon = 1
   theZones(1)%zarea = areas(1)
   IF ( zone_heights(1) > surf ) THEN
      theZones(1)%z_pc_wet = surf / zone_heights(1)
   ELSE
      theZones(1)%z_pc_wet = 1.0
   ENDIF
   DO lev=2, wlev
      IF (lheights(lev) <= zone_heights(zon)) THEN
         theZones(zon)%zarea = theZones(zon)%zarea + areas(lev) - areas(lev-1)
      ELSEIF (lheights(lev) > zone_heights(zon) .AND. lheights(lev-1) < zone_heights(zon)) THEN
         scale = (zone_heights(zon) - lheights(lev-1)) / (lheights(lev) - lheights(lev-1))
         theZones(zon)%zarea = theZones(zon)%zarea + (areas(lev) - areas(lev-1)) * scale
         zon = zon + 1
         theZones(zon)%zarea = theZones(zon)%zarea + (areas(lev) - areas(lev-1)) * (1-scale)
      ENDIF

      IF ( zone_heights(zon) > surf ) THEN
         IF (.NOT. w_zones) THEN
            w_zones = .TRUE.
            IF ( zon > 1 ) THEN
               theZones(zon)%z_pc_wet = surf / (zone_heights(zon) - zone_heights(zon-1))
            ELSE
               theZones(zon)%z_pc_wet = surf / zone_heights(zon)
            ENDIF
         ENDIF
      ELSE
         TheZones(zon)%z_pc_wet = 1.0
      ENDIF
   ENDDO

   theZones(1:n_zones)%zpres = -zone_heights(1:n_zones)
END SUBROUTINE calc_zone_areas
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE copy_to_zone(x_cc, x_diag, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,DIMENSION(:,:),INTENT(in) :: x_cc
   AED_REAL,DIMENSION(:,:),INTENT(in) :: x_diag
   INTEGER,INTENT(in) :: wlev
!
!LOCALS
   INTEGER  :: zon, lev, a_zones
   AED_REAL :: surf
   INTEGER  :: zcount(n_zones)
   LOGICAL  :: w_zones
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Reset zone data structure to zero, for copying in information from main
   ! lake data structures. Note that z_cc(:,:,nvars+1:nvars_ben) is not zeroed
   ! as this is the benthic zone data we need to preserve
   z_cc(1:nvars,:,:) = 0.
   z_diag(:,:,:) = 0.
   z_diag_hz(:,:) = 0.

   ! Initialise the zone environment
   theZones%zrad = 0.
   theZones%zsalt = 0.
   theZones%ztemp = 0.
   theZones%zrho = 0.
   theZones%zextc = 0.
   theZones%zlayer_stress = 0.
   theZones%ztss = 0.
   theZones%zpar = 0.
   theZones%znir = 0.
   theZones%zuva = 0.
   theZones%zuvb = 0.
   theZones(1)%z_sed_zones = 1.
   theZones%zvel = 0.

   ! Populate the 1st layer in each zone structure, with the zone-averaged quantity
   a_zones = 1
   zcount = 0
   w_zones = .FALSE.
   zon = 1
   DO lev=1,wlev
      IF ( lev > 1 .AND. lheights(lev) > zone_heights(zon) ) THEN
         zon = zon + 1
         IF (zon > n_zones) STOP 'Water level height is higher than highest zone height'
         theZones(zon)%z_sed_zones = zon
      ENDIF

      ! Populate the 1st layer in each zone structure, with the zone-averaged quantity,
      ! That is pelagic ("cc") variables need zonifying, in case benthic modules want
      ! to refer to them.
      ! Note that this is not be done on benthic vars (nvars+1->) because it will
      ! introduce errors in z_cc (split layers)
      ! Ideally this average would be based on volume weighting, rather than count-based

      z_cc(1:nvars,1,zon) = z_cc(1:nvars,1,zon) + x_cc(1:nvars,lev)  ! 1:nvars is the water coumn "cc". not benthic
      z_diag(:,zon,zon)   = z_diag(:,zon,zon)   + x_diag(:,lev)      ! store at level=zon so cell(zon) is correct with layer_idx=zon
!     z_diag_hz(:,zon)    = z_diag_hz(:,zon)    + x_diag_hz(:)       !

      theZones(zon)%ztemp         = theZones(zon)%ztemp + theLake(lev)%Temp
      theZones(zon)%zsalt         = theZones(zon)%zsalt + theLake(lev)%Salinity
      theZones(zon)%zrho          = theZones(zon)%zrho  + theLake(lev)%Density
      theZones(zon)%zrad          = theZones(zon)%zrad  + theLake(lev)%Light
      theZones(zon)%zpar          = theZones(zon)%zpar  + theLake(lev)%Light*0.45
      theZones(zon)%zvel          = theZones(zon)%zvel  + theLake(lev)%Umean
      theZones(zon)%zextc         = theZones(zon)%zextc + theLake(lev)%ExtcCoefSW
      theZones(zon)%zlayer_stress = theZones(zon)%zlayer_stress + theLake(lev)%LayerStress

      zcount(zon) = zcount(zon) + 1
   ENDDO
   a_zones = zon ! Available zones is the last zone we found, at the current water level

   DO zon=1,a_zones
      ! Finalise zone averaged information (1st layer in each zone structure reserved for zavg)
      z_cc(1:nvars,1,zon) = z_cc(1:nvars,1,zon)/zcount(zon)  ! water column state vars
      z_diag(:,zon,zon)   = z_diag(:,zon,zon)/zcount(zon)    ! water column diag vars (level=zon to match layer_idx=zon)
      z_diag_hz(:,zon)    = z_diag_hz(:,zon)/zcount(zon)     ! benthic diag vars

      ! Set the water column above a zone, to the respective water layer values
!     z_cc(1:nvars,2:wlev,zon) = x_cc(1:nvars,2:wlev)        ! water column state vars
!     z_diag(:,2:wlev,zon)     = x_diag(:,2:wlev)            ! water column diag vars
   ENDDO

   WHERE (zcount /= 0)
      theZones%ztemp         = theZones%ztemp / zcount
      theZones%zsalt         = theZones%zsalt / zcount
      theZones%zrho          = theZones%zrho  / zcount
      theZones%zrad          = theZones%zrad  / zcount
      theZones%zpar          = theZones%zpar  / zcount
      theZones%zvel          = theZones%zvel  / zcount
      theZones%zextc         = theZones%zextc / zcount
      theZones%zlayer_stress = theZones%zlayer_stress / zcount
   ELSEWHERE
      theZones%ztemp         = 0.
      theZones%zsalt         = 0.
      theZones%zrho          = 0.
      theZones%zrad          = 0.
      theZones%zpar          = 0.
      theZones%zvel          = 0.
      theZones%zextc         = 0.
      theZones%zlayer_stress = 0.
   ENDWHERE

   theZones%zdz = 0.
   surf = lheights(wlev)
   IF ( surf > zone_heights(1) ) THEN
      theZones(1)%zdepth = zone_heights(1)
      theZones(1)%zdz = zone_heights(1)
   ELSE
      theZones(1)%zdepth = surf
      theZones(1)%zdz = surf
      w_zones = .TRUE.
   ENDIF

   DO zon=2,a_zones
      theZones(zon)%z_sed_zones = zon
      IF ( .NOT. w_zones ) THEN
          IF ( surf > zone_heights(zon) ) THEN
             theZones(zon)%zdepth = zone_heights(zon)
             theZones(zon)%zdz =  zone_heights(zon) - zone_heights(zon-1)
          ELSE
             theZones(zon)%zdepth = surf
             theZones(zon)%zdz = surf - zone_heights(zon-1)
             w_zones = .TRUE.
          ENDIF
      ELSE
         theZones(zon)%zdepth = surf
         theZones(zon)%zdz = surf - zone_heights(zon-1)
      ENDIF
   ENDDO

END SUBROUTINE copy_to_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE copy_from_zone(n_aed_vars, x_cc, x_diag, x_diag_hz, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,DIMENSION(:,:),INTENT(inout) :: x_cc
   AED_REAL,DIMENSION(:,:),INTENT(inout) :: x_diag
   AED_REAL,DIMENSION(:),INTENT(inout) :: x_diag_hz
   INTEGER,INTENT(in) :: wlev, n_aed_vars
!
!LOCALS
   INTEGER  :: zon, lev, v_start, v_end
   AED_REAL :: scale, area
   LOGICAL  :: splitZone
!  LOGICAL  :: column_benthic_var_averaging = .false.
!  INTEGER  :: water_column_zone = 1
   INTEGER  :: i, j
   TYPE(aed_variable_t),POINTER :: tvar
!
!-------------------------------------------------------------------------------
!BEGIN
   v_start = nvars+1 ; v_end = nvars+nbenv

!print*,"size x_cc 1", size(x_cc, 1)
!print*,"size x_cc 2", size(x_cc, 2)
!print*,"size x_cc_hz 1", size(x_cc_hz, 1)
!
!print*,"size z_cc 1", size(z_cc, 1)
!print*,"size z_cc 2", size(z_cc, 2)
!print*,"size z_cc 3", size(z_cc, 3)
!print*,"size z_cc_hz 1", size(z_cc_hz, 1)
!print*,"size z_cc_hz 2", size(z_cc_hz, 2)
!
!print*,"Z1 cc(:,1) = ", x_cc(:,1), " xch ",x_cc_hz(:), " zcc ",z_cc(:,1,1), " zch ",z_cc_hz(:,1), nvars

   zon = n_zones
   ! Loop down through water layers
   DO lev=wlev,1,-1
      ! Check if zone boundary is in this water layer range
      IF ( zon .GT. 1 ) THEN
         IF (lev .GT. 1) THEN
            splitZone = lheights(lev-1) < zone_heights(zon-1)
         ELSE
            splitZone = 0.0 < zone_heights(zon-1)
         ENDIF
      ELSE
         splitZone = .FALSE.
      ENDIF

      ! Set water layer variables, based on zone infomration, where variable is flagged for zavg
      IF (splitZone) THEN
         ! Compute layer fraction
         IF (lev .GT. 1) THEN
            scale = (zone_heights(zon-1) - lheights(lev-1)) / (lheights(lev) - lheights(lev-1))
         ELSE
            scale = (zone_heights(zon-1) - 0.0) / (lheights(lev) - 0.0)
         ENDIF

         ! Select the diag vars that have zavg == true, and assign to layer
         ! (for first zone in the layer)
         j = 0
         DO i=1,n_aed_vars
            IF ( aed_get_var(i, tvar) ) THEN
               !print *,'zav', i, TRIM(tvar%name),tvar%zavg
               IF ( tvar%var_type == V_DIAGNOSTIC ) THEN
                  IF ( .NOT.  tvar%sheet ) THEN
                     j = j + 1
                     IF ( tvar%zavg ) THEN
                        x_diag(j, lev) = z_diag(j, zon, zon) * scale   ! was (lev,j) — swapped indices; level=zon matches copy_to_zone
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         ! Select the benthic vars (implicitly all zavg == true), and assign to layer
         ! (for first zone in the layer)
         x_cc(v_start:v_end,lev) = z_cc(v_start:v_end,1,zon) * scale ! KK only happening on benthic

         ! Select the diag vars that have zavg == true, and assign to layer
         ! (for second zone in the layer)
         zon = zon - 1
         j = 0
         DO i=1,n_aed_vars
            IF ( aed_get_var(i, tvar) ) THEN
               IF ( tvar%var_type == V_DIAGNOSTIC ) THEN
                  IF ( .NOT.  tvar%sheet ) THEN
                     j = j + 1
                     IF ( tvar%zavg ) THEN
                        x_diag(j, lev) = x_diag(j, lev) + (z_diag(j, zon, zon) * (1.0 - scale))
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         ! Select the benthic vars (implicitly all zavg == true), and assign to layer
         ! (for first zone in the layer)
         x_cc(v_start:v_end, lev) = x_cc(v_start:v_end, lev) + &
                                    z_cc(v_start:v_end, 1, zon) * (1.0 - scale)

      ELSE  ! Not a split layer - layer bounds are fully within zone
         ! Select the diag vars that have zavg == true, and assign to layer
         j = 0
         DO i=1,n_aed_vars
            IF ( aed_get_var(i, tvar) ) THEN
               IF ( tvar%var_type == V_DIAGNOSTIC ) THEN
                  IF ( .NOT.  tvar%sheet ) THEN
                     j = j + 1
                     IF ( tvar%zavg ) THEN
                        x_diag(j, lev) = z_diag(j, zon, zon)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         ! Select the benthic vars (implicitly all zavg == true), and assign to layer
         x_cc(v_start:v_end, lev) = z_cc(v_start:v_end, 1, zon)
      ENDIF
   ENDDO

   ! IF column_benthic_var_averaging, set single-value to the mean, weighted by area
   area = SUM(theZones(1:n_zones)%zarea)
   DO zon=1, n_zones
      x_diag_hz = x_diag_hz + (z_diag_hz(:, zon) * (theZones(zon)%zarea/area))
   ENDDO
END SUBROUTINE copy_from_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GInitialTemp(m,depth,wv,topTemp,botTemp,nSPinUpDays,tNew) BIND(C, name="InitialTemp")
!-------------------------------------------------------------------------------
   USE aed_util, ONLY : InitialTemp
!ARGUMENTS
   INTEGER,intent(in)   :: m
   AED_REAL,intent(in)  :: wv,depth(0:m+1)
   AED_REAL,intent(in)  :: topTemp,botTemp,nSPinUpDays
   AED_REAL,intent(out) :: tNew(0:m+1)
!-------------------------------------------------------------------------------
!BEGIN
   CALL InitialTemp(m, depth, wv, topTemp, botTemp, nSPinUpDays, tNew)
END SUBROUTINE GInitialTemp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GInitialTempFV(n,depth,vwc_init,surfTemp,deepTemp,spinup_days,dt,tNew) &
                                                  BIND(C, name="InitialTempFV")
!-------------------------------------------------------------------------------
!# C-callable spin-up for the cell-centred soil-temperature model.  Reads the
!# soil thermal properties from the (Fortran-owned) glm_types globals and
!# returns the n interior-cell temperatures in tNew(1:n).
!-------------------------------------------------------------------------------
   USE aed_util, ONLY : InitialTempFV, soil_props_t
!ARGUMENTS
   INTEGER,intent(in)   :: n
   AED_REAL,intent(in)  :: depth(0:n+1), vwc_init, surfTemp, deepTemp, spinup_days, dt
   AED_REAL,intent(out) :: tNew(n)
!LOCALS
   TYPE(soil_props_t) :: props
!-------------------------------------------------------------------------------
!BEGIN
   props%k_mineral = sed_k_mineral ; props%k_water = sed_k_water
   props%c_mineral = sed_c_mineral ; props%c_water = sed_c_water ; props%c_air = sed_c_air
   props%bulk_density = sed_bulk_density ; props%mineral_density = sed_mineral_density
   props%porosity = sed_porosity
   CALL InitialTempFV(n, depth, vwc_init, surfTemp, deepTemp, spinup_days, dt, props, tNew)
END SUBROUTINE GInitialTempFV
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE ZSoilTemp(izone) BIND(C, name="zZSoilTemp")
!-------------------------------------------------------------------------------
   USE aed_util, ONLY : SoilTempFV, soil_props_t
!ARGUMENTS
   !# izone is a C pointer passed BY VALUE from glm_surface.c (ZSoilTemp(&theZones[z]));
   !# it MUST carry the VALUE attribute or Fortran treats it as pass-by-reference and
   !# dereferences the struct's first word as the address.
   TYPE(C_PTR),VALUE :: izone
!LOCALS
   TYPE(ZoneType),POINTER :: zone
   TYPE(SedLayerType),DIMENSION(:),POINTER :: layers
   INTEGER :: nsl, n, kk
   AED_REAL :: hf, deepTemp
   AED_REAL,DIMENSION(:),ALLOCATABLE :: sdepth, cellvwc, cellt
   TYPE(soil_props_t) :: props
!-------------------------------------------------------------------------------
!BEGIN
   CALL C_F_POINTER(izone, zone);
   nsl = zone%n_sed_layers
   n   = nsl - 2
   CALL C_F_POINTER(zone%c_layers, layers, [nsl]);
   !# Cell-centred layout: node depths d(0:n+1) = layers(1:nsl)%depth; the n
   !# interior cells (the prognostic state) live in slots 2..nsl-1; layers(1) is
   !# the surface boundary node and layers(nsl) the deep boundary node.
   ALLOCATE(sdepth(0:n+1), cellvwc(n), cellt(n))
   DO kk = 1, nsl
      sdepth(kk-1) = layers(kk)%depth
   ENDDO
   DO kk = 1, n
      cellvwc(kk) = layers(kk+1)%vwc
      cellt(kk)   = layers(kk+1)%temp
   ENDDO
   deepTemp = layers(nsl)%temp
   props%k_mineral = sed_k_mineral ; props%k_water = sed_k_water
   props%c_mineral = sed_c_mineral ; props%c_water = sed_c_water ; props%c_air = sed_c_air
   props%bulk_density = sed_bulk_density ; props%mineral_density = sed_mineral_density
   props%porosity = sed_porosity
   CALL SoilTempFV(n, sdepth, cellvwc, zone%ztemp, deepTemp, soil_dt, props, cellt, hf)
   zone%heatflux = hf
   DO kk = 1, n
      layers(kk+1)%temp = cellt(kk)
   ENDDO
   layers(1)%temp = zone%ztemp
   DEALLOCATE(sdepth, cellvwc, cellt)
END SUBROUTINE ZSoilTemp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GSoilTemp(m,depth,wv,topTemp,temp,heatflux) BIND(C, name="SoilTemp")
!-------------------------------------------------------------------------------
   USE aed_util, ONLY : SoilTemp
!ARGUMENTS
   INTEGER,intent(in) :: m
   AED_REAL,intent(in) :: depth(0:m+1), wv(m), topTemp
   AED_REAL,intent(inout) :: temp(m+1)
   AED_REAL, intent(out) :: heatflux

!-------------------------------------------------------------------------------
!BEGIN
   CALL SoilTemp(m, depth, wv, topTemp, temp, heatflux)
END SUBROUTINE GSoilTemp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE glm_zones
