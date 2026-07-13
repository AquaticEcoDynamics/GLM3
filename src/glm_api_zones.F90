!###############################################################################
!#                                                                             #
!# glm_api_zones.F90                                                           #
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

#undef MISVAL

#include "glm.h"

#define _UNUSED(x) if (.FALSE.) print*,shape(x)

MODULE glm_api_zones

   USE ISO_C_BINDING

   USE glm_types
   USE aed_common
   USE aed_zones

   IMPLICIT NONE

   PRIVATE ! By default, make everything private

   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: z_cc      ! (n_vars, n_levs, nsed_zones)
   AED_REAL,DIMENSION(:,:),  ALLOCATABLE,TARGET :: z_cc_hz   ! (n_vars, nsed_zones)
   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: z_diag    ! (n_vars, n_levs, nsed_zones)
   AED_REAL,DIMENSION(:,:),  ALLOCATABLE,TARGET :: z_diag_hz ! (n_vars, nsed_zones)

   AED_REAL,DIMENSION(:),POINTER :: lheights

   AED_REAL,DIMENSION(:),POINTER :: zone_heights
   INTEGER :: n_aed_vars
   INTEGER :: nvars, nbenv, nvdiag, nvdiag_hz

   PUBLIC n_zones, zone_heights, lheights
   PUBLIC api_set_glm_zones, api_copy_from_zone, api_copy_to_zone, api_calc_zone_areas

   PUBLIC z_cc, z_cc_hz, z_diag, z_diag_hz, theZones, theLake

CONTAINS


!###############################################################################
SUBROUTINE api_set_glm_zones(numVars, numBenV, numDiagV, numDiagHzV, nAEDvars) &
                                               BIND(C, name="api_set_glm_zones")
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: numVars, numBenV, numDiagV, numDiagHzV, nAEDvars
!
!LOCALS
   INTEGER :: zon

   PROCEDURE(copy_to_zone_t),POINTER    :: copy_to
   PROCEDURE(copy_from_zone_t),POINTER  :: copy_from
   PROCEDURE(calc_zone_areas_t),POINTER :: calc_areas
!
!-------------------------------------------------------------------------------
!BEGIN
   nvars = numVars
   nbenv = numBenV
   nvdiag = numDiagV
   nvdiag_hz = numDiagHzV
   n_aed_vars = nAEDvars

   lheights => theLake%Height

   ALLOCATE(z_cc(numVars, MaxLayers, n_zones))    ; z_cc = 0.
   ALLOCATE(z_cc_hz(numBenV, n_zones))            ; z_cc_hz = 0.
   ALLOCATE(z_diag(numDiagV, MaxLayers, n_zones)) ; z_diag = 0.
   ALLOCATE(z_diag_hz(numDiagHzV, n_zones+1))     ; z_diag_hz = 0.

   CALL aed_init_zones(n_zones, 1, z_cc, z_cc_hz, z_diag, z_diag_hz)

   copy_to => api_copy_to_zone
   copy_from => api_copy_from_zone
   calc_areas => api_calc_zone_areas
   CALL api_set_zone_funcs(copy_to, copy_from, calc_areas)

   DO zon=1,n_zones
      aedZones(zon)%z_env%z_area = 0.
      aedZones(zon)%z_env%z_height = theZones(zon)%zheight
   ENDDO
   theZones%zarea = 0.
   zone_heights => theZones%zheight
END SUBROUTINE api_set_glm_zones
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_calc_zone_areas(aedZones, n_zones, areas, wheights, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_zone_t),DIMENSION(:),INTENT(inout) :: aedZones
   INTEGER,INTENT(in) :: n_zones
!
   AED_REAL,DIMENSION(:),POINTER,INTENT(in) :: areas
   AED_REAL,DIMENSION(:),POINTER,INTENT(in) :: wheights
   INTEGER,INTENT(in) :: wlev
!
!LOCALS
   INTEGER  :: lev, zon
   LOGICAL  :: w_zones
   AED_REAL :: surf, scale
!
!-------------------------------------------------------------------------------
!BEGIN
   DO zon=1, n_zones
      aedZones(zon)%z_env%z_area = 0.  ; aedZones(zon)%z_env%z_pc_wet = 0.
   ENDDO
   w_zones = .FALSE.
   surf = wheights(wlev)

   zon = 1
   aedZones(1)%z_env%z_area = areas(1)
   IF ( aedZones(1)%z_env%z_height > surf ) THEN
      aedZones(1)%z_env%z_pc_wet = surf / aedZones(1)%z_env%z_height
   ELSE
      aedZones(1)%z_env%z_pc_wet = 1.0
   ENDIF
   DO lev=2, wlev
      ! areas() here is the INCREMENTAL layer area, so summing it gives the zone band area
      IF (wheights(lev) <= aedZones(zon)%z_env%z_height) THEN
         aedZones(zon)%z_env%z_area = aedZones(zon)%z_env%z_area + areas(lev)
      ELSEIF (wheights(lev) > aedZones(zon)%z_env%z_height .AND. wheights(lev-1) < aedZones(zon)%z_env%z_height) THEN
         scale = (aedZones(zon)%z_env%z_height - wheights(lev-1)) / (wheights(lev) - wheights(lev-1))
         aedZones(zon)%z_env%z_area = aedZones(zon)%z_env%z_area + areas(lev) * scale
         zon = zon + 1
         aedZones(zon)%z_env%z_area = aedZones(zon)%z_env%z_area + areas(lev) * (1-scale)
      ENDIF

      IF ( aedZones(zon)%z_env%z_height > surf ) THEN
         IF (.NOT. w_zones) THEN
            w_zones = .TRUE.
            IF ( zon > 1 ) THEN
               aedZones(zon)%z_env%z_pc_wet = surf / (aedZones(zon)%z_env%z_height - aedZones(zon-1)%z_env%z_height)
            ELSE
               aedZones(zon)%z_env%z_pc_wet = surf / aedZones(zon)%z_env%z_height
            ENDIF
         ENDIF
      ELSE
         aedZones(zon)%z_env%z_pc_wet = 1.0
      ENDIF
   ENDDO

   DO zon=1, n_zones
      aedZones(zon)%z_env%z_pres = -aedZones(zon)%z_env%z_height
      aedZones(zon)%z_env%z_col_area = aedZones(zon)%z_env%z_area
      theZones(zon)%zarea = aedZones(zon)%z_env%z_area     ! for netCDF 'zarea' output
   ENDDO
END SUBROUTINE api_calc_zone_areas
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_copy_to_zone(aedZones, n_zones, wheights, x_cc, x_cc_hz, x_diag, x_diag_hz, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_zone_t),DIMENSION(:),INTENT(inout) :: aedZones
   INTEGER,INTENT(in) :: n_zones
!
   AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: wheights
   AED_REAL,DIMENSION(:,:),POINTER,INTENT(in) :: x_cc
   AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: x_cc_hz
   AED_REAL,DIMENSION(:,:),POINTER,INTENT(in) :: x_diag
   AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: x_diag_hz
!
   INTEGER,INTENT(in) :: wlev
!
!LOCALS
   INTEGER  :: zon, lev, a_zones
   AED_REAL :: surf
   INTEGER  :: zcount(n_zones)
   LOGICAL  :: w_zones
   INTEGER  :: i, sd
   TYPE(aed_variable_t),POINTER :: tv
!
!-------------------------------------------------------------------------------
!BEGIN
   _UNUSED(x_cc_hz)
   _UNUSED(x_diag)
   _UNUSED(x_diag_hz)

   DO zon=1,n_zones
      z_cc(1:nvars,:,zon) = 0.
      z_diag(:,:,zon) = 0.
      !# z_diag_hz is the HOME of the per-zone benthic diagnostics (bound to
      !# aedZones(zon)%z_cc_diag_hz), so clear it PER-VARIABLE: diagnostics
      !# flagged rezero=.FALSE. (e.g. the CGM running averages cgm_tavg/lavg/savg
      !# and slough counters) must persist across timesteps; the rest are cleared.
      sd = 0
      DO i=1,n_aed_vars
         IF ( aed_get_var(i, tv) ) THEN
            IF ( tv%var_type == V_DIAGNOSTIC .AND. tv%sheet ) THEN
               sd = sd + 1
               IF ( tv%rezero ) z_diag_hz(sd,zon) = 0.
            ENDIF
         ENDIF
      ENDDO

      aedZones(zon)%z_env%z_temp = 0.
      aedZones(zon)%z_env%z_salt = 0.
      aedZones(zon)%z_env%z_rho = 0.
      aedZones(zon)%z_env%z_rad = 0.
      aedZones(zon)%z_env%z_extc = 0.
      aedZones(zon)%z_env%z_layer_stress = 0.
      aedZones(zon)%z_env%z_tss = 0.
      aedZones(zon)%z_env%z_par = 0.
      aedZones(zon)%z_env%z_nir = 0.
      aedZones(zon)%z_env%z_uva = 0.
      aedZones(zon)%z_env%z_uvb = 0.
      aedZones(zon)%z_env%z_sed_zones = zon
      aedZones(zon)%z_env%z_vel = 0.
   ENDDO

   a_zones = 1
   zcount = 0
   w_zones = .FALSE.
   zon = 1
   DO lev=1,wlev
      IF ( lev > 1 .AND. wheights(lev) > aedZones(zon)%z_env%z_height ) THEN
         zon = zon + 1
         IF (zon > n_zones) STOP 'Water level height is higher than highest zone height'
!        aedZones(zon)%z_env(1)%z_sed_zones = zon
!        z_diag_hz(:,zon) = z_diag_hz(:,zon) + x_diag_hz(:)
      ENDIF

      ! Pelagic variables need zonifying, in case benthic people want them
      ! Note that this should not be done on benthic vars because it will
      ! introduce errors in z_cc (split layers)
      ! Ideally this average would be based on volume weighting

      !z_cc(1:nvars,lev,zon) = z_cc(1:nvars,lev,zon) + x_cc(1:nvars,lev)
!RQT
      z_cc(1:nvars,1,zon) = z_cc(1:nvars,1,zon) + x_cc(1:nvars,lev)  ! level=1 to match aedZones(zon)%z_cc => z_cc(:,1,:) pointer

     !z_diag(:,zon,zon) = z_diag(:,zon,zon) + x_diag(:,lev)      ! level=zon so cell(zon) correct with layer_idx=zon
      z_diag(:,1,zon) = z_diag(:,1,zon) + x_diag(:,lev)          ! level=1 matches agnostic library wiring z_cc_diag => z_diag(:,1,:) (cell(zon) -> z_diag(:,1,zon))

      aedZones(zon)%z_env%z_temp         = aedZones(zon)%z_env%z_temp + theLake(lev)%Temp
      aedZones(zon)%z_env%z_salt         = aedZones(zon)%z_env%z_salt + theLake(lev)%Salinity
      aedZones(zon)%z_env%z_rho          = aedZones(zon)%z_env%z_rho  + theLake(lev)%Density
      aedZones(zon)%z_env%z_rad          = aedZones(zon)%z_env%z_rad  + theLake(lev)%Light
      aedZones(zon)%z_env%z_par          = aedZones(zon)%z_env%z_par  + theLake(lev)%Light*0.45
      aedZones(zon)%z_env%z_vel          = aedZones(zon)%z_env%z_vel  + theLake(lev)%Umean
      aedZones(zon)%z_env%z_extc         = aedZones(zon)%z_env%z_extc + theLake(lev)%ExtcCoefSW
      aedZones(zon)%z_env%z_layer_stress = aedZones(zon)%z_env%z_layer_stress + theLake(lev)%LayerStress

      zcount(zon) = zcount(zon) + 1
   ENDDO
   a_zones = zon

   DO zon=1,a_zones
      z_cc(1:nvars,1,zon) = z_cc(1:nvars,1,zon)/zcount(zon)   ! level=1 to match pointer
     !z_diag(:,zon,zon)   = z_diag(:,zon,zon)/zcount(zon)    ! level=zon matches accumulation and layer_idx=zon
      z_diag(:,1,zon)     = z_diag(:,1,zon)/zcount(zon)      ! level=1 matches agnostic library wiring z_cc_diag => z_diag(:,1,:)
!     z_diag_hz(:,zon)    = z_diag_hz(:,zon)/zcount(zon)

      ! Set the water column above a zone, to the respective water layer values
!     z_cc(1:nvars,2:wlev,zon) = x_cc(1:nvars,2:wlev) ! water column state vars
!     z_diag(:,2:wlev,zon)     = x_diag(:,2:wlev)     ! water column diag vars
   ENDDO

   DO zon=1,a_zones
      IF (zcount(zon) /= 0) THEN
         aedZones(zon)%z_env%z_temp         = aedZones(zon)%z_env%z_temp / zcount(zon)
         aedZones(zon)%z_env%z_salt         = aedZones(zon)%z_env%z_salt / zcount(zon)
         aedZones(zon)%z_env%z_rho          = aedZones(zon)%z_env%z_rho  / zcount(zon)
         aedZones(zon)%z_env%z_rad          = aedZones(zon)%z_env%z_rad  / zcount(zon)
         aedZones(zon)%z_env%z_par          = aedZones(zon)%z_env%z_par  / zcount(zon)
         aedZones(zon)%z_env%z_vel          = aedZones(zon)%z_env%z_vel  / zcount(zon)
         aedZones(zon)%z_env%z_extc         = aedZones(zon)%z_env%z_extc / zcount(zon)
         aedZones(zon)%z_env%z_layer_stress = aedZones(zon)%z_env%z_layer_stress / zcount(zon)
      ELSE
         aedZones(zon)%z_env%z_temp         = 0.
         aedZones(zon)%z_env%z_salt         = 0.
         aedZones(zon)%z_env%z_rho          = 0.
         aedZones(zon)%z_env%z_rad          = 0.
         aedZones(zon)%z_env%z_par          = 0.
         aedZones(zon)%z_env%z_vel          = 0.
         aedZones(zon)%z_env%z_extc         = 0.
         aedZones(zon)%z_env%z_layer_stress = 0.
      ENDIF

      aedZones(zon)%z_env%z_dz = 0.
   ENDDO

   surf = wheights(wlev)

   ! Populate lake-wide sheet globals into each zone's env struct.
   ! In glm_aed.F90 these are wired as direct pointer-to-pointer; here the zone
   ! struct holds scalar fields so we copy the live values each timestep.
   DO zon=1,n_zones
      aedZones(zon)%z_env%z_col_depth = surf
      aedZones(zon)%z_env%z_air_pres  = MetData%AirPres
      aedZones(zon)%z_env%z_air_temp  = MetData%AirTemp
      aedZones(zon)%z_env%z_wind      = MetData%WindSpeed
      aedZones(zon)%z_env%z_humidity  = MetData%RelHum
      aedZones(zon)%z_env%z_rain      = MetData%Rain
      aedZones(zon)%z_env%z_evap      = SurfData%Evap
      aedZones(zon)%z_env%z_I_0       = MetData%ShortWave
      aedZones(zon)%z_env%z_longwave  = MetData%LongWave
   ENDDO

   IF ( surf > aedZones(1)%z_env%z_height ) THEN
      aedZones(1)%z_env%z_depth = aedZones(1)%z_env%z_height
      aedZones(1)%z_env%z_dz = aedZones(1)%z_env%z_height
   ELSE
      aedZones(1)%z_env%z_depth = surf
      aedZones(1)%z_env%z_dz = surf
      w_zones = .TRUE.
   ENDIF

   DO zon=2,a_zones
      aedZones(zon)%z_env%z_sed_zones = zon
      IF ( .NOT. w_zones ) THEN
          IF ( surf > aedZones(zon)%z_env%z_height ) THEN
             aedZones(zon)%z_env%z_depth = aedZones(zon)%z_env%z_height
             aedZones(zon)%z_env%z_dz =  aedZones(zon)%z_env%z_height - aedZones(zon-1)%z_env%z_height
          ELSE
             aedZones(zon)%z_env%z_depth = surf
             aedZones(zon)%z_env%z_dz = surf - aedZones(zon-1)%z_env%z_height
             w_zones = .TRUE.
          ENDIF
      ELSE
         aedZones(zon)%z_env%z_depth = surf
         aedZones(zon)%z_env%z_dz = surf - aedZones(zon-1)%z_env%z_height
      ENDIF
   ENDDO

END SUBROUTINE api_copy_to_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_copy_from_zone(aedZones, n_zones, wheights, x_cc, x_cc_hz, x_diag, x_diag_hz, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_zone_t),DIMENSION(:),INTENT(in) :: aedZones
   INTEGER,INTENT(in) :: n_zones
!
   AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: wheights
   AED_REAL,DIMENSION(:,:),POINTER,INTENT(inout) :: x_cc
   AED_REAL,DIMENSION(:),  POINTER,INTENT(inout) :: x_cc_hz
   AED_REAL,DIMENSION(:,:),POINTER,INTENT(inout) :: x_diag
   AED_REAL,DIMENSION(:),  POINTER,INTENT(inout) :: x_diag_hz
!
   INTEGER,INTENT(in) :: wlev
!
!LOCALS
   INTEGER  :: zon, lev, i, j
   AED_REAL :: scale, area
   LOGICAL  :: splitZone
   TYPE(aed_variable_t),POINTER :: tvar
!
!-------------------------------------------------------------------------------
!BEGIN
   _UNUSED(x_cc)
   zon = n_zones

   ! Loop down through water layers
   DO lev=wlev,1,-1
      ! Check if zone boundary is in this water layer range
      IF ( zon .GT. 1 ) THEN
         IF (lev .GT. 1) THEN
            splitZone = wheights(lev-1) < aedZones(zon-1)%z_env%z_height
         ELSE
            splitZone = 0.0 < aedZones(zon-1)%z_env%z_height
         ENDIF
      ELSE
         splitZone = .FALSE.
      ENDIF

      ! Set water layer variables, based on zone infomration, where variable is flagged for zavg
      IF (splitZone) THEN
         ! Compute layer fraction
         IF (lev .GT. 1) THEN
            scale = (aedZones(zon-1)%z_env%z_height - wheights(lev-1)) / (wheights(lev) - wheights(lev-1))
         ELSE
            scale = (aedZones(zon-1)%z_env%z_height - 0.0) / (wheights(lev) - 0.0)
         ENDIF

         ! Select the diag vars that have zavg == true, and assign to layer
         ! (for first zone in the layer)
         j = 0
         DO i=1,n_aed_vars
            IF ( aed_get_var(i, tvar) ) THEN
               IF ( tvar%var_type == V_DIAGNOSTIC ) THEN
                  IF ( .NOT.  tvar%sheet ) THEN
                     j = j + 1
                     IF ( tvar%zavg ) THEN
                       !x_diag(j,lev) = z_diag(j,zon,zon) * scale
                        x_diag(j,lev) = z_diag(j,1,zon) * scale   ! level=1 matches library z_cc_diag => z_diag(:,1,:)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         ! Select the benthic vars (implicitly all zavg == true), and assign to layer
         ! (for first zone in the layer)
         x_cc_hz(:) = z_cc_hz(:,zon) * scale

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
                       !x_diag(j,lev) = x_diag(j,lev) + (z_diag(j,zon,zon) * (1.0 - scale))
                        x_diag(j,lev) = x_diag(j,lev) + (z_diag(j,1,zon) * (1.0 - scale))   ! level=1 matches library z_cc_diag => z_diag(:,1,:)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         ! Select the benthic vars (implicitly all zavg == true), and assign to layer
         ! (for first zone in the layer)
         x_cc_hz(:) = x_cc_hz(:) + z_cc_hz(:,zon) * (1.0 - scale)
      ELSE  ! Not a split layer - layer bounds are fully within zone
         ! Select the diag vars that have zavg == true, and assign to layer
         j = 0
         DO i=1,n_aed_vars
            IF ( aed_get_var(i, tvar) ) THEN
               IF ( tvar%var_type == V_DIAGNOSTIC ) THEN
                  IF ( .NOT.  tvar%sheet ) THEN
                     j = j + 1
                     IF ( tvar%zavg ) THEN
                       !x_diag(j,lev) = z_diag(j,zon,zon)
                        x_diag(j,lev) = z_diag(j,1,zon)   ! level=1 matches library z_cc_diag => z_diag(:,1,:)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         ! Select the benthic vars (implicitly all zavg == true), and assign to layer
         x_cc_hz(:) = z_cc_hz(:,zon)
      ENDIF
   ENDDO

   ! Set the normal sheet diagnostics to the mean of the zone, weighted by area
   area = 0.
   DO zon=1,n_zones
      area = area + aedZones(zon)%z_env%z_area
   ENDDO
   DO zon=1,n_zones
      x_diag_hz = x_diag_hz + (z_diag_hz(:,zon) * (aedZones(zon)%z_env%z_area/area))
   ENDDO

END SUBROUTINE api_copy_from_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE glm_api_zones
