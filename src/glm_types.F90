!###############################################################################
!#                                                                             #
!#  glm_types.F90                                                              #
!#                                                                             #
!#  A module to define constants and types.                                    #
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
#include "glm.h"


!************************* Important Note **************************************
!* The order of entries in these structures MUST match those in glm_types.h    *
!************************* Important Note **************************************

!*******************************************************************************
MODULE glm_types

   USE ISO_C_BINDING

   IMPLICIT NONE

!===============================================================================
!GLOBAL CONSTANTS

   AED_REAL,PARAMETER :: missing = MISVAL

!===============================================================================
!TYPE DECLARATIONS

   !#===========================================================#!
   TYPE,BIND(C) :: StringT
      CINTEGER :: Len
      CCHARACTER :: S(40)
   END TYPE StringT

   !#===========================================================#!
   !# Structured type for key global lake environmental vars
   !# A Lake will be an allocated array of MaxLayers of these
   TYPE,BIND(C) :: LakeDataType
      AED_REAL :: Density          !# density kg/m3
      AED_REAL :: Temp             !# temperature
      AED_REAL :: Salinity         !# salinity
      AED_REAL :: Height           !# layer heights above the bottom
      AED_REAL :: MeanHeight       !# mean height of a layer
      AED_REAL :: LayerVol         !# volume of layer
      AED_REAL :: LayerArea        !# area of layer

      AED_REAL :: Light            !# PAR, photosynthetically active radiation
      AED_REAL :: ExtcCoefSW       !# Kd, light extinction coefficient

      AED_REAL :: Vol1             !# Cumulative volume to this layer top
      AED_REAL :: Epsilon          !# Diffusivity

      AED_REAL :: Umean            !# Mean velocity
      AED_REAL :: Uorb             !# Orbital velocity
      AED_REAL :: LayerStress      !# Layer Stress
   END TYPE LakeDataType

   !#===========================================================#!
   !# Structured type for Met vars
   TYPE,BIND(C) :: MetDataType
      AED_REAL :: Rain             !# rainfall
      AED_REAL :: RelHum           !# relative humidty
      AED_REAL :: SatVapDef        !# vapour pressure
      AED_REAL :: LongWave         !# longwave radiation
      AED_REAL :: ShortWave        !# shortwave radiation
      AED_REAL :: AirTemp          !# temperature
      AED_REAL :: AirPres          !# air pressure
      AED_REAL :: WindSpeed        !# windspeed
      AED_REAL :: Snow             !# snowfall
      AED_REAL :: RainConcPO4      !# Concentration of PO4 in rain
      AED_REAL :: RainConcTP       !# Concentration of TP in rain
      AED_REAL :: RainConcNO3      !# Concentration of NO3 in rain
      AED_REAL :: RainConcNH4      !# Concentration of NH4 in rain
      AED_REAL :: RainConcTN       !# Concentration of TN in rain
      AED_REAL :: RainConcSi       !# Concentration of SI in rain
      AED_REAL :: WindDir          !# Wind direction
      AED_REAL :: As               !# Area of sheltering
   END TYPE MetDataType

   !#===========================================================#!
   !# Structured type for Surface Data vars
   TYPE,BIND(C) :: SurfaceDataType
      AED_REAL :: Evap             !# Evaporation
      AED_REAL :: delzBlueIce      !# Thickness of blue ice layer
      AED_REAL :: delzWhiteIce     !# Thickness of white ice layer
      AED_REAL :: delzSnow         !# Thickness of snow layer
      AED_REAL :: dHt              !# Change in thickness of snow / ice layer
      AED_REAL :: RhoSnow          !# Density of snow layer (kg/m^3)
      AED_REAL :: dailyEvap        !# Daily Evaporation (m3/day)
      AED_REAL :: dailyRain        !# Daily Rain (m3/day)
      AED_REAL :: dailyRunoff      !# Daily Rain (m3/day)
      AED_REAL :: dailySnow        !# Daily Snow (m3/day)
      AED_REAL :: Qsw              !# Shortwave radiation (W/m2)
      AED_REAL :: Qe               !# Latent heat flux (W/m2)
      AED_REAL :: Qh               !# Sensible heat flux (W/m2)
      AED_REAL :: Qlw              !# Net longwave radiation (W/m2)
      AED_REAL :: dailyQsw         !# Daily Short Wave Radiation (J/day)
      AED_REAL :: dailyQe          !# Daily Latent Heat(J/day)
      AED_REAL :: dailyQh          !# Daily Sensible Heat (J/day)
      AED_REAL :: dailyQlw         !# Daily Long Wave Radiation (J/day)
      AED_REAL :: dailyInflow      !# Total Daily Inflow (m3/day)
      AED_REAL :: dailyOutflow     !# Total Daily Outflow (m3/day)
      AED_REAL :: dailyOverflow    !# Total Daily Overflow (m3/day)
      AED_REAL :: dailySeepage     !# Total Daily Seepage (m3/day)
      AED_REAL :: albedo           !# Daily surface albedo
      AED_REAL :: dailyzonL        !# Average z/L value, daily atmospheric stability
      AED_REAL :: u_star           !# Wind friction velocity (m/s)
      AED_REAL :: Q_net            !# Net non-penetrative heat flux (W/m2)
   END TYPE SurfaceDataType

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !#===========================================================#!
   !# Structured type for Sediment Layer
   TYPE,BIND(C) :: SedLayerType
      AED_REAL :: depth            !# Layer depth
      AED_REAL :: temp             !# Layer temperature
      AED_REAL :: vwc
      AED_REAL :: wq
   END TYPE SedLayerType

   !#===========================================================#!
   !# Structured type for iSediment Zones
   TYPE,BIND(C) :: ZoneType
      AED_REAL :: zheight
      AED_REAL :: zrad
      AED_REAL :: zsalt
      AED_REAL :: ztemp
      AED_REAL :: zrho
      AED_REAL :: zarea
      AED_REAL :: zextc
      AED_REAL :: zlayer_stress
      AED_REAL :: ztss
      AED_REAL :: zdz
      AED_REAL :: zvel
      AED_REAL :: zpar
      AED_REAL :: znir
      AED_REAL :: zuva
      AED_REAL :: zuvb
      AED_REAL :: zpres
      AED_REAL :: zdepth
      AED_REAL :: z_sed_zones
      AED_REAL :: z_pc_wet
      AED_REAL :: heatflux
      CINTEGER :: n_sed_layers;    !# number of sediment layers
      TYPE(C_PTR) :: c_layers      !# array of sed layers
   END TYPE ZoneType

#if 0

!# These are not used in the fortran part - will need to check against .h if they are to be used
   !#===========================================================#!
   !# Structured type for Particle Transport Model (PTM)
   TYPE,BIND(C) :: ParticleDataType
       INTEGER  :: Status         ! indivdual particle status
       INTEGER  :: Flag           ! indivdual particle flag indicating if BED (1) or SCUM (2), or neither (0)
       AED_REAL :: Height
       AED_REAL :: Mass
       AED_REAL :: Diam
       AED_REAL :: Density
       AED_REAL :: Velocity
       AED_REAL :: vvel
       CINTEGER :: Layer
   END TYPE ParticleDataType

   !#===========================================================#!
   !# NEW structured type for Particle Transport Model (PTM), following AED API
   TYPE,BIND(C) :: partgroup
      INTEGER(KIND=4) :: NP                                ! Number of Particles
      INTEGER(KIND=4) :: id_stat, id_i2, id_i3, id_layer   ! Particle ISTAT Index Values
      INTEGER(KIND=4) :: id_bed_layer, id_motility         ! Particle ISTAT Index Values
      INTEGER(KIND=4) :: id_uvw0, id_uvw, id_nu, id_wnd    ! Particle PROP Index Values
      INTEGER(KIND=4) :: id_wsel, id_watd, id_partd        ! Particle PROP Index Values
      INTEGER(KIND=4) :: id_age, id_state                  ! Particle TSTAT Index Values
      INTEGER(KIND=4) :: i_next                            ! next particle index
      INTEGER(KIND=4),POINTER,DIMENSION(:,:) :: istat      ! Particle Integer Status/Cell-index variables (4,NPart)
      REAL(KIND=8),POINTER,DIMENSION(:,:) :: tstat         ! Particle Time/Age Vector (2,Npart)
      REAL(KIND=8),POINTER,DIMENSION(:,:) :: xyz           ! particle position vector
      REAL(KIND=4),POINTER,DIMENSION(:,:) :: prop          ! Particle Property Vector (12,Npart)
      REAL(KIND=4),POINTER,DIMENSION(:,:) :: U             ! Particle Conserved Variable Vector (NU,NP)
   END TYPE partgroup

   TYPE,BIND(C) :: partgroup_p
      INTEGER :: idx, grp
   END TYPE partgroup_p

   TYPE,BIND(C) :: partgroup_cell
      INTEGER :: count, n
      TYPE(partgroup_p),ALLOCATABLE,DIMENSION(:) :: prt
   END TYPE partgroup_cell
#endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !================================================================
   !# variables in C code of GLM
   !----------------------------------------------------------------
   TYPE(CINTEGER),BIND(C, name="MaxLayers")   :: MaxLayers
   TYPE(C_PTR),BIND(C, name="Lake")           :: cLake
   TYPE(CINTEGER),BIND(C, name="n_zones")     :: n_zones
   TYPE(C_PTR),BIND(C, name="theZones")       :: cZones

   TYPE(C_PTR),BIND(C, name="pMetData")       :: cMetData
   TYPE(C_PTR),BIND(C, name="pSurfData")      :: cSurfData

   TYPE(LakeDataType),   DIMENSION(:),POINTER :: theLake
   TYPE(MetDataType),                 POINTER :: MetData   !# Meteorological data
   TYPE(SurfaceDataType),             POINTER :: SurfData  !# Surface Data
   TYPE(ZoneType),       DIMENSION(:),POINTER :: theZones

   TYPE(FLOGICAL),BIND(C, name="mobility_off")       :: mobility_off      = .FALSE.
   TYPE(FLOGICAL),BIND(C, name="bioshade_feedback")  :: bioshade_feedback = .TRUE.
   TYPE(FLOGICAL),BIND(C, name="repair_state")       :: repair_state      = .TRUE.
   TYPE(FLOGICAL),BIND(C, name="do_plots")           :: do_plots          = .FALSE.
   TYPE(FLOGICAL),BIND(C, name="link_rain_loss")     :: link_rain_loss    = .FALSE.
   TYPE(FLOGICAL),BIND(C, name="link_solar_shade")   :: link_solar_shade  = .FALSE.
   TYPE(FLOGICAL),BIND(C, name="link_bottom_drag")   :: link_bottom_drag  = .FALSE.
   TYPE(FLOGICAL),BIND(C, name="ice")                :: ice               = .FALSE.

   TYPE(CINTEGER),BIND(C, name="split_factor")       :: split_factor = 1
   TYPE(CINTEGER),BIND(C, name="ode_method")         :: ode_method
   TYPE(CINTEGER),BIND(C, name="benthic_mode")       :: benthic_mode

   TYPE(AED_REAL),TARGET,BIND(C, name="rain_factor") :: rain_factor
   TYPE(AED_REAL),TARGET,BIND(C, name="sw_factor")   :: sw_factor
   TYPE(AED_REAL),TARGET,BIND(C, name="friction")    :: friction

   TYPE(AED_REAL),TARGET,BIND(C, name="Kw")          :: Kw
   TYPE(AED_REAL),TARGET,BIND(C, name="dt")          :: dt = 0.

   TYPE(AED_REAL),TARGET,BIND(C, name="yearday")     :: yearday   = 0.
   TYPE(AED_REAL),TARGET,BIND(C, name="timestep")    :: timestep  = 0.
   TYPE(AED_REAL),TARGET,BIND(C, name="Longitude")   :: longitude = 0.
   TYPE(AED_REAL),TARGET,BIND(C, name="Latitude")    :: latitude  = 0.

   TYPE(FLOGICAL),BIND(C, name="do_particle_bgc")    :: do_particle_bgc = .FALSE.
   TYPE(FLOGICAL),BIND(C, name="link_ext_par")       :: link_ext_par = .FALSE.

   !# Dynamic sediment soil-temperature model (sed_heat_model = 2): the per-step
   !# timestep (set from noSecs by glm_surface.c) and the soil thermal properties
   !# (mirror intertidal-soil SoilParams; optional &sediment overrides).
   !# Fortran owns these symbols; the C side extern's them (glm_globals.h).
   TYPE(AED_REAL),TARGET,BIND(C, name="soil_dt")             :: soil_dt = 3600.
   TYPE(AED_REAL),TARGET,BIND(C, name="sed_k_mineral")       :: sed_k_mineral = 2.5
   TYPE(AED_REAL),TARGET,BIND(C, name="sed_k_water")         :: sed_k_water = 0.57
   TYPE(AED_REAL),TARGET,BIND(C, name="sed_k_air")           :: sed_k_air = 0.025
   TYPE(AED_REAL),TARGET,BIND(C, name="sed_c_mineral")       :: sed_c_mineral = 2.0e6
   TYPE(AED_REAL),TARGET,BIND(C, name="sed_c_water")         :: sed_c_water = 4.18e6
   TYPE(AED_REAL),TARGET,BIND(C, name="sed_c_air")           :: sed_c_air = 1.25e3
   TYPE(AED_REAL),TARGET,BIND(C, name="sed_bulk_density")    :: sed_bulk_density = 1.5
   TYPE(AED_REAL),TARGET,BIND(C, name="sed_mineral_density") :: sed_mineral_density = 2.6
   TYPE(AED_REAL),TARGET,BIND(C, name="sed_porosity")        :: sed_porosity = -1.0
   TYPE(AED_REAL),TARGET,BIND(C, name="sed_deep_temp")       :: sed_deep_temp = -9999.0

CONTAINS

!#
!# These are 2 useful routines for converting between fortran and C strings
!# They are in here because, well, I guess they are sort of type conversions
!#

!###############################################################################
SUBROUTINE make_string(s1,s2,len)
   CHARACTER(len=*),INTENT(out) :: s1
   CCHARACTER,INTENT(in) :: s2(*)
   CSIZET,INTENT(in)    :: len
!LOCALS
   CHARACTER(len=len) :: s3
!
!-------------------------------------------------------------------------------
!BEGIN
   s1 = trim(transfer(s2(1:len),s3))
END SUBROUTINE make_string
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION make_c_string(s1,s2) RESULT(len)
   CCHARACTER,INTENT(out) :: s1(*)
   CHARACTER(len=*),INTENT(in) :: s2
!LOCALS
   INTEGER :: i
   INTEGER :: len
!
!-------------------------------------------------------------------------------
!BEGIN
   len = len_trim(s2)
   DO i=1,len
      s1(i) = s2(i:i)
   ENDDO
   s1(len+1) = ACHAR(0)
END FUNCTION make_c_string
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE glm_init_fortran_support() BIND(C, name="glm_init_fortran_support")
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL C_F_POINTER(cLake, theLake, [MaxLayers]);
   CALL C_F_POINTER(cMetData, MetData)
   CALL C_F_POINTER(cSurfData, SurfData)
   CALL C_F_POINTER(cZones, theZones, [n_zones]);
END SUBROUTINE glm_init_fortran_support
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE glm_types
