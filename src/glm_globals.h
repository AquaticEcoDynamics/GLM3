/******************************************************************************
 *                                                                            *
 * glm_globals.h                                                              *
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
#ifndef _GLM_GLOBALS_H_
#define _GLM_GLOBALS_H_

#ifndef __STDC__
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 INTERFACE

     SUBROUTINE set_c_wqvars_ptr(iwqvars) BIND(C, name="set_c_wqvars_ptr")
        USE ISO_C_BINDING
#       if defined( _WIN32 ) && USE_DL_LOADER
        !DEC$ ATTRIBUTES DLLIMPORT :: set_c_wqvars_ptr
#       endif
        AED_REAL,INTENT(in) :: iwqvars(*)
     END SUBROUTINE set_c_wqvars_ptr

     SUBROUTINE set_c_wqsvars_ptr(iwqsvars) BIND(C, name="set_c_wqsvars_ptr")
        USE ISO_C_BINDING
        AED_REAL,INTENT(in) :: iwqsvars(*)
     END SUBROUTINE set_c_wqsvars_ptr

     SUBROUTINE set_c_wqdvars_ptr(iwqdvars,iwqdsvars,nwqd,nwqds) BIND(C, name="set_c_wqdvars_ptr")
        USE ISO_C_BINDING
#       if defined( _WIN32 ) && USE_DL_LOADER
        !DEC$ ATTRIBUTES DLLIMPORT :: set_c_wqvars_ptr
#       endif
        AED_REAL,INTENT(in) :: iwqdvars(*)
        AED_REAL,INTENT(in) :: iwqdsvars(*)
        CINTEGER,INTENT(in) :: nwqd
        CINTEGER,INTENT(in) :: nwqds
     END SUBROUTINE set_c_wqdvars_ptr

     SUBROUTINE set_c_ptmstat_ptr(iptms) BIND(C, name="set_c_ptmstat_ptr")
        USE ISO_C_BINDING
#       if defined( _WIN32 ) && USE_DL_LOADER
        !DEC$ ATTRIBUTES DLLIMPORT :: set_c_ptmstat_ptr
#       endif
        CINTEGER,INTENT(in) :: iptms(*)
     END SUBROUTINE set_c_ptmstat_ptr

     SUBROUTINE set_c_ptmenv_ptr(iptmv) BIND(C, name="set_c_ptmenv_ptr")
        USE ISO_C_BINDING
#       if defined( _WIN32 ) && USE_DL_LOADER
        !DEC$ ATTRIBUTES DLLIMPORT :: set_c_ptmenv_ptr
#       endif
        AED_REAL,INTENT(in) :: iptmv(*)
     END SUBROUTINE set_c_ptmenv_ptr

     SUBROUTINE set_c_num_ptm_vars(n) BIND(C, name="set_c_num_ptm_vars")
        USE ISO_C_BINDING
        CINTEGER,INTENT(in) :: n
     END SUBROUTINE set_c_num_ptm_vars

# if DEBUG
     SUBROUTINE debug_print_lake() BIND(C, name="debug_print_lake")
     END SUBROUTINE debug_print_lake

     SUBROUTINE debug_initialisation(which) BIND(C, name="debug_initialisation_")
        USE ISO_C_BINDING
        CINTEGER,INTENT(in) :: which
     END SUBROUTINE debug_initialisation
# endif

 END INTERFACE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#else
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "glm.h"
#include "glm_types.h"

/* from glm_ncdf */
extern int ncid;

/* from glm_surf.F90 */
extern AED_REAL AvgSurfTemp;
/*----------------------------------------------------------------------------*/
extern CINTEGER NumLayers;   //# current number of layers

extern AED_REAL DMin;    //# minimum layer thickness
extern AED_REAL DMax;    //# maximum layer thickness
extern AED_REAL VMin;    //# minimum layer volume
extern AED_REAL VMax;    //# maximum layer volume

extern int wq_calc;         //# are we doing water quality calcs
extern int Num_WQ_Vars;     //# number of water quality variables
extern int Num_WQ_Ben;      //# number of benthic water quality variables
extern int Tot_WQ_Vars;     //# nVars + nBen
extern int Num_WQD_Vars;    //# number of water quality diagnostic variables
extern int Num_WQDS_Vars;   //# number of water quality diagnostic "sheet" variables

extern AED_REAL *WQ_Vars;   //# water quality array : [nvars, nlayers]
extern AED_REAL *WQS_Vars;  //# water quality benthics array : [nvars]
extern AED_REAL *WQD_Vars;  //# water quality diagnostics array : [nvars, nlayers]
extern AED_REAL *WQDS_Vars; //# water quality sheet diagnostics array : [nvars]

extern int *PTM_Stat;       //# water quality diagnostics array : [nlayers, nvars]
extern AED_REAL *PTM_Vars;  //# water quality sheet diagnostics array : [nvars]

/*----------------------------------------------------------------------------*/

extern AED_REAL MaxHeight;   //# maxmimum height of reservoir
extern AED_REAL CrestHeight; //# crest height of reservoir
extern AED_REAL LenAtCrest;  //# length of reservoir at crest
extern AED_REAL WidAtCrest;  //# width of reservoir at crest
extern AED_REAL VolAtCrest;  //# volume at crest level
extern AED_REAL MaxVol;      //# volume at max level
extern AED_REAL Base;        //# bottom elevation of reservoir

extern AED_REAL Benthic_Light_pcArea;
extern AED_REAL Benthic_Imin;
extern AED_REAL MaxArea;

/*----------------------------------------------------------------------------*/
// INFLOWS & OUTFLOWS
extern int NumInf;                 //# number of inflows
extern InflowDataType Inflows[];   //# Array of Inflows

extern int NumOut;                 //# Number of outflows
extern OutflowDataType Outflows[]; //# Array of Outflows
extern AED_REAL crit_val;
extern AED_REAL crit_dep;
extern int crit_days;
extern CLOGICAL CRITabove;
extern CLOGICAL MIXwithdraw;
extern CLOGICAL COUPLoxy;
extern AED_REAL WithdrawalTemp;
extern AED_REAL fac_range_upper, fac_range_lower;
extern AED_REAL MINlaketemp;

extern AED_REAL crest_width;
extern AED_REAL crest_factor;

extern CLOGICAL single_layer_draw;
extern AED_REAL outflow_thick_limit;

extern CLOGICAL evap_from_file;
extern AED_REAL f_evap_ts_prop;

/*----------------------------------------------------------------------------*/
//
extern int NumDif;
extern AED_REAL mol_diffusivity[];

/*----------------------------------------------------------------------------*/
// SURFACE
extern SurfaceDataType SurfData; //# Surface Data
extern MetDataType MetData;      //# Meteorological data

// Set when a restart file has been loaded this run. Gates the restart-resume
// handling: skipping the init_model resize_internals/check_layer_thickness
// (the saved derived layer fields are used verbatim), and the mid-day-resume
// overrides in do_model / do_model_non_avg.
extern int Restart_loaded;

// Sub-daily-level SWold at the end of the last sub-iter executed in the
// previous run. On a mid-day resume it seeds the first sub-iter's SWold so it
// matches what a continuous run had at iclock=startTOD. Only consumed in
// non-subdaily mode (subdaily=.false.), where calculate_qsw interpolates the
// timestep's shortwave from yesterday's (SWold) and today's daily totals; in
// subdaily mode calculate_qsw ignores SWold so this has no effect.
extern AED_REAL Restart_SWold;
extern AED_REAL runoff_coef;
extern AED_REAL rain_threshold;
extern CLOGICAL catchrain;
extern CLOGICAL atm_stab;         //# Account for non-neutral atmospheric stability
extern AED_REAL coef_wind_drag;   //# = 0.0013;
extern AED_REAL coef_wind_chwn;   //# = 0.0013;
extern AED_REAL CD;               //# = 0.0013;
extern AED_REAL CE;               //# = 0.0013;
extern AED_REAL CH;               //# = 0.0013;
extern int      subdaily;         //# = FALSE;
extern int      rad_mode;
extern int      albedo_mode;
extern int      cloud_mode;
extern int      light_mode;
extern int      n_bands;
extern AED_REAL *light_extc;
extern AED_REAL *energy_frac;

extern CLOGICAL use_met_atm_pres;
extern AED_REAL biodrag;

extern AED_REAL salt_fall;
/*----------------------------------------------------------------------------*/
// MORPHOMETRY
extern int Nmorph;                //# Number of data points

extern AED_REAL  MphInc;
extern AED_REAL *MphLevelArea;    //# area of each layer determined by linear interpolation
extern AED_REAL *dMphLevelArea;   //# gradients of area between 0.1m layers
extern AED_REAL *dMphLevelVol;    //# gradients of volume between 0.1m layers
extern AED_REAL *dMphLevelVolda;  //#
extern AED_REAL *MphLevelVol;     //# volume of each layer determined by linear interpolation
extern AED_REAL *MphLevelVoldash; //#

/*----------------------------------------------------------------------------*/
// MIXING
extern AED_REAL vel;
extern AED_REAL WaveNumSquared;
extern AED_REAL XMoment1;
extern AED_REAL einff;          //# change in potential energy (see do_inflows)
extern AED_REAL coef_mix_KH;    //# Kelvin-Helmholtz billows
extern AED_REAL coef_mix_conv;  //# convective overturn
extern AED_REAL coef_mix_shear; //# shear efficiency
extern AED_REAL coef_mix_turb;  //# unsteady effects
extern AED_REAL coef_wind_stir; //# wind stirring
extern AED_REAL coef_mix_hyp;   //# efficiency of hypolimnetic mixing
extern AED_REAL coef_mix_shreq; //# unsteady effects

extern CLOGICAL non_avg;
extern int      deep_mixing;          //# = 0 => off > 0 => on
extern int      surface_mixing;

extern int      density_model;

extern AED_REAL albedo_mean;          //# mean albedo
extern AED_REAL albedo_amplitude;     //#  albedo seasonal amplitude
extern AED_REAL lw_factor ;
extern AED_REAL lw_offset ;

//# DepMX is the layer height of the meta top on the previous timestep.
extern AED_REAL DepMX;

extern AED_REAL PrevThick;            //# mixed layer thickness from previous time step

extern AED_REAL gPrimeTwoLayer;       //# Reduced gravity for int wave estimate

extern AED_REAL Energy_AvailableMix;  //# Total available energy to mix (carries over from previous timesteps)

extern AED_REAL Mass_Epi;             //# Sigma mass of Epilimnion (surface layer after Kelvin-Helmholtz) kg

extern AED_REAL OldSlope;
extern AED_REAL Time_end_shear;       //# Time left before shear cut off [hours]
extern AED_REAL Time_start_shear;     //# Time count since start of sim for shear period start [hours]
extern AED_REAL Time_count_end_shear; //# Time count since start of sim for shear period end [hours]
extern AED_REAL Time_count_sim;       //# Time count since start of simulation [hours]

extern AED_REAL Half_Seiche_Period;   //# One half the seiche period
extern AED_REAL Thermocline_Height;   //# Height at the top of the metalimnion [m]
extern AED_REAL FO;
extern AED_REAL FSUM;
extern AED_REAL u_f;
extern AED_REAL u0;
extern AED_REAL u_avg;
extern int Mixer_Count;



/*----------------------------------------------------------------------------*/
// SNOWICE
extern AED_REAL snow_albedo_factor;
extern AED_REAL snow_rho_max;
extern AED_REAL snow_rho_min;
extern AED_REAL snow_water_equivalent;
extern AED_REAL snow_rain_compact;
extern AED_REAL K_ice_white;
extern AED_REAL K_ice_blue;
extern AED_REAL K_water;
extern AED_REAL f_sw_wl1;
extern AED_REAL f_sw_wl2;
extern AED_REAL attn_ice_blue_wl1;
extern AED_REAL attn_ice_blue_wl2;
extern AED_REAL attn_ice_white_wl1;
extern AED_REAL attn_ice_white_wl2;
extern AED_REAL attn_snow_wl1;
extern AED_REAL attn_snow_wl2;
extern AED_REAL rho_ice_blue;
extern AED_REAL rho_ice_white;
extern AED_REAL min_ice_thickness;
extern AED_REAL dt_iceon_avg;
extern AED_REAL avg_surf_temp_thres;

/*----------------------------------------------------------------------------*/
// SEDIMENT
extern CLOGICAL sed_heat_sw;
extern int      sed_heat_model;
extern AED_REAL sed_heat_Ksoil;
extern AED_REAL sed_temp_depth;
extern AED_REAL *sed_temp_mean;
extern AED_REAL *sed_temp_amplitude;
extern AED_REAL *sed_temp_peak_doy;
extern AED_REAL *sed_reflectivity;
extern AED_REAL *sed_roughness;
// sed_heat_model == 2 (dynamic soil/sediment temperature model)
extern int      n_sed_layers;
extern AED_REAL *sed_layer_depth;
extern AED_REAL *sed_vwc;
extern AED_REAL sed_spinup_days;
// dynamic soil-temperature model thermal properties + per-step dt.
// These symbols are OWNED by Fortran (glm_types.F90 BIND(C)); here we only extern them.
extern AED_REAL soil_dt;
extern AED_REAL sed_k_mineral, sed_k_water, sed_k_air;
extern AED_REAL sed_c_mineral, sed_c_water, sed_c_air;
extern AED_REAL sed_bulk_density, sed_mineral_density, sed_porosity, sed_deep_temp;

/*----------------------------------------------------------------------------*/
// GROUNDWATER
extern int   gw_mode;    //# mode
extern char *gw_file;    //# name of gw file
extern AED_REAL *K_gw;   //# turn off evaporation
extern AED_REAL *L_gw;   //# turn off evaporation

/*----------------------------------------------------------------------------*/
// FETCH
extern CLOGICAL  fetch_sw;
extern int       fetch_ndirs;
extern AED_REAL *fetch_dirs;
extern AED_REAL *fetch_scale;
extern AED_REAL  fetch_height;
extern AED_REAL  fetch_porosity;

extern int       fetch_mode;
extern AED_REAL  fetch_aws;
extern AED_REAL  fetch_xws;
extern char     *fetch_fws;

/*----------------------------------------------------------------------------*/
// LITTORAL
extern CLOGICAL littoral_sw;

//------------------------------------------------------------------------------
// PARTICLE TRANSPORT MODEL
extern CLOGICAL ptm_sw;
extern int max_particle_num;   //# number of particles
extern ParticleDataType *Particle;
extern AED_REAL particle_density;
extern AED_REAL particle_diameter;
extern AED_REAL settling_velocity;
extern int init_particle_num;
extern AED_REAL settling_efficiency;
extern AED_REAL *inflow_conc;  //# concentration of particles per ?? in the inflow

extern partgroup *Particles;
extern int Num_PTM_Vars;  //# number of AED particle-tracked WQ variables (n_ptm_vars)

/*----------------------------------------------------------------------------*/
// TIME
extern AED_REAL timezone_r, timezone_m, timezone_i, timezone_o;
extern int nDays;          //# number of days to simulate
extern AED_REAL timestep;
extern int noSecs;
extern AED_REAL yearday;   //# day of year
/*----------------------------------------------------------------------------*/
// HEAT PUMP SYSTEM
extern int heat_pump_switch;                 //# enable/disable heat pump (0=off, 1=on, 2=heat flux mode) - add other heatflux type after
extern int heat_pump_inflow_idx;             //# index of the heat pump inflow to be linked to the outflow
extern int heat_pump_outflow_idx;            //# index of the heat pump outflow to be linked with the inflow
extern AED_REAL heat_pump_temp_change;       //# temperature increase (°C)
extern AED_REAL heat_pump_heat_flux;         //# heat flux input (W)
extern AED_REAL heat_pump_dynamic_heat_flux; //# dynamic heat flux from CSV (W)
extern AED_REAL heat_pump_current_heat_flux; //# current dynamic heat flux from CSV (W)

/*----------------------------------------------------------------------------*/
// OXYGENATION SYSTEM (artificial aeration / oxygenator)
extern int      oxygenation_mode;         //# 0=off, 1=direct(constant), 2=direct(CSV), 3=recirculation
extern int      oxy_num;                  //# number of direct-addition devices (modes 1 & 2)
extern char    *oxy_name;                 //# AED dissolved-oxygen variable name (default 'OXY_oxy')
extern AED_REAL oxy_max;                  //# optional O2 concentration cap (<=0 disables)
extern int      oxy_o2_idx;               //# resolved WQ index of the O2 variable (set at init)
extern int      oxy_input_type[];               //# per device: 1=mass rate, 2=flow*conc
extern AED_REAL oxy_height[];             //# per device: height above bottom for addition (m)
extern AED_REAL oxy_load[];               //# per device: mode 1 O2 mass loading rate (mass/day)
extern AED_REAL oxy_flow[];               //# per device: mode 2 nominal flow (m3/day)
extern AED_REAL oxy_conc[];               //# per device: mode 2 O2 concentration
//# Approach 3 (withdraw + re-inject with added O2) - self-contained (oxygenation_mode==3)
extern AED_REAL oxy_recirc_withdraw_height; //# height above bottom to withdraw from (m)
extern AED_REAL oxy_recirc_return_height;   //# height above bottom to return to (m)
extern AED_REAL oxy_recirc_flow;            //# recirculation rate (m3/s, or from CSV)
extern AED_REAL oxy_recirc_add;             //# O2 mass loading rate (mass/day, or from CSV)

/*----------------------------------------------------------------------------*/
// DEBUGGING
extern CLOGICAL dbg_mix;   //# debug output from mixer
extern CLOGICAL no_evap;   //# turn off evaporation
extern int      quiet;     //# turn down output messages

/*----------------------------------------------------------------------------*/
// C-Fortran shared (see glm_types.F90
extern CINTEGER MaxLayers;   //# Maximum number of layers in this sim
extern LakeDataType *Lake;
extern CINTEGER  n_zones;    //# number of sediment zones
extern ZoneType *theZones;

//extern C_PTR pMetData;
//extern C_PTR pSurfData;

//extern LakeDataType theLake;
//extern MetDataType MetData;
//extern SurfaceDataType SurfData;
//extern ZoneType theZones;

extern FLOGICAL mobility_off;
extern FLOGICAL bioshade_feedback;
extern FLOGICAL repair_state;
extern FLOGICAL do_plots;
extern FLOGICAL link_rain_loss;
extern FLOGICAL link_solar_shade;
extern FLOGICAL link_bottom_drag;
extern FLOGICAL ice;
extern FLOGICAL do_particle_bgc;
extern FLOGICAL link_ext_par;

extern CINTEGER split_factor;
extern CINTEGER ode_method;
extern CINTEGER benthic_mode;

extern AED_REAL rain_factor;
extern AED_REAL sw_factor;
extern AED_REAL friction;

extern AED_REAL Kw;
extern AED_REAL dt;

extern AED_REAL Longitude;
extern AED_REAL Latitude;

/******************************************************************************/
void allocate_storage(void);
void set_c_wqvars_ptr(AED_REAL *iwqvars);
void set_c_wqsvars_ptr(AED_REAL *iwqsvars);
void set_c_wqdvars_ptr(AED_REAL *iwqd, AED_REAL *iwqds, int *nwqd, int *nwqds);
void set_c_ptmstat_ptr(int *iptms);
void set_c_ptmenv_ptr(AED_REAL *iptmv);
void set_c_num_ptm_vars(int *n);
void debug_print_lake(void);
void debug_initialisation(int which);
void debug_initialisation_(int *which);

//# NB: The order of array indices has been reversed as of V4
// #  define _WQ_Vars(var,lyr) WQ_Vars[_IDX_2d(MaxLayers,Tot_WQ_Vars,lyr,var)]
// //#define _WQS_Vars(var,lyr) WQ_Vars[_IDX_2d(MaxLayers,Tot_WQ_Vars,lyr,var)]
// #  define _WQD_Vars(var,lyr) WQD_Vars[_IDX_2d(MaxLayers,Num_WQD_Vars,lyr,var)]
// #  define _WQDS_Vars(var,lyr) WQDS_Vars[_IDX_2d(Num_WQDS_Vars,lyr,var)]

#  define _WQ_Vars(var,lyr) WQ_Vars[_IDX_2d(Tot_WQ_Vars,MaxLayers,var,lyr)]
//#define _WQS_Vars(var,lyr) WQ_Vars[_IDX_2d(Tot_WQ_Vars,MaxLayers,var,lyr)]
#  define _WQD_Vars(var,lyr) WQD_Vars[_IDX_2d(Num_WQD_Vars,MaxLayers,var,lyr)]
#  define _WQDS_Vars(var,lyr) WQDS_Vars[_IDX_2d(Num_WQDS_Vars,var,lyr)]

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif

/*============================================================================*/
#endif
