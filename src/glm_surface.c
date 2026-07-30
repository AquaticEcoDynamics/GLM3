/******************************************************************************
 *                                                                            *
 * glm_surface.c                                                              *
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glm.h"

#include "glm_types.h"
#include "glm_const.h"
#include "glm_globals.h"

#include "glm_util.h"
#include "aed_time.h"
#include "glm_mixu.h"
#include "glm_layers.h"

#include "glm_wqual.h"  // for sediment heating - move?

#include "glm_debug.h"

#if DEBUG
#define dbgprt(...) fprintf(stderr, __VA_ARGS__)
//#define dbgprt(...) /* __VA_ARGS__ */
#else
#  define dbgprt(...) /* __VA_ARGS__ */
#endif

// define the function type for reflection/refraction/...
typedef double (*ref_f)(double, double, double, double, int);

void solpond(int nband, int npoint,
             double depth, double *delz, double rb, double hdir, double anglei, double hdif,
             double *energy, double *absorb, double *gx);

static double expint(double ri, double cr, double h, double cmu, int n);
static double fct(double ri, double cr, double h, double cmu, int n);
static double fdif(double rindex, double critw, double h, double cmu, int n);

static double ref(double ri, double cr, double cmu);
static void   gaus10(double ri, double cr, double c, double d, ref_f f, int n, double h, double *v);

static AED_REAL grischenko_albedo(AED_REAL Latitude, AED_REAL mDays);

int  atmos_stability(     AED_REAL *Q_latentheat,
                          AED_REAL *Q_sensible,
                          AED_REAL  wind_speed,
                          AED_REAL  temp_water,
                          AED_REAL  temp_air,
                          AED_REAL  humidity_surface,
                          AED_REAL  humidity_altitude,
                          AED_REAL  rho_air,
                          AED_REAL  rho_o,
                          AED_REAL  p_atm,
                          AED_REAL  latent_heat_vap,
                          AED_REAL *coef_wind_drag,
                          AED_REAL *coef_wind_chwn,
                          AED_REAL *zonL );

/******************************************************************************
 * Module variables                                                           *
 ******************************************************************************/

AED_REAL  AvgSurfTemp= 6.;  // Recent average of surface temp, for ice-on.

// Heat fluxes; these are made available for the lake.csv output
AED_REAL zonL;                 // average z/L - MO atmospheric stability
AED_REAL Q_shortwave;          // Solar radiation at water surface
AED_REAL Q_sensibleheat;       // Sensible heat flux
AED_REAL Q_latentheat;         // Evaporative heat flux
AED_REAL Q_longwave;           // Net longwave heat flux
static AED_REAL  Q_iceout;     // Heat flux out of snow/ice surface
static AED_REAL  Q_icemet;     // Heat flux at the surface due to meteorological forcing
static AED_REAL  Q_watermol;   // Heat flux through water due to molecular conductivity
static AED_REAL  Q_icewater;   // Heat flux across water/ice interface
static AED_REAL  Q_surflayer;  // Heat flux through the water surface
static AED_REAL  Q_underflow;  // Heat flux through water due to flow under the ice
//static AED_REAL  U_flow;     // Velocity estimate of the underflow

//static AED_REAL  snow_rain_compact = 1. ; //update based on timestep and scaling

void recalc_surface_salt(void);

AED_REAL calculate_qsw(int kDays, int mDays, int iclock,
                       AED_REAL Latitude, AED_REAL SWOld,
                       AED_REAL ShortWave, AED_REAL WindSp);


/******************************************************************************
 *                                                                            *
 *  atm_density  [kg m^3]                                                     *
 *     = 0.348 * (1 + r)/(1 + 1.61*r) * p/T   of TVA, p 5.12                  *
 *                                                                            *
 *  Returns the density of the atmosphere given total atmospheric pressure    *
 *  and water vapour pressure.                                                *
 *                                                                            *
 *  Ref: TVA, [1972,p5.12]                                                    *
 *                                                                            *
 ******************************************************************************/
AED_REAL atm_density(AED_REAL atmosPressure, // (total) atmospheric pressure    [Pa]
                     AED_REAL vapPressure,   // water vapour (partial) pressure [Pa]
                     AED_REAL AirTemp)       // dry bulb air temperature        [Cel]
{
/*
    // Dry air
    return (p_atm)/(287.058 * (AirTemp + Kelvin));
*/
    //r_o = humidity_surface/(1-humidity_surface/c_gas);

    // Moist air
    // mixing ratio: r = Mwater/(Mwater+Mdry_air)
    // (the grams of vapor per kg of dry air)
    AED_REAL r = mwrw2a * vapPressure/(atmosPressure - vapPressure);

    //rho_a = 100.*p_atm / (c_gas*(AirTemp+Kelvin));

    return 1.0/c_gas * (1 + r)/(1 + r/mwrw2a) * atmosPressure/(AirTemp+Kelvin);
}


/******************************************************************************
 * Performs thermal transfers across the lake surface (water and ice)         *
 ******************************************************************************/
void do_surface_thermodynamics(int jday, int iclock, int LWModel,
                          AED_REAL Latitude, AED_REAL SWOld, AED_REAL ShortWave)

   // LWModel   // type of longwave radiation
   // SWOld     // Total solar radiation at the surface for yesterday
   // ShortWave // Total solar radiation at the surface for today
{
/*----------------------------------------------------------------------------*/
    const AED_REAL  Temp_melt = 0.0;               //# temperature at which ice melts ~= temperature at which water freezes = 0.0oC
    const AED_REAL  Latent_Heat_Fusion = 334000.0; //# latent heat of fusion J/kg both ice & snow
    const AED_REAL  eps_water = 0.985;             //# emissivity of the water surface
/*----------------------------------------------------------------------------*/

//  int nband = 4;
    int npoint = 10;

    AED_REAL *LayerThickness;
    AED_REAL *heat;
    int *layer_zone;

    AED_REAL *gx;

    AED_REAL p_atm;          //# Atmospheric pressure in hectopascals, eg. 101300 Pa
    AED_REAL rho_air, rho_o; //# atm_density
    AED_REAL SatVap_surface; //# Saturated vapour pressure at surface layer or top ice layer
    AED_REAL altitude;
    AED_REAL latent_heat_vap;

    AED_REAL WindSp;         //# Wind speed corrected by multiplicative factor
    AED_REAL CloudCover;     //# Cloud cover as a fraction of 1 used to calculate incident long wave radiation

    AED_REAL dTemp;          //# Change in temperature due to atmospheric heat exchange
    AED_REAL dHt_WhiteIce;   //# Change in height of white ice when snow weight exceeds buoyancy potential

    AED_REAL Q_lw_in;        //# Long wave radiation from atmosphere W/m2
    AED_REAL Q_lw_out;       //# Long wave radiation emission W/m2
    AED_REAL Q_whiteice;     //# heat flux due to snow turning into white ice
    AED_REAL Q_snowice;      //# heat flux due to snow turning into white ice
    AED_REAL Q_rain;         //# heat flux due to rain either melting ice or freezing into ice

    AED_REAL AAA, BBB, CCC, DDD, EEE, FFF, GGG;
    AED_REAL Temp_ice;        //# temperature of the top layer of the ice cover model oC
    AED_REAL T01_OLD,T01_NEW; //# Not sure how these are used as Temp_ice always == 0 >= melting
    AED_REAL rho_snow;        //# Density of snow layer kg/m3
    AED_REAL rho_snow_old;    //# Density of snow layer compacted by rain kg/m3
    AED_REAL compact_snow;    //# Fraction of snow compacted by rain []
    AED_REAL K_snow;          //# thermal conductivity of snow [W/m2/oC]
    AED_REAL Q_latent_ice;    //# Latent heat of melting or freezing of ice [W/m2]

    AED_REAL SUMPO4, SUMTP, SUMNO3, SUMNH4, SUMTN, SUMSI;

    AED_REAL flankArea, NotPARLight_s, NotPARLight_sm1;


    //# NB These are set to 1.0 to remove a "possibly uninitialised" warning
    //#     They are initialised if litoral_sw is true, and only used if litoral_sw
    //#     so it should be OK, but should check from time to time.
    AED_REAL onshoreDensity = 1.0, offshoreDensity = 1.0;
    AED_REAL onshoreVol = 1.0, offshoreVol = 1.0;

    //# New parameters for sediment heat flux estimate
    AED_REAL KSED, TYEAR, ZSED;

    int i, z, j, wqidx;
    int kDays;
//  int non_neutral_converged;

    //int nband, npoint;
    AED_REAL depth, rb, anglei, hdir, hdif;


//  int sed_layers = 12;
//  AED_REAL ztemp = 20.0;
    AED_REAL soil_heat_flux = 0.0;
//  AED_REAL *sed_depths=NULL;
//  AED_REAL *sed_vwc=NULL;
//  AED_REAL *sed_temps=NULL;


/*----------------------------------------------------------------------------*/


    /**********************************************************************
    * Prior to surface heating, check for local runoff inputs from rain
    * and merge with the surface layer properties
    ***********************************************************************/

    AED_REAL catch_runoff = zero;


    if ( catchrain && MetData.Rain>rain_threshold ) {

        // Compute runoff in m3 for this time step
        catch_runoff = ( MaxArea - Lake[surfLayer].LayerArea ) *
                       ( MetData.Rain - rain_threshold ) * (noSecs/SecsPerDay) * runoff_coef;

        // Catchment runoff as a special inflow to the surface layer
        AED_REAL RunoffFlowRate = catch_runoff;
        Lake[surfLayer].Temp = combine(Lake[surfLayer].Temp,
                                       Lake[surfLayer].LayerVol,
                                       Lake[surfLayer].Density,
                                       MetData.AirTemp, RunoffFlowRate,
                                       calculate_density(MetData.AirTemp,zero+0.01));
        Lake[surfLayer].Salinity = combine(Lake[surfLayer].Salinity,
                                       Lake[surfLayer].LayerVol,
                                       Lake[surfLayer].Density,
                                       zero+0.01, RunoffFlowRate,
                                       calculate_density(MetData.AirTemp,zero+0.01));

        for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
            _WQ_Vars(wqidx, surfLayer) = combine_vol(_WQ_Vars(wqidx, surfLayer),
                                       Lake[surfLayer].LayerVol,
                                       zero, RunoffFlowRate);

        Lake[surfLayer].Density = calculate_density(Lake[surfLayer].Temp,
                                                    Lake[surfLayer].Salinity);
        Lake[surfLayer].LayerVol = Lake[surfLayer].LayerVol+RunoffFlowRate;

        Lake[botmLayer].Vol1 = Lake[botmLayer].LayerVol;
        if (surfLayer != botmLayer) {
            for (j = (botmLayer+1); j <= surfLayer; j++)
                    Lake[j].Vol1 = Lake[j-1].Vol1 + Lake[j].LayerVol;
        }

        //# Make adjustments to correct layer heights following vol addition
        resize_internals(2,botmLayer);
        check_layer_thickness();

        SurfData.dailyRunoff += catch_runoff;
    }

    /**********************************************************************
    * Now get ready for surface heating
    ***********************************************************************/

    LayerThickness = calloc(MaxLayers, sizeof(AED_REAL));
    heat = calloc(MaxLayers, sizeof(AED_REAL));
    layer_zone = calloc(MaxLayers, sizeof(int));
    gx = calloc(MaxLayers, sizeof(AED_REAL));
    for (i=0; i < MaxLayers; i++) {
        LayerThickness[i] = 0.0; heat[i] = 0.0;
        layer_zone[i] = 0; gx[i] = 0.0;
    }

    SUMPO4 = zero; SUMTP = zero; SUMNO3 = zero;
    SUMNH4 = zero; SUMTN = zero; SUMSI = zero;

    T01_NEW = zero; T01_OLD = zero;

    SurfData.Evap = zero;
    Q_longwave = zero;
    Q_latent_ice = zero;
    Q_whiteice = zero;
    Q_snowice = zero;
    Q_rain = zero;
    rho_snow = SurfData.RhoSnow;

    //# Snow conductivity (eqn. 12 Rogers et al. 1995)
    //  To avoid division by zero set to 0.31 when no snow
    if (SurfData.delzSnow  ==  0.)
        K_snow = 0.31;
    else
        K_snow = 0.021+(0.0042*rho_snow)+(2.2E-09*pow(rho_snow, 3));

    // Initialize the ice temperature to the air temperature
    Temp_ice = MetData.AirTemp;

    // Modify wind speed experienced by the surafce layer, if ice is present
    if (ice) WindSp = 0.0001;
    else     WindSp = MetData.WindSpeed;

    /**********************************************************************
     * ATMOSPHERIC CONDITIONS
     * Set atmospheric properties
     *********************************************************************/

    // Get atmospheric pressure and density

    // get altitude
    altitude = Base + Lake[surfLayer].Height + 5.;
    // update air pressure at altitude
//MetData.AirPres = atm_pressure_sl;
    p_atm = ((100.*MetData.AirPres) * pow((1 - 2.25577e-5*altitude),5.25588))/100.; //> hPa
////fprintf(stderr, "MetData.AirPres %f p_atm %f altitude %f\n", MetData.AirPres, p_atm, altitude);
    // gte the saturation vapour pressure at the water surface
    SatVap_surface = saturated_vapour(Lake[surfLayer].Temp); //> hPa
    // calculate gas constant for moist air
// This is not used anywhere
//  c_gas_m = c_gas*(1 + 0.608*(mwrw2a*MetData.SatVapDef/p_atm)); //> J/kg/K
    // update the latent heat of vaporisation, based on temperature
    latent_heat_vap = 2.501e6 - 2370*Lake[surfLayer].Temp; //> J/kg
    // update air density, right above the lake surface and at the ref height
    rho_air = atm_density(p_atm*100.0,MetData.SatVapDef*100.0,MetData.AirTemp); //>kg/m3
    rho_o = atm_density(p_atm*100.0,SatVap_surface*100.0,Lake[surfLayer].Temp); //>kg/m3

/*      printf(">>>>");
        printf(">altitude = %10.5f\n",altitude);
        printf(">p_atm = %10.5f\n",p_atm);
        printf(">SatVap_surface = %10.5f\n",SatVap_surface);
        printf(">c_gas_m = %10.5f\n",c_gas_m);
        printf(">latent_heat_vap = %10.5f\n",latent_heat_vap);
        printf(">rho_air = %10.5f\n",rho_air);
        printf(">rho_o = %10.5f\n",rho_o);
        printf(">MetData.SatVapDef = %10.5f\n",MetData.SatVapDef);
*/

    /**********************************************************************
     * SHORTWAVE INPUT
     * Assess surface shortwave flux hitting the water or ice
     * and add to the surface layer, plus work out heating
     *********************************************************************/

    //# Get shortwave intensity (W/m2) - either from provided sub-daily data,
    //  or approximate from daily total
    //  note, convert "julian days" to the "days since the start of the year"
    Q_shortwave = calculate_qsw(day_of_year(jday - 1), day_of_year(jday),
                                    iclock, Latitude, SWOld, ShortWave, WindSp);


    //# Into the surfLayer goes Qsw(surfLayer)-Qsw(surfLayer-1) over the
    //  area common to layers surfLayer and surfLayer-1, and Qsw(surfLayer)
    //  over the rest of the area. The units of heat[surfLayer] are
    //  joules/sec; units of area(surfLayer) are m**2
    for (i = (botmLayer+1); i <= surfLayer; i++)
        LayerThickness[i] = Lake[i].Height-Lake[i-1].Height;

    LayerThickness[botmLayer] = Lake[botmLayer].Height;

    //# PAR entering the upper most (surface) layer
    if (!ice)
        // Assume PAR only (45% of the incident shortwave radiation
        // penetrates beyond the surface layer)
        Lake[surfLayer].Light = 0.45 * Q_shortwave;
    else
        // If there is ice cover the heating due to shortwave radiation
        // is attenuated across the three ice layers
        Lake[surfLayer].Light =
            (0.45 * Q_shortwave *
              ((f_sw_wl1 * exp(-attn_snow_wl1 * SurfData.delzSnow)) +
                        (f_sw_wl2 * exp(-attn_snow_wl2 * SurfData.delzSnow)))) *
              ((f_sw_wl1 * exp(-attn_ice_white_wl1 * SurfData.delzWhiteIce)) +
                        (f_sw_wl2 * exp(-attn_ice_white_wl2 * SurfData.delzWhiteIce))) *
              ((f_sw_wl1 * exp(-attn_ice_blue_wl1  * SurfData.delzBlueIce)) +
                        (f_sw_wl2 * exp(-attn_ice_blue_wl2  * SurfData.delzBlueIce)));

    //# PAR flux entering the layer below the surface layer
    if (surfLayer > 0)
        Lake[surfLayer-1].Light = Lake[surfLayer].Light *
                     exp(-Lake[surfLayer].ExtcCoefSW*LayerThickness[surfLayer]);

    //# NIR/UV flux (Not using Q_shortwave as it hasnt been through ice)
    NotPARLight_s = (1-0.45)*(Lake[surfLayer].Light/0.45) ;
    NotPARLight_sm1 = NotPARLight_s *
                  exp(-2.*Lake[surfLayer].ExtcCoefSW*LayerThickness[surfLayer]);

    //# Heating of the surface layer due to PAR & NIR/UV:
    if (surfLayer > 0) {
        flankArea = Lake[surfLayer].LayerArea - Lake[surfLayer-1].LayerArea;
        //  1. PAR amount in "deep water" area
        heat[surfLayer] = Lake[surfLayer-1].LayerArea * (Lake[surfLayer].Light - Lake[surfLayer-1].Light);
        //  2. PAR amount in "littoral water" area
        heat[surfLayer] = heat[surfLayer] + (flankArea * Lake[surfLayer].Light)
              - (flankArea * (1.-MAX(1.,sed_reflectivity[0])) * Lake[surfLayer-1].Light);
        //  3. NIR/UV amount in "deep water" area
        heat[surfLayer] = heat[surfLayer] + Lake[surfLayer-1].LayerArea * (NotPARLight_s - NotPARLight_sm1);
    } else {
        flankArea = Lake[surfLayer].LayerArea;
        //  1. PAR amount in "deep water" area
        heat[surfLayer] = Lake[surfLayer].Light;
        //  2. PAR amount in "littoral water" area
        heat[surfLayer] = heat[surfLayer] + (flankArea * Lake[surfLayer].Light);
    }
    //  4. NIR/UV amount in "littoral water" area
    heat[surfLayer] = heat[surfLayer] + (flankArea * NotPARLight_s)
          - (flankArea * (1.-MAX(1.,sed_reflectivity[0])) * NotPARLight_sm1);

//    heat[surfLayer] = heat[surfLayer] + flankArea * (NotPARLight_s - (1.-MAX(1.,sed_reflectivity[0]))*NotPARLight_sm1);
//    heat[surfLayer] = Lake[surfLayer-1].LayerArea * (Lake[surfLayer].Light - (1.-MAX(1.,sed_reflectivity[0]))*Lake[surfLayer-1].Light) ;
//    heat[surfLayer] = heat[surfLayer] + flankArea * (Lake[surfLayer].Light - Lake[surfLayer-1].Light*0.9);
//    heat[surfLayer] = heat[surfLayer] + (Lake[surfLayer].LayerArea - Lake[surfLayer-1].LayerArea) * Lake[surfLayer].Light * 0.1);
//    //# Heating due to NIR/UV/
//    //heat[surfLayer] = heat[surfLayer] + 0.55 * Q_shortwave * Lake[surfLayer].LayerArea;
//    heat[surfLayer] = heat[surfLayer] + Lake[surfLayer-1].LayerArea * (NotPARLight_s - NotPARLight_sm1);
//    heat[surfLayer] = heat[surfLayer] + flankArea * (NotPARLight_s - (1.-MAX(1.,sed_reflectivity[0]))*NotPARLight_sm1);
//
//    if (littoral_sw) {
//      onshoreDensity = calculate_density(Lake[onshoreLayer].Temp,Lake[surfLayer].Salinity);
//      offshoreDensity = calculate_density(Lake[offshoreLayer].Temp,Lake[surfLayer].Salinity);
//      offshoreVol = Lake[surfLayer-1].LayerArea * LayerThickness[surfLayer];
//      onshoreVol = Lake[surfLayer].LayerVol - offshoreVol;/
//
//      heat[offshoreLayer] = Lake[surfLayer-1].LayerArea * (Lake[surfLayer].Light - Lake[surfLayer-1].Light)
//                          + Lake[surfLayer-1].LayerArea * (NotPARLight_s - NotPARLight_sm1);
//      heat[onshoreLayer]  = flankArea * (Lake[surfLayer].Light)
//                          + flankArea * (NotPARLight_s - NotPARLight_sm1*0.95);
//    }

    //# Summarise total daily short wave radiation (J/day), stored for lake.csv
    SurfData.dailyQsw += Lake[surfLayer].LayerArea * Q_shortwave * noSecs;
    SurfData.Qsw = Q_shortwave;


    // ---- MH TEST SOLPOND IN PROGRESS ---- //
    if (light_mode == 2){

        Lake[surfLayer].Light = Q_shortwave;

        //# Advanced option - compute the light penetration suing the integral of light adsorption
        depth = Lake[surfLayer].Height;
        rb = 0.3;
        anglei = 10;
        hdir = Q_shortwave * 0.9;
        hdif = Q_shortwave * 0.1;
        npoint = NumLayers+1;
        memset(gx, 0, sizeof(AED_REAL)*MaxLayers);

        //  printf(">solpond = %10.1f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",energy_frac[0],energy_frac[1],energy_frac[2],energy_frac[3],light_extc[0],light_extc[1],light_extc[2],light_extc[3]);

        //solpond(n_bands, npoint, depth, rb, hdir, anglei, hdif, energy, absorb, gx);

        solpond(n_bands, npoint, depth, LayerThickness, rb, hdir, anglei, hdif, energy_frac, light_extc, gx);

        printf(">solpond = %10.1f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
                           gx[0],gx[1],gx[2],gx[3],gx[4],gx[5],gx[6],gx[7],gx[8]);


        printf("0light  = %10.1f\n", Lake[surfLayer].Light);
        for (i = surfLayer-1; i >= botmLayer; i-- ) {
            Lake[i].Light = Lake[i+1].Light - (gx[i]-gx[i+1])*LayerThickness[i+1];
            printf(">light  = %10.1f %6.2f \n",Lake[i].Light,LayerThickness[i+1]);
        }
    }
    // MH

    /**********************************************************************
     * SURFACE MASS FLUXES ON ICE
     * Assess surface mass fluxes of snow and rain on the ice
     * and check ice buoyancy for crackign and white ice formation
     *********************************************************************/

    if (iclock == 0 && ice) {
        AED_REAL BuoyantPotential;

        if (SurfData.delzSnow > 0.0) {
          // Snow cover as well as ice cover
          if (MetData.Snow > 0.0 && MetData.Rain >= 0.0) {
                // Snowfall on snow
                if (MetData.Rain == 0.0) MetData.Rain = MetData.Snow*snow_water_equivalent; //0.10;
                if (MetData.AirTemp > 0.0)
                    compact_snow = 0.166+0.834*(1.-exp(-1.*MetData.Rain/snow_rain_compact));
                else
                    compact_snow = 0.088+0.912*(1.-exp(-1.*MetData.Rain/snow_rain_compact));

                // Compact snow first
                rho_snow_old = rho_snow+(snow_rho_max-rho_snow)*compact_snow;
                SurfData.delzSnow = SurfData.delzSnow*rho_snow/rho_snow_old;

                // Determine the avg density of combined snow
                rho_snow = 1000.0*((rho_snow_old*SurfData.delzSnow/1000.0)+MetData.Rain)/(SurfData.delzSnow+(MetData.Snow));

                SurfData.delzSnow = SurfData.delzSnow+(MetData.Snow);

                // This should be a good estimate of water equivalent in from snow
                // Use rain because it should be volumetric equivalent of water (converted above, or supplied by user)
                SurfData.dailySnow += MetData.Rain * Lake[surfLayer].LayerArea;

                if (rho_snow > snow_rho_max) rho_snow = snow_rho_max;
                if (rho_snow < snow_rho_min) rho_snow = snow_rho_min;
                Q_rain = 0.0;

            } else if (MetData.Snow == 0.0 && MetData.Rain == 0.0) {
                // No snowfall or rainfall, compaction only of snow
                if (MetData.AirTemp > 0.0) compact_snow = 0.166;
                else                       compact_snow = 0.088;

                rho_snow_old = rho_snow+(snow_rho_max-rho_snow)*compact_snow;
                SurfData.delzSnow = SurfData.delzSnow*rho_snow/rho_snow_old;
                rho_snow = rho_snow_old;
                Q_rain = 0.0;

            } else if (MetData.Snow == 0.0 && MetData.Rain > 0.0) {
                // Rainfall on Snow.
                if (MetData.AirTemp > 0.0) {
                    //Check the air temperature and if AirTemp > 0 then
                    //  add the rain into the lake;
                    compact_snow = 0.166+0.834*(1.-exp(-1.*MetData.Rain/snow_rain_compact));
                    rho_snow_old = rho_snow+(snow_rho_max-rho_snow)*compact_snow;
                    SurfData.delzSnow = SurfData.delzSnow*rho_snow/rho_snow_old;
                    rho_snow = rho_snow_old;

                    Lake[surfLayer].Height = Lake[surfLayer].Height+MetData.Rain;
                    SurfData.dailyRain += MetData.Rain * Lake[surfLayer].LayerArea;
                    recalc_surface_salt();

                    if (Temp_ice == Temp_melt)
                        Q_rain = SPHEAT*(MetData.AirTemp-Temp_ice)*(MetData.Rain)/noSecs;
                    else
                        Q_rain = Latent_Heat_Fusion*MetData.Rain/noSecs;

                } else {
                    // If AirTemp < 0 then add the rainfall to snow
                    //MetData.Snow = MetData.Rain/snow_rho_max;
                    MetData.Snow = MetData.Rain/(snow_rho_max/rho0);
                    rho_snow_old = snow_rho_max;
                    SurfData.delzSnow = SurfData.delzSnow*rho_snow/rho_snow_old;
                    rho_snow = snow_rho_max;
                    SurfData.delzSnow = SurfData.delzSnow+MetData.Snow;
                    Q_rain = 0.0;
                    SurfData.dailyRain += MetData.Rain * Lake[surfLayer].LayerArea;
                }
            }
        } else {

            // No snow cover on the ice
            if (MetData.Snow > 0.0 && MetData.Rain >= 0.0) {

                // Snowfall on ice
                if (MetData.Rain == 0.0)MetData.Rain = MetData.Snow*0.10;
                SurfData.delzSnow = MetData.Snow;

                rho_snow = rho0*MetData.Rain/MetData.Snow;
                SurfData.dailySnow += MetData.Rain * Lake[surfLayer].LayerArea;

                if (rho_snow > snow_rho_max)rho_snow = snow_rho_max;
                if (rho_snow < snow_rho_min)rho_snow = snow_rho_min;
                Q_rain = 0.0;

            } else if (MetData.Snow == 0.0 && MetData.Rain == 0.0) {

                // No snowfall or rainfall
                Q_rain = 0.0;

            } else if (MetData.Snow == 0.0 && MetData.Rain > 0.0) {

                // Rainfall on ice - need to know whether it will contribute to water
                if (MetData.AirTemp > 0.0) {
                    Lake[surfLayer].Height = Lake[surfLayer].Height+MetData.Rain;
                    SurfData.dailyRain += MetData.Rain * Lake[surfLayer].LayerArea;
                    recalc_surface_salt();

                    if (Temp_ice == Temp_melt)
                        Q_rain = SPHEAT*(MetData.AirTemp-Temp_ice)*(MetData.Rain)/noSecs;
                    else
                        Q_rain = Latent_Heat_Fusion*MetData.Rain/noSecs;

                } else {
                    // If AirTemp < 0
                    MetData.Snow = MetData.Rain/snow_rho_max;
                    rho_snow = snow_rho_max;
                    SurfData.delzSnow = MetData.Snow;
                    SurfData.dailySnow += MetData.Snow * Lake[surfLayer].LayerArea * rho_snow/1000.0;
                    Q_rain = 0.0;
                }
            }
        }

        //# Archmides principle for the weight of Snow sitting on Ice -
        //  when the weight of snow exceeds the buoyancy of the ice
        //  cover, the ice will crack, and surface water will seep into snow
        //  layer leading to formation of White Ice
        if (rho_snow == 0.) rho_snow = snow_rho_min;

        BuoyantPotential = ((SurfData.delzBlueIce*(rho0-rho_ice_blue) +
                             SurfData.delzWhiteIce*(rho0-rho_ice_white))/rho_snow);

        if (SurfData.delzSnow > BuoyantPotential) {

            dHt_WhiteIce = SurfData.delzSnow-BuoyantPotential;

            Q_whiteice = ((Lake[surfLayer].Temp*SPHEAT+Latent_Heat_Fusion)*Lake[surfLayer].Density
                 *dHt_WhiteIce*(1.-(rho_snow/Lake[surfLayer].Density)))/SecsPerDay;

            SurfData.delzWhiteIce = SurfData.delzWhiteIce + SurfData.delzSnow - BuoyantPotential;

            //Adjust surface layer height down based on water moving into white ice
            Lake[surfLayer].Height -= dHt_WhiteIce * (rho_ice_white - rho_snow) / Lake[surfLayer].Density;

            Q_snowice = (Q_whiteice*SurfData.delzSnow)/(2.0*K_snow);
            SurfData.delzSnow = BuoyantPotential;

            recalc_surface_salt();

        } else {
            dHt_WhiteIce = 0.0;
            Q_whiteice = 0.0;
            Q_snowice = 0.0;
        }

        // rainfall/snowfall and dry deposition effects on nutrients
        SUMPO4 = SUMPO4 + 0.145 + MetData.RainConcPO4 * MetData.Rain ;
        SUMTP  = SUMTP  + 0.368 + MetData.RainConcTp  * MetData.Rain ;
        SUMNO3 = SUMNO3 + 0.507 + MetData.RainConcNO3 * MetData.Rain ;
        SUMNH4 = SUMNH4 + 0.507 + MetData.RainConcNH4 * MetData.Rain ;
        SUMTN  = SUMTN  + 2.576 + MetData.RainConcTn  * MetData.Rain ;
        SUMSI  = SUMSI  + 0.0   + MetData.RainConcSi  * MetData.Rain ;

    }  // end iclock == 0 && ice

    /**********************************************************************
     * NON-PENETRATIVE HEAT FLUXES
     * Now do non-pentrative heat fluxes, depending on ice presence
     *********************************************************************/
    if (!ice) {
        //#  Evaporative heat flux affects top layer only
        SatVap_surface = saturated_vapour(Lake[surfLayer].Temp);

        //# Q_latentheat [W/m2] = CE * rho_air * latent heat * psychro const / air_presssure * windspeed * VPD
        Q_latentheat = -CE * rho_air * latent_heat_vap * (mwrw2a/p_atm) * WindSp * (SatVap_surface - MetData.SatVapDef);

        if (Q_latentheat > 0.0) Q_latentheat = 0.0; // no condensation

        //# Conductive/sensible heat gain only affects the top layer.
        //  Q_sensibleheat [W/m2] = CH * rho_air * specific heat * windspeed * temp diff
        Q_sensibleheat = -CH * rho_air * cp_air * WindSp * (Lake[surfLayer].Temp - MetData.AirTemp);

        //# If chosen by user, do atmospheric stability correction routine
        coef_wind_chwn = CH;
        coef_wind_drag = CD;

        if ( atm_stab > 0 ) {
//           non_neutral_converged =
               atmos_stability(&Q_latentheat,
                               &Q_sensibleheat,
                                WindSp,
                                Lake[surfLayer].Temp,
                                MetData.AirTemp,
                                mwrw2a*SatVap_surface/p_atm,
                                mwrw2a*MetData.SatVapDef/p_atm,
                                rho_air,
                                rho_o,
                                p_atm,
                                latent_heat_vap,
                                &coef_wind_drag,
                                &coef_wind_chwn,
                                &zonL );

             SurfData.dailyzonL += zonL;
             //fprintf(stderr, " and now z/L %e = ...%e Q_l %e Q_s %e\n", zonL, coef_wind_drag, Q_latentheat, Q_sensibleheat);
        }

        //# Compute evaporative mass flux, in m/s, based on heat flux
        if ( no_evap )
            SurfData.Evap = 0.0;
        else {
            if (evap_from_file)
                SurfData.Evap = f_evap_ts_prop;
            else
                SurfData.Evap = Q_latentheat / (latent_heat_vap * Lake[surfLayer].Density);
        }

        //# Longwave emission (ie. longwave out), affecting only the top layer
        Q_lw_out = -Stefan_Boltzman * eps_water * pow((Kelvin+Lake[surfLayer].Temp), 4.0);

        //# Longwave absorption (ie. longwave in) also affects only top layer
        //  see Henderson-Sellers (1986) or Flerchinger (2009) for a good summary
        if (LWModel  ==  LW_CC) {
            // Cloud data is available
            AED_REAL eps_star = 0.8;  // default in case of a duff cloudmode value
            CloudCover = MetData.LongWave;
            switch (cloud_mode) {
                case 1:
                    // Idso and Jackson (1969) style
                    // eps_star = (1.0 + 0.275*CloudCover)*(1.0 - 0.261 * exp(-0.000777 * pow(-MetData.AirTemp, 2.0))); //
                    eps_star = (1.0 + 0.275*CloudCover)*(1.0 - 0.261 * exp(-0.000777 * pow(MetData.AirTemp, 2.0)));
                    break;
                case  2:
                    // Swinbank (1963) style
                    eps_star = (1.0 + 0.17*CloudCover*CloudCover) * (9.365e-6 * pow(MetData.AirTemp + Kelvin, 2.0));
                    break;
                case 3:
                    // Brutsaert (1975) style
                    eps_star = (1.0 + 0.275*CloudCover) * 1.24 * pow(MetData.SatVapDef/(MetData.AirTemp+Kelvin), 1/7) ;
                    break;
                case 4:
                    // Yajima 2014 - Tono Dam
                    eps_star = (1.0 - pow(CloudCover, 2.796) ) * 1.24
                                   * pow(MetData.SatVapDef/(MetData.AirTemp+Kelvin), 1/7)
                                          + 0.955 * pow(CloudCover, 2.796) ;
                    break;
            }

            Q_lw_in = (1 - 0.03) * eps_star * Stefan_Boltzman * pow((Kelvin+MetData.AirTemp), 4.0);
            Q_longwave = Q_lw_out + Q_lw_in * lw_factor;

        } else if (LWModel  ==  LW_IN) {
            // Incoming long-wave data is provided (assumes not corrected for lw albedo)
            Q_lw_in = (1 - 0.03) * MetData.LongWave;
            Q_longwave = Q_lw_out + Q_lw_in * lw_factor;

        } else if (LWModel  ==  LW_NET)
            // Net long-wave data is provided (ignore internally computed Q_lw_out)
            Q_longwave = MetData.LongWave * lw_factor;

        //# Apply longwave factor to net longwave rather than input data
        Q_longwave = Q_longwave + lw_offset ;

        //# Update surface layer heat (J/s)
        heat[surfLayer] = heat[surfLayer]+(Q_latentheat+Q_sensibleheat+Q_longwave)*Lake[surfLayer].LayerArea;

        if (littoral_sw) {
          heat[offshoreLayer] = heat[offshoreLayer] +
                                (Q_latentheat+Q_sensibleheat+Q_longwave)*Lake[surfLayer-1].LayerArea;
          heat[onshoreLayer] = heat[onshoreLayer] +
                                (Q_latentheat+Q_sensibleheat+Q_longwave)*flankArea;
        }

        //# Daily heat budget (J/day)
        SurfData.dailyQe += Q_latentheat * Lake[surfLayer].LayerArea * noSecs;
        SurfData.dailyQh += Q_sensibleheat * Lake[surfLayer].LayerArea * noSecs;
        SurfData.dailyQlw += Q_longwave * Lake[surfLayer].LayerArea * noSecs;
        SurfData.Qe = Q_latentheat;
        SurfData.Qh = Q_sensibleheat;
        SurfData.Qlw = Q_longwave;
    } else {
        //---------------------------------------------------------------------+
        //# Ice cover specific computation:
        //  The various atmospheric fluxes - evaporative, sensible heat, longwave
        //  for ice cover all depend on the ice surface temperature.  Note that
        //  wind here is not set to zero (USE MetData.WindSpeed, NOT WindSp].
        //  Emissivity for longwave OUT is 0.985 (eps_water)
        //---------------------------------------------------------------------+
        T01_NEW =  50.;   T01_OLD = -50.;
        Temp_ice = zero;
        while (1) {
            //# Saturated vapor pressure above snow and ice is different than above
            //  water. Algoirthm taken from Jeong (2009), ultimately from
            //  Mellor (1964) "Properties of Snow"
            SatVap_surface = (1.+(0.00972*Temp_ice)+(0.000042*pow(Temp_ice, 2.)))
                             *saturated_vapour(Temp_ice);

            // I think this might be wrong, resulting value seems way too small, even for ice
            // Q_latentheat = -3.9 * MetData.WindSpeed * (SatVap_surface - MetData.SatVapDef);
            //# Q_latentheat [W/m2] = CE * rho_air * latent heat * psychro const / air_presssure * windspeed * VPD
            Q_latentheat = -CE * rho_air * latent_heat_vap * (mwrw2a/p_atm) * MetData.WindSpeed
                           * (SatVap_surface - MetData.SatVapDef);
            if (Q_latentheat > zero) Q_latentheat = zero;

            //Q_sensibleheat = -CH * MetData.WindSpeed * (Temp_ice - MetData.AirTemp);
            //# Q_sensibleheat [W/m2] = CH * rho_air * specific heat * windspeed * temp diff
            Q_sensibleheat = -CH * rho_air * cp_air * WindSp * (Temp_ice - MetData.AirTemp);

            Q_lw_out = -Stefan_Boltzman * eps_water * pow((Kelvin+Temp_ice), 4.0);

            if (LWModel  ==  LW_CC) {
                CloudCover = MetData.LongWave;
                Q_lw_in = (1.0 + 0.275 * CloudCover) * Stefan_Boltzman * pow((Kelvin+MetData.AirTemp), 4.0) *
                          (1.0 - 0.261 * exp(-0.000777E-4 * pow(MetData.AirTemp, 2.0)));
                Q_longwave = Q_lw_out+Q_lw_in;
            } else if (LWModel  ==  LW_IN) {
                Q_lw_in = MetData.LongWave;
                Q_longwave = Q_lw_out+Q_lw_in;
            } else if (LWModel  ==  LW_NET)
                Q_longwave = MetData.LongWave;

            //# This is the net meteorological flux that drives the ice algorithm
            //  flux due to precipitation on ice or snow
            Q_icemet = Q_latentheat+Q_sensibleheat+Q_longwave+Q_rain;

            //# Now determine the new ice/snow surface temperature based on the balance
            //  fluxes h_ice and the expression for the upward heat flux
            AAA = (1. - exp(-attn_snow_wl1 * SurfData.delzSnow)) / (K_snow * attn_snow_wl1);
            BBB = exp(-attn_snow_wl1 * SurfData.delzSnow)
                * (1.-exp(-attn_ice_white_wl1 * SurfData.delzWhiteIce))/(attn_ice_white_wl1 * K_ice_white);
            FFF = exp((-attn_snow_wl1*SurfData.delzSnow)-(attn_ice_white_wl1*SurfData.delzWhiteIce))
                * (1.-exp(-attn_ice_blue_wl1*SurfData.delzBlueIce))/(attn_ice_blue_wl1*K_ice_blue);

            CCC = (1.-exp(-attn_snow_wl2*SurfData.delzSnow))/(K_snow*attn_snow_wl2);
            DDD = exp(-attn_snow_wl2*SurfData.delzSnow)
                *(1.-exp(-attn_ice_white_wl2*SurfData.delzWhiteIce))/(attn_ice_white_wl2*K_ice_white);
            GGG = exp((-attn_snow_wl2*SurfData.delzSnow)-(attn_ice_white_wl2*SurfData.delzWhiteIce))
                *(1.-exp(-attn_ice_blue_wl2*SurfData.delzBlueIce))/(attn_ice_blue_wl2*K_ice_blue);
            EEE = (K_snow*K_ice_white*K_ice_blue)
                / ((SurfData.delzSnow*K_ice_white*K_ice_blue)
                +(SurfData.delzBlueIce*K_snow*K_ice_blue)
                +(SurfData.delzWhiteIce*K_snow*K_ice_white));

            Q_iceout = ((Temp_melt-Temp_ice-(Q_shortwave*(f_sw_wl1*(AAA+BBB+FFF)
                           + f_sw_wl2*(CCC+DDD+GGG)))+Q_snowice)*EEE)+Q_shortwave+Q_whiteice;

            //# Now compare the balance between the surface fluxes - and
            //  iterate by bisection (note the loop back in the iteration)

            if (fabs(Q_iceout+Q_icemet) > 1.0 && fabs(T01_NEW-T01_OLD) > 0.001) {
                if ((Q_iceout+Q_icemet) < 0.0) T01_NEW = Temp_ice;
                if ((Q_iceout+Q_icemet) > 0.0) T01_OLD = Temp_ice;
                Temp_ice = (T01_NEW+T01_OLD)/2.0;
            } else {
                break;  // start melting
            }
        } // end while

        //# Reset
        T01_NEW =  50.0;
        T01_OLD = -50.0;

        //# Evaporation in water equivalents (so use rho0)
        if ( evap_from_file )
            SurfData.Evap = f_evap_ts_prop;
        else
            SurfData.Evap = Q_latentheat / Latent_Heat_Evap / rho0;

        //# Increment daily summaries
        SurfData.dailyEvap += (SurfData.Evap * noSecs * Lake[surfLayer].LayerArea);
        if (SurfData.delzSnow > 0.){
            SurfData.delzSnow += Q_latentheat/Latent_Heat_Evap/rho_snow*noSecs;
        } else if (SurfData.delzWhiteIce > 0.) {
            SurfData.delzWhiteIce += Q_latentheat/Latent_Heat_Evap/rho_ice_white*noSecs;
        } else {
            SurfData.delzBlueIce += Q_latentheat/Latent_Heat_Evap/rho_ice_blue*noSecs;
        }

        //---------------------------------------------------------------------+
        // ICE MELTING & FREEZING @ SURFACE
        // Now compare the ice/snow surface temperature with the melting
        // temperature - if it is above the melting temp, then adjust
        // the thickness of the surface ice/snow layer accordingly
        //---------------------------------------------------------------------+
        if (Temp_ice >= Temp_melt) {

            Temp_ice = Temp_melt;
            SatVap_surface = (1+(0.00972*Temp_ice)+(0.000042*pow((Temp_ice), 2)))*saturated_vapour(Temp_ice);
            //# Q_latentheat [W/m2] = CE * rho_air * latent heat * psychro const / air_presssure * windspeed * VPD
            Q_latentheat = -CE * rho_air * Latent_Heat_Evap * (0.622/p_atm) * WindSp * (SatVap_surface - MetData.SatVapDef);
            if (Q_latentheat > 0.0) Q_latentheat = 0.0;

            //# Q_sensibleheat [W/m2] = CH * rho_air * specific heat * windspeed * temp diff
            Q_sensibleheat = -CH * (rho_air * 1005.) * WindSp * (Temp_ice - MetData.AirTemp);

            Q_lw_out = -Stefan_Boltzman * eps_water * pow((Kelvin+Temp_ice), 4.0);

            if (LWModel  ==  LW_CC) {
                CloudCover = MetData.LongWave;
                Q_lw_in = (1.0 + 0.275 * CloudCover) * Stefan_Boltzman * pow((Kelvin+MetData.AirTemp), 4.0) *
                          (1.0 - 0.261 * exp(-0.000777E-4 * pow(MetData.AirTemp, 2.0)));
                Q_longwave = Q_lw_out + Q_lw_in;
            } else if (LWModel  ==  LW_IN) {
                Q_lw_in = MetData.LongWave;
                Q_longwave = Q_lw_out + Q_lw_in;
            } else if (LWModel  ==  LW_NET)
                Q_longwave = MetData.LongWave;

            Q_icemet = (Q_latentheat+Q_sensibleheat+Q_longwave+Q_rain);
            Q_iceout = Q_shortwave+Q_whiteice+((Temp_melt-Temp_ice-(Q_shortwave
               *(f_sw_wl1*(AAA+BBB+FFF)+f_sw_wl2*(CCC+DDD+GGG)))+Q_snowice)*EEE);

            //-----------------------------------------------------------------+
            // Now determine the new ice/snow/water surface temperature and
            // thickness based on the balance between the surface fluxes, h_ice
            // and the expression of the upward heat flux as given by
            // Patterson and Hamblin (1988). Note that melting can occur at the
            // surface and therefore water level will increase. Note also the
            // assumption that ice thickness won't change >1 cm in a timestep
            //-----------------------------------------------------------------+

            if (SurfData.delzSnow > 0.0) {
            //if (SurfData.delzSnow != 0.0)  //!MH Getting -ve Snow Height in lake.csv

                // If there is snow, melt that first
                if (rho_snow == 0.0) rho_snow = snow_rho_max;
                SurfData.dHt = (1/(Latent_Heat_Fusion*rho_snow))*(Q_icemet+Q_iceout)*noSecs;
                if ((SurfData.delzSnow-SurfData.dHt) < 0.0)SurfData.dHt = SurfData.delzSnow;
                SurfData.delzSnow = SurfData.delzSnow-SurfData.dHt;
                Lake[surfLayer].Height += SurfData.dHt*(rho_snow/Lake[surfLayer].Density);
                recalc_surface_salt();

            } else if (SurfData.delzWhiteIce > 0.){
                // Otherwise melt the white ice
                SurfData.dHt = (1/(Latent_Heat_Fusion*rho_ice_white))*(Q_icemet+Q_iceout)*noSecs;

                if ((SurfData.delzWhiteIce-SurfData.dHt) < 0.0)
                    SurfData.dHt = SurfData.delzWhiteIce;

                SurfData.delzWhiteIce -= SurfData.dHt;

                Lake[surfLayer].Height += SurfData.dHt*(rho_ice_white/Lake[surfLayer].Density);
                recalc_surface_salt();

            } else {
                // Lastly, melt the blue ice
                SurfData.dHt = (1/(Latent_Heat_Fusion*rho_ice_blue))*(Q_icemet+Q_iceout)*noSecs;
                if ((SurfData.delzBlueIce-SurfData.dHt) < 0.) {
                    SurfData.dHt = SurfData.delzBlueIce;
                }
                SurfData.delzBlueIce = SurfData.delzBlueIce-SurfData.dHt;

                Lake[surfLayer].Height += SurfData.dHt*(rho_ice_blue/Lake[surfLayer].Density);
                recalc_surface_salt();
            }   // end melting snow/white/blue ice

        } // end melting if

       // Increment daily heat budget (MJ/day)
       SurfData.dailyQe += Q_latentheat * Lake[surfLayer].LayerArea * noSecs;
       SurfData.dailyQh += Q_sensibleheat * Lake[surfLayer].LayerArea * noSecs;
       SurfData.dailyQlw += Q_longwave * Lake[surfLayer].LayerArea * noSecs;
       SurfData.Qe = Q_latentheat;
       SurfData.Qh = Q_sensibleheat;
       SurfData.Qlw = Q_longwave;
       SurfData.Q_net = Q_latentheat + Q_sensibleheat + Q_longwave;
    }

    /***************************************************************************
     * APPLY HEATING TO WATER LAYERS
     * Now look at the ice or water interface
     **************************************************************************/
    if (surfLayer > botmLayer) {
        AED_REAL npl_p1, npli;

        for (i = surfLayer-1; i >= botmLayer; i-- )
            Lake[i].Light = Lake[i+1].Light * exp(-Lake[i+1].ExtcCoefSW*LayerThickness[i+1]);

        //---------------------------------------------------------------------+
        // Into layer i goes QSW[i]-QSW(i-1) over the area common to layers
        // i and i-1 and QSW[i] over the rest of AREA[i]. The units
        //  of heat[i) are joules/sec; units of AREA[i] are m**2
        //---------------------------------------------------------------------+
        npl_p1 = NotPARLight_sm1;
        for (i = surfLayer-1; i >= (botmLayer+1); i--){
            flankArea = Lake[i].LayerArea - Lake[i-1].LayerArea;

            heat[i] = Lake[i-1].LayerArea * (Lake[i].Light - Lake[i-1].Light);
            heat[i] = heat[i] + flankArea * (Lake[i].Light - Lake[i-1].Light*(1.-sed_reflectivity[0]));

            //MH non PAR temporary hack:
            npli = npl_p1 * exp(-Lake[i+1].ExtcCoefSW*2*LayerThickness[i+1]) ;
            heat[i] = heat[i] + (npl_p1-npli)*Lake[i].LayerArea;
            npl_p1 = npli;
        }
        heat[botmLayer] = Lake[botmLayer].Light * (1-(1.-sed_reflectivity[0]) *
          exp(-Lake[botmLayer].ExtcCoefSW*LayerThickness[botmLayer])) * Lake[botmLayer].LayerArea;
        //heat[botmLayer] = heat[botmLayer] + npl_p1 * Lake[botmLayer].LayerArea;
    }

    //# Compute the temperature increase in non-surface layers over noSecs
    for (i = botmLayer; i <= surfLayer; i++) {
        AED_REAL layerHeat = (SPHEAT*Lake[i].Density*Lake[i].LayerVol*Lake[i].Temp);
        if ( (heat[i] < 0.) && ( layerHeat < (fabs(heat[i])*noSecs) ) )
            heat[i] = -(layerHeat/noSecs);

        if (fabs(heat[i]) >= 1E-20 && Lake[i].Density != 0.0 && Lake[i].LayerVol != 0.0)
            dTemp = heat[i]*noSecs/(SPHEAT*Lake[i].Density*Lake[i].LayerVol);
        else
            dTemp = 0.;

        Lake[i].Temp += dTemp;
    }

    if (littoral_sw) {
        dTemp = heat[offshoreLayer]*noSecs/(SPHEAT*offshoreDensity*offshoreVol);
        Lake[offshoreLayer].Temp += dTemp;

        dTemp = heat[onshoreLayer]*noSecs/(SPHEAT*onshoreDensity*onshoreVol);
        Lake[onshoreLayer].Temp += dTemp;
    }

    /**************************************************************************
     * Check if the ice melted in code above
     *************************************************************************/
    if ((SurfData.delzBlueIce+SurfData.delzWhiteIce) < min_ice_thickness && ice) {
      Lake[surfLayer].Height = Lake[surfLayer].Height
      + SurfData.delzBlueIce  * (rho_ice_blue/Lake[surfLayer].Density)
      + SurfData.delzWhiteIce * (rho_ice_white/Lake[surfLayer].Density)
      + SurfData.delzSnow     * (rho_snow/Lake[surfLayer].Density);

      recalc_surface_salt();

      ice = FALSE;
      SurfData.delzBlueIce  = 0.0;
      SurfData.delzWhiteIce = 0.0;
      SurfData.delzSnow    = 0.0;
    }


    /***************************************************************************
     * ICE MELTING & FREEZING @ BOTTOM
     * The change in ice thickness at the bottom can now be determined
     * since the temperature of the lake water is readjusted for the surface
     * heat exchange, including the influence of the Ice/Snow
     **************************************************************************/
    if (ice) {
        //---------------------------------------------------------------------+
        //# Both ablation and accretion of the ice can occur at the ice-water
        //  interface; The thickness change dht will be governed by the flux
        //  coming down from the surface to the water below the ice.
        //  see Patterson & Hamblin (1988)
        //---------------------------------------------------------------------+

        //# First calculate the flux through the ice
        AAA = (1. - exp( -attn_snow_wl1*SurfData.delzSnow ))/(K_snow*attn_snow_wl1);
        BBB = exp(-attn_snow_wl1*SurfData.delzSnow)*(1.-exp(-attn_ice_white_wl1*SurfData.delzWhiteIce))/ (attn_ice_white_wl1*K_ice_white);
        FFF = exp((-attn_snow_wl1*SurfData.delzSnow)-(attn_ice_white_wl1*SurfData.delzWhiteIce))*
                                      (1.-exp(-attn_ice_blue_wl1*SurfData.delzBlueIce))/(attn_ice_blue_wl1*K_ice_blue);
        CCC = (1.-exp(-attn_snow_wl2*SurfData.delzSnow))/(K_snow*attn_snow_wl2);
        DDD = exp(-attn_snow_wl2*SurfData.delzSnow)*(1.-exp(-attn_ice_white_wl2*SurfData.delzWhiteIce))/ (attn_ice_white_wl2*K_ice_white);
        GGG = exp((-attn_snow_wl2*SurfData.delzSnow)-(attn_ice_white_wl2*SurfData.delzWhiteIce))*
                                      (1.-exp(-attn_ice_blue_wl2*SurfData.delzBlueIce))/(attn_ice_blue_wl2*K_ice_blue);
        EEE = (K_snow*K_ice_white*K_ice_blue)/((SurfData.delzSnow*K_ice_white*K_ice_blue)+(SurfData.delzBlueIce*K_snow*K_ice_blue)+
                                       (SurfData.delzWhiteIce*K_snow*K_ice_white));

        Q_iceout = ((Temp_melt-Temp_ice-(Q_shortwave*(f_sw_wl1*(AAA+BBB+FFF)+f_sw_wl2*(CCC+DDD+GGG))) +Q_snowice)*EEE)+Q_shortwave+Q_whiteice;

        Q_icewater = Q_iceout-(Q_shortwave*f_sw_wl1*(1.-exp(-(attn_snow_wl1*SurfData.delzSnow+attn_ice_blue_wl1*SurfData.delzBlueIce+attn_ice_white_wl1*SurfData.delzWhiteIce))))
                             -(Q_shortwave*f_sw_wl2*(1.-exp(-(attn_snow_wl2*SurfData.delzSnow+attn_ice_blue_wl2*SurfData.delzBlueIce+attn_ice_white_wl2*SurfData.delzWhiteIce))));

        //---------------------------------------------------------------------+
        //# Now determine the flux through the lake water below the ice. For
        //  the temperature flux at the ice-water interface an exponential
        //  decrease is used for fickian diffusion so use a gaussian
        // distribution and adjust the diffuse accordingly - see Farmer (1978)
        //  and eq 22 of Rogers et al, noting that dz has been adjusted
        //---------------------------------------------------------------------+
        Q_watermol = -K_water*(Temp_melt-Lake[surfLayer].Temp) / 0.039;

        //# Heat transfer due to ice underflow
        //  if (underFlow) Q_underflow = CSEN*(Lake[surfLayer].Density - rho0)*SPHEAT*U_FLOW*(Lake[surfLayer].Temp-Temp_melt);
        //  else
        Q_underflow = 0.0; //  Currently not being used as underFlow == .FALSE.

        //---------------------------------------------------------------------+
        //# Heat transfer to the water
        //---------------------------------------------------------------------+
        Q_surflayer = Q_watermol + Q_underflow + Q_latent_ice;

        //---------------------------------------------------------------------+
        //# Now we determine the amount of ablation or accretion of ice as
        //  given by qf-qw.  Once the ice has melted or formed, assume fluxes
        //  are in equilibrium. Correction for thermal contraction given
        //  by ratios of d is constant and area can be considered as constant
        //  as change in depth and time step is small
        //---------------------------------------------------------------------+
        SurfData.dHt = (Q_icewater-Q_surflayer)*noSecs/(Latent_Heat_Fusion*rho_ice_blue);
        if (SurfData.dHt < -1.*SurfData.delzBlueIce) {
            if (SurfData.dHt < -1. * (SurfData.delzBlueIce+SurfData.delzWhiteIce) )
                SurfData.dHt = -1.*(SurfData.delzBlueIce+SurfData.delzWhiteIce);
            SurfData.delzWhiteIce = SurfData.delzWhiteIce+SurfData.delzBlueIce+SurfData.dHt;
            SurfData.delzBlueIce = 0.0;
        } else
            SurfData.delzBlueIce = SurfData.delzBlueIce+SurfData.dHt;

        //# Adjust water temperature for ice/water exchange, determine latent heat
        Q_latent_ice = Latent_Heat_Fusion*rho_ice_blue*SurfData.dHt/noSecs;

        //Lake[surfLayer].Temp = Lake[surfLayer].Temp+(((K_water*((Temp_melt-Lake[surfLayer].Temp)/0.039))*Lake[surfLayer].LayerArea*
        //                         noSecs)+Q_latent_ice)/(SPHEAT*Lake[surfLayer].Density*Lake[surfLayer].LayerVol);

        Lake[surfLayer].Temp = Lake[surfLayer].Temp
                +((-Q_watermol + Q_latent_ice)*Lake[surfLayer].LayerArea*noSecs)
                /(SPHEAT*Lake[surfLayer].Density*Lake[surfLayer].LayerVol);

        Lake[surfLayer].Height = Lake[surfLayer].Height
                               -SurfData.dHt*(rho_ice_blue/Lake[surfLayer].Density);
        recalc_surface_salt();
    }

    /**************************************************************************
     * SEDIMENT HEATING
     * Sediment "heating" factor now applied to any layer based on which zone
     * it sits above.
     **************************************************************************/
    if (sed_heat_sw) {
        //# Input of heat from the sediments.
        kDays = day_of_year(jday);
        ZSED = sed_temp_depth;
        KSED = sed_heat_Ksoil;

        if (benthic_mode == 1) {
            //# Apply the same sediment heating parameters across flanks of all layers
            kDays = day_of_year(jday);
            TYEAR = sed_temp_mean[0] + sed_temp_amplitude[0] * cos(((kDays-sed_temp_peak_doy[0])*2.*Pi)/365.);
            for (i = botmLayer+1; i <= surfLayer; i++) {
                Lake[i].Temp += ((KSED * (TYEAR - Lake[i].Temp) / ZSED) *
                             (Lake[i].LayerArea - Lake[i-1].LayerArea) * noSecs) /
                                              (SPHEAT*Lake[i].Density * Lake[i].LayerVol);
            }
            Lake[botmLayer].Temp += ((KSED * (TYEAR - Lake[botmLayer].Temp) / ZSED) *
                                   Lake[botmLayer].LayerArea * noSecs) /
                                          (SPHEAT*Lake[botmLayer].Density * Lake[botmLayer].LayerVol);
        } else  if (benthic_mode == 2) {
            //# Apply the sediment zone specific heating parameters to overlying
            //  layers. First find which layers correspond to which zone
            for (i = botmLayer; i <= surfLayer; i++) {
                layer_zone[i] = 0;
                for (z = 1; z < n_zones; z++) {
                    if ((Lake[i].Height < theZones[z].zheight) && (Lake[i].Height > theZones[z-1].zheight))
                        layer_zone[i] = z;
                }
            }
            //# Now compute layer-specifc sed heating and increment temperature
            if ( sed_heat_model == 2 ){
              soil_dt = noSecs;   /* per-step timestep for the dynamic soil solver (s) */
//            memset(sed_depths, 0, sizeof(AED_REAL)*sed_layers);
//            memset(sed_vwc, 0, sizeof(AED_REAL)*sed_layers);
//            memset(sed_temps, 0, sizeof(AED_REAL)*sed_layers);

              //# Volume-weighted zone-average water temperature: the soil model's
              //# surface boundary condition.  Computed here rather than in the WQ
              //# coupler because only the legacy 'aed' coupler's copy_to_zone fills
              //# ztemp; the 'api' coupler and wq_calc=.false. runs would leave it 0.
              //# Zones with no wet layers keep their previous ztemp.
              for (z = 0; z < n_zones; z++) {
                  AED_REAL tsum = 0., tvol = 0.;
                  for (i = botmLayer; i <= surfLayer; i++) {
                      if (layer_zone[i] == z) {
                          tsum += Lake[i].Temp * Lake[i].LayerVol;
                          tvol += Lake[i].LayerVol;
                      }
                  }
                  if (tvol > 0.) theZones[z].ztemp = tsum / tvol;
              }

              for (z = 0; z < n_zones; z++) {
                  // call the dynamic soil/sediment temperature model
                  /*
                  SoilTemp( &theZones[z].n_sed_layers,
                             sed_depths,
                             sed_vwc,
                             theZones[z].ztemp,
                             sed_temps );
                            // &soil_heat_flux );
                  */
#ifdef AED
                  // ZSoilTemp advances this zone's sediment temperature profile
                  // and stores the sediment-water interface heat flux (W/m2) in
                  // theZones[z].heatflux.
                  if ( p_wq_ZSoilTemp != NULL ) ZSoilTemp(&theZones[z]);
#endif
                  // flux heat from the soil into the water, if the layer is over z
                  for (i = botmLayer+1; i <= surfLayer; i++) {
                    if (layer_zone[i] == z){
                      Lake[i].Temp += theZones[z].heatflux
                              * ((Lake[i].LayerArea - Lake[i-1].LayerArea) * noSecs)
                              / (SPHEAT * Lake[i].Density * Lake[i].LayerVol);
                    }
                  }
               }
            } else if ( sed_heat_model == 1 ){
              for (i = botmLayer+1; i <= surfLayer; i++) {


                                if(i == 5){
                    //fprintf(stdout, "Lake[i].Temp %f soil_heat_flux, %f Lake[i].LayerArea, %f Lake[i-1].LayerArea, %f Lake[i].Density, %f Lake[i].LayerVol %f\n",  Lake[i].Temp, soil_heat_flux, Lake[i].LayerArea,Lake[i-1].LayerArea, Lake[i].Density, Lake[i].LayerVol);
                  }
                TYEAR = sed_temp_mean[layer_zone[i]]
                        + sed_temp_amplitude[layer_zone[i]]
                        * cos(((kDays-sed_temp_peak_doy[layer_zone[i]])*2.*Pi)/365.);
                soil_heat_flux = KSED * (TYEAR - Lake[i].Temp) / ZSED;
                Lake[i].Temp += soil_heat_flux
                              * ((Lake[i].LayerArea - Lake[i-1].LayerArea) * noSecs)
                              / (SPHEAT * Lake[i].Density * Lake[i].LayerVol);
              }

              TYEAR = sed_temp_mean[0] + sed_temp_amplitude[0] * cos(((kDays-sed_temp_peak_doy[0])*2.*Pi)/365.);
              Lake[botmLayer].Temp += ((KSED * (TYEAR - Lake[botmLayer].Temp) / ZSED) *
                                      Lake[botmLayer].LayerArea * noSecs) /
                                      (SPHEAT * Lake[botmLayer].Density*Lake[botmLayer].LayerVol);
            }
        }
//        if (littoral_sw) {
//            TYEAR = sed_temp_mean[2] + sed_temp_amplitude[2] * cos(((kDays-sed_temp_peak_doy[2])*2.*Pi)/365.);
//            Lake[onshoreLayer].Temp += ((KSED * (TYEAR - Lake[onshoreLayer].Temp) / ZSED) * onshoreVol * noSecs) /
//                                            (SPHEAT * onshoreDensity * onshoreVol);
    }

    /**************************************************************************
     * SURFACE MASS FLUXES (NO ICE COVER PRESENT)
     * Precipitation, evaporation in the absence of ice
     **************************************************************************/
    if (!ice) {
        AED_REAL evap = MAX( SurfData.Evap*noSecs,-0.9*Lake[surfLayer].Height );

        AED_REAL evapvol = MIN( evap,zero ) * Lake[surfLayer].LayerArea;

        AED_REAL rainvol = MAX( MetData.Rain,zero )
                             * (noSecs / SecsPerDay) * Lake[surfLayer].LayerArea;

        AED_REAL snowvol = MAX( MetData.Snow,zero ) * (1./10.)
                             * (noSecs / SecsPerDay) * Lake[surfLayer].LayerArea;

        SurfData.dailyEvap += evapvol;
        SurfData.dailyRain += rainvol;
        SurfData.dailySnow += snowvol;

        //---------------------------------------------------------------------+
        //# Add rain directly to surface layer height
        //---------------------------------------------------------------------+
        Lake[surfLayer].Height += rainvol / Lake[surfLayer].LayerArea;

        //---------------------------------------------------------------------+
        //# Add snow directly to surface layer height, if there is no ice.
        //  If there is ice, snow will be handled in the next block
        //  Use 1:10 rule for snow water equivalent (SWE)
        //---------------------------------------------------------------------+
        Lake[surfLayer].Height += snowvol / Lake[surfLayer].LayerArea;

        //---------------------------------------------------------------------+
        //# Dilution by rainfall/snowfall. Note that rain can deposit
        //  consituents at a rate dependent on the composition of rainfall. # TBC
        //---------------------------------------------------------------------+
        Lake[surfLayer].Temp = combine(Lake[surfLayer].Temp,
                                       Lake[surfLayer].LayerVol,
                                       Lake[surfLayer].Density,
                                       MetData.AirTemp, rainvol+snowvol,
                                       calculate_density(MetData.AirTemp, MAX(zero,salt_fall)));
        Lake[surfLayer].Salinity = combine(Lake[surfLayer].Salinity,
                                       Lake[surfLayer].LayerVol,
                                       Lake[surfLayer].Density,
                                       MAX(zero,salt_fall), rainvol+snowvol,
                                       calculate_density(MetData.AirTemp, MAX(zero,salt_fall)));
        for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
            _WQ_Vars(wqidx, surfLayer) = combine_vol(_WQ_Vars(wqidx, surfLayer),
                                                     Lake[surfLayer].LayerVol,
                                                     zero, rainvol+snowvol);

//      resize_internals(1, surfLayer);  // recompute surflayer volume

        //---------------------------------------------------------------------+
        //# Evaporation and evapo-concentration, as evaporation leaves consituents
        //---------------------------------------------------------------------+
        Lake[surfLayer].Height += evapvol / Lake[surfLayer].LayerArea;

        recalc_surface_salt();
        //recalc_surface_wq();
        resize_internals(1, surfLayer);  // recompute surflayer volume
    }

    //# Recalculate densities
    for (i = botmLayer; i <= surfLayer; i++)
        Lake[i].Density = calculate_density(Lake[i].Temp,Lake[i].Salinity);

    //if (littoral_sw) {
    //    onshoreDensity = calculate_density(Lake[onshoreLayer].Temp,Lake[surfLayer].Salinity);
    //    offshoreDensity = calculate_density(Lake[onshoreLayer].Temp,Lake[surfLayer].Salinity);
    //  }


    //printf(">min_ice_thickness = %10.5f\n",min_ice_thickness);

    //# Check and set ice cover flag. To ensure ice is due a moving average
    //  of surface temperature is computed and used to assess ice-on event
    //AED_REAL dt_ice_avg = MAX(0.5,noSecs/SecsPerDay);
    AED_REAL dt_ice_avg = MAX(dt_iceon_avg,noSecs/SecsPerDay);
    AvgSurfTemp = AvgSurfTemp * (1 - (noSecs/SecsPerDay)/dt_ice_avg)
                       + Lake[surfLayer].Temp * (noSecs/SecsPerDay)/dt_ice_avg ;

    if (AvgSurfTemp <= avg_surf_temp_thres && SurfData.delzBlueIce == 0.0 && Lake[surfLayer].Height>0.1) {
        // Start a new blue ice layer
        ice                     = TRUE;
        SurfData.delzBlueIce    = min_ice_thickness; //0.05;
        SurfData.delzWhiteIce   = 0.0;
        SurfData.delzSnow       = 0.0;
        Lake[surfLayer].Height -= min_ice_thickness * (rho_ice_blue/Lake[surfLayer].Density);

        recalc_surface_salt();
    }
    if ((SurfData.delzBlueIce+SurfData.delzWhiteIce) < min_ice_thickness && ice) {
        Lake[surfLayer].Height = Lake[surfLayer].Height
                + SurfData.delzBlueIce  * (rho_ice_blue/Lake[surfLayer].Density)
                + SurfData.delzWhiteIce * (rho_ice_white/Lake[surfLayer].Density)
                + SurfData.delzSnow     * (rho_snow/Lake[surfLayer].Density);

        recalc_surface_salt();

        ice = FALSE;
        SurfData.delzBlueIce  = 0.0;
        SurfData.delzWhiteIce = 0.0;
        SurfData.delzSnow    = 0.0;
    }
    SurfData.RhoSnow = rho_snow;

    free(LayerThickness);  free(heat);  free(layer_zone);

    free(gx);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
AED_REAL calculate_qsw(int kDays,          // Days since start of year for yesterday
                       int mDays,          // Days since start of year for today
                       int iclock,         // No. seconds into the day
                       AED_REAL Latitude,  // Latitude
                       AED_REAL SWOld,     // Total solar rad at the surface for yesterday
                       AED_REAL ShortWave, // Total solar rad at the surface for today
                       AED_REAL WindSp)    // Current wind speed (adjusted for possible ice)
{
    AED_REAL RADTIM0;  //# Start of time step as fraction of day in radians
    AED_REAL RADTIM1;  //# End of time step as fraction of day in radians
    AED_REAL SOLAR0;   //# Solar angle yesterday
    AED_REAL SOLAR1;   //# Solar angle today
    AED_REAL SOLARW;   //# Sun's arc from sunrise to sunset for yesterday
    AED_REAL SOLARY;   //# Sun's arc from sunrise to sunset for today
    AED_REAL SW0;      //# Total solar radiation at the surface for yesterday after ALBEDO
    AED_REAL SW1;      //# Total solar radiation at the surface for today after ALBEDO
    AED_REAL TD0;      //# Day length for yesterday
    AED_REAL TD1;      //# Day length for today
    int jClock;        //# End of the current time step

    AED_REAL Albedo0;  //# Photosynthetically active radiation scattered previous day
    AED_REAL Albedo1;  //# Photosynthetically active radiation scattered present day
    AED_REAL Temp_ice; //# Temperature of the top layer of the ice cover model [0C]
    AED_REAL lQSW;     //# Shortwave entering the uppermost water layer

/*----------------------------------------------------------------------------*/

    lQSW = 0.;
    RADTIM1 = 0.;
    RADTIM0 = 0.;
    Albedo0 = 0.;
    Albedo1 = 0.;

    //# Initialize Temp_ice to the Air Temperature
    Temp_ice = MetData.AirTemp;

    //# Get current Solar Zenith Angle
    AED_REAL sza  = MIN(90.0,zenith_angle(Longitude, Latitude, mDays, iclock, timezone_r));  //degrees
    AED_REAL csza = cos(Pi * sza/180);

    //# Check ice flag, determine resultant albedo
    if (ice) {
        //# Albedo of ice cover depends on ice thickness, cover type & temperature
        //  This algotihm is based on Table 1 of Vavrus et al (1996).
        if ((SurfData.delzBlueIce+SurfData.delzWhiteIce) > 0.55) {
            if (Temp_ice <= -5.)                       Albedo0 = 0.6;
            else if (Temp_ice > -5.0 && Temp_ice < 0.) Albedo0 = 0.44 - 0.032*Temp_ice;
            else if (Temp_ice >= 0.)                   Albedo0 = 0.44;
        } else{
            Albedo0 = (albedo_mean + 0.44 *
                  pow((SurfData.delzBlueIce+SurfData.delzWhiteIce-min_ice_thickness), 0.28));
        }
        // Adjust albedo based on multiplicative factor, limited to 1
        Albedo0 = MIN( 1.0, snow_albedo_factor * Albedo0 );

        if (SurfData.delzSnow > 0.0) {
            if (Temp_ice <= -5.)                      Albedo1 = 0.7;
            else if (Temp_ice > -5. && Temp_ice < 0.) Albedo1 = 0.5 - 0.04*Temp_ice;
            else if (Temp_ice >= 0.)                  Albedo1 = 0.5;

            if (SurfData.delzSnow < 0.1)
                Albedo0 = Albedo1-(((0.1-SurfData.delzSnow)/0.1)*(Albedo1-Albedo0));
            else
                Albedo0 = Albedo1;

            // Adjust albedo based on multiplicative factor, limited to 1
            Albedo0 = MIN( 1.0, snow_albedo_factor * Albedo0 );

        }
        Albedo1 = Albedo0;

    } else {
        //# Open water conditions
        Albedo0 = albedo_mean;
        Albedo1 = albedo_mean;

        switch (albedo_mode) {
            case 1: // simple daily albedo calc as in Hamilton and Schladow (1997)
                if (Latitude > 5*deg2rad) { //# Lake is in the northern hemisphere
                    Albedo0 = albedo_mean -
                         albedo_amplitude * sin( (two_Pi * kDays/365) - halfPi);
                    Albedo1 = albedo_mean -
                         albedo_amplitude * sin( (two_Pi * mDays/365) - halfPi);
                } else if(Latitude < -5*deg2rad) { //# Lake is in the southern hemisphere
                    Albedo0 = albedo_mean -
                         albedo_amplitude * sin( (two_Pi * kDays/365) + halfPi);
                    Albedo1 = albedo_mean -
                         albedo_amplitude * sin( (two_Pi * mDays/365) + halfPi);
                }
                break;
            case 2: { // Briegleb et al. (1986), B scheme in Scinocca et al 2006
                    Albedo1 = ((2.6/(1*pow(csza, 1.7)+0.065)) + (15*(csza-0.1)*(csza-0.5)*(csza-1)))/100;
                    Albedo0 = Albedo1;
                }
                break;
            case 3: { // Yajima and Yamamoto (2014)
                    Albedo1 = MAX(0.02, 0.001 * MetData.RelHum * pow(1-csza, 0.33) -
                                        0.001 * WindSp * pow(1-csza, -0.57) -
                                        0.001 * 6 * pow(1-csza, 0.829));
                    Albedo0 = Albedo1;
                }
                break;
            case 4: // simple daily albedo calc from Cogley (1979)
                Albedo1 = grischenko_albedo(Latitude,mDays);
                Albedo0 = Albedo1;
                break;
            default : break;
        }
    }

    //# Increment the daily albedo for this time-step, if within daylight hours
    //if (sza<89.) SurfData.albedo += Albedo1;
    if (Albedo1<SurfData.albedo) SurfData.albedo = Albedo1; // noon albedo

    //# Now either assign hourly insolation flux from BC value if sub-daily,
    //  or compute if BC is daily only, based on time of day
    if ( subdaily )
        lQSW = ShortWave * (1.0-Albedo1);
    else {
        //# Determine the daylength (hours) from the latitude (entered as -ve
        //  for and the number of days since the start of the year
        SOLAR0 = -23.45 * sin((284 + kDays) * 2.0 * Pi / 365.0);
        SOLAR1 = -23.45 * sin((284 + mDays) * 2.0 * Pi / 365.0);

        //# Determine the angle of the sun's arc from sunrise to sunset after
        //  converting to radians, first determine solar declination
        SOLAR0 = SOLAR0 * Pi / 180;
        SOLAR1 = SOLAR1 * Pi / 180;
        SOLARW = 2.0 * acos(-tan(Latitude) * tan(SOLAR0));
        SOLARY = 2.0 * acos(-tan(Latitude) * tan(SOLAR1));

        //# Calculate the daylength (Secs) using angular velocity 15 degrees/hr
        TD0 = 3600 * SOLARW / (15 * Pi/180);
        TD1 = 3600 * SOLARY / (15 * Pi/180);

        //# Set the end of the GLM time step
        jClock = iclock + noSecs;

        //# Convert the start and end of the day to radians
        if (iclock < 0.5*TD0)
            RADTIM0 = iclock*Pi/TD0+Pi/2;
        else if (iclock >= 0.5*TD0 && iclock <= SecsPerDay-0.5*TD0)
            RADTIM0 = Pi;
        else if (iclock > SecsPerDay-0.5*TD0)
            RADTIM0 = (0.5*TD0+iclock-SecsPerDay)*Pi/TD0;

        if (jClock < 0.5*TD1)
            RADTIM1 = jClock*Pi/TD1+Pi/2;
        else if (jClock >= 0.5*TD1 && jClock <= SecsPerDay-0.5*TD1)
            RADTIM1 = Pi;
        else if (jClock > SecsPerDay-0.5*TD1)
            RADTIM1 = (0.5*TD1+jClock-SecsPerDay)*Pi/TD1;

        SW0 = SWOld*(1.0-Albedo0) * 86.4;
        SW1 = ShortWave*(1.0-Albedo1) * 86.4;

        //# Determine the area under a sinusoidal curve and convert to a rate
        if (iclock <= (SecsPerDay - 0.5 * TD0) && jClock > (SecsPerDay - 0.5 * TD1) )
            lQSW = (0.5*SW0*1000.0*(cos(-1.0*RADTIM0)+1.0) + 0.5*SW1*1000.0*(1.0-cos(-1.0*RADTIM1)))/noSecs;
        else if (iclock <= (SecsPerDay - 0.5 * TD0) && jClock <= (SecsPerDay - 0.5 * TD1))
            lQSW = (0.5*SW0*1000.0*(cos(-1.0*RADTIM0)-cos(-1.0*RADTIM1)))/noSecs;
        else if (iclock > (SecsPerDay - 0.5 * TD0) && jClock > (SecsPerDay - 0.5 * TD1))
            lQSW = (0.5*SW1*1000.0*(cos(-1.0*RADTIM0)-cos(-1.0*RADTIM1)))/noSecs;
    }

    if (lQSW < 0.1) lQSW = 0.0;

    return lQSW;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * MELTWATER DILUTION
 ******************************************************************************/
void recalc_surface_salt()
{
    AED_REAL OldVol, AddDensity, WaterMass;

    OldVol = Lake[surfLayer].LayerVol;

    resize_internals(1, surfLayer);

    AddDensity = calculate_density(Lake[surfLayer].Temp, zero);

    WaterMass = OldVol * (Lake[surfLayer].Density - AddDensity) +
                          Lake[surfLayer].LayerVol * AddDensity;

    Lake[surfLayer].Salinity = Lake[surfLayer].Salinity *
                              (Lake[surfLayer].Density / WaterMass) * OldVol;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/




/******************************************************************************
 * ATMOSPHERIC STABILITY
 ******************************************************************************/

#define K_air    0.0280
#define D_air    0.0000214

#define WIND_HEIGHT 10.0
#define HUMIDITY_HEIGHT 10.0

#define SIGN(a,b)  ( ((b) >= 0.) ? fabs(a) : -fabs(a) )

static AED_REAL psi_m(AED_REAL zL);
static AED_REAL psi_hw(AED_REAL zL);

/******************************************************************************/
int atmos_stability(      AED_REAL *Q_latentheat,
                          AED_REAL *Q_sensible,
                          AED_REAL  wind_speed,
                          AED_REAL  temp_water,
                          AED_REAL  temp_air,
                          AED_REAL  humidity_surface,
                          AED_REAL  humidity_altitude,
                          AED_REAL  rho_air,
                          AED_REAL  rho_o,
                          AED_REAL  p_atm,
                          AED_REAL  latent_heat_vap,
                          AED_REAL *coef_wind_drag,
                          AED_REAL *coef_wind_chwn,
                          AED_REAL *zonL                )
{

    AED_REAL U10, U_sensM, U_sensH, Ux;
    AED_REAL zL, L, zL0, z0, zS, G1, G2, G3, G5, G6, Ldenom;
    AED_REAL CDN10, CHWN10, CDN4, CDN3, CHWN, rCDN, CD4, CHW;
    AED_REAL P1, P2, P4;
    AED_REAL T_virt, dT, dq;
    AED_REAL SH, LH, momn_flux;
    AED_REAL alpha_e, alpha_h, visc_k_air;
    AED_REAL Q_latentheat_still, Q_sensible_still;

    int atmos_count, atmos_status;

    AED_REAL vonK = 0.4100;    // von Karman's constant
    AED_REAL c_z0 = 0.0001;    // Default surafce roughness
    AED_REAL zL_MAX = -15.0;   // Bound the iteration (eg. 15 for 10m, 3 for 2m)

/*----------------------------------------------------------------------------*/
    atmos_status = 0;

    //# Some initial windspeed checks
    U_sensM = wind_speed;
    if (fabs(WIND_HEIGHT-10.0) > 0.5)
        U10 = wind_speed * (log(10.0/c_z0)/log(WIND_HEIGHT/c_z0));
    else
        U10 = wind_speed;

    CHWN10 = CH;

    //# Surface temperature and humidity gradients
    dT = temp_water - temp_air;
    dq = humidity_surface - humidity_altitude;

    visc_k_air = (1./rho_air)*(4.94e-8*temp_air + 1.7184e-5); //>m2/s

    //# First calculate still air approximations (free convenction) and use
    //  this as a minimum limit to compare with forced convection value
    //  For the fluxes to still air, see TVA Section 5.311 and 5.314

    //# Still air flux computation
    if (rho_air - rho_o > zero) {
        alpha_h = 0.137*0.5*K_air * pow( (g* fabs(rho_air-rho_o)/(rho_air*visc_k_air*D_air)), (1/3.0));
        alpha_e = alpha_h/cp_air;
        Q_sensible_still = -alpha_h * dT;
        Q_latentheat_still = -alpha_e * latent_heat_vap * dq ;
        // printf(">alpha_e = %10.5f\n",alpha_e);
    } else {
        Q_sensible_still = zero;
        Q_latentheat_still = zero;
    }

    if ( atm_stab == 2 ) {
      // Assign free vs forced
      if (Q_sensible_still < *Q_sensible)
         *Q_sensible = Q_sensible_still;
      if (Q_latentheat_still < *Q_latentheat)
         *Q_latentheat = Q_latentheat_still;
      *zonL = zero;
      atmos_status = -2;

      return atmos_status;
    }

    /////////////

    printf("top = %10.5f\n",0.137*0.5*K_air);
    printf("bit = %10.5f\n",pow( (g* fabs(rho_air-rho_o)/(rho_air*visc_k_air*D_air)), (1/3.0)));
    printf("visc_k_air = %10.5f\n",visc_k_air);
    printf("D_air = %10.5f\n",D_air);
    printf("K_air = %10.5f\n",K_air);
    printf("dq = %10.5f\n",dq);
    printf("*Q_sensible_still = %10.5f\n",Q_sensible_still);
    printf("*Q_latentheat_still = %10.5f\n",Q_latentheat_still);

    /////////////



    //# Now check windspeed
    CDN10 = 0.001;
    if (U_sensM<0.01) {
        *Q_sensible = Q_sensible_still;
        *Q_latentheat = Q_latentheat_still;
    } else {
        // Neutral Drag Coefficient is a function of windspeed @ 10m
        //Option1
          //  if (U10 > 5.0)
          //      CDN10 = (1.0 + 0.07*(U10-5.0))/1000.0;
        //Option2
        CDN10 = 1.92E-7 * U10*U10*U10 + 0.00096;

        //Check
        if (CDN10>0.0025) CDN10 = 0.0025;
    }


    //# Charnock computation of roughness, from Ux estimate.
    //  0.00001568 is kinematic viscosity ?
    //  0.012 is Charnock constant (alpha)
    Ux = sqrt(CDN10  * U_sensM * U_sensM);
    z0 = (0.012*Ux*Ux/g) + 0.11*visc_k_air/Ux;
    CDN10 = pow(vonK/log(10./z0),2.0);

    //# Estimate surface roughness lengths
    z0 = 10.0/(exp(vonK/sqrt(CDN10)));
    zS = 10.0/(exp(vonK*vonK/(CHWN10*log(10.0/z0))));

    //# Height correction factors
    G1 = log(10.0/z0);
    G2 = log(10.0/zS);
    G3 = log(HUMIDITY_HEIGHT/zS);
    G5 = log(HUMIDITY_HEIGHT/z0);
    G6 = log(WIND_HEIGHT/z0);

    CDN4 = CDN10*(G1*G1)/(G6*G6);    // Scale down to sensor heights
    CDN3 = CDN10*(G1*G1)/(G5*G5);
    CHWN = CHWN10*(G1*G2)/(G5*G3);
    CD4  = CDN4;                     // Initialize
    CHW  = CHWN;

    //# Windspeed at the humidity sensor height
    U_sensH = U_sensM*(G5/G6);

    //# Virtual air temperature
    T_virt = (temp_air+Kelvin) * (1.0 + 0.61*humidity_altitude);

    //# Heat fluxes based on bulk transfer (forced convection)
    SH = CHW * rho_air * cp_air * U_sensH  * dT;  //> W/m2
    LH = CHW * rho_air * U_sensH * dq;            //>

    //# Friction velocity
    momn_flux = CD4 * rho_air * U_sensM*U_sensM;
    Ux = sqrt(momn_flux/rho_air);

    //# Monin-Obukhov Length
    Ldenom = (vonK * g * ((SH/cp_air) + 0.61*(temp_air+Kelvin)*LH));
    if (fabs(Ldenom) < 1e-5) {
        zL = SIGN(zL_MAX,dT);
        L = HUMIDITY_HEIGHT/zL;
    } else {
        L = -(rho_air*Ux*Ux*Ux*T_virt) / Ldenom;
        zL = HUMIDITY_HEIGHT/L;
    }

    printf("U_sensM = %10.5f\n",U_sensM);
    printf("L = %10.5f\n",L);
    printf("zL = %10.5f\n",zL);

    //# Start iterative sequence for heat flux calculations
    atmos_count = 1;
    atmos_status = 1;
    zL0 = zero;
    while ((fabs(zL - zL0) >= 0.0001) ){ // && (fabs(zL) <= fabs(zL_MAX)+1)) {
        zL0 = zL;
        zL = WIND_HEIGHT/L;

        if (++atmos_count>=100){
            atmos_status = -1;
            break;
        }

        // Calculate drag coefficient, CD
        P4 = psi_m(zL);
        rCDN = sqrt(CDN4);
        CD4 = CDN4/(1.0+CDN4*(P4*P4 - 2.0*vonK*P4/rCDN)/(vonK*vonK));

        // Calculate Humdity/Temp coefficient, CHW
        zL = HUMIDITY_HEIGHT/L;

        P1 = psi_m(zL);
        P2 = psi_hw(zL);
        rCDN = sqrt(CDN3);
        CHW = CHWN/(1.0 + CHWN*(P1*P2 - (vonK*P2/rCDN)
                            - vonK*P1*rCDN/CHWN)/(vonK*vonK));

        // Recalculate heat and momn fluxes
        SH = CHW * rho_air * cp_air * U_sensH * dT;
        LH = CHW * rho_air * U_sensH  * dq;
        momn_flux = CD4 * rho_air * U_sensM*U_sensM;

        // Recalculate friction velocity
        Ux = sqrt(momn_flux/rho_air);

        // Recalculate Monin - Obukhov length
        L = -rho_air *Ux*Ux*Ux * T_virt / (vonK * g
                                * ((SH/cp_air) + 0.61*(temp_air+Kelvin)*LH));

        //printf("L = %10.5f\n",L);
        //printf("CHW = %10.5f\n",CHW);
        //printf("CD4 = %10.5f\n",CD4);

      //  if (fabs(L) < 0.5)
      //      L = SIGN(1.0e-20,dT);
        zL = HUMIDITY_HEIGHT/L;
    } // enddo

    if (atmos_status==-1)
       return atmos_status;


    //# Last calculation ... but 1st, check for high values
    if (fabs(zL)>fabs(zL_MAX))
        zL = SIGN(fabs(zL_MAX),zL);
    else
        zL = WIND_HEIGHT/L;

    P4 = psi_m(zL);
    rCDN = sqrt(CDN4);
    CD4 = CDN4/(1.0+CDN4*(P4*P4 - 2.0*vonK*P4/rCDN)/(vonK*vonK));
    zL = zL*HUMIDITY_HEIGHT/WIND_HEIGHT;
    P1 = psi_m(zL);
    P2 = psi_hw(zL);
    rCDN = sqrt(CDN3);
    CHW = CHWN/(1.0 + CHWN*(P1*P2 - (vonK*P2/rCDN)
            - vonK*P1*rCDN/CHWN)/(vonK*vonK));

    *coef_wind_drag = CD4;
    *coef_wind_chwn = CHW;

    *Q_sensible = -CHW * rho_air * cp_air * U_sensH * dT;
    *Q_latentheat = -CHW * rho_air * U_sensH * dq * latent_heat_vap;

    printf("*atmos_count = %10d\n",atmos_count);
    printf("*CHW = %10.6f\n",CHW);
    printf("*CD4 = %10.6f\n",CD4);
    printf("zL = %10.5f\n",zL);
    printf("dT = %10.5f\n",dT);
    printf("dq = %10.5f\n",dq);
    printf("*Q_sensible = %10.5f\n",*Q_sensible);
    printf("*Q_latentheat = %10.5f\n",*Q_latentheat);

    *zonL = zL;
    if ( atm_stab == 3 )
       return atmos_status;

    //# Limit minimum to still air value
    if (Q_sensible_still < *Q_sensible)
        *Q_sensible = Q_sensible_still;
    if (Q_latentheat_still < *Q_latentheat)
        *Q_latentheat = Q_latentheat_still;

    return atmos_status;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
static AED_REAL psi_m(AED_REAL zL)
{
    AED_REAL X;

    AED_REAL AA = 5.0;
    AED_REAL PIE = 3.14159;

/*----------------------------------------------------------------------------*/

    if (zL < 0.0) {
        X = pow(1.0 - 16.0*zL, 0.25);
        return 2.0*log((1.0+X)/2.0)+log((1.0+X*X)/2.0) - 2.0*atan(X) + PIE/2.0;
    } else if (zL > 0.0) {
        if (zL > 0.5) {
            if (zL > 10.0)
                return log(1.0*zL) - 0.76*zL - 12.093;
            return (0.5/(zL*zL)) -4.25/zL -7.0*log(zL)-0.852;
        }
        return -AA * zL;
    }
    return 0.0;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
static AED_REAL psi_hw(AED_REAL zL)
{
    AED_REAL X;

    AED_REAL AA = 5.0;

/*----------------------------------------------------------------------------*/

    if (zL < 0.0) {
        X = pow(1.0 - 16.0*zL, 0.25);
        return 2.0*log((1.0+X*X)/2.0);
    } else if (zL > 0.0) {
        if (zL > 0.5) {
            if (zL > 10.0)
                return log(zL) - 0.76*zL - 12.093;
            return (0.5/(zL*zL)) -4.25/zL -7.0*log(zL) -0.852;
        }
        return -AA * zL;
    }

    return 0.0;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/



/******************************************************************************
 *                                                                            *
 * solpond watfiv ncs. ez. mae,cengel, time = (, 10), p = 20                  *
 *                                                                            *
 *       This program is to determine the local rate of solar energy          *
 *       absorption by water in a water column over specified number of       *
 *       positions.                                                           *
 *                                                                            *
 *       energy(k) = fraction of energy in wavelength band k                  *
 *       absorb(k) = volume absorption coefficient of water in band k (l/m)   *
 *       hdir   = direct solar rad. incident on a harizontal surface (w/m2)   *
 *       hdif   = diffuse solar rad. incident on a horizontal surface (w/m2)  *
 *       nband  = number of wavelength bands                                  *
 *       depth  = solar pond depth (m)                                        *
 *       gx(j)  = rate of solar energy absorption in water at level d (w/m3)  *
 *       npoint = number of positions at which gx i s to be calculated        *
 *       rindex = relative refractive index of water                          *
 *       rinver = i / rindex                                                  *
 *       critw  = cosine of critical angle of incidence for water side        *
 *       crita  = cosine of critical angle of incidence for air side (=0)     *
 *       anqlei = angle of incidence of direct radiation                      *
 *       canqle = cos(anqlei)                                                 *
 *       rmu    = cosine of the angle of refraction of direct solar rad.      *
 *       rb     = reflectivity of bottom surface (diffuse)                    *
 *                                                                            *
 ******************************************************************************/




/******************************************************************************
 *                                                                            *
 * This is the main routine with the inputs passed as parameters and the      *
 *  result returned in the array "gx"                                         *
 *                                                                            *
 ******************************************************************************/
void solpond(int nband, int npoint,
    double depth, double *delz, double rb, double hdir, double anglei, double hdif,
    double *energy, double *absorb, double *gx)
{
// implicit real*8(a-h,o-z),integer*4(i,n)

    const double crita = 0.0;
    const double rindex = 1.333;
    const double pi = 3.14159265358979;

    /* local variables */
    int i, j, l;
//  double rinver, critw, cangle, rmu, del, hkdif, hkdir;
    double rinver, critw, cangle, rmu, hkdif, hkdir;

    double a, a2, vdif, vint1, vint2, em2, vdifx, vdirx, vdir1, vdir2;
    double alpha, gbk, x, xpa, xm, gxk;

/*
    // These were originally read in:
    read, nband, depth, rb, hdir, anglei, hdif, npoint
    read, (energy(i), i=1, nband), (absorb(j), j=1, nband)
*/

    for (i = 0; i < npoint; i++)
        gx[i] = zero;

    rinver = 1.0 / rindex;
    critw  = sqrt(1.0 - pow(rinver, 2));
    cangle = cos(anglei * pi/180.0);
    rmu    = sqrt(pow(rindex, 2) + pow(cangle, 2) - 1.0) / rindex;

    //printf(">solpond = %10.1f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
    //                         gx[0],gx[1],gx[2],gx[3],gx[4],gx[5],gx[6],gx[7],gx[8]);


    rb = 0.0;

    for (l = 0; l < nband; l++) {
        hkdir = energy[l] * hdir;
        hkdif = energy[l] * hdif;
        a = depth * absorb[l];

        gaus10(rindex, critw, critw, 1.0, fdif, 1, a, &vdif);
        a2 = a * 2.0;
        gaus10(rinver, critw, 0.0, critw, fct, 1, a2, &vint1);
        gaus10(rinver, critw, critw, 1.0, fct, 1, a2, &vint2);
        alpha = vint1 + vint2;

        gbk = ((1.0-ref(rindex, crita, cangle)) * hkdir * exp(-a/rmu) + 2.0 * pow(rindex, 2) * hkdif*vdif) / (1.0-2.0*rb*alpha);

        x = 0.0;
        //del = a/(npoint-1);
        //del = zero;


        gx[0] += hkdir + hkdif;

        for (j = 0; j < npoint; j++) {
            xpa = a+x;
            xm  = a-x;

            gaus10(rindex, critw, 0.0, 1.0, expint, 2, xm, &em2);
            gaus10(rindex, critw, critw, 1.0, fdif, 0, x, &vdifx);
            gaus10(rinver, critw, 0.0, critw, fct, 1, xpa, &vdir1);
            gaus10(rinver, critw, critw, 1.0, fct, 1, xpa, &vdir2);
            vdirx = vdir1 + vdir2;

            gxk = absorb[l]*((1.0-ref(rindex, crita, cangle))*exp(-x/rmu)*hkdir/rmu + 2.0 * pow(rindex, 2) * hkdif*vdifx+2.0*rb*gbk*(vdirx+em2));

            gx[j+1] += gxk;

            printf(">solpond = %10.1f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
                                   gx[0],x,gxk,gbk,a,xpa,xm,rmu);

            //del += delz[((npoint-1)-j)-1];
            //x = (j+1) * del;
            x += delz[((npoint-1)-j)-1];
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    This subroutine performs numerical integration                          *
 *    using a 10 - point gausian quadrature                                   *
 *                                                                            *
 *    c = lower bound on integration                                          *
 *    d = upper bound on integration                                          *
 *    v = value of integral                                                   *
 *    f = function to be integrated                                           *
 *                                                                            *
 ******************************************************************************/
void gaus10(double ri, double cr, double c, double d, ref_f f, int n, double h, double *v)
{
//  implicit real*8(a-h,o-z),integer*4(i,n)

    double x[10], w[10];
    double a, b;
    int i;

    *v = 0.0;

    if (c == d) return;

    x[0] = 0.148874338981631;
    x[1] = 0.433395394129247;
    x[2] = 0.679409568299024;
    x[3] = 0.865063366688985;
    x[4] = 0.973906528517172;
    w[0] = 0.295524224714753;
    w[1] = 0.269266719309996;
    w[2] = 0.219086362515982;
    w[3] = 0.149451349150581;
    w[4] = 0.066671344308688;

    for (i = 0; i < 5; i++) {
        x[i+5] = -x[i];
        w[i+5] =  w[i];
    }
    for (i = 0; i < 10; i++) {
        a = ((d-c) * x[i] + d + c) / 2.0;
        b = (d-c) * w[i] / 2.0;
        *v += (*f)(ri, cr, h, a, n) * b;
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
static double fct(double ri, double cr, double h, double cmu, int n)
{
//  implicit real*8(a-h,o-z), integer*4(i,n)
    return ref(ri, cr, cmu) * pow(cmu, n) * exp(-h/cmu);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
static double expint(double ri, double cr, double h, double cmu, int n)
//    this function is to determine exponential
{
//  implicit real*8(a-h,o-z),integer*4(i,n)
    return pow(cmu, (n-2)) * exp(-h/cmu);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/******************************************************************************
 *                                                                            *
 *    This function calculates the fresnel reflectivity                       *
 *    cmu = cosine of the angle of incidence                                  *
 *    ri = relative refractive index of the refracting medium                 *
 *    cr = critical cmu value below which reflectivity is unity               *
 *                                                                            *
 ******************************************************************************/
//    integral functions
static double ref(double ri, double cr, double cmu)
{
//  implicit real*8(a-h,o-z),integer*4(i,n)
    double ssine, rcos;

    if (cmu < cr) return 1.0;

    ssine = sqrt(pow(ri, 2) - 1.0+ pow(cmu, 2));
    rcos = pow(ri, 2) * cmu;
    return 0.50 * pow(((ssine-cmu)/(ssine+cmu)), 2) +
           0.50 * pow(((rcos-ssine)/(rcos+ssine)), 2);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/******************************************************************************
 *                                                                            *
 *    This function is to determine the refracted diffuse radiation           *
 *                                                                            *
 ******************************************************************************/
static double fdif(double rindex, double critw, double h, double cmu, int n)
{
//  implicit real*8(a-h,o-z),integer*4(i,n)
    double cmup;

    if (cmu < critw) return 0.0;

    cmup = sqrt(1.0 - pow(rindex, 2) * (1.0 - pow(cmu, 2)));
    return (1.0-ref(rindex, 0.0, cmup)) * pow(cmu, n) * exp(-h/cmu);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/



/******************************************************************************
 *                                                                            *
 *    This function has a lookup table for albedo, based on Cogley (1975)     *
 *                                                                            *
 ******************************************************************************/
static AED_REAL grischenko_albedo(AED_REAL Latitude, AED_REAL mDays)
{
    AED_REAL albedo, lat, day;

    lat = rad2deg * (Latitude - two_Pi); //# Convert latitude from rad to deg

    day = mDays;
    //# For southern hemishpere, trick the algorithm by shifting 183 days
    if ( Latitude <0 )  {
      day = mDays - 183;
      if (day <0) day = day + 365;
    }

    albedo = 7.5/100;

    if (lat >= -5 && lat < 5 ) {
      if      ( mDays >=0  && mDays <= 31 ){albedo = 6.6; }
      else if ( mDays >31  && mDays <= 59 ){albedo = 6.4; }
      else if ( mDays >59  && mDays <= 90 ){albedo = 6.3; }
      else if ( mDays >90  && mDays <= 120){albedo = 6.4; }
      else if ( mDays >120 && mDays <= 151){albedo = 6.6; }
      else if ( mDays >151 && mDays <= 181){albedo = 6.8; }
      else if ( mDays >181 && mDays <= 212){albedo = 6.7; }
      else if ( mDays >212 && mDays <= 243){albedo = 6.4; }
      else if ( mDays >243 && mDays <= 273){albedo = 6.3; }
      else if ( mDays >273 && mDays <= 304){albedo = 6.4; }
      else if ( mDays >304 && mDays <= 334){albedo = 6.6; }
      else if ( mDays >334 && mDays <= 367){albedo = 6.8; }
    }
    else if (lat >= 5 && lat < 15 ) {
      if      ( mDays >=0  && mDays <= 31 ){albedo = 7.2; }
      else if ( mDays >31  && mDays <= 59 ){albedo = 6.7; }
      else if ( mDays >59  && mDays <= 90 ){albedo = 6.4; }
      else if ( mDays >90  && mDays <= 120){albedo = 6.3; }
      else if ( mDays >120 && mDays <= 151){albedo = 6.4; }
      else if ( mDays >151 && mDays <= 181){albedo = 6.4; }
      else if ( mDays >181 && mDays <= 212){albedo = 6.4; }
      else if ( mDays >212 && mDays <= 243){albedo = 6.3; }
      else if ( mDays >243 && mDays <= 273){albedo = 6.3; }
      else if ( mDays >273 && mDays <= 304){albedo = 6.6; }
      else if ( mDays >304 && mDays <= 334){albedo = 7.1; }
      else if ( mDays >334 && mDays <= 367){albedo = 7.4; }
    }
    else if (lat >= 15 && lat < 25 ) {
      if      ( mDays >=0  && mDays <= 31 ){albedo = 8.3; }
      else if ( mDays >31  && mDays <= 59 ){albedo = 7.4; }
      else if ( mDays >59  && mDays <= 90 ){albedo = 6.7; }
      else if ( mDays >90  && mDays <= 120){albedo = 6.4; }
      else if ( mDays >120 && mDays <= 151){albedo = 6.3; }
      else if ( mDays >151 && mDays <= 181){albedo = 6.3; }
      else if ( mDays >181 && mDays <= 212){albedo = 6.3; }
      else if ( mDays >212 && mDays <= 243){albedo = 6.4; }
      else if ( mDays >243 && mDays <= 273){albedo = 6.6; }
      else if ( mDays >273 && mDays <= 304){albedo = 7.2; }
      else if ( mDays >304 && mDays <= 334){albedo = 8.1; }
      else if ( mDays >334 && mDays <= 367){albedo = 8.7; }
    }
    else if (lat >= 25 && lat < 35 ) {
      if      ( mDays >=0  && mDays <= 31 ){albedo = 10.3; }
      else if ( mDays >31  && mDays <= 59 ){albedo = 8.6; }
      else if ( mDays >59  && mDays <= 90 ){albedo = 7.3; }
      else if ( mDays >90  && mDays <= 120){albedo = 6.7; }
      else if ( mDays >120 && mDays <= 151){albedo = 6.5; }
      else if ( mDays >151 && mDays <= 181){albedo = 6.4; }
      else if ( mDays >181 && mDays <= 212){albedo = 6.4; }
      else if ( mDays >212 && mDays <= 243){albedo = 6.6; }
      else if ( mDays >243 && mDays <= 273){albedo = 7.1; }
      else if ( mDays >273 && mDays <= 304){albedo = 8.2; }
      else if ( mDays >304 && mDays <= 334){albedo = 10.0; }
      else if ( mDays >334 && mDays <= 367){albedo = 11.1; }
    }
    else if (lat >= 35 && lat < 45 ) {
      if      ( mDays >=0  && mDays <= 31 ){albedo = 14.5; }
      else if ( mDays >31  && mDays <= 59 ){albedo = 11.1; }
      else if ( mDays >59  && mDays <= 90 ){albedo = 8.5; }
      else if ( mDays >90  && mDays <= 120){albedo = 7.3; }
      else if ( mDays >120 && mDays <= 151){albedo = 6.8; }
      else if ( mDays >151 && mDays <= 181){albedo = 6.7; }
      else if ( mDays >181 && mDays <= 212){albedo = 6.8; }
      else if ( mDays >212 && mDays <= 243){albedo = 7.1; }
      else if ( mDays >243 && mDays <= 273){albedo = 8.0; }
      else if ( mDays >273 && mDays <= 304){albedo = 10.3; }
      else if ( mDays >304 && mDays <= 334){albedo = 13.8; }
      else if ( mDays >334 && mDays <= 367){albedo = 16.1; }
    }
    else if (lat >= 45 && lat < 55 ) {
      if      ( mDays >=0  && mDays <= 31 ){albedo = 22.0; }
      else if ( mDays >31  && mDays <= 59 ){albedo = 16.1; }
      else if ( mDays >59  && mDays <= 90 ){albedo = 10.8; }
      else if ( mDays >90  && mDays <= 120){albedo = 8.4; }
      else if ( mDays >120 && mDays <= 151){albedo = 7.5; }
      else if ( mDays >151 && mDays <= 181){albedo = 7.3; }
      else if ( mDays >181 && mDays <= 212){albedo = 7.4; }
      else if ( mDays >212 && mDays <= 243){albedo = 8.0; }
      else if ( mDays >243 && mDays <= 273){albedo = 9.0; }
      else if ( mDays >273 && mDays <= 304){albedo = 14.4; }
      else if ( mDays >304 && mDays <= 334){albedo = 21.0; }
      else if ( mDays >334 && mDays <= 367){albedo = 24.1; }
    }
    else if (lat >= 55 && lat < 65 ) {
      if      ( mDays >=0  && mDays <= 31 ){albedo = 33.9; }
      else if ( mDays >31  && mDays <= 59 ){albedo = 24.0; }
      else if ( mDays >59  && mDays <= 90 ){albedo = 15.5; }
      else if ( mDays >90  && mDays <= 120){albedo = 10.5; }
      else if ( mDays >120 && mDays <= 151){albedo = 8.8; }
      else if ( mDays >151 && mDays <= 181){albedo = 8.4; }
      else if ( mDays >181 && mDays <= 212){albedo = 8.6; }
      else if ( mDays >212 && mDays <= 243){albedo = 9.8; }
      else if ( mDays >243 && mDays <= 273){albedo = 13.6; }
      else if ( mDays >273 && mDays <= 304){albedo = 21.6; }
      else if ( mDays >304 && mDays <= 334){albedo = 32.1; }
      else if ( mDays >334 && mDays <= 367){albedo = 35.5; }
    }
    else if (lat >= 65 && lat < 75 ) {
      if      ( mDays >=0  && mDays <= 31 ){albedo = 30.1; }
      else if ( mDays >31  && mDays <= 59 ){albedo = 33.8; }
      else if ( mDays >59  && mDays <= 90 ){albedo = 22.9; }
      else if ( mDays >90  && mDays <= 120){albedo = 14.8; }
      else if ( mDays >120 && mDays <= 151){albedo = 11.6; }
      else if ( mDays >151 && mDays <= 181){albedo = 11.2; }
      else if ( mDays >181 && mDays <= 212){albedo = 11.4; }
      else if ( mDays >212 && mDays <= 243){albedo = 13.4; }
      else if ( mDays >243 && mDays <= 273){albedo = 20.2; }
      else if ( mDays >273 && mDays <= 304){albedo = 31.3; }
      else if ( mDays >304 && mDays <= 334){albedo = 30.1; }
      else if ( mDays >334 && mDays <= 367){albedo = 40.0; }
    }
    else if (lat >= 75 && lat < 85 ) {
      if      ( mDays >=0  && mDays <= 31 ){albedo = 40.0; }
      else if ( mDays >31  && mDays <= 59 ){albedo = 30.1; }
      else if ( mDays >59  && mDays <= 90 ){albedo = 31.9; }
      else if ( mDays >90  && mDays <= 120){albedo = 22.5; }
      else if ( mDays >120 && mDays <= 151){albedo = 16.0; }
      else if ( mDays >151 && mDays <= 181){albedo = 13.1; }
      else if ( mDays >181 && mDays <= 212){albedo = 14.5; }
      else if ( mDays >212 && mDays <= 243){albedo = 20.6; }
      else if ( mDays >243 && mDays <= 273){albedo = 29.4; }
      else if ( mDays >273 && mDays <= 304){albedo = 30.5; }
      else if ( mDays >304 && mDays <= 334){albedo = 40.0; }
      else if ( mDays >334 && mDays <= 367){albedo = 40.0; }
    }
    else if (lat > 85 ) {
      if      ( mDays >=0  && mDays <= 31 ){albedo = 40.0; }
      else if ( mDays >31  && mDays <= 59 ){albedo = 40.0; }
      else if ( mDays >59  && mDays <= 90 ){albedo = 30.1; }
      else if ( mDays >90  && mDays <= 120){albedo = 29.3; }
      else if ( mDays >120 && mDays <= 151){albedo = 17.1; }
      else if ( mDays >151 && mDays <= 181){albedo = 14.8; }
      else if ( mDays >181 && mDays <= 212){albedo = 16.0; }
      else if ( mDays >212 && mDays <= 243){albedo = 24.6; }
      else if ( mDays >243 && mDays <= 273){albedo = 34.2; }
      else if ( mDays >273 && mDays <= 304){albedo = 40.0; }
      else if ( mDays >304 && mDays <= 334){albedo = 40.0; }
      else if ( mDays >334 && mDays <= 367){albedo = 40.0; }
    }
  return albedo/100;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/








/******************************************************************************
 *                                                                            *
 ******************************************************************************/
#if 0

static void soiltemp(int m, double wv)
{
    AED_REAL *w, *t, *tn, *k, *cp, *a, *b, *c, *d, *z;

    int ta = 20, am = 15, bd = 1.3, tb = 20;

    AED_REAL ti, dt, da, f, g, mc, c1, c2, c3, c4;

    int i;

    w =  calloc((m+1), sizeof(AED_REAL) );
    t =  calloc((m+1), sizeof(AED_REAL) );
    tn = calloc((m+1), sizeof(AED_REAL) );
    k =  calloc((m+1), sizeof(AED_REAL) );
    cp = calloc(m,     sizeof(AED_REAL) );
    a =  calloc((m+1), sizeof(AED_REAL) );
    b =  calloc(m,     sizeof(AED_REAL) );
    c =  calloc(m,     sizeof(AED_REAL) );
    d =  calloc(m,     sizeof(AED_REAL) );
    z =  calloc((m+1), sizeof(AED_REAL) );

    k[0] = 20;    // boundary layer conductance in w/(m^2 k)

    for (i = 1; i <= m; i++) {
        z[i+1] = z[i]+.005* pow(1.5, (i-1));
        t[i] = tb;
    }

    t[m+1] = tb ; tn[m+1] = t[m+1] ; t[0] = tb;

    // ti is time of day;dt is time step (sec);da is day number
    ti = 0 ; dt = 3600 ; da = 0;
    f = .6 ; g = 1-f;

    mc = .12;   // clay fraction
    c1 = .65-.78*bd+.6*bd*bd ; c2 = 1.06*bd ; c3 = 1+2.6/sqrt(mc) ; c4 = .3+.1*bd*bd;

    for (i = 1; i <= m; i++) {
        cp[i] = (2400000*bd/2.65+4180000*wv)*(z[i+1]-z[i-1])/(2*dt);
        k[i] = (c1+c2*wv-(c1-c4)*exp(-pow((c3*wv), 4)))/(z[i+1]-z[i]);
    }

    do {
        ti = ti+dt/3600 ; if ( ti > 24 ) { ti = ti-24 ; da = da+1; }
        tn[0] = ta+am*sin(.261799*(ti-6));
        for (i = 1; i <= m; i++) {
            c[i] = -k[i]*f  ;  a[i+1] = c[i];
            b[i] = f*(k[i]+k[i-1])+cp[i];
            d[i] = g*k[i-1]*t[i-1]+(cp[i]-g*(k[i]+k[i-1]))*t[i]+g*k[i]*t[i+1];
        }
        d[1] = d[1]+k[0]*tn[0]*f;
        d[m] = d[m]+k[m]*f*tn[m+1];
        for (i = 1; i <= m-1; i++) {
            c[i] = c[i]/b[i];
            d[i] = d[i]/b[i];
            b[i+1] = b[i+1]-a[i+1]*c[i];
            d[i+1] = d[i+1]-a[i+1]*d[i];
        }
        tn[m] = d[m]/b[m];
        for (i = m-1; i>= 1; i--) {
            tn[i] = d[i]-c[i]*tn[i+1];
        }
//      print "day =";da,"hour =";ti
//      print "heat flux =";k[0]*(g*(t[0]-t[1])+f*(tn[0]-tn[1]));"w/m2"
//      print "depth","temperature","k[i]"
        for (i = 0; i <= m+1; i++) {
//          print z[i],tn[i],k[i];
            t[i] = tn[i];
        }
    } while (da < 5);

    free(w); free(t); free(tn); free(k); free(cp);
    free(a); free(b); free(c);  free(d); free(z);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif
