/******************************************************************************
 *                                                                            *
 * glm_ncdf.c                                                                 *
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

#include <netcdf.h>

#include "glm.h"
#include "glm_const.h"
#include "glm_types.h"
#include "glm_globals.h"
#include "glm_ncdf.h"

int ncid=-1;

static int set_no = -1;

//# dimension sizes
int x_dim, y_dim, z_dim, zone_dim, time_dim, ptm_dim;

//# dimension lengths
static int lon_len=1;
static int lat_len=1;
static int height_len;
static size_t start[4],edges[4];

//# variable ids
static int lon_id, lat_id, z_id, H_id, V_id, A_id, TV_id, Taub_id, NS_id, time_id;
static int SL_id, HICE_id, HSNOW_id, HWICE_id;
static int precip_id, evap_id, rho_id, rad_id, extc_id, i0_id, wnd_id;
static int temp_id, salt_id, umean_id, uorb_id, epsilon_id;

//# from lake.csv
static int SA_id, VSnow_id,VBIce_id,VWIce_id, TotInVol_id, TotOutVol_id;
static int OverFVol_id, Evap_id, Rain_id, LocRunoff_id, SnowF_id, LakeLvl_id, Seep_id;
static int SnowDns_id, Alb_id, MaxT_id, MinT_id, SfT_id, DQsw_id, DQe_id;
static int DQh_id, DQlw_id, Light_id, BenLight_id, SWH_id, SWL_id, SWP_id;
static int Qsw_id, Qe_id, Qh_id, Qlw_id;
static int LkNum_id, maxdtz_id, CD_id, CHE_id, zL_id;
static int sedlayer_dim, sed_temp_id, sed_hf_id, sed_zt_id, sed_ld_id;
static int zarea_id;

#ifdef _WIN32
    char *strndup(const char *s, size_t len);
#endif

/*============================================================================*/
void check_nc_error(int err);
static void check_nc_error_x(int err, int ncid, int id);
static void xyt_store_nc_scalar(int ncid, int id, AED_REAL scalar);

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
int init_glm_ncdf(const char *fn, const char *title, AED_REAL lat,
                                AED_REAL lon, int nlev, const char *start_time)
{
    int ncid;
    char time_str[128], history[128];
    extern char wq_lib[];

    int dims[4];
/*----------------------------------------------------------------------------*/

    if ( fn == NULL || strcmp(fn, "")==0 ) return -1;

    check_nc_error(nc_create(fn, NC_CLOBBER, &ncid));

    snprintf(time_str, 128, "hours since %s", start_time);
    snprintf(history,  128, "Created by glm3/%s v. %s", wq_lib, GLM_VERSION);

    height_len = nlev;

    //# define dimensions
    check_nc_error(nc_def_dim(ncid, "lon", 1, &x_dim));
    check_nc_error(nc_def_dim(ncid, "lat", 1, &y_dim));
    check_nc_error(nc_def_dim(ncid, "z", nlev, &z_dim));
    check_nc_error(nc_def_dim(ncid, "particle", 1000000, &ptm_dim));

    if ( n_zones > 0 )
        check_nc_error(nc_def_dim(ncid, "nzones", n_zones+1, &zone_dim));
    check_nc_error(nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim));

    //# define coordinates
    dims[0] = x_dim;
    check_nc_error(nc_def_var(ncid, "lon", NC_REALTYPE, 1, dims, &lon_id));

    dims[0] = y_dim;
    check_nc_error(nc_def_var(ncid, "lat", NC_REALTYPE, 1, dims, &lat_id));

    dims[0] = time_dim;
    check_nc_error(nc_def_var(ncid, "NS",    NC_INT,      1, dims, &NS_id));
    check_nc_error(nc_def_var(ncid, "time",  NC_REALTYPE, 1, dims, &time_id));
    check_nc_error(nc_def_var(ncid, "blue_ice_thickness",  NC_REALTYPE, 1, dims, &HICE_id));
    check_nc_error(nc_def_var(ncid, "snow_thickness",      NC_REALTYPE, 1, dims, &HSNOW_id));
    check_nc_error(nc_def_var(ncid, "white_ice_thickness", NC_REALTYPE, 1, dims, &HWICE_id));
    check_nc_error(nc_def_var(ncid, "surface_layer",       NC_INT,      1, dims, &SL_id));

    /**************************************************************************
     * define 2D variables                                                    *
     **************************************************************************/

    //# x,y,t
    dims[2] = x_dim;
    dims[1] = y_dim;
    dims[0] = time_dim;

    check_nc_error(nc_def_var(ncid, "precipitation",  NC_REALTYPE, 3, dims, &precip_id));
    check_nc_error(nc_def_var(ncid, "evap_mass_flux", NC_REALTYPE, 3, dims, &evap_id));
    check_nc_error(nc_def_var(ncid, "solar",          NC_REALTYPE, 3, dims, &i0_id));
    check_nc_error(nc_def_var(ncid, "wind",           NC_REALTYPE, 3, dims, &wnd_id));
    check_nc_error(nc_def_var(ncid, "lake_volume",    NC_REALTYPE, 3, dims, &TV_id));

    /**************************************************************************
     * Define 2D variables from lake.csv                                      *
     **************************************************************************/
    check_nc_error(nc_def_var(ncid, "surface_area",        NC_REALTYPE, 3, dims, &SA_id));
    check_nc_error(nc_def_var(ncid, "vol_snow",            NC_REALTYPE, 3, dims, &VSnow_id));
    check_nc_error(nc_def_var(ncid, "vol_blue_ice",        NC_REALTYPE, 3, dims, &VBIce_id));
    check_nc_error(nc_def_var(ncid, "vol_white_ice",       NC_REALTYPE, 3, dims, &VWIce_id));
    check_nc_error(nc_def_var(ncid, "tot_inflow_vol",      NC_REALTYPE, 3, dims, &TotInVol_id));
    check_nc_error(nc_def_var(ncid, "tot_outflow_vol",     NC_REALTYPE, 3, dims, &TotOutVol_id));
    check_nc_error(nc_def_var(ncid, "overflow_vol",        NC_REALTYPE, 3, dims, &OverFVol_id));
    check_nc_error(nc_def_var(ncid, "evaporation",         NC_REALTYPE, 3, dims, &Evap_id));
    check_nc_error(nc_def_var(ncid, "rain",                NC_REALTYPE, 3, dims, &Rain_id));
    check_nc_error(nc_def_var(ncid, "local_runoff",        NC_REALTYPE, 3, dims, &LocRunoff_id));
    check_nc_error(nc_def_var(ncid, "snowfall",            NC_REALTYPE, 3, dims, &SnowF_id));
    check_nc_error(nc_def_var(ncid, "seepage_vol",         NC_REALTYPE, 3, dims, &Seep_id));
    check_nc_error(nc_def_var(ncid, "lake_level",          NC_REALTYPE, 3, dims, &LakeLvl_id));
    check_nc_error(nc_def_var(ncid, "snow_density",        NC_REALTYPE, 3, dims, &SnowDns_id));
    check_nc_error(nc_def_var(ncid, "albedo",              NC_REALTYPE, 3, dims, &Alb_id));
    check_nc_error(nc_def_var(ncid, "max_temp",            NC_REALTYPE, 3, dims, &MaxT_id));
    check_nc_error(nc_def_var(ncid, "min_temp",            NC_REALTYPE, 3, dims, &MinT_id));
    check_nc_error(nc_def_var(ncid, "surface_temp",        NC_REALTYPE, 3, dims, &SfT_id));
    check_nc_error(nc_def_var(ncid, "daily_qsw",           NC_REALTYPE, 3, dims, &DQsw_id));
    check_nc_error(nc_def_var(ncid, "daily_qe",            NC_REALTYPE, 3, dims, &DQe_id));
    check_nc_error(nc_def_var(ncid, "daily_qh",            NC_REALTYPE, 3, dims, &DQh_id));
    check_nc_error(nc_def_var(ncid, "daily_qlw",           NC_REALTYPE, 3, dims, &DQlw_id));
    check_nc_error(nc_def_var(ncid, "Qsw",                 NC_REALTYPE, 3, dims, &Qsw_id));
    check_nc_error(nc_def_var(ncid, "Qe",                  NC_REALTYPE, 3, dims, &Qe_id));
    check_nc_error(nc_def_var(ncid, "Qh",                  NC_REALTYPE, 3, dims, &Qh_id));
    check_nc_error(nc_def_var(ncid, "Qlw",                 NC_REALTYPE, 3, dims, &Qlw_id));
    check_nc_error(nc_def_var(ncid, "light",               NC_REALTYPE, 3, dims, &Light_id));
    check_nc_error(nc_def_var(ncid, "benthic_light",       NC_REALTYPE, 3, dims, &BenLight_id));
    check_nc_error(nc_def_var(ncid, "surface_wave_height", NC_REALTYPE, 3, dims, &SWH_id));
    check_nc_error(nc_def_var(ncid, "surface_wave_length", NC_REALTYPE, 3, dims, &SWL_id));
    check_nc_error(nc_def_var(ncid, "surface_wave_period", NC_REALTYPE, 3, dims, &SWP_id));
    check_nc_error(nc_def_var(ncid, "lake_number",         NC_REALTYPE, 3, dims, &LkNum_id));
    check_nc_error(nc_def_var(ncid, "max_dT_dz",           NC_REALTYPE, 3, dims, &maxdtz_id));
    check_nc_error(nc_def_var(ncid, "CD",                  NC_REALTYPE, 3, dims, &CD_id));
    check_nc_error(nc_def_var(ncid, "CHE",                 NC_REALTYPE, 3, dims, &CHE_id));
    check_nc_error(nc_def_var(ncid, "z_L",                 NC_REALTYPE, 3, dims, &zL_id));

    /**************************************************************************
     * define 3D variables                                                    *
     **************************************************************************/

    //# x,y,z,t
    dims[3] = x_dim;
    dims[2] = y_dim;
    dims[1] = z_dim;
    dims[0] = time_dim;

    check_nc_error(nc_def_var(ncid, "z",       NC_REALTYPE, 4, dims, &z_id));
    check_nc_error(nc_def_var(ncid, "H",       NC_REALTYPE, 4, dims, &H_id));
    check_nc_error(nc_def_var(ncid, "V",       NC_REALTYPE, 4, dims, &V_id));
    check_nc_error(nc_def_var(ncid, "A",       NC_REALTYPE, 4, dims, &A_id));
    check_nc_error(nc_def_var(ncid, "salt",    NC_REALTYPE, 4, dims, &salt_id));
    check_nc_error(nc_def_var(ncid, "temp",    NC_REALTYPE, 4, dims, &temp_id));

    check_nc_error(nc_def_var(ncid, "dens",    NC_REALTYPE, 4, dims, &rho_id));
    check_nc_error(nc_def_var(ncid, "radn",    NC_REALTYPE, 4, dims, &rad_id));
    check_nc_error(nc_def_var(ncid, "extc",    NC_REALTYPE, 4, dims, &extc_id));

    check_nc_error(nc_def_var(ncid, "umean",   NC_REALTYPE, 4, dims, &umean_id));
    check_nc_error(nc_def_var(ncid, "uorb",    NC_REALTYPE, 4, dims, &uorb_id));
    check_nc_error(nc_def_var(ncid, "taub",    NC_REALTYPE, 4, dims, &Taub_id));
    check_nc_error(nc_def_var(ncid, "epsilon", NC_REALTYPE, 4, dims, &epsilon_id));

    dims[3] = x_dim;
    dims[2] = y_dim;
    dims[1] = zone_dim;
    dims[0] = time_dim;
    check_nc_error(nc_def_var(ncid, "zarea", NC_REALTYPE, 4, dims, &zarea_id));

    /**************************************************************************
     * assign attributes                                                      *
     **************************************************************************/

    //# coordinates
    set_nc_attributes(ncid, lon_id,    "degrees_east",  NULL        PARAM_FILLVALUE);
    set_nc_attributes(ncid, lat_id,    "degrees_north", NULL        PARAM_FILLVALUE);

    set_nc_attributes(ncid, time_id,   time_str,        NULL        PARAM_FILLVALUE);

    nc_put_att(ncid, NS_id, "long_name", NC_CHAR, 16, "Number of Layers");
    nc_put_att(ncid, SL_id, "long_name", NC_CHAR, 13, "Surface Layer");

    set_nc_attributes(ncid, HICE_id,   "meters",  "Height of Ice"   PARAM_FILLVALUE);
    set_nc_attributes(ncid, HSNOW_id,  "meters",  "Height of Snow"  PARAM_FILLVALUE);
    set_nc_attributes(ncid, HWICE_id,  "meters",  "Height of WhiteIce" PARAM_FILLVALUE);

    //# x,y,t
    set_nc_attributes(ncid, precip_id, "m/s",     "precipitation"   PARAM_FILLVALUE);
    set_nc_attributes(ncid, evap_id,   "m/s",     "evaporation"     PARAM_FILLVALUE);
    set_nc_attributes(ncid, i0_id,     "10E-6m",  "Shortwave"       PARAM_FILLVALUE);
    set_nc_attributes(ncid, wnd_id,    "m/s",     "wind"            PARAM_FILLVALUE);
    set_nc_attributes(ncid, TV_id,     "m3",      "lake volume"     PARAM_FILLVALUE);

    // 2D variables from lake.csv
    set_nc_attributes(ncid, SA_id,        "m2",      "Surface Area"         PARAM_FILLVALUE);
    set_nc_attributes(ncid, VSnow_id,     "m3",      "Vol Snow"             PARAM_FILLVALUE);
    set_nc_attributes(ncid, VBIce_id,     "m3",      "Vol Blue Ice"         PARAM_FILLVALUE);
    set_nc_attributes(ncid, VWIce_id,     "m3",      "Vol White Ice"        PARAM_FILLVALUE);
    set_nc_attributes(ncid, TotInVol_id,  "m3",      "Total Inflow Volume"  PARAM_FILLVALUE);
    set_nc_attributes(ncid, TotOutVol_id, "m3",      "Total Outflow Volume" PARAM_FILLVALUE);
    set_nc_attributes(ncid, OverFVol_id,  "m3",      "Overflow Volume"      PARAM_FILLVALUE);
    set_nc_attributes(ncid, Evap_id,      "m3",      "Evaporation"          PARAM_FILLVALUE);
    set_nc_attributes(ncid, Rain_id,      "m3",      "Rain"                 PARAM_FILLVALUE);
    set_nc_attributes(ncid, LocRunoff_id, "m3",      "Local Runoff"         PARAM_FILLVALUE);
    set_nc_attributes(ncid, SnowF_id,     "m3",      "Snowfall"             PARAM_FILLVALUE);
    set_nc_attributes(ncid, Seep_id,      "m3",      "Seepage Volume"       PARAM_FILLVALUE);
    set_nc_attributes(ncid, LakeLvl_id,   "meters",  "Lake Level"           PARAM_FILLVALUE);
    set_nc_attributes(ncid, SnowDns_id,   "unknown", "Snow Density"         PARAM_FILLVALUE);
    set_nc_attributes(ncid, Alb_id,       "unknown", "Surface Albedo"       PARAM_FILLVALUE);
    set_nc_attributes(ncid, MaxT_id,      "celsius", "Maximum Temperature"  PARAM_FILLVALUE);
    set_nc_attributes(ncid, MinT_id,      "celsius", "Minimum Temperature"  PARAM_FILLVALUE);
    set_nc_attributes(ncid, SfT_id,       "celsius", "Surface Temperature"  PARAM_FILLVALUE);
    set_nc_attributes(ncid, DQsw_id,      "unknown", "Daily Qsw"            PARAM_FILLVALUE);
    set_nc_attributes(ncid, DQe_id,       "unknown", "Daily Qe"             PARAM_FILLVALUE);
    set_nc_attributes(ncid, DQh_id,       "unknown", "Daily Qh"             PARAM_FILLVALUE);
    set_nc_attributes(ncid, DQlw_id,      "unknown", "Daily Qlw"            PARAM_FILLVALUE);
    set_nc_attributes(ncid, Qsw_id,       "W/m2",   "Shortwave Radiation"  PARAM_FILLVALUE);
    set_nc_attributes(ncid, Qe_id,        "W/m2",   "Latent Heat Flux"     PARAM_FILLVALUE);
    set_nc_attributes(ncid, Qh_id,        "W/m2",   "Sensible Heat Flux"   PARAM_FILLVALUE);
    set_nc_attributes(ncid, Qlw_id,       "W/m2",   "Net Longwave Flux"    PARAM_FILLVALUE);
    set_nc_attributes(ncid, Light_id,     "unknown", "Light"                PARAM_FILLVALUE);
    set_nc_attributes(ncid, BenLight_id,  "unknown", "Benthic Light"        PARAM_FILLVALUE);
    set_nc_attributes(ncid, SWH_id,       "meters",  "Surface Wave Height"  PARAM_FILLVALUE);
    set_nc_attributes(ncid, SWL_id,       "meters",  "Surface Wave Length"  PARAM_FILLVALUE);
    set_nc_attributes(ncid, SWP_id,       "meters",  "Surface Wave Period"  PARAM_FILLVALUE);
    set_nc_attributes(ncid, LkNum_id,     "unknown", "Lake Number"          PARAM_FILLVALUE);
    set_nc_attributes(ncid, maxdtz_id,    "m/s",     "Max dT_dz"            PARAM_FILLVALUE);
    set_nc_attributes(ncid, CD_id,        "m/s",     "Coeff. of Wind Drag"  PARAM_FILLVALUE);
    set_nc_attributes(ncid, CHE_id,       "unknown", "CHE"                  PARAM_FILLVALUE);
    set_nc_attributes(ncid, zL_id,        "unknown", "z_L"                  PARAM_FILLVALUE);

    //# x,y,z,t
    set_nc_attributes(ncid, z_id,       "meters",  "layer heights"    PARAM_FILLVALUE);
    set_nc_attributes(ncid, H_id,       "meters",  "layer heights"    PARAM_FILLVALUE);
    set_nc_attributes(ncid, V_id,       "m3",      "layer volume"     PARAM_FILLVALUE);
    set_nc_attributes(ncid, A_id,       "m2",      "layer area"       PARAM_FILLVALUE);
    set_nc_attributes(ncid, salt_id,    "g/kg",    "salinity"         PARAM_FILLVALUE);
    set_nc_attributes(ncid, temp_id,    "celsius", "temperature"      PARAM_FILLVALUE);

    set_nc_attributes(ncid, extc_id,    "unknown", "extc_coef"        PARAM_FILLVALUE);
    set_nc_attributes(ncid, rho_id,     "unknown", "density"          PARAM_FILLVALUE);
    set_nc_attributes(ncid, rad_id,     "unknown", "solar radiation"  PARAM_FILLVALUE);

    set_nc_attributes(ncid, umean_id,   "m/s",     "mean velocity"    PARAM_FILLVALUE);
    set_nc_attributes(ncid, uorb_id,    "m/s",     "orbital velocity" PARAM_FILLVALUE);

    set_nc_attributes(ncid, Taub_id,    "N/m2",    "layer stress"     PARAM_FILLVALUE);
    set_nc_attributes(ncid, epsilon_id, "m2/s",    "diffusivity"      PARAM_FILLVALUE);
    set_nc_attributes(ncid, zarea_id,   "m2",      "sediment zone area" PARAM_FILLVALUE);

    //# global attributes
    nc_put_att(ncid, NC_GLOBAL, "Title", NC_CHAR, strlen(title), title);

    nc_put_att(ncid, NC_GLOBAL, "history", NC_CHAR, strlen(history), history);
    check_nc_error(nc_put_att(ncid, NC_GLOBAL, "Conventions", NC_CHAR, 6, "COARDS"));

    nc_put_att(ncid, NC_GLOBAL, "start_time", NC_CHAR, strlen(start_time), start_time);

    /**************************************************************************
     * Save some static data                                                  *
     **************************************************************************/

    /**************************************************************************
     * sediment temperature model (sed_heat_model == 2)                       *
     **************************************************************************/
    if ( sed_heat_model == 2 && n_zones > 0 && theZones != NULL ) {
        int nsl = theZones[0].n_sed_layers;
        check_nc_error(nc_def_dim(ncid, "sed_layers", nsl, &sedlayer_dim));

        dims[0] = time_dim; dims[1] = zone_dim; dims[2] = sedlayer_dim;
        check_nc_error(nc_def_var(ncid, "sed_temp", NC_REALTYPE, 3, dims, &sed_temp_id));

        dims[0] = time_dim; dims[1] = zone_dim;
        check_nc_error(nc_def_var(ncid, "sed_heatflux", NC_REALTYPE, 2, dims, &sed_hf_id));
        check_nc_error(nc_def_var(ncid, "zone_ztemp",   NC_REALTYPE, 2, dims, &sed_zt_id));

        dims[0] = sedlayer_dim;
        check_nc_error(nc_def_var(ncid, "sed_layer_depth", NC_REALTYPE, 1, dims, &sed_ld_id));

        set_nc_attributes(ncid, sed_temp_id, "celsius", "sediment temperature"         PARAM_FILLVALUE);
        set_nc_attributes(ncid, sed_hf_id,   "W/m2",    "sediment-water heat flux"      PARAM_FILLVALUE);
        set_nc_attributes(ncid, sed_zt_id,   "celsius", "zone water temp (sed top BC)"  PARAM_FILLVALUE);
        set_nc_attributes(ncid, sed_ld_id,   "meters",  "sediment node depth"           PARAM_FILLVALUE);
    }

    //# leave define mode
    check_nc_error(nc_enddef(ncid));

    //# save latitude and longitude
    store_nc_scalar(ncid, lon_id, POINT, lon);
    store_nc_scalar(ncid, lat_id, POINT, lat);

    //# save static sediment node depths
    if ( sed_heat_model == 2 && n_zones > 0 && theZones != NULL ) {
        int nsl = theZones[0].n_sed_layers, k;
        AED_REAL *ld = malloc(nsl * sizeof(AED_REAL));
        for (k = 0; k < nsl; k++) ld[k] = theZones[0].layers[k].depth;
        check_nc_error(nc_put_var_double(ncid, sed_ld_id, ld));
        free(ld);
    }

    check_nc_error(nc_sync(ncid));

    return ncid;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void write_glm_ncdf(int ncid, int wlev, int nlev, int stepnum, AED_REAL timestep)
{
    AED_REAL temp_time, LakeVolume;
    AED_REAL *heights, *vols, *area, *salts, *temps, *dens, *qsw, *extc_coef;
    AED_REAL *u_mean, *u_orb, *taub, *epsilon;
    int i, littoralLayer = 0;

    if (ncid == -1) return;

    set_no++;

    store_nc_integer(ncid, NS_id, T_SHAPE, wlev);
    store_nc_integer(ncid, SL_id, T_SHAPE, wlev);

    temp_time = (stepnum * timestep) / 3600.;
    LakeVolume = 0.0;
    for (i = 0; i < wlev; i++)
        LakeVolume += Lake[i].LayerVol;

    //# Time varying data : t
    /*------------------------------------------------------------------------*/
    store_nc_scalar(ncid,  time_id, T_SHAPE, temp_time);
    store_nc_scalar(ncid,  HICE_id, T_SHAPE, SurfData.delzBlueIce);
    store_nc_scalar(ncid, HWICE_id, T_SHAPE, SurfData.delzWhiteIce);
    store_nc_scalar(ncid, HSNOW_id, T_SHAPE, SurfData.delzSnow);

    store_nc_scalar(ncid, precip_id, XYT_SHAPE, MetData.Rain);
    store_nc_scalar(ncid,   evap_id, XYT_SHAPE, SurfData.Evap);

    store_nc_scalar(ncid,  i0_id, XYT_SHAPE, MetData.ShortWave);
    store_nc_scalar(ncid, wnd_id, XYT_SHAPE, MetData.WindSpeed);

    store_nc_scalar(ncid, Qsw_id, XYT_SHAPE, SurfData.Qsw);
    store_nc_scalar(ncid,  Qe_id, XYT_SHAPE, SurfData.Qe);
    store_nc_scalar(ncid,  Qh_id, XYT_SHAPE, SurfData.Qh);
    store_nc_scalar(ncid, Qlw_id, XYT_SHAPE, SurfData.Qlw);

    store_nc_scalar(ncid,  TV_id, XYT_SHAPE, LakeVolume);

    //# Time varying profile data : z,t
    /*------------------------------------------------------------------------*/
    start[3] = 0;      edges[3] = lon_len;
    start[2] = 0;      edges[2] = lat_len;
    start[1] = 0;      edges[1] = height_len;
    start[0] = set_no; edges[0] = 1;

    heights   = malloc(nlev*sizeof(AED_REAL));
    vols      = malloc(nlev*sizeof(AED_REAL));
    area      = malloc(nlev*sizeof(AED_REAL));
    salts     = malloc(nlev*sizeof(AED_REAL));
    temps     = malloc(nlev*sizeof(AED_REAL));
    dens      = malloc(nlev*sizeof(AED_REAL));
    qsw       = malloc(nlev*sizeof(AED_REAL));
    extc_coef = malloc(nlev*sizeof(AED_REAL));
    u_mean    = malloc(nlev*sizeof(AED_REAL));
    u_orb     = malloc(nlev*sizeof(AED_REAL));
    taub      = malloc(nlev*sizeof(AED_REAL));
    epsilon   = malloc(nlev*sizeof(AED_REAL));

    for (i = 0; i < wlev; i++) {
        heights[i] = Lake[i].Height;
        area[i] = Lake[i].LayerArea;
        vols[i] = Lake[i].LayerVol;
        salts[i] = Lake[i].Salinity;
        temps[i] = Lake[i].Temp;
        dens[i] = Lake[i].Density;
        qsw[i] = Lake[i].Light;
        extc_coef[i] = Lake[i].ExtcCoefSW;
        u_mean[i] = Lake[i].Umean;
        u_orb[i] = Lake[i].Uorb;
        taub[i] = Lake[i].LayerStress;
        epsilon[i] = Lake[i].Epsilon;
    }
    for (i = wlev; i < nlev; i++) {
        heights[i] = NC_FILLER;
        vols[i] = NC_FILLER;
        area[i] = NC_FILLER;
        salts[i] = NC_FILLER;
        temps[i] = NC_FILLER;
        dens[i] = NC_FILLER;
        qsw[i] = NC_FILLER;
        extc_coef[i] = NC_FILLER;
        u_mean[i] = NC_FILLER;
        u_orb[i] = NC_FILLER;
        taub[i] = NC_FILLER;
        epsilon[i] = NC_FILLER;
    }
    if ( littoral_sw ) {
        littoralLayer = wlev;
        heights[littoralLayer] = Lake[littoralLayer].Height;
        vols[littoralLayer] = Lake[littoralLayer].LayerVol;
        salts[littoralLayer] = Lake[littoralLayer].Salinity;
        temps[littoralLayer] = Lake[littoralLayer].Temp;
        dens[littoralLayer] = NC_FILLER;
        qsw[littoralLayer] = NC_FILLER;
        extc_coef[littoralLayer] = NC_FILLER;
        u_mean[littoralLayer] = NC_FILLER;
        u_orb[littoralLayer] = NC_FILLER;
        taub[littoralLayer] = NC_FILLER;
        epsilon[littoralLayer] = NC_FILLER;
    }

    check_nc_error(nc_put_vara(ncid,       z_id, start, edges, heights));
    check_nc_error(nc_put_vara(ncid,       H_id, start, edges, heights));
    check_nc_error(nc_put_vara(ncid,       V_id, start, edges, vols));
    check_nc_error(nc_put_vara(ncid,       A_id, start, edges, area));
    check_nc_error(nc_put_vara(ncid,    salt_id, start, edges, salts));
    check_nc_error(nc_put_vara(ncid,    temp_id, start, edges, temps));
    check_nc_error(nc_put_vara(ncid,     rho_id, start, edges, dens));
    check_nc_error(nc_put_vara(ncid,     rad_id, start, edges, qsw));
    check_nc_error(nc_put_vara(ncid,    extc_id, start, edges, extc_coef));
    check_nc_error(nc_put_vara(ncid,   umean_id, start, edges, u_mean));
    check_nc_error(nc_put_vara(ncid,    uorb_id, start, edges, u_orb));
    check_nc_error(nc_put_vara(ncid,    Taub_id, start, edges, taub));
    check_nc_error(nc_put_vara(ncid, epsilon_id, start, edges, epsilon));

    free(heights); free(vols); free(area); free(salts);
    free(temps);  free(dens); free(qsw);
    free(extc_coef); free(u_mean); free(u_orb);
    free(taub); free(epsilon);

    /* sediment temperature model time series (sed_heat_model == 2) */
    if ( sed_heat_model == 2 && n_zones > 0 && theZones != NULL ) {
        int nsl = theZones[0].n_sed_layers, z, k;
        AED_REAL *st = malloc(n_zones * nsl * sizeof(AED_REAL));
        AED_REAL *hf = malloc(n_zones * sizeof(AED_REAL));
        AED_REAL *zt = malloc(n_zones * sizeof(AED_REAL));
        for (z = 0; z < n_zones; z++) {
            hf[z] = theZones[z].heatflux;
            zt[z] = theZones[z].ztemp;
            for (k = 0; k < nsl; k++) st[z*nsl + k] = theZones[z].layers[k].temp;
        }
        start[0] = set_no; edges[0] = 1;
        start[1] = 0;      edges[1] = n_zones;
        start[2] = 0;      edges[2] = nsl;
        check_nc_error(nc_put_vara(ncid, sed_temp_id, start, edges, st));
        check_nc_error(nc_put_vara(ncid, sed_hf_id,   start, edges, hf));
        check_nc_error(nc_put_vara(ncid, sed_zt_id,   start, edges, zt));
        free(st); free(hf); free(zt);
    }

    if ( n_zones > 0 ) {
        AED_REAL *zarea = malloc(n_zones * sizeof(AED_REAL));
        for (i = 0; i < n_zones; i++) zarea[i] = theZones[i].zarea;
        store_nc_array(ncid, zarea_id, XYNT_SHAPE, n_zones, n_zones, zarea);
        free(zarea);
    }

    check_nc_error(nc_sync(ncid));
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * Write the various output files and/or plots.                               *
 ******************************************************************************/
void write_glm_diag_ncdf(int ncid, AED_REAL LakeNum,
                    AED_REAL max_temp, AED_REAL min_temp, AED_REAL max_dtdz_at)
{
    extern AED_REAL Hs, L, T;
    AED_REAL lake_level = Lake[surfLayer].Height;

    if (ncid == -1) return;

    if (lake_level<0.011) lake_level=zero;

    /*------------------------------------------------------------------------*/
    start[2] = 0;      edges[2] = lon_len;
    start[1] = 0;      edges[1] = lat_len;
    start[0] = set_no; edges[0] = 1;

    xyt_store_nc_scalar(ncid, SA_id, Lake[surfLayer].LayerArea);
    xyt_store_nc_scalar(ncid, VSnow_id, SurfData.delzSnow*Lake[surfLayer].LayerArea*SurfData.RhoSnow/1e3);
    xyt_store_nc_scalar(ncid, VBIce_id, SurfData.delzBlueIce*Lake[surfLayer].LayerArea*917.0/1e3);
    xyt_store_nc_scalar(ncid, VWIce_id, SurfData.delzWhiteIce*Lake[surfLayer].LayerArea*890.0/1e3);
    xyt_store_nc_scalar(ncid, TotInVol_id, SurfData.dailyInflow);
    xyt_store_nc_scalar(ncid, TotOutVol_id, SurfData.dailyOutflow);
    xyt_store_nc_scalar(ncid, OverFVol_id, SurfData.dailyOverflow);
    xyt_store_nc_scalar(ncid, Evap_id, SurfData.dailyEvap);
    xyt_store_nc_scalar(ncid, Rain_id, SurfData.dailyRain);
    xyt_store_nc_scalar(ncid, LocRunoff_id, SurfData.dailyRunoff);
    xyt_store_nc_scalar(ncid, SnowF_id, SurfData.dailySnow);
    xyt_store_nc_scalar(ncid, Seep_id, SurfData.dailySeepage);
    xyt_store_nc_scalar(ncid, LakeLvl_id, lake_level);

//  write_csv_lake("Blue Ice Thickness", SurfData.delzBlueIce);
//  write_csv_lake("White Ice Thickness",SurfData.delzWhiteIce);
//  write_csv_lake("Snow Thickness",  SurfData.delzSnow);

    xyt_store_nc_scalar(ncid, SnowDns_id, SurfData.RhoSnow);
    xyt_store_nc_scalar(ncid, Alb_id, SurfData.albedo);
    xyt_store_nc_scalar(ncid, MaxT_id, max_temp);
    xyt_store_nc_scalar(ncid, MinT_id, min_temp);
    xyt_store_nc_scalar(ncid, SfT_id, Lake[surfLayer].Temp);
    xyt_store_nc_scalar(ncid, DQsw_id, SurfData.dailyQsw / Lake[surfLayer].LayerArea/SecsPerDay);
    xyt_store_nc_scalar(ncid, DQe_id, SurfData.dailyQe / Lake[surfLayer].LayerArea/SecsPerDay);
    xyt_store_nc_scalar(ncid, DQh_id, SurfData.dailyQh / Lake[surfLayer].LayerArea/SecsPerDay);
    xyt_store_nc_scalar(ncid, DQlw_id, SurfData.dailyQlw / Lake[surfLayer].LayerArea/SecsPerDay);
    xyt_store_nc_scalar(ncid, Light_id, Lake[surfLayer].Light);
    xyt_store_nc_scalar(ncid, BenLight_id, Benthic_Light_pcArea);
    xyt_store_nc_scalar(ncid, SWH_id, Hs);
    xyt_store_nc_scalar(ncid, SWL_id, L);
    xyt_store_nc_scalar(ncid, SWP_id, T);
    xyt_store_nc_scalar(ncid, LkNum_id, LakeNum);
    xyt_store_nc_scalar(ncid, maxdtz_id, max_dtdz_at);
    xyt_store_nc_scalar(ncid, CD_id, coef_wind_drag);
    xyt_store_nc_scalar(ncid, CHE_id, coef_wind_chwn);
    xyt_store_nc_scalar(ncid, zL_id, SurfData.dailyzonL*(noSecs/SecsPerDay));
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void close_glm_ncdf(int ncid)
{
   if (ncid != -1)
       check_nc_error(nc_close(ncid));

   set_no = -1;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
int new_nc_variable(int ncid, const char *name, int data_type,
                                                      int ndim, const int *dims)
{
   int id;

    if (ncid == -1) return -1;
    check_nc_error(nc_def_var(ncid,name,data_type,ndim,dims,&id));
    return id;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void set_nc_attributes(int ncid, int id, const char *units,
                                      const char *long_name, AED_REAL FillValue)
{
    if (ncid == -1) return;
    nc_put_att(ncid, id, "units", NC_CHAR, strlen(units), units);
    if ( long_name != NULL ) {
        nc_put_att(ncid, id, "long_name", NC_CHAR, strlen(long_name), long_name);
        nc_put_att(ncid, id, "_FillValue", NC_REALTYPE, 1, &FillValue);
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void store_nc_integer(int ncid, int id, int var_shape, int iscalar)
{
    int iret;

    if (ncid == -1) return;
    /*------------------------------------------------------------------------*/
    if (var_shape == POINT)
        iret = nc_put_var(ncid, id, &iscalar);
    else if (var_shape == T_SHAPE) {
        start[0] = set_no; edges[0] = 1;
        iret = nc_put_vara(ncid, id, start, edges, &iscalar);
    } else {
        fprintf(stderr, "store_nc_integer : non valid shape %d\n", var_shape);
        exit(1);
    }

    if ( iret != NC_NOERR ) check_nc_error_x(iret, ncid, id);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void store_nc_scalar(int ncid, int id, int var_shape, AED_REAL scalar)
{
    int iret;

    if (ncid == -1) return;
    /*------------------------------------------------------------------------*/
    if (var_shape == POINT)
        iret = nc_put_var(ncid, id, &scalar);
    else {
        if (var_shape == T_SHAPE) {
            start[0] = set_no; edges[0] = 1;
        } else if (var_shape == XYT_SHAPE) {
            start[2] = 0;      edges[2] = lon_len;
            start[1] = 0;      edges[1] = lat_len;
            start[0] = set_no; edges[0] = 1;
        } else {
            fprintf(stderr, "store_nc_scalar : non valid shape %d\n", var_shape);
            exit(1);
        }
        iret = nc_put_vara(ncid, id, start, edges, &scalar);
    }

    if ( iret != NC_NOERR ) check_nc_error_x(iret, ncid, id);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void xyt_store_nc_scalar(int ncid, int id, AED_REAL scalar)
{
    int iret;
    iret = nc_put_vara(ncid, id, start, edges, &scalar);

    if ( iret != NC_NOERR ) check_nc_error_x(iret, ncid, id);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void store_nc_array(int ncid, int id, int var_shape, int nvals,
                                                   int maxvals, AED_REAL *array)
{
    int iret, i;
    AED_REAL *tarr;

    if (ncid == -1) return;
    /*------------------------------------------------------------------------*/
    if (var_shape == Z_SHAPE) {
        start[0] = 0;      edges[0] = height_len;
    } else if (var_shape == XYZT_SHAPE || var_shape == XYNT_SHAPE) {
        start[3] = 0;      edges[3] = lon_len;
        start[2] = 0;      edges[2] = lat_len;
        if ( var_shape == XYZT_SHAPE ) {
            start[1] = 0;      edges[1] = height_len;
        } else {
            start[1] = 0;      edges[1] = n_zones;
        }
        start[0] = set_no; edges[0] = 1;
    } else {
        fprintf(stderr, "store_nc_array : non valid shape %d\n", var_shape);
        exit(1);
    }
    tarr = malloc(maxvals*sizeof(AED_REAL));
    for (i = 0; i < nvals; i++) tarr[i] = array[i];
    for (i = nvals; i < maxvals; i++) tarr[i] = NC_FILLER;
    iret = nc_put_vara(ncid, id, start, edges, tarr);
    free(tarr);

    if ( iret != NC_NOERR ) check_nc_error_x(iret, ncid, id);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


//#undef DEBUG
/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void check_nc_error(int err)
{
    if (err != NC_NOERR) {
        fprintf(stderr, "Error : %s (%d)\n", nc_strerror(err), err);
#if DEBUG
        CRASH("check_nc_error");
#endif
    }
}
/*----------------------------------------------------------------------------*/
/* This variant should only be called if err is not NC_NOERR                  */
void check_nc_error_x(int err, int ncid, int id)
{
    char name[NC_MAX_NAME+1];
    nc_inq_varname(ncid, id, name);
    fprintf(stderr, "Error : %s (%d) on variable %3d : \"%s\"\n", nc_strerror(err), err, id, name);
#if DEBUG
    CRASH("check_nc_error");
#endif
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 * for wq                                                                     *
 *                                                                            *
 ******************************************************************************/
void define_mode_on(int *ncid) { nc_redef(*ncid); }
void define_mode_off(int *ncid) { nc_enddef(*ncid); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int new_nc_variable_(int *ncid, const char *name, int *len, int *data_type,
                                                     int *ndim, const int *dims)
{
   int i, n = *ndim, ret;
   char *s = strndup(name,*len);
   int dim[6];

   for (i = 0; i < n; i++) dim[n-i-1] = dims[i];
   ret = new_nc_variable(*ncid, s, *data_type, n, dim);
   free(s);

   return ret;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void set_nc_attributes_(int *ncid, int *id, const char *units,
                                     const char *long_name, AED_REAL *FillValue)
{ set_nc_attributes(*ncid, *id, units, long_name, *FillValue); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void store_nc_scalar_(int *ncid, int *id, int *var_shape, AED_REAL *array)
{ store_nc_scalar(*ncid, *id, *var_shape, *array); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void store_nc_array_(int *ncid, int *id, int *var_shape, int *nvals,
                                                  int *maxvals, AED_REAL *array)
{ store_nc_array(*ncid, *id, *var_shape, *nvals, *maxvals, array); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
