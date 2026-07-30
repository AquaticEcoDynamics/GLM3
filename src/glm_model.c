/******************************************************************************
 *                                                                            *
 * glm_model.c                                                                *
 *                                                                            *
 *        run_model   This subroutine sets up the conditions for a            *
 *                    GLM simulation run, and controls the running            *
 *                    of the simulation.                                      *
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
#ifndef _WIN32
#include <unistd.h>
#endif
#include <string.h>
#include <math.h>
#include <time.h>

#include "glm.h"

#include "glm_types.h"
#include "glm_const.h"
#include "glm_globals.h"

#include "glm_util.h"
#include "glm_input.h"
#include "glm_output.h"
#include "glm_deep.h"
#include "glm_mixu.h"
#include "glm_mixer.h"
#include "glm_layers.h"
#include "glm_flow.h"
#include "glm_surface.h"
#include "glm_init.h"
#include "glm_lnum.h"
#include "glm_wqual.h"
#include "glm_ptm.h"
#include "glm_stress.h"
#include "glm_balance.h"
#include "glm_heatexchange.h"
#include "glm_oxygenation.h"
#include "glm_bubbler.h"
#include "glm_restart.h"
#ifdef PLOTS
#include <libplot.h>
#include "glm_plot.h"
#endif

#include "aed_time.h"
#include "aed_csv.h"

#include "glm_debug.h"

/******************************************************************************/
#if DEBUG
#  define dbgprt(...) fprintf(stderr, __VA_ARGS__)
//#  define dbgprt(...) /* __VA_ARGS__ */
#else
#  define dbgprt(...) /* __VA_ARGS__ */
#endif


extern int lw_ind;

static char EOLN = '\r';

#define mod(a,b) ((a) % (b))

/*----------------------------------------------------------------------------*/
#if defined(_WIN32)
__declspec(dllexport) void __cdecl init_model(int *jstart, int *nsave);
__declspec(dllexport) void __cdecl do_model_coupled(int step_start, int step_end,
        AED_REAL *FlowNew, AED_REAL *DrawNew, AED_REAL *elevation, int nsave);
__declspec(dllexport) void __cdecl end_model(void);
#else
void init_model(int *jstart, int *nsave);
void do_model_coupled(int step_start, int step_end,
        AED_REAL *FlowNew, AED_REAL *DrawNew, AED_REAL *elevation, int nsave);
void end_model(void);
#endif
void do_model(int jstart, int nsave);
void do_model_non_avg(int jstart, int nsave);
int do_subdaily_loop(int stepnum, int jday, int stoptime, int nsave, AED_REAL SWold, AED_REAL SWnew);

//int n_steps_done = 0;
//#define END_STEPS 30
int startTOD = 0;
int stopTOD = 0;
int nDates = 1;
int write_step;
int last_step;


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void run_model()
{
    int jstart, nsave, lth, ltm, lts, ltd;
    time_t begn, done;
    char buf[64];

    init_model(&jstart, &nsave);

    begn = time(NULL);
    if (quiet < 10) printf("\n     Wall clock start time :  %s", ctime_r(&begn, buf));

    if (non_avg)
        do_model_non_avg(jstart, nsave);
    else
        do_model(jstart, nsave);

    done = time(NULL);
    if (quiet < 10) printf("\n     Wall clock finish time : %s", ctime_r(&done, buf));
    ltd = difftime(done, begn);
    lth = ltd / 3600;
    ltm = (ltd - (lth * 3600)) / 60;
    lts = ltd - (lth * 3600) - (ltm * 60);
    if (quiet < 10)
        printf("     Wall clock runtime was %d seconds : %02d:%02d:%02d [hh:mm:ss]\n", ltd, lth, ltm, lts);

    end_model();
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void init_model(int *jstart, int *nsave)
{
    char out_dir[PATH_MAX];
    char out_fn[PATH_MAX];

/*----------------------------------------------------------------------------*/
#ifndef _WIN32
    // This is for the progress line :
    // If output is going to a file use a normal end of line - but if going to a
    // screen - use CR and overwrite the line.
    if (!isatty(fileno(stdout))) EOLN = '\n';
#endif

    init_glm(jstart, out_dir, out_fn, nsave);

#ifdef PLOTS
    psubday = timestep * (*nsave) / SecsPerDay;
    plotstep = 0;
#endif

    if ( bubbler_on ) init_bubbler();

    //# Create the output files.
    init_output(*jstart, out_dir, out_fn, MaxLayers, Longitude, Latitude);

    //# Calculate cumulative layer volumes from layer depths.
    //# Skip when state was loaded from a restart file: the restart now
    //# preserves LayerVol/LayerArea/Vol1/MeanHeight directly, and recomputing
    //# them here from Heights introduces ~1e-12 roundoff drift vs. the
    //# continuous run's stored values.
    if (!Restart_loaded) resize_internals(1, botmLayer);

    //# Check layers for vmax,vmin.
    //# Skip when state was loaded from a restart file: the saved state
    //# already had subdaily-level check_layer_thickness run; invoking it
    //# here on a partial-end restart would merge layers at an instant a
    //# continuous run never visits, breaking bit-for-bit resume.
    if (!Restart_loaded) check_layer_thickness();

    if(DepMX == 0.0) init_mixer();

    Latitude = two_Pi + Latitude * deg2rad; //# Convert latitude from degrees to radians
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
static void fix_radiation(AED_REAL Light_Surface)
{
    int i;

    Lake[surfLayer].Light = Light_Surface;
    for (i = surfLayer-1; i >= botmLayer; i-- )
        Lake[i].Light = Lake[i+1].Light *
               exp(-Lake[i+1].ExtcCoefSW * (Lake[i+1].Height - Lake[i].Height));
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
static AED_REAL calc_benthic_light()
{
    int i;
    AED_REAL Benthic_Light_Area, photic_depth, depth;


    // I/I0 = exp(-Kw*z);
    photic_depth = log(Benthic_Imin)/(-Kw);


    //# Calculate the percent benthic area where the light level is greater
    //# than the minimum level required for production
    Benthic_Light_Area = 0.;
    depth = 0.;
    for (i = surfLayer; i > botmLayer; i-- ) {
        depth = depth + (Lake[i].Height-Lake[i-1].Height);
        if (photic_depth>depth) {
          Benthic_Light_Area += (Lake[i].LayerArea - Lake[i-1].LayerArea);
      }
    }
  //  if (Lake[botmLayer].Light * exp(-Lake[botmLayer].ExtcCoefSW*Lake[botmLayer].Height) >= Benthic_Imin)
  //      Benthic_Light_Area = Benthic_Light_Area + Lake[botmLayer].LayerArea;

    return Benthic_Light_Area / Lake[surfLayer].LayerArea * 100./(SecsPerDay/noSecs);
  //return Benthic_Light_Area / Lake[surfLayer].LayerArea * 100.;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void calc_mass_temp(const char *msg)
{
    AED_REAL Lake_Mass; //# Total mass of lake [kg]
    AED_REAL Lake_Temp; //# Mass averaged lake temperature [oC]
    int i;

    Lake_Mass = zero;
    for (i = surfLayer; i >= botmLayer; i-- )
        Lake_Mass += Lake[i].Density * Lake[i].LayerVol;

    Lake_Temp = zero;
    for (i = surfLayer; i >= botmLayer; i-- )
        Lake_Temp += Lake[i].Temp * Lake[i].Density * Lake[i].LayerVol;
    Lake_Temp = Lake_Temp / Lake_Mass;

    if ( quiet < 5)
        printf("     %s Lake_Mass = %10.5f\t, Lake_Temp = %10.5f\n", msg, Lake_Mass/1e6, Lake_Temp);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void do_model(int jstart, int nsave)
{
    AED_REAL FlowNew[MaxInf], DrawNew[MaxOut], WithdrTempNew;
    AED_REAL FlowOld[MaxInf], DrawOld[MaxOut], WithdrTempOld;

    MetDataType MetOld, MetNew;
    AED_REAL    SWold, SWnew, DailyKw, DailyEvap;

   /***************************************************************************
    *CAB Note: these WQ arrays should be sized to Num_WQ_Vars not MaxVars,    *
    *           look into that later ....                                     *
    ***************************************************************************/
    AED_REAL SaltNew[MaxInf], TempNew[MaxInf], WQNew[MaxInf * MaxVars];
    AED_REAL SaltOld[MaxInf], TempOld[MaxInf], WQOld[MaxInf * MaxVars];
    // NOTE: despite the name, Elev[] holds inflow HEIGHT above the lake bottom
    // (same units as InflowDataType.SubmHeight / NML subm_height), NOT elevation
    // above sea level.  The name and the CSV column ("elev") are legacy from when
    // the field was called SubmElev.  The CSV column name is kept to avoid breaking
    // existing user inflow files.
    AED_REAL Elev[MaxInf];
    AED_REAL ElevOut[MaxOut];     // Outflow elevation (above sea level) for Type 6 submerged outflow
    AED_REAL HeatFluxOut[MaxOut]; // Outflow heat flux array for dynamic heat pump support
    int jday, ntot, stepnum, stoptime;
    int i, j;
    AED_REAL day_fraction;

    /*------------------------------------------------------------------------*/

    memset(WQNew, 0, sizeof(AED_REAL)*MaxInf*MaxVars);
    memset(WQOld, 0, sizeof(AED_REAL)*MaxInf*MaxVars);

    /* See note in do_model_non_avg: a mid-day restart resume needs the loaded
     * sub-daily state used verbatim; a midnight resume behaves like a normal
     * day boundary. startTOD is read before do_subdaily_loop resets it. */
    int resuming_midday = (Restart_loaded && startTOD > 0);

    /**************************** Start Simulation ****************************/
    fputs("\n     Simulation begins...\n", stdout);

    ntot = 0;
    stepnum = 0;
    stoptime = iSecsPerDay;

    read_daily_inflow(jstart, NumInf, FlowOld, TempOld, SaltOld, Elev, WQOld);
//  read_daily_gw(jstart, gw_mode, GWFlOld);
    read_daily_outflow(jstart, NumOut, DrawOld, ElevOut, HeatFluxOut, NULL); // No WQ for outflows yet
    read_daily_withdraw_temp(jstart, &WithdrTempOld);
    read_daily_met(jstart, &MetOld);
    MetData = MetOld;
    SWold = MetOld.ShortWave;
    /* On a mid-day resume, seed SWold from the previous run's stop instant so
     * iter-1's first sub-iter matches a continuous run (matters in non-subdaily
     * mode; calculate_qsw ignores SWold when subdaily). */
    if (resuming_midday) SWold = Restart_SWold;

    jday = jstart - 1;
#ifdef PLOTS
    today = jday;
#endif

//  write_output(jday, SecsPerDay, nsave, stepnum);
    /**************************************************************************
     * Loop over all days                                                     *
     **************************************************************************/
    while (ntot < nDates) {
        ntot++;
        jday++;
        yearday = day_of_year(jday);

        //# If it is the last day, adjust the stop time for the day if necessary
        if (ntot == nDates) stoptime = stopTOD;
        if (stoptime == 0) break;
        //# See note in do_model_non_avg: on a mid-day resume, apply the full
        //# day's end-of-iter lump (the previous run deferred the [0, startTOD]
        //# portion).
        if (resuming_midday && ntot == 1)
            day_fraction = stoptime / iSecsPerDay;
        else
            day_fraction = (stoptime - startTOD) / iSecsPerDay;

        //# Initialise daily values for volume & heat balance reporting (lake.csv)
        SurfData.dailyRain = 0.; SurfData.dailyEvap = 0.;
        SurfData.dailyQsw = 0.; SurfData.dailyQe = 0.;
        SurfData.dailyQh = 0.; SurfData.dailyQlw = 0.;
        SurfData.dailyInflow = 0.; SurfData.dailySnow = 0.;
        SurfData.dailyOutflow = 0.; SurfData.dailyOverflow = 0.;
        SurfData.dailySeepage = 0.; SurfData.dailyzonL = 0.;
        SurfData.dailyRunoff = 0.; SurfData.albedo = 1.;

        //# Read & set today's inflow properties
        read_daily_inflow(jday, NumInf, FlowNew, TempNew, SaltNew, Elev, WQNew);
//      read_daily_gw(jday, gw_mode, GWFlNew);
        read_daily_evap(jday, &DailyEvap);
        f_evap_ts_prop = DailyEvap / (noSecs * Lake[surfLayer].LayerArea);
        //# Averaging of flows
        //# To get daily inflow (i.e. m3/day) times by SecsPerDay
        //# (stoptime - startTOD) allow for partial dates at the the beginning and end of
        //# simulation
        for (i = 0; i < NumInf; i++) {
            Inflows[i].FlowRate = (FlowOld[i] + FlowNew[i]) / 2.0 * day_fraction * iSecsPerDay;
            Inflows[i].TemInf   = (TempOld[i] + TempNew[i]) / 2.0;
            Inflows[i].SalInf   = (SaltOld[i] + SaltNew[i]) / 2.0;
            Inflows[i].SubmHeight = Elev[i];
            for (j = 0; j < Num_WQ_Vars; j++) {
                Inflows[i].WQInf[j] = (WQ_INF_(WQOld,i, j) + WQ_INF_(WQNew, i, j)) / 2.0;
            }
            wq_inflow_update(Inflows[i].WQInf, &Num_WQ_Vars, &Inflows[i].TemInf, &Inflows[i].SalInf);
        }
        // Heat pump insert the captured values from outflow (AFTER averaging to avoid overwrite)
        heat_pump_insert_inflow();

        //# Read & set today's outflow properties
        read_daily_outflow(jday, NumOut, DrawNew, ElevOut, HeatFluxOut, NULL);
        //# Read today's oxygenation values (CSV-driven devices and recirc stream)
        read_daily_oxygenation(jday);
        //# To get daily outflow (i.e. m3/day) times by the seconds in the current day
        for (i = 0; i < NumOut; i++)
            Outflows[i].Draw = (DrawOld[i] + DrawNew[i]) / 2.0 * day_fraction * iSecsPerDay ;

        //# Update outflow elevations for Type 6 with dynamic elevation
        for (i = 0; i < NumOut; i++) {
            if (Outflows[i].Type == 6 && Outflows[i].SubmElevDynamic) {
                Outflows[i].SubmElev = ElevOut[i];
            }
        }
        //# Update heat pump heat flux if dynamic heat flux is available
        if (heat_pump_switch == 2 && heat_pump_outflow_idx >= 0 && heat_pump_outflow_idx < NumOut) {
            heat_pump_dynamic_heat_flux = HeatFluxOut[heat_pump_outflow_idx];
        }
        read_daily_withdraw_temp(jday, &WithdrTempNew);
        WithdrawalTemp = (WithdrTempOld + WithdrTempNew) / 2.0;

        read_daily_kw(jday, &DailyKw);
        /* On a mid-day resume's first iter, keep the loaded ExtcCoefSW (it
         * carries the bioshade_feedback contribution from the previous run);
         * see note in glm_restart.c "lake_extc". A midnight resume resets
         * normally, exactly as a continuous run does at a day boundary. */
        if (!(resuming_midday && ntot == 1))
            for (i = 0; i < MaxLayers; i++) Lake[i].ExtcCoefSW = DailyKw;

        //# Read & set today's meteorological data
        if (jday != jstart) read_daily_met(jday, &MetNew);
        if ( !subdaily ) {
            MetData.Rain =        (MetOld.Rain + MetNew.Rain) / 2.0;
            MetData.SatVapDef =   (MetOld.SatVapDef + MetNew.SatVapDef) / 2.0;
            MetData.LongWave =    (MetOld.LongWave + MetNew.LongWave) / 2.0;
            MetData.ShortWave =   (MetOld.ShortWave + MetNew.ShortWave) / 2.0;
            MetData.AirTemp =     (MetOld.AirTemp + MetNew.AirTemp) / 2.0;
            MetData.AirPres =     (MetOld.AirPres + MetNew.AirPres) / 2.0;
            MetData.WindSpeed =   (MetOld.WindSpeed + MetNew.WindSpeed) / 2.0;
            MetData.Snow =        (MetOld.Snow + MetNew.Snow) / 2.0;
            MetData.RainConcPO4 = (MetOld.RainConcPO4 + MetNew.RainConcPO4) / 2.0;
            MetData.RainConcTp =  (MetOld.RainConcTp + MetNew.RainConcTp) / 2.0;
            MetData.RainConcNO3 = (MetOld.RainConcNO3 + MetNew.RainConcNO3) / 2.0;
            MetData.RainConcNH4 = (MetOld.RainConcNH4 + MetNew.RainConcNH4) / 2.0;
            MetData.RainConcTn =  (MetOld.RainConcTn + MetNew.RainConcTn) / 2.0;
            MetData.RainConcSi =  (MetOld.RainConcSi + MetNew.RainConcSi) / 2.0;
        }
        SWnew = MetNew.ShortWave;

        //# Now enter into sub-daily calculations

        stepnum = do_subdaily_loop(stepnum, jday, stoptime, nsave, SWold, SWnew);

        //# End of forcing-mixing-diffusion loop

        /* See note in do_model_non_avg: skip end-of-iter state mutations on a
         * partial-end last iter when a restart will be written, so the writer
         * and the resumer produce identical pre-end-of-iter state. */
        int is_partial_end_with_restart =
            (ntot == nDates) && (stoptime != iSecsPerDay) && (restart_fname != NULL);

        if (!is_partial_end_with_restart) {
        //# Insert inflows for all streams
        SurfData.dailyInflow = do_inflows(); //# Do inflow for all streams

        //# Extract withdrawal from all offtakes
        SurfData.dailyOutflow = do_outflows(jday, day_fraction);

        //# Oxygenation: re-inject recirculated water (A3) and add direct O2 (A1/A2)
        oxy_do_recirculation(day_fraction);
        do_oxygenation(day_fraction);

        //# Take care of any overflow
        SurfData.dailyOverflow = do_overflow(jday, day_fraction);

        //# Enforce layer limits
        check_layer_thickness();
        }

        // # Write output on last time step within a day
        // # after including the inflow and outflows.
        // # Output is not written on the last time step in a daily in the subdaily loop
        if ( stepnum == write_step){

#ifdef PLOTS
            today = jday;
#endif
            write_output(jday, SecsPerDay, nsave, stepnum);
            write_diags(jday, calculate_lake_number());
            write_balance(jday);
            write_step += nsave;
#ifdef PLOTS
            plotstep++;
            today = -1;
#endif

        }

        if ( (ntot == nDates) && (stepnum < nsave)) {
            fprintf(stderr, "     ERROR: NO netcdf output generated because nsave is less total number of time steps in simuluation\n");
        }


        /**********************************************************************
         * End of daily calculations, Prepare for next day and return.        *
         **********************************************************************/
        for (i = 0; i < NumInf; i++) {
            FlowOld[i] = FlowNew[i];
            TempOld[i] = TempNew[i];
            SaltOld[i] = SaltNew[i];
            for (j = 0; j < Num_WQ_Vars; j++)
                WQ_INF_(WQOld, i, j) = WQ_INF_(WQNew, i, j);
        }
        for (i = 0; i < MaxOut; i++) DrawOld[i] = DrawNew[i];
        WithdrTempOld = WithdrTempNew;
        MetOld = MetNew;
        SWold = SWnew;

#ifdef XPLOTS
        if ( xdisp )
            flush_all_plots();
        else
#endif
          if (quiet < 2) {
            printf("     Running day %8d [DoY %3ld], %4.2f%% of days complete%c",
                   jday, lrint(trunc(yearday)), ntot*100./nDates, EOLN);
            fflush(stdout);
        }
    }   //# do while (ntot < nDates)
    if (quiet < 2) { printf("\n"); fflush(stdout); }
    /*----------########### End of main daily loop ################-----------*/
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void do_model_non_avg(int jstart, int nsave)
{
    AED_REAL FlowNew[MaxInf], DrawNew[MaxOut], WithdrTempNew;
    AED_REAL SWold, SWnew, DailyKw, DailyEvap;
    int jday, ntot, stepnum, stoptime;
    int i, j;
    AED_REAL day_fraction;

   /***************************************************************************
    *CAB Note: these WQ arrays should be sized to Num_WQ_Vars not MaxVars,    *
    *           look into that later ....                                     *
    ***************************************************************************/
    AED_REAL SaltNew[MaxInf], TempNew[MaxInf], WQNew[MaxInf * MaxVars];
    // NOTE: despite the name, Elev[] holds inflow HEIGHT above the lake bottom
    // (same units as InflowDataType.SubmHeight / NML subm_height), NOT elevation
    // above sea level.  The name and the CSV column ("elev") are legacy from when
    // the field was called SubmElev.  The CSV column name is kept to avoid breaking
    // existing user inflow files.
    AED_REAL Elev[MaxInf];
    AED_REAL ElevOut[MaxOut];     // Outflow elevation (above sea level) for Type 6 submerged outflow
    AED_REAL HeatFluxOut[MaxOut]; // Outflow heat flux array for dynamic heat pump support

    /*------------------------------------------------------------------------*/
    memset(WQNew, 0, sizeof(AED_REAL)*MaxInf*MaxVars);

    /* "Mid-day resume": a restart was loaded AND the run starts part-way
     * through a calendar day (startTOD>0). Only then does iter-1 continue a
     * calendar day the previous run began — requiring the loaded sub-daily
     * state to be used verbatim (no top-of-day resets) and the day's
     * end-of-iter lump applied in full. A midnight resume (startTOD==0) begins
     * a fresh calendar day and must behave exactly like a continuous run's
     * day boundary, so none of these overrides apply. startTOD is read here
     * before do_subdaily_loop resets it to 0. */
    int resuming_midday = (Restart_loaded && startTOD > 0);

    /**************************** Start Simulation ****************************/
    fputs("\n     Simulation begins...\n", stdout);

    ntot = 0;
    stepnum = 0;
    stoptime = iSecsPerDay;
    SWold = 0.;
    /* On a mid-day resume, seed SWold from the previous run's stop instant so
     * iter-1's first sub-iter matches a continuous run (matters in non-subdaily
     * mode; calculate_qsw ignores SWold when subdaily). */
    if (resuming_midday) SWold = Restart_SWold;

    jday = jstart - 1;

//  write_output(jday, SecsPerDay, nsave, stepnum);
    /**************************************************************************
     * Loop over all days                                                     *
     **************************************************************************/
    while (ntot < nDates) {
        ntot++;
        jday++;

        //# If it is the last day, adjust the stop time for the day if necessary
        if (ntot == nDates) stoptime = stopTOD;
        if (stoptime == 0) break;
        //# day_fraction scales the daily inflow/outflow/seepage/overflow lump
        //# applied at the end of this calendar-day iter. On a mid-day resume the
        //# previous run simulated the [0, startTOD] part of this day but deferred
        //# its end-of-iter lump (partial-end ops skipped), so apply the full day's
        //# lump here (as a continuous run does at this midnight) rather than only
        //# the post-resume fraction.
        if (resuming_midday && ntot == 1)
            day_fraction = stoptime / iSecsPerDay;
        else
            day_fraction = (stoptime - startTOD) / iSecsPerDay;

        //# Initialise daily values for volume & heat balance reporting (lake.csv)
        SurfData.dailyRain    = 0.; SurfData.dailyEvap     = 0.;
        SurfData.dailyQsw     = 0.; SurfData.dailyQe       = 0.;
        SurfData.dailyQh      = 0.; SurfData.dailyQlw      = 0.;
        SurfData.dailyInflow  = 0.; SurfData.dailySnow     = 0.;
        SurfData.dailyOutflow = 0.; SurfData.dailyOverflow = 0.;
        SurfData.dailySeepage = 0.; SurfData.dailyzonL     = 0.;
        SurfData.dailyRunoff  = 0.; SurfData.albedo        = 1.;

        //# Read & set today's inflow properties
        read_daily_inflow(jday, NumInf, FlowNew, TempNew, SaltNew, Elev, WQNew);
//      read_daily_gw(jday, gw_mode, GWFlNew);
        read_daily_evap(jday, &DailyEvap);
        f_evap_ts_prop = DailyEvap / (noSecs * Lake[surfLayer].LayerArea);

        //# To get daily inflow (i.e. m3/day) times by SecsPerDay
        for (i = 0; i < NumInf; i++) {
            Inflows[i].FlowRate = FlowNew[i] * day_fraction * iSecsPerDay;
            Inflows[i].TemInf   = TempNew[i];
            Inflows[i].SalInf   = SaltNew[i];
            Inflows[i].SubmHeight = Elev[i];
            for (j = 0; j < Num_WQ_Vars; j++) {
                Inflows[i].WQInf[j] = WQ_INF_(WQNew, i, j);
            }
            wq_inflow_update(Inflows[i].WQInf, &Num_WQ_Vars, &Inflows[i].TemInf, &Inflows[i].SalInf);
        }

        //# Read & set today's outflow properties
        read_daily_outflow(jday, NumOut, DrawNew, ElevOut, HeatFluxOut, NULL); // No WQ for outflows yet
        //# Read today's oxygenation values (CSV-driven devices and recirc stream)
        read_daily_oxygenation(jday);
        //# To get daily outflow (i.e. m3/day) times by the seconds in the current day
        //# (stoptime - startTOD) allow for partial dates at the the beginning and end of
        //# simulation
        for (i = 0; i < NumOut; i++)
            Outflows[i].Draw = DrawNew[i] * day_fraction * iSecsPerDay ;

        //# Update outflow elevations for Type 6 with dynamic elevation
        for (i = 0; i < NumOut; i++) {
            if (Outflows[i].Type == 6 && Outflows[i].SubmElevDynamic) {
                Outflows[i].SubmElev = ElevOut[i];
            }
        }

        //# Update heat pump heat flux if dynamic heat flux is available
        if (heat_pump_switch == 2 && heat_pump_outflow_idx >= 0 && heat_pump_outflow_idx < NumOut) {
            heat_pump_dynamic_heat_flux = HeatFluxOut[heat_pump_outflow_idx];
        }

        read_daily_withdraw_temp(jday, &WithdrTempNew);
        WithdrawalTemp = WithdrTempNew;

        //# Read & set today's Kw (if it is being read in)
        read_daily_kw(jday, &DailyKw);
        /* On a mid-day resume's first iter, keep the loaded ExtcCoefSW (it
         * carries the bioshade_feedback contribution from the previous run);
         * see note in glm_restart.c "lake_extc". A midnight resume resets
         * normally, exactly as a continuous run does at a day boundary. */
        if (!(resuming_midday && ntot == 1))
            for (i = 0; i < MaxLayers; i++) Lake[i].ExtcCoefSW = DailyKw;

        //# Read & set today's meteorological data
        read_daily_met(jday, &MetData);
        SWnew = MetData.ShortWave;

        //# Now enter into sub-daily calculations

        stepnum = do_subdaily_loop(stepnum, jday, stoptime, nsave, SWold, SWnew);

        //# End of forcing-mixing-diffusion loop

        /* Skip end-of-iter state-mutating ops on a partial-end last iter when
         * a restart will be written. Otherwise the post-subdaily check_layer_thickness
         * (and overflow) mutate state at a mid-day instant that a continuous run
         * never visits — breaking bit-for-bit restart. The same skip is applied
         * on the writer and the resumer so they produce identical restart files. */
        int is_partial_end_with_restart =
            (ntot == nDates) && (stoptime != iSecsPerDay) && (restart_fname != NULL);

        if (!is_partial_end_with_restart) {
         //# Insert inflows for all streams
        SurfData.dailyInflow = do_inflows();

        if (Lake[surfLayer].Vol1 > zero) {
           //# Extract withdrawal from all offtakes
           SurfData.dailyOutflow = do_outflows(jday, day_fraction);

           //# Oxygenation: re-inject recirculated water (A3) and add direct O2 (A1/A2)
           oxy_do_recirculation(day_fraction);
           do_oxygenation(day_fraction);

           //# Take care of any overflow
           SurfData.dailyOverflow = do_overflow(jday, day_fraction);
        }

        //# Enforce layer limits
        check_layer_thickness();
        }

        // # Write output on last time step within a day
        // # after including the inflow and outflows.
        // # Output is not written on the last time step in a daily in the subdaily loop
        if (stepnum == write_step) {
#ifdef PLOTS
            today = jday;
#endif
            write_output(jday, SecsPerDay, nsave, stepnum);
            write_diags(jday, calculate_lake_number());
            write_balance(jday);
            write_step += nsave;

            //if ( write_step > last_step ) write_step = last_step;
#ifdef PLOTS
            plotstep++;
            today = -1;
#endif
        }

        if ( (ntot == nDates) && (stepnum < nsave)) {
            fprintf(stderr, "     ERROR: NO netcdf output generated because nsave is less total number of time steps in simuluation\n");
        }

        /**********************************************************************
         * End of daily calculations, Prepare for next day and return.        *
         **********************************************************************/
        SWold = SWnew;

#ifdef XPLOTS
        if ( xdisp )
            flush_all_plots();
        else
#endif
          if (quiet < 2) {
            printf("     Running day %8d, %4.2f%% of days complete%c", jday, ntot*100./nDates, EOLN);
            fflush(stdout);
        }
    }   //# do while (ntot < nDates)
    if (quiet < 2) { printf("\n"); fflush(stdout); }
    /*----------########### End of main daily loop ################-----------*/
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void do_model_coupled(int step_start, int step_end,
           AED_REAL *FlowNew, AED_REAL *DrawNew, AED_REAL *elevation, int nsave )
{
    AED_REAL SWold, SWnew, DailyKw;

   /***************************************************************************
    *CAB Note: these WQ arrays should be sized to Num_WQ_Vars not MaxVars,    *
    *           look into that later ....                                     *
    ***************************************************************************/
    AED_REAL WQNew[MaxInf * MaxVars];
    int jday, ntot, stepnum, stoptime, cDays;
    int i, j;
    AED_REAL day_fraction;

    /*------------------------------------------------------------------------*/
    memset(WQNew, 0, sizeof(AED_REAL)*MaxInf*MaxVars);

    /**************************** Start Simulation ****************************/
    fputs("\n     Simulation begins...\n", stdout);

    ntot = 0;
    stepnum = 0;
    stoptime = iSecsPerDay;
    SWold = 0.;



    cDays = step_end - step_start + 1;
    jday = step_start - 1;

//   write_output(jday, SecsPerDay, nsave, stepnum);
    /**************************************************************************
     * Loop over all days                                                     *
     **************************************************************************/
    while (ntot < cDays) {
        ntot++;
        jday++;

        //# If it is the last day, adjust the stop time for the day if necessary
        if (ntot == nDates) stoptime = stopTOD;
        if (stoptime == 0) break;
        day_fraction = (stoptime - startTOD) / iSecsPerDay;

        //# Initialise daily values for volume & heat balance reporting (lake.csv)
        SurfData.dailyRain = 0.; SurfData.dailyEvap = 0.;
        SurfData.dailyQsw = 0.; SurfData.dailyQe = 0.;
        SurfData.dailyQh = 0.; SurfData.dailyQlw = 0.;
        SurfData.dailyInflow = 0.; SurfData.dailySnow = 0.;
        SurfData.dailyOutflow = 0.; SurfData.dailyOverflow = 0.;
        SurfData.dailySeepage = 0.; SurfData.dailyzonL = 0.;
        SurfData.dailyRunoff = 0.; SurfData.albedo = 1.;

        //# Read & set today's inflow properties
    //  read_daily_inflow(jday, NumInf, FlowNew, TempNew, SaltNew, WQNew);
        //# Set to today's inflow
        //# To get daily outflow (i.e. m3/day) times by the seconds in the current day
        //# (stoptime - startTOD) allow for partial dates at the the beginning and end of
        //# simulation
        for (i = 0; i < NumInf; i++) {
            Inflows[i].FlowRate = FlowNew[i] * day_fraction * iSecsPerDay;
//          Inflows[i].TemInf   = TempNew[i];
//          Inflows[i].SalInf   = SaltNew[i];
            for (j = 0; j < Num_WQ_Vars; j++) {
                Inflows[i].WQInf[j] = WQ_INF_(WQNew, i, j);
            }
            wq_inflow_update(Inflows[i].WQInf, &Num_WQ_Vars, &Inflows[i].TemInf, &Inflows[i].SalInf);
        }

    //  read_daily_outflow(jday, NumOut, DrawNew);
        //# To get daily outflow (i.e. m3/day) times by SecsPerDay
        for (i = 0; i < NumOut; i++)
            Outflows[i].Draw = DrawNew[i] * day_fraction * iSecsPerDay;

    //  read_daily_withdraw_temp(jday, &WithdrTempNew);
    //  WithdrawalTemp = WithdrTempNew;

        //# Read & set today's Kw (if it is being read in)
        read_daily_kw(jday, &DailyKw);
        for (i = 0; i < MaxLayers; i++) Lake[i].ExtcCoefSW = DailyKw;

        //# Read & set today's meteorological data
        read_daily_met(jday, &MetData);
        SWnew = MetData.ShortWave;

        //# Now enter into sub-daily calculations

        stepnum = do_subdaily_loop(stepnum, jday, stoptime, nsave, SWold, SWnew);

        //# End of forcing-mixing-diffusion loop

         //# Insert inflows for all streams
        SurfData.dailyInflow = do_inflows();

        //# Read today's oxygenation values (CSV-driven devices and recirc stream)
        read_daily_oxygenation(jday);

        if (Lake[surfLayer].Vol1 > zero) {
           //# Extract withdrawal from all offtakes
           SurfData.dailyOutflow = do_outflows(jday, day_fraction);

           //# Oxygenation: re-inject recirculated water (A3) and add direct O2 (A1/A2)
           oxy_do_recirculation(day_fraction);
           do_oxygenation(day_fraction);

           //# Take care of any overflow
           SurfData.dailyOverflow = do_overflow(jday, day_fraction);
        }

        //# Enforce layer limits
        check_layer_thickness();

        // # Write output on last time step within a day
        // # after including the inflow and outflows.
        // # Output is not written on the last time step in a daily in the subdaily loop
        if (stepnum == write_step) {
#ifdef PLOTS
            today = jday;
#endif
            write_output(jday, SecsPerDay, nsave, stepnum);
            write_diags(jday, calculate_lake_number());
            write_balance(jday);
            write_step += nsave;

            //if ( write_step > last_step ) write_step = last_step;
#ifdef PLOTS
            plotstep++;
            today = -1;


#endif
        }

        if ( (ntot == nDates) && (stepnum < nsave)) {
            fprintf(stderr, "     ERROR: NO netcdf output generated because nsave is less total number of time steps in simuluation\n");
        }

        /**********************************************************************
         * End of daily calculations, Prepare for next day and return.        *
         **********************************************************************/
        SWold = SWnew;

#ifdef XPLOTS
        if ( xdisp )
            flush_all_plots();
        else
#endif
        if (quiet < 2) {
            printf("     Running day %8d, %4.2f%% of days complete%c", jday, ntot*100./nDates, EOLN);
            fflush(stdout);
        }


    }   //# do while (ntot < nDates)
    if (quiet < 2) { printf("\n"); fflush(stdout); }
    /*----------########### End of main daily loop ################-----------*/

    *elevation = Lake[surfLayer].Height;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
int do_subdaily_loop(int stepnum, int jday, int stoptime, int nsave, AED_REAL SWold, AED_REAL SWnew)
{
    int iclock;  //# The seconds counter during a day
    AED_REAL Light_Surface; //# Light at the surface of the lake after do_surface
    AED_REAL part_day_per_step;

    yearday = day_of_year(jday);
    noSecs = timestep;
    coef_wind_drag = CD;
    part_day_per_step = timestep / SecsPerDay;
    MetData.AirPres = atm_pressure_sl;

    if ( bubbler_on ) read_bubbler(jday);

    /**************************************************************************
     *  Loop for each second in a day (86400 = #seconds in a day)             *
     **************************************************************************/
    iclock = startTOD;
    last_step = stepnum + ((stoptime - iclock) / noSecs);
    if(stepnum == 0) write_step = stepnum + nsave;

    //if ( iclock != 0 ) {
    //    int t = mod(((stoptime - iclock) / noSecs), nsave);
    //    if ( t != 0 ) write_step = stepnum + t;
    //}

//  printf("last step %d write_step %d\n", last_step, write_step);

    startTOD = 0; /* from now on start at the beginning of the day */
    Benthic_Light_pcArea = 0.;
    while (iclock < stoptime) { //# iclock = seconds counter
        if ( subdaily ) {
//fprintf(stderr, "1 MetData.AirPres %f\n", MetData.AirPres);
            read_sub_daily_met(jday, iclock, &MetData);
            SWnew = MetData.ShortWave;
        }
//fprintf(stderr, "2 MetData.AirPres %f\n", MetData.AirPres);

        stepnum++;
        _dbg_time(jday, iclock);

        //# Thermal transfers are done by do_surface_thermodynamics
        do_surface_thermodynamics(jday, iclock, lw_ind, Latitude, SWold, SWnew);

        //# Save surface light to use at end of sub-daily time loop
        Light_Surface = Lake[surfLayer].Light/0.45;

        //# Mixing is done by do_mixing
        if ( surface_mixing > 0 )
            do_mixing();

//      calc_mass_temp("After do_mixing");

        //# Mix out instabilities, combine/split  layers
        check_layer_thickness();

        fix_radiation(Light_Surface);

        if ( surface_mixing > -1 ){
             check_layer_stability();
             fix_radiation(Light_Surface);
        }

        // flag in &glm_setup (int deep_mixing 0 = off, >0 = on)
        if ( deep_mixing > 0 ) {
            //# Estimate dissipation from energy inputs, buoyancy frequency etc
            if (NumLayers > 3) do_dissipation();

//          calc_mass_temp("After do_dissipation");

            //# Do deep mixing integrations
            //# If reservoir is mixed (NumLayers<3) then skip deep mixing
            if (NumLayers > 3) do_deep_mixing();
        }

            //# Check mixed layers for volume
            check_layer_thickness();
        fix_radiation(Light_Surface);

        //# Calculate the percent benthic area where the light level is greater
        //# than the minimum level required for production
        Benthic_Light_pcArea += calc_benthic_light();

        calc_layer_stress(MetData.WindSpeed,
                      sqrt( (Lake[surfLayer].LayerArea)/Pi ) * 2 );

        /**********************************************************************
         *## Start PTM calls                                                  *
         **********************************************************************/
        if (ptm_sw) do_ptm_update();

        /**********************************************************************
         *## Start Water Quality calls                                        *
         **********************************************************************/
        if (wq_calc) wq_do_glm(&NumLayers);

        //# If an output write is requested for the last time step of the day
        //# then do not output in the subdaily.  Output writing is moved to the
        //# daily loop so it occurs after the inflow and output calculations.

        if ( stepnum == write_step && (iclock +  noSecs) != SecsPerDay && stepnum != last_step) {

#ifdef PLOTS
            today = jday;
#endif
            write_output(jday, iclock, nsave, stepnum);
            write_step += nsave;
#ifdef PLOTS
//if (++n_steps_done > END_STEPS) { int i; for (i = 0; i < NumLayers; i++) show_l_line(2, Lake[i].Height); flush_all_plots(); }
            plotstep++;
            today = -1;
#endif
        }

//      calc_mass_temp("End sub_daily");

        //#If sub-daily re-set SWold
        if ( subdaily ) SWold = SWnew;

        if ( bubbler_on ) do_bubbler(jday, iclock);

//      printf("stepnum %d\n", stepnum);
        iclock += noSecs;
        yearday += part_day_per_step;
    }   //# do while (iclock < iSecsPerDay)

    //if ( write_step > last_step ) write_step = last_step;
    /**************************************************************************
     * End of sub-daily loop                                                  *
     **************************************************************************/
//  printf("end subdaily loop\n");

    /* Capture the SWold value at loop exit — the value a continuing simulation
     * would use as SWold on the next sub-iter. Saved to the restart so a
     * non-subdaily mid-day resume matches a continuous run. (In subdaily mode
     * SWold is unchanged inside the loop, so this just preserves it.) */
    Restart_SWold = SWold;

#ifdef PLOTS
    plotstep = 0;
#endif

    if ( bubbler_on ) write_bubbler(jday);

    return stepnum;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void end_model()
{
    fputc('\n', stdout);

    if ( bubbler_on ) print_bubbler();

    close_kw_files();
    close_met_files();
    close_inflow_files();
    close_outflow_files();
    close_withdrtemp_files();

    /* Write the final NetCDF restart file BEFORE any WQ state is freed.
     * write_glm_restart() reads WQ_Vars / WQS_Vars / WQD_Vars, which
     * wq_clean_glm() deallocates, so it must run first (otherwise the restart
     * write is a use-after-free — harmless for the small pelagic array, but a
     * hard segfault once the larger cc_diag/WQD_Vars array is read). */
    if (restart_fname != NULL)
        write_glm_restart(restart_fname);

    if (wq_calc) wq_clean_glm();    //# deallocataes wq stuff

    if ( bubbler_on ) done_bubbler();

    close_output();
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
