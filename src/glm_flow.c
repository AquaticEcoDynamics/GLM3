/******************************************************************************
 *                                                                            *
 * glm_flow.c                                                                 *
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
 *  Algorithms implmented in the below functions are adapted from those       *
 *  originally developed by the Centre for Water Research, University of      *
 *  Western Australia. Readers are referred to Imberger and Patterson (1981)  *
 *  for a comprehensive overview of original algorithm approaches of the model*
 *                                                                            *
 *  Imberger, J. and Patterson, J.C., 1981. A dynamic reservoir simulation    *
 *    model-DYRESM:5. In: H.B. Fisher (ed.), Transport Models for Inland      *
 *    and Coastal Waters. Academic Press, New York: 310-361.                  *
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
#include "glm_wqual.h"

#include "glm_util.h"
#include "glm_mixu.h"
#include "glm_layers.h"
#include "glm_output.h"
#include "glm_heatexchange.h"

#include "glm_balance.h"
#include "glm_ptm.h"

#include "glm_debug.h"

//#define _WQ_VarsTmp(i,j,k)  WQ_VarsTmp[_IDX_3d(Num_WQ_Vars,NumInf,MaxPar,i,j,k)]

#if DEBUG
#  define dbgprt(...) fprintf(stderr, __VA_ARGS__)
//#  define dbgprt(...) /* __VA_ARGS__ */
#else
#  define dbgprt(...) /* __VA_ARGS__ */
#endif

typedef AED_REAL wq_vars_t[MaxVars];
typedef wq_vars_t wq_partic_t[MaxPar];

/******************************************************************************/
static int flowing_depth(AED_REAL *HF, AED_REAL *DL, AED_REAL VOL,
                                  AED_REAL DEPTH, AED_REAL Phi, AED_REAL Alpha);
static AED_REAL extra_volume(AED_REAL d1, AED_REAL d2,
                                 AED_REAL SurfaceHeight, InflowDataType *Inflow);
static void do_single_outflow(AED_REAL HeightOfOutflow, AED_REAL flow, OutflowDataType *outf);
static int new_storage(int *iRiv);
AED_REAL delta_volume(AED_REAL z1, AED_REAL z2, AED_REAL da, AED_REAL avdel,
                 AED_REAL hh, AED_REAL delt, AED_REAL delb);
/*----------------------------------------------------------------------------*/

static AED_REAL *WQ_VarsS = NULL;
static wq_partic_t *WQ_VarsTmp = NULL;

static AED_REAL *Delta_V = NULL; //# The delta V from each layer taken by outflow

static int checkjday = -1;

LOGICAL seepage = FALSE;
AED_REAL seepage_rate = 0.0;
AED_REAL hBot;           //# Height of bottom of withdrawal layer
AED_REAL hTop;           //# Height of top of withdrawal layer

static FILE **sc_files = NULL;
static AED_REAL SpecialConditionDraw(int jday, int i);

/******************************************************************************
 * Removes the outflow at level Outflow_LayerNum                              *
 ******************************************************************************/
void do_single_outflow(AED_REAL HeightOfOutflow, AED_REAL flow, OutflowDataType *outf)
{
    const AED_REAL min_flow = 0.0001;

    AED_REAL LenAtOutflow;   //# Reservoir length at height of withdrawal
    AED_REAL WidthAtOutflow; //# Reservoir width at height of withdrawal
    AED_REAL ThickOutflow;   //# Thickness of withdrawal layer
    AED_REAL DeltaBot;       //# Lower half withdrawal layer thickness
    AED_REAL DeltaTop;       //# Upper half withdrawal layer thickness
    AED_REAL DeltaTotal;
    AED_REAL DeltaAvg;       //#
    AED_REAL Delta_rho;      //# Density difference between outflow layer i above or below
    AED_REAL Delta_z;        //# Width over which density gradient is calculated
    AED_REAL Q_outf;         //# Flow rate of the outflow, converted to m3/s
    AED_REAL F2,F3;          //* Froude number powers
    AED_REAL Grashof;        //# Grashof number
    AED_REAL R;              //# Outflow parameter
    AED_REAL DeltaSinkDepth; //# Level above W.L. from which fluid is drawn; length scale water above will drop
    AED_REAL HeightLayerBot; //# Height of GLM water layer immediately below the withdrawal thckness
    AED_REAL HeightLayerTop; //# Height of GLM water layer at the top of the withdrawal thickness
    AED_REAL Q_outf_star;    //# Flow rate of the outflow estimated from summing dV's
    AED_REAL Viscos;         //# Vertical diffusivity of momentum averaged over the withdrawal thickness
    AED_REAL Nsqrd_outflow;  //# Brunt-Vaisala frequency squared

    int Layer_i;             //# Layer count either above or below outflow height
    int i,iBot,iTop;
    int Outflow_half;        //# Which half of outflow calculating +1 = top half, -1 = bottom half
    int Outflow_LayerNum;    //# Layer number of outflow height
    int ILP;
    int flag;                //# Used to describe the nature of the withdraw

/*----------------------------------------------------------------------------*/
//BEGIN
    if (outf != NULL) outf->DrawnFrom = -1;       // layer index where last drawn

    if (HeightOfOutflow < 0 && outf->Type != 0 ) {
        fprintf(stderr, "HeightOfOutflow < 0; Outflow type is %d\n", outf->Type);
        exit(1);
    }

    //# Find number of layer (Outflow_LayerNum) opposite offtake
    for (i = botmLayer; i <= surfLayer; i++)
        if (Lake[i].Height >=  HeightOfOutflow) break;

    //# Return if reservoir surface is below outlet level
    if (i > surfLayer) return;

    Outflow_LayerNum = i;
//  printf("HeightOfOutflow is %10.5f; Outflow type is %d in layer %d\n",HeightOfOutflow, outf->Type, Outflow_LayerNum);

    if (outf != NULL) {                                                     // Check if outflow structure pointer is valid (not NULL)
        outf-> DrawnFrom = i;                                               // Store the layer index (i) where water was last drawn from for this outflow

        // Collect WQ variables from the layer being drawn from Type 6 outflow
        if (outf->Type == 6 && WQ_Vars != NULL && Num_WQ_Vars > 0) {
            int wqidx;
            for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++) {
                if (outf->WQ_Outflow != NULL && i >= 0 && i < MaxLayers) {
                    outf->WQ_Outflow[wqidx] = _WQ_Vars(wqidx, i);
                }
            }
        }
    }

    WidthAtOutflow = 0.;
    LenAtOutflow = 0.;
    ThickOutflow = 0.;
    DeltaTop = 0.;
    Delta_z = 0.;
    ILP = FALSE;

    HeightLayerBot = zero;   HeightLayerTop = zero;

    //# Reset Delta_V (layer specific withdrawal amount) for all layers as zero
    for (i = botmLayer; i < MaxLayers; i++) Delta_V[i]=0.0;

    /***********************************************************************
     * Determine offtake level, length, width and flow                     *
     *                                                                     *
     * withdrawal height selection                                         *
     * switch control by FloatOff                                          *
     * Assume:                                                             *
     *      1) that lake approximates as an ellipse                        *
     *      2) Area = pi/4 * Length * Width                                *
     *      3) ratio Length:Width at outflow is same as crest              *
     *                                                                     *
     ***********************************************************************/
    if (outf == NULL) {
        if ( HeightOfOutflow > zero) { //# Is an overflow so use crest height and width
            LenAtOutflow = LenAtCrest;
            WidthAtOutflow = WidAtCrest;
        } else { //# Seepage use bottom layer (Eq. x & x in GLM manual)
            LenAtOutflow = sqrt(Lake[Outflow_LayerNum].LayerArea*4.0/Pi*(LenAtCrest/WidAtCrest));
            WidthAtOutflow = LenAtOutflow * WidAtCrest/LenAtCrest;
        }
    } else {
        if (outf->FloatOff) { // Floating offtake
            //# Assume offtake level close to surface
            LenAtOutflow = Lake[surfLayer].LayerArea / WidAtCrest;
            WidthAtOutflow = WidAtCrest;
        } else { //# Fixed Offtake (Eq. x & x in GLM manual)
            LenAtOutflow = sqrt(Lake[Outflow_LayerNum].LayerArea*4.0/Pi*(LenAtCrest/WidAtCrest));
            WidthAtOutflow = LenAtOutflow * WidAtCrest/LenAtCrest;
        }
    }

    Q_outf = flow / SecsPerDay;

    if (flow <= min_flow) return;

    //# Calculate withdrawal layer thickness and volume removed from each layer
    if ((outf == NULL && HeightOfOutflow == zero) || surfLayer == botmLayer) {
        //# Is a seepage or single layer lake, therefore take all from bottom layer
        Delta_V[0] = flow;

        //# Correction if bottom layer should be emptied. Note that for a single layer
        //# lake, less than the specified outflow may be removed and this is what is
        //# noted in the lake.csv diagnostic file
        if (surfLayer == botmLayer) {
            if (Delta_V[0] >= Lake[0].LayerVol)
                Delta_V[0] = Lake[0].LayerVol - 1.0;
            Delta_V[0] = 0.9 * Lake[0].LayerVol;
        } else {
            for (i = botmLayer; i < surfLayer; i++) {
                if (Delta_V[i] >= Lake[i].LayerVol) {
                    Delta_V[i+1] = Delta_V[i+1] + Delta_V[i] - 0.9*(Lake[i].LayerVol);
                    Delta_V[i] = 0.9 * (Lake[i].LayerVol);
                }
                //printf("Delta_V[i]p = %d, %10.5f,%10.5f,%10.5f\n",i, Lake[i].LayerVol,Delta_V[i],Delta_V[i+1]);
             }
        }
    } else {
        Outflow_half = 1;      //# Start with withdraw from layers above outlet
        Layer_i = Outflow_LayerNum + 1;
        while (1) {
            if (Layer_i > surfLayer) Layer_i = surfLayer;
            if (Layer_i < botmLayer) Layer_i = botmLayer;
            Delta_rho = Lake[Outflow_LayerNum].Density - Lake[Layer_i].Density;

            flag = FALSE;
            if (Delta_rho * Outflow_half <= zero)
                //# Density of outflow layer less than layer i above or greater than layer i below
                //# i.e. unstable density
                flag = TRUE;
            else {
                Delta_z = Lake[Layer_i].MeanHeight - HeightOfOutflow;
                if ((Delta_z * Outflow_half) <= zero)
                    //# Counting up and Delta_z is negative or counting down and Delta_z is positive
                    flag = TRUE;
                else { //# Calculate the Brunt-Vaisala frequency of outflow (Eqn x in GLM Manual)
                    Nsqrd_outflow = Delta_rho * g / (Lake[Outflow_LayerNum].Density * Delta_z);
                    if (Nsqrd_outflow <= zero) flag = TRUE;
                }
            }

            if (flag) {
                //# check for zero gradient
                if (Outflow_half == 1) //# Taking from layers above only
                    ThickOutflow = Lake[surfLayer].Height - HeightOfOutflow;
                else //# Taking from layers below only
                    ThickOutflow = HeightOfOutflow;
            } else {
                //# vertical diffusivity of momentum averaged over withdraw layer
                Viscos = Lake[Outflow_LayerNum].Epsilon * 20.0;
                //# limit to molecular diffusivity
                if (Viscos <= zero) Viscos = Visc;
                //# Grashof number (Eq. x in GLM manual)
                Grashof = Nsqrd_outflow * sqr(LenAtOutflow) * sqr(LenAtOutflow) / sqr(Viscos);
                flag = TRUE;

                if (!ILP) {
                    //# Point sink case
                    F3 = Q_outf/(sqrt(Nsqrd_outflow)*LenAtOutflow*sqr(LenAtOutflow));
                    ThickOutflow = LenAtOutflow* pow(F3, (1.0/3.0));
                    if ((2.0 * ThickOutflow <= WidthAtOutflow) || (Outflow_half != 1)) {
                        R = F3*sqrt(Grashof);
                        if (R <= 1.0)
                            ThickOutflow = 2.9 * (pow((Viscos*Visc), (1.0/6.0)) *
                                             (pow((LenAtOutflow/sqrt(Nsqrd_outflow)), (1.0/3.0))));
                        flag = FALSE;
                    }
                }

                if ( flag ) {
                    //# Line sink case
                    //# Froude number (Eq. x in GLM Manual)
                    F2 = Q_outf/WidthAtOutflow/sqrt(Nsqrd_outflow)/sqr(LenAtOutflow);
                    ILP = TRUE;
                    ThickOutflow = 2.0 * LenAtOutflow * sqrt(F2);
                    //# R parameter (Eq. x in GLM Manual)
                    R = F2 * pow(Grashof, (1.0/3.0));
                    //# Thickness of outflow (Eq. x in GLM manual)
                    if (R <= 1.0) ThickOutflow = 2.0*LenAtOutflow/pow(Grashof, (1.0/6.0));
                }
            }

            /******************************************************************
             * Check that the withdrawal layer does not intersect the top or  *
             * bottom of the reservoir, and that ThickOutflow is less than    *
             * the range over which the density gradient is calculated.       *
             ******************************************************************/
            if (Outflow_half == 1) {
                //# check if upper bound of withdrawal layer is reached
                if (ThickOutflow <= fabs(Delta_z) || Layer_i == surfLayer) {
                    DeltaTop = ThickOutflow;
                    if (DeltaTop > (Lake[surfLayer].Height - HeightOfOutflow) )
                        DeltaTop = Lake[surfLayer].Height - HeightOfOutflow;
                    Outflow_half = -1;
                    Layer_i = Outflow_LayerNum-1;
                } else {
                    ILP = FALSE;
                    Layer_i++;
                }
            } else if (ThickOutflow <= fabs(Delta_z) || Layer_i == botmLayer) {
                DeltaBot = ThickOutflow;
                if (DeltaBot > HeightOfOutflow) DeltaBot = HeightOfOutflow;
                break;
            } else
                Layer_i--;
        }


        /**********************************************************************
         * Force a hard check of layer thickness                              *
         **********************************************************************/
        if (DeltaTop>outflow_thick_limit) DeltaTop=outflow_thick_limit;
        if (DeltaBot>outflow_thick_limit) DeltaBot=outflow_thick_limit;


        /**********************************************************************
         * Calculate top and bottom of withdrawal layer and distance above    *
         * withdrawal layer from which fluid is drawn.                        *
         **********************************************************************/
        hTop = HeightOfOutflow + DeltaTop;
        hBot = HeightOfOutflow - DeltaBot;

        DeltaTotal = DeltaTop + DeltaBot;
        DeltaAvg = DeltaTotal / 2.0;
        DeltaSinkDepth = flow / (WidthAtOutflow * LenAtOutflow);
        hTop += DeltaSinkDepth;
        if (hTop >  Lake[surfLayer].Height) hTop = Lake[surfLayer].Height;

        //# Find the indices of the layers containing hTop and hBot
        // first, locate iBot
        for (i = botmLayer; i <= Outflow_LayerNum; i++)
            if (Lake[i].Height  >=  hBot) break;
        // if hBot higher than the bottom of the sink, note an error.
        if (i > Outflow_LayerNum)
            fprintf(stderr,"ERROR: do_outflows - bottom of withdrawal layer above outlet bottom\n");
        // otherwise set layer index
        iBot = i;

        // now, locate iTop
        for (i = Outflow_LayerNum; i <= surfLayer; i++)
            if (Lake[i].Height >= hTop) break;

        // and set layer index
        iTop = i;

        //# Assign the volumes to be taken, after checking if all drawn from one layer
        if (iBot == iTop || single_layer_draw )
            // entire flow from the chosen layer matching the draw height
            Delta_V[Outflow_LayerNum] = flow;
        else {
            // calculate Delta_V[i], the portion of water taken from the ith layer
            for (i = iBot; i <= iTop; i++) {
                if (i != botmLayer) HeightLayerBot = Lake[i-1].Height;
                if (i == iBot) HeightLayerBot = hBot;
                HeightLayerTop = Lake[i].Height;
                if (i == iTop) HeightLayerTop = hTop;
                Delta_V[i] = delta_volume(HeightLayerBot, HeightLayerTop,
                                DeltaSinkDepth, DeltaAvg, HeightOfOutflow, DeltaTop, DeltaBot);
            }

            //# Proportion drawn from each layer is known. Now match the
            //# total volume (flow) to the actual volume drawn out.

            Q_outf_star = zero;
            for (i = botmLayer; i <= surfLayer; i++)
                Q_outf_star += Delta_V[i];

            for (i = botmLayer; i <= surfLayer; i++)
                Delta_V[i] = (flow / Q_outf_star) * Delta_V[i];
        }

        //# Correction if any layer should be emptied.
        for (i = botmLayer; i < surfLayer; i++) {
            if (Delta_V[i] >= Lake[i].LayerVol) {
                Delta_V[i+1] = Delta_V[i+1] + (Delta_V[i] - 0.9*Lake[i].LayerVol);
                Delta_V[i] = 0.9*(Lake[i].LayerVol) ;
            }
        }
    }

    /**********************************************************************
     * Now we have Delta_V[i] for all layers we can remove it             *
     **********************************************************************/
    for (i = botmLayer; i <= surfLayer; i++){
        if (Delta_V[i] > zero) {
/*
            printf("%d DeltaV %8.4f; flow %10.4f;%10.4f %d %d %d %10.1f %10.1f \n",
               i,Delta_V[i],flow,Q_outf_star,Outflow_LayerNum,iBot,iTop,hBot,hTop);
*/
            if ( ptm_sw ) {
                ptm_removeparticles(i, Delta_V[i], Lake[i].LayerVol, max_particle_num);
            }
            Lake[i].LayerVol -= Delta_V[i];
        }
        mb_sub_outflows(i, Delta_V[i]);
    }

    /**********************************************************************
     * Recompute volumes                                                  *
     **********************************************************************/
    Lake[botmLayer].Vol1 = Lake[botmLayer].LayerVol;
    for (i = (botmLayer+1); i <= surfLayer; i++)
        Lake[i].Vol1 = Lake[i-1].Vol1 + Lake[i].LayerVol;

    /**********************************************************************
    * Update layer heights                                                *
    **********************************************************************/
    resize_internals(2, botmLayer);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * Withdraw a volume of water at a given height using GLM's selective-         *
 * withdrawal physics, for the oxygenation recirculation (Approach 3).         *
 * Samples the withdrawal-layer T/S/WQ into the output args and returns the     *
 * actual volume removed (which the caller re-injects elsewhere with added O2).*
 ******************************************************************************/
AED_REAL oxy_recirc_extract(AED_REAL height, AED_REAL want_vol,
                            AED_REAL *temp, AED_REAL *salt, AED_REAL *wq)
{
    int lvl, wqidx;
    AED_REAL vol_before, drawn;
    OutflowDataType tmp;

    if (want_vol <= zero) return zero;

    //# Find the layer at the withdrawal height and sample its properties as
    //# representative of the withdrawn water (same approximation the heat pump uses).
    for (lvl = botmLayer; lvl <= surfLayer; lvl++)
        if (Lake[lvl].Height >= height) break;
    if (lvl > surfLayer) lvl = surfLayer;

    if (temp != NULL) *temp = Lake[lvl].Temp;
    if (salt != NULL) *salt = Lake[lvl].Salinity;
    if (wq != NULL)
        for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
            wq[wqidx] = _WQ_Vars(wqidx, lvl);

    //# Run GLM's selective-withdrawal physics with a temporary fixed outlet;
    //# it recomputes Vol1 and resizes layers internally.
    memset(&tmp, 0, sizeof(tmp));
    tmp.Type = 1;
    tmp.FloatOff = FALSE;
    tmp.WQ_Outflow = NULL;

    vol_before = Lake[surfLayer].Vol1;
    do_single_outflow(height, want_vol, &tmp);
    drawn = vol_before - Lake[surfLayer].Vol1;
    if (drawn < zero) drawn = zero;
    return drawn;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * Loop through all outflows and process - return the difference between      *
 *  total volume before and after.                                            *
 ******************************************************************************/
AED_REAL do_outflows(int jday, AED_REAL day_fraction)
{
    int i;
    AED_REAL DrawHeight = -1; //# Height of withdraw [m from bottom]
    AED_REAL VolSum = Lake[surfLayer].Vol1; //# Lake volume prior to outflows [m3]
    AED_REAL SeepDraw = 0.0; //# Seepage volume [m3]

    if ( Delta_V == NULL ) Delta_V = calloc(MaxLayers, sizeof(AED_REAL));

    /**************************************************************************
     * Do withdrawal for each offtake                                         *
     * Include water quality and particles                                    *
     * offtake type is switch controlled by FloatOff                          *
     **************************************************************************/
    for (i = 0; i < NumOut; i++) {
        AED_REAL tVolSum = Lake[surfLayer].Vol1;

        /**********************************************************************
         * Outlet type     .                                                  *
         * Type 1 is fixed outlet heights                                     *
         * Type 2 is floating offtake                                         *
         * Type 3 is fixed outlet heights + check for crit. hypol. oxygen     *
         * Type 4 is variable outlet heights for ISOTHERM                     *
         * Type 5 is variable outlet heights for TEMPERATURE TIME SERIES      *
         **********************************************************************/

        /**********************************************************************
         * Outlet type     .                                                  *
         * Type 1 is fixed outlet heights                                     *
         **********************************************************************/
        if (Outflows[i].Type == 1) {
            DrawHeight = Outflows[i].OLev;  //# Fixed height offtake

        /**********************************************************************
         * Floating offtake.                                                  *
         * OLev(i) is distance below surface level for the offtake            *
         **********************************************************************/
        } else if ( Outflows[i].Type == 2 ) {
            DrawHeight = Lake[surfLayer].Height - Outflows[i].OLev;
            //# Let it sit on the bottom if the lake has drained
            if (DrawHeight < 0.) DrawHeight = 0.11;
        /***********************************************************************
         * Type 6 is submerged outflow (fixed or dynamic elevation from bottom)*
         * Special additions by Taynara Fernandes                              *
         **********************************************************************/
        } else if (Outflows[i].Type == 6) {
            DrawHeight = Outflows[i].SubmElev;  // Use submerged elevation
        /**********************************************************************
         * 3, 4 and 5 are the special additions by Michael Weber              *
         **********************************************************************/
        } else
            DrawHeight = SpecialConditionDraw(jday, i);

        Outflows[i].Draw *= Outflows[i].Factor;

        do_single_outflow(DrawHeight, Outflows[i].Draw, &Outflows[i]);
        // Set LastDrawn for diagnostic calculations
        Outflows[i].LastDrawn = Outflows[i].Draw;
        // DrawHeight is layer where particles are; if know # of particles in

        //Heat_pump captures the outflow
        if (heat_pump_switch > 0 && i == heat_pump_outflow_idx) {
            // substep discharge
            AED_REAL qout = Outflows[i].Draw;   /* m3/day */
            // find the layer that has the withdrawal elevation
            int lvl;
            for (lvl = botmLayer; lvl <= surfLayer; ++lvl) {   //Look up and finds withdrawl elevation
                if (Lake[lvl].Height >= DrawHeight) break;
            }
            // clamp in case DrawHeight is above the current surface
            if (lvl > surfLayer) lvl = surfLayer;

            AED_REAL Tin = Lake[lvl].Temp;
            AED_REAL Sin = Lake[lvl].Salinity;

            // Use WQ variables already stored in outflow structure if available (Type 6), otherwise collect fresh
            if (Outflows[i].Type == 6 && Outflows[i].WQ_Outflow != NULL && Num_WQ_Vars > 0) {
                // Use pre-collected/stored WQ data from Type 6 outflow structure
                heat_pump_capture_outflow(jday, DrawHeight, qout, Tin, Sin, Outflows[i].WQ_Outflow);
            } else if (WQ_Vars != NULL && Num_WQ_Vars > 0) {
                // Collect WQ variables fresh (for non-Type 6 outflows or if WQ_Outflow is NULL)
                AED_REAL wq_vars[MaxVars];
                for (int wqidx = 0; wqidx < Num_WQ_Vars; wqidx++) {
                    wq_vars[wqidx] = _WQ_Vars(wqidx, lvl);
                }
                heat_pump_capture_outflow(jday, DrawHeight, qout, Tin, Sin, wq_vars);
            } else {
                // No WQ variables available
                heat_pump_capture_outflow(jday, DrawHeight, qout, Tin, Sin, NULL);
            }
        }

        write_outflow(i, jday, DrawHeight, tVolSum-Lake[surfLayer].Vol1, Outflows[i].Draw, hBot, hTop);
    }
    if (seepage) {
        if (seepage_rate>zero) {
            // Darcy's Law used, so input rate is hydraulic conductivity (m/day) x hydraulic head
            SeepDraw = seepage_rate * Lake[surfLayer].Height * Lake[surfLayer].LayerArea * 0.95 * day_fraction;
        } else {
            // Constant seepage assumed, so input rate is dh (m/day)
            // 0.95 added since the effective area of seeping is probably
            // a bit less than max area of water innundation???
            SeepDraw = -seepage_rate * Lake[surfLayer].LayerArea * 0.95 * day_fraction;
        }
        SurfData.dailySeepage += SeepDraw;
        do_single_outflow(0., SeepDraw, NULL);
    }

    return VolSum - Lake[surfLayer].Vol1;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * Calculate overflow.  Overflow occurs when the current reservoir            *
 * volume is greater than the volume of the reservoir at the crest            *
 * level.  Add in stack volumes when calculating overflow.                    *
 * After overflow, delete stack volumes from the structure.                   *
 ******************************************************************************/
AED_REAL do_overflow(int jday, AED_REAL day_fraction)
{
    AED_REAL VolSum = Lake[surfLayer].Vol1;
    AED_REAL DrawHeight = 0.;
    AED_REAL overflow = 0.;

    // Too much water for the entire lake domain, remove this first


    if (VolSum > MaxVol){
        do_single_outflow((CrestHeight+(MaxHeight-CrestHeight)*0.9), (VolSum - MaxVol), NULL);
        overflow = VolSum - Lake[surfLayer].Vol1;
    }
    VolSum = Lake[surfLayer].Vol1;
    // Water above the crest, which will overflow based on a weir equation
    if (VolSum > VolAtCrest){
        AED_REAL ovfl_Q, ovfl_dz;

        ovfl_dz = MAX( Lake[surfLayer].Height - CrestHeight, zero );
        ovfl_Q = 2./3. * crest_factor * pow(2*g,0.5) * crest_width * pow(ovfl_dz,1.5) * day_fraction * iSecsPerDay;
        ovfl_Q  = MIN( (VolSum - VolAtCrest) , ovfl_Q );

        do_single_outflow(CrestHeight, ovfl_Q , NULL);
        overflow += ovfl_Q ;
    }

    write_outflow(MaxOut, jday, DrawHeight, overflow, zero, CrestHeight, MaxHeight);

    return overflow;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
static int insert_inflow(int k, //# Inflow parcel counter
    int iRiver,                 //# River inflow number
    AED_REAL Height_t0,         //# Height of lake before inflows [m]
    AED_REAL Alpha,             //# Stream half angle [radians]
    AED_REAL Phi,               //# Stream slope [radians]
    AED_REAL EntrainmentCoeff,  //# Entrainment coefficient [-]
    AED_REAL Ri)                //# Bulk Richardson number of the inflow [-]
{
    int kk;
    int flag;
    int wqidx;                 //# WQ variable column id's
    int Layer;                 //# Layer number adjacent to inflow parcel, entrainment layer
    int Layer_t0;              //# Layer number at which inflow parcel starts the day
    AED_REAL Delta_Q;          //# Delta Q the rate of entrainment of inflow aliquot [ML/day]- MHcheck unit
    AED_REAL Delta_t;          //# Delta t the time taken for downflow [days]
    AED_REAL Downflow_Depth;   //# Depth of downflow [m]
    AED_REAL Inflow_Energy;    //# Energy of inflow
    AED_REAL Inflow_time;      //# Time taken for inflow to reach neutral buoyancy or bottom of lake [days]
    AED_REAL Inflow_dx;        //# Distance travelled by inflow [m]
    AED_REAL Inflow_gprime;    //# Reduced gravity (gprime) between inflow parcel and adjacent layer [-]
    AED_REAL Inflow_height_prev; //# Thickness of inflow from previous layer [m]
    AED_REAL Inflow_height;    //# Thickness of inflow, post entrainment [m]
    AED_REAL Inflow_Flowrate;  //# Inflow flowrate [ML/day] - MHcheck unit
    AED_REAL Inflow_Velocity;  //# Inflow velocity [m/s]
    AED_REAL Inflow_Temp;      //# Inflow temperature [oC]
    AED_REAL Inflow_Salinity;  //# Inflow salinity [psu]
    AED_REAL Inflow_Density;   //# Inflwo density [kg/m3]

/*----------------------------------------------------------------------------*/
//BEGIN

    //# Check to see if have reached the number of inflow parcels for this day.
    if (k >= Inflows[iRiver].iCnt) return FALSE;

    Inflow_Energy = zero;

    //# If fresh inflow for that day set depth as lake height
    if (Inflows[iRiver].DDown[k] == MISVAL) {
        Inflows[iRiver].TotIn += Inflows[iRiver].QDown[k];
        Inflows[iRiver].DDown[k] = Height_t0;
    }

    //# Initialise inflow time
    Inflow_time = zero;

    //# Set the depth the inflow aliquot reached on previous step/day (DOld)
    Inflows[iRiver].DOld[k] = Inflows[iRiver].DDown[k];

    //# Calculate inflow density for this inflow aliquot
    Inflow_Density = calculate_density(Inflows[iRiver].TDown[k], Inflows[iRiver].SDown[k]);

    //# Calculate the layer in which this element starts the day.
    for (Layer_t0 = botmLayer; Layer_t0 <= surfLayer; Layer_t0++)
        if (Lake[Layer_t0].Height >= Inflows[iRiver].DDown[k]) break;

    if (Layer_t0 > surfLayer) Layer_t0 = surfLayer; //# Limit to surface layer
    Layer = Layer_t0;

    //# Loop for progression of the parcel through the next layer down
    Inflow_Flowrate = Inflows[iRiver].QDown[k];
    Inflow_Temp = Inflows[iRiver].TDown[k];
    Inflow_Salinity = Inflows[iRiver].SDown[k];
    for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
        WQ_VarsS[wqidx] = Inflows[iRiver].WQDown[k][wqidx];

    mb_add_inflows(Inflows[iRiver].QDown[k], Inflow_Temp, Inflow_Salinity, WQ_VarsS);

    Downflow_Depth = Inflows[iRiver].DDown[k];

    //# Check if this parcel lies below level of neutral buoyancy.
    //# If so, then add to insertion queue, and renumber the stack of parcels
    if (Inflow_Density <= Lake[Layer].Density){

        Inflows[iRiver].InPar[Inflows[iRiver].NoIns] = k;
        Inflows[iRiver].QIns[Inflows[iRiver].NoIns]  = Inflow_Flowrate;
        Inflows[iRiver].TIns[Inflows[iRiver].NoIns]  = Inflow_Temp;
        Inflows[iRiver].SIns[Inflows[iRiver].NoIns]  = Inflow_Salinity;
        Inflows[iRiver].DIIns[Inflows[iRiver].NoIns] = Inflow_Density;
        for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
            Inflows[iRiver].WQIns[Inflows[iRiver].NoIns][wqidx] = WQ_VarsS[wqidx];

        Inflows[iRiver].NoIns++;

        return TRUE;
    }

    //# Otherwise keep moving the parcel down into the lake
    //# Calculate the height and velocity of the inflow and then the entrainment
    while (1) {
        //# Reduced gravity of downflow (Eq. x in GLM manual)
        Inflow_gprime = g*(Inflow_Density-Lake[Layer].Density)/Lake[Layer].Density;

        //# Initial height of the inflow as it plunges (Eq. x in GLM manual)
        Inflow_height_prev = pow((2.0*Ri*pow((Inflow_Flowrate/SecsPerDay*cos(Alpha)/sin(Alpha)), 2) / Inflow_gprime), 0.2);

        //# Distance travelled by inflow aliquot (Eq. x in GLM manual)
        if (Layer == botmLayer) Inflow_dx = Downflow_Depth / sin(Phi);
        else                    Inflow_dx = (Downflow_Depth - Lake[Layer-1].Height) / sin(Phi);

        //# New inflow thickness due to entrainment (Eq. x in GLM manual)
        Inflow_height = 1.2 * EntrainmentCoeff * Inflow_dx  +  Inflow_height_prev;

        //# Inflow velocity in m/day (Eq. x in GLM manual)
        Inflow_Velocity = Inflow_Flowrate * cos(Alpha) / (pow(Inflow_height, 2) * sin(Alpha));

        //# Time taken for inflow aliquot to travel in this day
        Delta_t = Inflow_dx/Inflow_Velocity;

        //# If time is greater than one day then carry on to next day and only increment this day's entrainment
        //# to calculation of new inflow layer thickness
        if ((Inflow_time + Delta_t) > 1.0) {
            Inflow_dx = Inflow_dx * (1.0-Inflow_time) / Delta_t;
            Delta_t = 1.0-Inflow_time;
            Inflow_height = 1.2 * EntrainmentCoeff * Inflow_dx  +  Inflow_height_prev;
        }

        Inflow_time += Delta_t;

        //# Estimate increase in flow rate due to entrainment for this day (Eq. x in GLM manual)
        //Delta_Q = 0.2 * Inflow_Flowrate * (pow((Inflow_height/Inflow_height_prev), (5.0/3.0)) - 1.0);
        Delta_Q = Inflow_Flowrate * (pow((Inflow_height/Inflow_height_prev), (5.0/3.0)) - 1.0);

        //# Check for negative inflow layer volume
        if (Lake[Layer].LayerVol < 0.0)
            fprintf(stderr, "Vol[%d] is negative = %f, surfLayer = %d\n", Layer, Lake[Layer].LayerVol, surfLayer);

        //# Check that the entrainment is less than 90% of the layer volume
        if (Delta_Q > 0.9*Lake[Layer].LayerVol) Delta_Q = 0.9 * Lake[Layer].LayerVol;

        //# Determine physical properties of the inflow aliquot once water
        //  from adjacent layer has been entrained
        Inflow_Salinity = combine(Inflow_Salinity, Inflow_Flowrate,
                             Inflow_Density, Lake[Layer].Salinity, Delta_Q, Lake[Layer].Density);
        Inflow_Temp = combine(Inflow_Temp,
                             Inflow_Flowrate, Inflow_Density, Lake[Layer].Temp, Delta_Q, Lake[Layer].Density);
        Inflow_Density = calculate_density(Inflow_Temp, Inflow_Salinity);

        //# Add entrained water to inflow aliquot and take from adjacent layer
        Inflow_Flowrate += Delta_Q;
        Lake[Layer].LayerVol -= Delta_Q;

        //# Update other water attributes
        for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
            WQ_VarsS[wqidx] = combine_vol(WQ_VarsS[wqidx], Inflow_Flowrate, _WQ_Vars(wqidx, Layer), Delta_Q);

        //# Calculate energy of inflowing streams.
        if (Layer == botmLayer)
            Inflow_Energy += (Inflow_Density - Lake[Layer].Density) * Downflow_Depth * g * Inflow_Flowrate / SecsPerDay;
        else
            Inflow_Energy += (Inflow_Density - Lake[Layer].Density) * (Downflow_Depth - Lake[Layer-1].Height) * g * Inflow_Flowrate / SecsPerDay;

        //# Update the downflowing parcel with the new attributes
        Inflows[iRiver].QDown[k] = Inflow_Flowrate;
        Inflows[iRiver].TDown[k] = Inflow_Temp;
        Inflows[iRiver].SDown[k] = Inflow_Salinity;
        for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
            Inflows[iRiver].WQDown[k][wqidx] = WQ_VarsS[wqidx];

        Inflows[iRiver].TotIn = Inflows[iRiver].TotIn + Delta_Q;
        Inflows[iRiver].DDown[k] = Downflow_Depth - Inflow_dx * sin(Phi);
        Downflow_Depth = Inflows[iRiver].DDown[k];

        if (Inflows[iRiver].DDown[k] < Inflows[iRiver].Dlwst)
            Inflows[iRiver].Dlwst = Inflows[iRiver].DDown[k];

        //# If the inflow parcel has ended it's days travel, reached the
        //# level of neutral buoyancy or has reached the bottom, put it in
        //# the insertion queue
        flag = TRUE;
        if (Layer != botmLayer)
            if (Inflow_Density > Lake[Layer-1].Density) flag = FALSE;

        if (Inflow_time >= 1.0)  flag = TRUE ;

        if ( flag ) {
            Inflows[iRiver].InPar[Inflows[iRiver].NoIns] = k;
            Inflows[iRiver].TIns[Inflows[iRiver].NoIns]  = Inflow_Temp;
            Inflows[iRiver].SIns[Inflows[iRiver].NoIns]  = Inflow_Salinity;
            Inflows[iRiver].DIIns[Inflows[iRiver].NoIns] = Inflow_Density;
            Inflows[iRiver].QIns[Inflows[iRiver].NoIns]  = Inflow_Flowrate;
            for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
                Inflows[iRiver].WQIns[Inflows[iRiver].NoIns][wqidx] = WQ_VarsS[wqidx];

            Inflows[iRiver].NoIns++;
        }

        //# If the inflow parcel has ended it's days travel, reached its level of
        //# neutral buoyancy or the bottom of the reservoir, go to the next parcel
        flag = FALSE;
        if (Layer == botmLayer)
            flag = TRUE;
        else if (Inflow_time >= 1.0)
            flag = TRUE;
        else if (Inflow_Density <= Lake[Layer-1].Density)
            flag = TRUE;

        if ( flag ) { //#Insert inflow parcel and re-calculate cumulative volumes, then return to do_inflows
            einff += Inflow_Energy/Inflow_time;
            if (Layer == botmLayer) Lake[botmLayer].Vol1 = Lake[botmLayer].LayerVol;
            if (Layer != botmLayer) Lake[Layer].Vol1 = Lake[Layer-1].Vol1 + Lake[Layer].LayerVol;
            for (kk = Layer+1; kk <= surfLayer; kk++)
                 Lake[kk].Vol1 = Lake[kk-1].Vol1 + Lake[kk].LayerVol;

            return TRUE;
        }
        Layer--;
    }

    return FALSE;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * Inserts the inflow at the level of neutral buoyancy, after calculating the *
 * entrainment                                                                *
 * Return the difference between total volume before and after.               *
 ******************************************************************************/
AED_REAL do_inflows()
{
    const AED_REAL cwnsq2 = 39.48;

    AED_REAL Alpha;             //# Stream half angle in radians
    AED_REAL Phi;               //# Stream slope in radians
    AED_REAL DragCoeff;         //# Stream drag coefficient
    AED_REAL EntrainmentCoeff;  //# Entrainment coefficient
    AED_REAL Ri;                //# Bulk Richardson number of the inflow [-]
    AED_REAL Inflow_Density;    //# Density of the inflow
    AED_REAL Inflow_Height_t0;  //# Initial inflow plunge point [m]
    AED_REAL Inflow_gprime_t0;  //# Reduced gravity (gprime) given surface and inflow densities

    int j, k, ll, jk;
    int iRiver;
    int Layer_subm;             //# Layer number to insert submerged inflow

    int wqidx;

    AED_REAL Inflow_width;      //# Width of inflow [m]
    AED_REAL VolSum = Lake[surfLayer].Vol1; //# Total lake volume before inflows

    int new_particles;
    AED_REAL height_start = Lake[surfLayer].Height; //# Total lake height before inflows
    AED_REAL double_particles;
    AED_REAL upper_height, lower_height;

/*----------------------------------------------------------------------------*/
//BEGIN

    if ( NumInf <= 0 ) return zero; // nothing more to do

    // Array allocation for WQ in the inflow parcels
    if ( WQ_VarsS == NULL && Num_WQ_Vars > 0)
         WQ_VarsS = calloc(Num_WQ_Vars, sizeof(AED_REAL));
    if ( WQ_VarsTmp == NULL ) WQ_VarsTmp = calloc(NumInf, sizeof(wq_partic_t));


    /**************************************************************************
     * Initially:                                                             *
     * 1) Calculate initial reduced gravity and plunge point of the inflow.   *
     * 2) Calculate WaveNumSquared and vel for use in the calculation of deep *
     * water mixing. See do_deep_mixing where they're used to calculate the   *
     * dispersion coefficient, Epsilon.                                       *
     * Count down from smallest river, overwriting each time, in case river   *
     * has zero inflow.                                                       *
     **************************************************************************/
     WaveNumSquared = zero;
     vel = zero;

     for (iRiver = NumInf-1; iRiver >= 0; iRiver--) {
        if (Inflows[iRiver].FlowRate*Inflows[iRiver].Factor != zero && !Inflows[iRiver].SubmFlag) {
            Alpha = Inflows[iRiver].Alpha;
            Phi = Inflows[iRiver].Phi;
            Inflow_Density = calculate_density(Inflows[iRiver].TemInf,Inflows[iRiver].SalInf);
            //# Reduced gravity of inflow (Eq. X in GLM manual)
            Inflow_gprime_t0 = gprime(Lake[surfLayer].Density,Inflow_Density);

            if (Inflow_gprime_t0 > zero) //#inflow density > surface
                //# Initial plunge point (Eq. x in GLM manual)
                Inflow_Height_t0 = pow(Inflows[iRiver].FlowRate * Inflows[iRiver].Factor / SecsPerDay, 0.4)/
                                   pow(Inflow_gprime_t0, 0.2);
            else { //#inflow density < surface inflow thickness so plunge to depth of surface layer
                if (surfLayer > 0)
                    Inflow_Height_t0 = Lake[surfLayer].Height-Lake[surfLayer-1].Height;
                else
                    Inflow_Height_t0 = Lake[surfLayer].Height;
            }

            if (Inflow_Height_t0 > Lake[surfLayer].Height) //#Minimum limit = depth of the surface layer
                Inflow_Height_t0 = Lake[surfLayer].Height;

            //# Wave number squared (Eq. y in GLM manual)
            WaveNumSquared = cwnsq2/sqr(Inflow_Height_t0);
            //# Velocity of the plunging inflow (Eq. y in GLM manual)
            vel = 0.1 * Inflows[iRiver].FlowRate*Inflows[iRiver].Factor / SecsPerDay /
                                          (sqr(Inflow_Height_t0)*sin(Alpha)/cos(Alpha));
        }
     }

    /**************************************************************************
     * Here the integer iCnt is an inflow count of the number of inflows      *
     * associated with iRiver for that day including any inflows from         *
     * previous days.                                                         *
     * iCnt is initialised at the beginning of the simulation as zero and     *
     * incremented  each day that the inflow > 0.                             *
     * Once inflow has been inserted, iCnt = iCnt - 1.                        *
     **************************************************************************/

    for (iRiver = 0; iRiver < NumInf; iRiver++) {
        if (Inflows[iRiver].FlowRate * Inflows[iRiver].Factor > zero && !Inflows[iRiver].SubmFlag) {
            Inflows[iRiver].QDown[Inflows[iRiver].iCnt] = Inflows[iRiver].FlowRate * Inflows[iRiver].Factor;
            Inflows[iRiver].TDown[Inflows[iRiver].iCnt] = Inflows[iRiver].TemInf;
            Inflows[iRiver].SDown[Inflows[iRiver].iCnt] = Inflows[iRiver].SalInf;
            Inflows[iRiver].DDown[Inflows[iRiver].iCnt] = MISVAL;

            for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
                Inflows[iRiver].WQDown[Inflows[iRiver].iCnt][wqidx] = Inflows[iRiver].WQInf[wqidx];

            Inflows[iRiver].iCnt++;
        }
    }

    einff = zero; //# At the start of the day initialise the deltaPE to zero
    iRiver = 0;   //# Start with first inflow
    while (1) {
       /**************************************************************************
        * Work through each element in the downflow stacks and calculate the     *
        * travel distance and entrainment for the present day, and whether or    *
        * not it reaches its level of neutral buoyancy and hence can be inserted.*
        **************************************************************************/
        while (iRiver < NumInf) {
            if  (!Inflows[iRiver].SubmFlag) {
               Alpha = Inflows[iRiver].Alpha;            // stream cross-section half-angle
               Phi = Inflows[iRiver].Phi;                // river thalweg slope
               DragCoeff = Inflows[iRiver].DragCoeff;    // bottom roughness of thalweg

               //# Bulk Richardson's number of the inflow (Eq. x in GLM manual)
               Ri = DragCoeff * ( 1.0 + 0.21 * sqrt(DragCoeff) * sin(Alpha) ) / (sin(Alpha) * sin(Phi) / cos(Phi));
               //# ...or... Imberger and Patterson 1981 Eq 56
               //Ri = (DragCoeff / (sin(Alpha) * sin(Phi) / cos(Phi))) * 1/( 1.0 - 0.85 * sqrt(DragCoeff) * sin(Alpha) ) ;

               //# Entrainment coefficient
               EntrainmentCoeff = 1.6 * pow(DragCoeff, 1.5) / Ri;

               Inflows[iRiver].NoIns = 0;  //# Initialise number of insertions as zero

               k = -1;
               //# Loop for calculation of each separate inflow element.
               do  {
                   k++;
               } while (insert_inflow(k, iRiver, Lake[surfLayer].Height,
                                           Alpha, Phi, EntrainmentCoeff, Ri));
            }
            iRiver++;
        }

        /**********************************************************************
         * Insert all of the parcels which reached their level of NB on this  *
         * day.  Adjust the stacking to note the removal.                     *
         **********************************************************************/

        Inflow_width = 0.0;
        for (iRiver = 0; iRiver < NumInf; iRiver++) {
            if (Inflows[iRiver].SubmFlag) {
                //# Submerged inflows have no momentum, simply insert at specified level and
                //# adjust layer properties including water quality variables
                //# Get layer number at submerged inflow level
                for (Layer_subm = botmLayer; Layer_subm <= surfLayer; Layer_subm++)
                    if (Lake[Layer_subm].Height >= Inflows[iRiver].SubmHeight) break;

                Lake[Layer_subm].Temp = combine(Lake[Layer_subm].Temp, Lake[Layer_subm].LayerVol, Lake[Layer_subm].Density,
                                               Inflows[iRiver].TemInf, (Inflows[iRiver].FlowRate*Inflows[iRiver].Factor),
                                               calculate_density(Inflows[iRiver].TemInf, Inflows[iRiver].SalInf));
                Lake[Layer_subm].Salinity = combine(Lake[Layer_subm].Salinity, Lake[Layer_subm].LayerVol, Lake[Layer_subm].Density,
                                               Inflows[iRiver].SalInf, (Inflows[iRiver].FlowRate*Inflows[iRiver].Factor),
                                               calculate_density(Inflows[iRiver].TemInf, Inflows[iRiver].SalInf));

                for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
                    _WQ_Vars(wqidx, Layer_subm) = combine_vol(_WQ_Vars(wqidx, Layer_subm), Lake[Layer_subm].LayerVol,
                                               Inflows[iRiver].WQInf[wqidx], (Inflows[iRiver].FlowRate*Inflows[iRiver].Factor));

                Lake[Layer_subm].Density = calculate_density(Lake[Layer_subm].Temp, Lake[Layer_subm].Salinity);
                Lake[Layer_subm].LayerVol = Lake[Layer_subm].LayerVol+(Inflows[iRiver].FlowRate*Inflows[iRiver].Factor);

                if ( ptm_sw ) {
                    // insert particles ---
                    upper_height = Lake[Layer_subm].Height;
                    lower_height = 0.0; if (Layer_subm>botmLayer) lower_height = Lake[Layer_subm-1].Height;
                    double_particles = floor(Inflows[iRiver].ParticleConc * (Inflows[iRiver].FlowRate*Inflows[iRiver].Factor)); //@MEL implement this
                    new_particles = (int) floor(double_particles);
                    ptm_addparticles(new_particles, max_particle_num, upper_height, lower_height);
                    // insert particles ---
                }

                Lake[botmLayer].Vol1 = Lake[botmLayer].LayerVol;
                if (surfLayer != botmLayer) {
                    for (j = (botmLayer+1); j <= surfLayer; j++)
                        Lake[j].Vol1 = Lake[j-1].Vol1 + Lake[j].LayerVol;
                }
            } else {
                for (j = 0; j < Inflows[iRiver].NoIns; j++) {
                    for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
                        WQ_VarsTmp[iRiver][j][wqidx] = Inflows[iRiver].WQIns[j][wqidx];

                    insert(Inflows[iRiver].QIns[j], Inflows[iRiver].DIIns[j], Inflows[iRiver].Phi,
                              Inflows[iRiver].TIns[j], Inflows[iRiver].SIns[j],
                                         WQ_VarsTmp[iRiver][j], SecsPerDay, &Inflow_width, &ll);
                                         // QIns is flow; D depth, WQ vars defined above; need to pass particle
                                         // fraction to insert routine

                    for (jk = Inflows[iRiver].InPar[j]-1; jk < Inflows[iRiver].iCnt-1; jk++) {
                        Inflows[iRiver].QDown[jk] = Inflows[iRiver].QDown[jk+1];
                        Inflows[iRiver].TDown[jk] = Inflows[iRiver].TDown[jk+1];
                        Inflows[iRiver].SDown[jk] = Inflows[iRiver].SDown[jk+1];

                        for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
                            Inflows[iRiver].WQDown[jk][wqidx] = Inflows[iRiver].WQDown[jk+1][wqidx];
                            // jk is the inflow parcel; wqidx is wq attribute; would need NDown for particles

                        Inflows[iRiver].DDown[jk] = Inflows[iRiver].DDown[jk+1];
                        Inflows[iRiver].DOld[jk] = Inflows[iRiver].DOld[jk+1];
                    }
                    Inflows[iRiver].TotIn -= Inflows[iRiver].QIns[j];
                    Inflows[iRiver].iCnt--;
                    if (Inflows[iRiver].iCnt == 0) {
                        Inflows[iRiver].TotIn = zero;
                        Inflows[iRiver].Dlwst = MISVAL;
                    }
                    for (k = j; k < Inflows[iRiver].NoIns; k++)
                        Inflows[iRiver].InPar[k]--;
                }
            }
        }

        for (iRiver = 0; iRiver < NumInf; iRiver++) {
            //# Reset the number of insertions per river to be zero.
            Inflows[iRiver].NoIns = 0;

            //# Calculate the front of the downflow for each river.
            Inflows[iRiver].Dlwst = MISVAL;
            for (j = 0; j < Inflows[iRiver].iCnt; j++) {
                if (Inflows[iRiver].DDown[j] < Inflows[iRiver].Dlwst)
                    Inflows[iRiver].Dlwst = Inflows[iRiver].DDown[j];
            }
        }

        iRiver = -1;
        /**********************************************************************
         * If a flow fit error has occured, go back and reroute first inflow  *
         * for the river of concern - iRiver.                                   *
         **********************************************************************/
        if ( new_storage(&iRiver) ) break;
    }

    //# Make adjustments to update the layer heights, based on these vol changes
    resize_internals(2, botmLayer);
    check_layer_thickness();

    //# Update particle vertical position due to inflow insertion
    if ( ptm_sw )
        ptm_layershift(0.0, Lake[surfLayer].Height - height_start);  // !!!! ASSUMING SHIFT IS ALL LAYERS

    return Lake[surfLayer].Vol1 - VolSum;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * Solve the cubic for flowing depth in a triangular river valley.            *
 ******************************************************************************/
static int flowing_depth(AED_REAL *HF, AED_REAL *DL, AED_REAL Volume,
                                   AED_REAL DEPTH, AED_REAL Phi, AED_REAL Alpha)
{
    AED_REAL Theta, H, TAG, CON;
    int all_ok = TRUE;

/*----------------------------------------------------------------------------*/

    //# Set the coefficients of the cubic.
    H = zero;
    if (Volume > zero) {
        TAG = tan(Alpha) / tan(Phi);

        while(1) {
            CON = DEPTH - *DL;
            if (fabs(1.0 - 6.0 * Volume  / TAG / pow(CON, 3)) <= 1.0) {
                Theta = acos(1.0 - 6.0 * Volume / TAG / pow(CON, 3));
                H = (2.0 * cos(Theta / 3.0) + 1.0) * CON / 2.0;
                if (H > zero && H < CON) break;
                H = (2.0 * cos(Theta / 3.0 + 2.0 * Pi / 3.0) + 1.0) * CON / 2.0;
                if (H > zero && H < CON) break;
                H = (2.0 * cos(Theta / 3.0 + 4.0 * Pi / 3.0) + 1.0) * CON / 2.0;
                if (H > zero && H < CON) break;
            }
            *DL -= 0.25;

            if (*DL <= zero) {
                all_ok = FALSE;
                break;
            }
        }
    }

    //# Then the flowing depth for this river is H.
    *HF = H;
    return all_ok;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * Calculate the new temporary storage table for use in resize_internals      *
 ******************************************************************************/
static int new_storage(int *iRiv)
{
   int iTop,j,k;

   AED_REAL D;
   AED_REAL DLOW;
   AED_REAL DOL;
   AED_REAL EXTRA;
   AED_REAL LAYPRO;
   AED_REAL LayerVol;
   AED_REAL VlSum;
   AED_REAL VOLSUM;
   AED_REAL MphLevelVoltemp;

/*----------------------------------------------------------------------------*/

    //# Must begin by calculating the new depth of the entire reservoir.
    VlSum = zero;
    for (k = 0; k < NumInf; k++)
        VlSum += Inflows[k].TotIn;

    VOLSUM = VlSum + Lake[surfLayer].Vol1;
    j = 0;
    while (j < Nmorph) {
        if (VOLSUM <= MphLevelVol[j]) {
            j--;
            break;
        }
        j++;
    }
    if (j >= Nmorph) j = Nmorph - 1;
    if (j < 0) j = 0;
    Lake[surfLayer].Height = ((j+1) + (VOLSUM - MphLevelVol[j])/dMphLevelVol[j]) * 0.1;

    //# Calculate the flowing depths of the inflowing streams.
    //#
    DLOW = MISVAL;
    for (k = 0; k < NumInf; k++) {
        if (!flowing_depth(&Inflows[k].HFlow, &Inflows[k].Dlwst, Inflows[k].TotIn,
                         Lake[surfLayer].Height, Inflows[k].Phi, Inflows[k].Alpha) ) {
             //# If flow error then go back to INFLOW
             *iRiv = k;
             return FALSE;
        }
        if (DLOW > Inflows[k].Dlwst) DLOW = Inflows[k].Dlwst;
        if (Inflows[k].HFlow > Lake[surfLayer].Height) {
            fprintf(stderr, "Error in height of flow\n");
            exit(1);
        }
    }

    //# Account for the extra volumes of the downflows.
    iTop = round(Lake[surfLayer].Height * MphInc) - 1;
    if (iTop < botmLayer) iTop = botmLayer;
    if ( iTop >= Nmorph ) iTop = Nmorph - 1;
    EXTRA = zero;
    DOL = zero;

    for (j = botmLayer; j <= iTop; j++) {
        D = (j+1) / MphInc;
        EXTRA = extra_volume(DOL, D, Lake[surfLayer].Height, Inflows);
        if (j == botmLayer)
            LayerVol = MphLevelVol[j];
        else
            LayerVol = MphLevelVol[j] - MphLevelVol[j-1];

        LAYPRO = (D-DLOW) / D;
        if (LAYPRO < zero) LAYPRO = zero;
        MphLevelVoltemp = zero;
        if (j != botmLayer) MphLevelVoltemp = MphLevelVoldash[j-1];
        MphLevelVoldash[j] = MphLevelVol[j] - EXTRA;
        if (MphLevelVoldash[j] < (MphLevelVoltemp + (1.0 - LAYPRO) * LayerVol))
            MphLevelVoldash[j] =  MphLevelVoltemp + (1.0 - LAYPRO) * LayerVol;
    }

    for (j = iTop+1; j < Nmorph; j++)
        MphLevelVoldash[j] = MphLevelVol[j] - VlSum;

    for (j = 0; j < Nmorph-1; j++)
        dMphLevelVolda[j] = MphLevelVoldash[j+1] - MphLevelVoldash[j];

    dMphLevelVolda[Nmorph-1] = dMphLevelVolda[Nmorph-2];

    return TRUE;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * Compute the volume in inflow stacks which lies between depths d1 and d2.   *
 ******************************************************************************/
static AED_REAL extra_volume(AED_REAL d1, AED_REAL d2,
                                  AED_REAL SurfaceHeight, InflowDataType *Inflow)
{
    AED_REAL DTOP, DBOT, SUM, TAG, ZB;
    int k;

    SUM = zero;
    for (k = 0; k < NumInf; k++) {
        ZB = Inflow[k].Dlwst;
        if (! ((d2 <= ZB) || (SurfaceHeight <= Inflow[k].Dlwst)) ) {
            TAG = tan(Inflow[k].Alpha)/tan(Inflow[k].Phi);
            DBOT = d1;
            if (d1 < ZB) DBOT = ZB;
            DTOP = d2;
            if (d2 > SurfaceHeight) DTOP=SurfaceHeight;

            if (DTOP < Inflow[k].Dlwst+Inflow[k].HFlow)
                SUM += TAG*(pow(DTOP-ZB,3)-pow(DBOT-ZB,3))/3.0;
            else if (DBOT > Inflow[k].Dlwst + Inflow[k].HFlow)
                SUM += TAG*pow(Inflow[k].HFlow,2)*(DTOP-DBOT);
            else
                SUM += TAG*(pow(Inflow[k].HFlow,3)-pow(DBOT-ZB,3))/3.0 +
                                TAG*pow(Inflow[k].HFlow,2) *
                                   (DTOP-Inflow[k].Dlwst-Inflow[k].HFlow);
        }
    }
    return SUM;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 * Type 3 is fixed outlet heights + check for crit. hypol. oxygen             *
 * Type 4 is variable outlet heights for ISOTHERM                             *
 * Type 5 is variable outlet heights for TEMPERATURE TIME SERIES              *
 *                                                                            *
 * These conditions were provided by :                                        *
 *                                                                            *
 *     Dipl.-Hydrol. Michael Weber / PhD Student                              *
 *                                                                            *
 *     Helmholtz-Zentrum für Umweltforschung (UFZ)                            *
 *     Helmholtz Centre for Environmental Research (UFZ)                      *
 *                                                                            *
 *     Department Seenforschung                                               *
 *     Department of Lake Research                                            *
 *                                                                            *
 ******************************************************************************/

int find_layer_by_height(AED_REAL height)
{
    int i;
    for (i = 0; i < NumLayers; i++) {
        if (Lake[i].Height >=  height)
             return i;
    }
    return -1;
}


/******************************************************************************/
static AED_REAL NewDrawHeight(int jday, int i, int idx_dep, AED_REAL height,
                    AED_REAL lWithdrawalTemp, AED_REAL maxtemp, AED_REAL mintemp,
                    int within_temp_range, int within_facility_range,
                    int upper_bound, int lower_bound, FILE *sc_file)
{
    if (height > Lake[surfLayer].Height)
        height = Lake[surfLayer].Height;
    if (idx_dep < 0) {
        fprintf(sc_file,"%8d,%8d,%8d,%12.4lf,%12.4lf,%12.4lf,%8d,%8d,%8d,%8d,%12.4lf,%8d,%8d,%8d,%8d,%8d\n",
                         jday,Outflows[i].Type,0,lWithdrawalTemp,maxtemp,mintemp,
                         within_temp_range,within_facility_range,upper_bound,lower_bound,height,
                         999,999,999,999,999);
    } else {
        fprintf(sc_file,"%8d,%8d,%8d,%12.4lf,%12.4lf,%12.4lf,%8d,%8d,%8d,%8d,%12.4lf,%12.4lf,%8d,%8d,%8d,%8d\n",
                         jday,Outflows[i].Type,0,lWithdrawalTemp,maxtemp,mintemp,
                         within_temp_range,within_facility_range,upper_bound,lower_bound,height,
                         _WQ_Vars(Outflows[i].crit_idx,idx_dep),999,999,999,999);
    }
    return height;
}

/******************************************************************************/
static AED_REAL NewDrawHeightTempdep(int jday, int i, int idx_dep, AED_REAL height,
                    AED_REAL lWithdrawalTemp, AED_REAL maxtemp, AED_REAL mintemp,
                    int within_temp_range, int within_facility_range,
                    int upper_bound, int lower_bound, FILE *sc_file)
{
    AED_REAL MASSdepvarwith = 0, MASSbotout = 0, TEMPbotout = 0;
    AED_REAL TEMPdepvarwith_real = -1, Tmix = -1;

    int idx_lay;

    if (height > Lake[surfLayer].Height)
    height = Lake[surfLayer].Height;

    idx_lay = find_layer_by_height(height);

    TEMPdepvarwith_real = Lake[idx_lay].Temp;
    Tmix = ((TEMPdepvarwith_real*MASSdepvarwith)+(TEMPbotout*MASSbotout))/(MASSdepvarwith+MASSbotout);
    if (idx_dep < 0) {
        fprintf(sc_file,"%8d,%8d,%8d,%12.4lf,%12.4lf,%12.4lf,%8d,%8d,%8d,%8d,%12.4lf,%8d,%8d,%12.4lf,%12.4lf,%12.4lf\n",
                            jday,Outflows[i].Type,0,lWithdrawalTemp,maxtemp,mintemp,
                            within_temp_range,within_facility_range,upper_bound,lower_bound,height,
                            999,1,(Outflows[i].Draw * Outflows[i].Factor),
                            (Outflows[i+1].Draw * Outflows[i+1].Factor),Tmix);
    } else {
        fprintf(sc_file,"%8d,%8d,%8d,%12.4lf,%12.4lf,%12.4lf,%8d,%8d,%8d,%8d,%12.4lf,%12.4lf,%8d,%12.4lf,%12.4lf,%12.4lf\n",
                         jday,Outflows[i].Type,0,lWithdrawalTemp,maxtemp,mintemp,
                         within_temp_range,within_facility_range,upper_bound,lower_bound,height,
                         _WQ_Vars(Outflows[i].crit_idx,idx_dep),1,(Outflows[i].Draw * Outflows[i].Factor),
                           (Outflows[i+1].Draw * Outflows[i+1].Factor),Tmix);
    }
    return height;
}

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
AED_REAL SpecialConditionDraw(int jday, int i)
{
    int j, doit;
    AED_REAL DrawHeight = 0.;
    int idx_dep = -1, idx_lay = -1, l_crit = 1;

    if ( sc_files == NULL ) {
        int j;
        sc_files = (FILE**)malloc(NumOut * sizeof(FILE*));
        for (j = 0; j < NumOut; j++) sc_files[j] = NULL;
    }

    /**********************************************************************
     * Type 3 is fixed outlet heights + check for critical variable       *
     **********************************************************************/
    if (Outflows[i].Type == 3) {
        if (sc_files[i] == NULL) {
            sc_files[i] = fopen("outlet_values_type_3.txt","w");
            fprintf(sc_files[i],"JDay,OutletType,CritFlag,DrawHeight,ActVal\n");
        }

        idx_dep = find_layer_by_height(crit_dep-Base);

        if (jday < checkjday) {
            DrawHeight = Outflows[i].Hcrit-Base;       //# withdraw at another fixed height (e.g. bottom outlet)
            if (DrawHeight > Lake[surfLayer].Height)   //get TargetLayerTemp
                DrawHeight = Lake[surfLayer].Height;
        } else {
            AED_REAL ActVal = _WQ_Vars(Outflows[i].crit_idx,idx_dep);
            checkjday = -1;
            if ( CRITabove ? (ActVal >= crit_val) : (ActVal <= crit_val) ) {
                if (checkjday < 0)
                    checkjday = jday+crit_days;
                DrawHeight = Outflows[i].Hcrit-Base;              //# withdraw at another fixed height (e.g. bottom outlet)
                if (DrawHeight > Lake[surfLayer].Height)
                    DrawHeight = Lake[surfLayer].Height;
                else if (Outflows[i].Hcrit < Base) {
                    DrawHeight = 1; //1 m above bottom
                    fprintf(stderr,"outlet_crit %12.4lf < base_elev %12.4lf - set to %12.4lf\n",Outflows[i].Hcrit,Base,Base+1);
                }
            } else {
                DrawHeight = Outflows[i].OLev;             //# Fixed height offtake
                if (DrawHeight > Lake[surfLayer].Height)   //get TargetLayerTemp
                    DrawHeight = Lake[surfLayer].Height;
                l_crit = 0;
            }
        }
        fprintf(sc_files[i],"%8d,%8d,%8d,%12.4lf,%12.4lf\n",
                   jday,Outflows[i].Type,l_crit,DrawHeight,_WQ_Vars(Outflows[i].crit_idx,idx_dep));

    /**********************************************************************
     * Type 4 is variable outlet heights for ISOTHERM                     *
     * Type 5 is variable outlet heights for TEMPERATURE TIME SERIES      *
     **********************************************************************/
    } else if (Outflows[i].Type == 4 || Outflows[i].Type == 5) {
        int targetlyr = 0;
        AED_REAL mintemp = 100., maxtemp = -100., mindifftemp;
        AED_REAL lWithdrawalTemp = WithdrawalTemp;
        AED_REAL laketemp;

        if (Outflows[i].Type == 4)
            lWithdrawalTemp = Outflows[i].TARGETtemp;


        if (sc_files[i] == NULL) {
            if (Outflows[i].Type == 4)
                sc_files[i] = fopen("outlet_values_type_4.txt","w");
            else
                sc_files[i] = fopen("outlet_values_type_5.txt","w");
            fprintf(sc_files[i],"JDay,OutletType,CritOXY,TargetTemp,LakeMAXTemp,LakeMINTemp");
            fprintf(sc_files[i],",within_temp_range,within_facility_range,upper_bound,lower_bound");
            fprintf(sc_files[i],",DrawHeight,ActOXY,mix_withdraw,DISdepvarwith,DISbotout,Tmix\n");
        }

        doit = FALSE;
        if (COUPLoxy) {
            // oxygen coupling — always checks OXY_oxy regardless of crit_varname
            static int oxy_idx = -1;
            if (oxy_idx < 0) {
                size_t tl = strlen("OXY_oxy");
                oxy_idx = wq_var_index_c("OXY_oxy", &tl);
            }
            idx_dep = find_layer_by_height(crit_dep-Base);

            if (jday < checkjday) {
                DrawHeight = Outflows[i].Hcrit-Base;         //# withdraw at another fixed height (e.g. bottom outlet)
                if (DrawHeight > Lake[surfLayer].Height)   //get TargetLayerTemp
                    DrawHeight = Lake[surfLayer].Height;
                fprintf(sc_files[i],"%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d,%12.4lf,%12.4lf,%8d,%8d,%8d,%8d\n",
                              jday,Outflows[i].Type,1,999,999,999,999,999,999,999,DrawHeight,
                              (oxy_idx >= 0 ? _WQ_Vars(oxy_idx,idx_dep) : -1.),999,999,999,999);
            } else {
                AED_REAL ActOxy = (oxy_idx >= 0) ? _WQ_Vars(oxy_idx,idx_dep) : 0.0;
                checkjday = -1;
                if (ActOxy <= crit_val) {  // oxygen coupling always checks below threshold
                    if (checkjday < 0)
                        checkjday = jday+crit_days;
                    DrawHeight = Outflows[i].Hcrit-Base;               //# withdraw at another fixed height (e.g. bottom outlet)
                    if (DrawHeight > Lake[surfLayer].Height)
                        DrawHeight = Lake[surfLayer].Height;
                    else if (Outflows[i].Hcrit < Base) {
                        DrawHeight = 1; //1 m above bottom
                        fprintf(stderr,"outlet_crit %12.4lf < base_elev %12.4lf - set to %12.4lf\n",Outflows[i].Hcrit,Base,Base+1);
                    }
                    fprintf(sc_files[i],"%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d,%12.4lf,%12.4lf,%8d,%8d,%8d,%8d\n",
                                  jday,Outflows[i].Type,1,999,999,999,999,999,999,999,DrawHeight,
                                  ActOxy,999,999,999,999);
                } else
                    doit = TRUE;
            }
        } else {
            // without oxygen coupling
            //fprintf(stderr,"De-coupled - Num_WQ_Vars %8d\n",Num_WQ_Vars);
            idx_dep = -1;
            doit = TRUE;
        }

        if ( doit ) {
            for (j = 0; j < NumLayers; j++) {
                if (maxtemp < Lake[j].Temp) maxtemp = Lake[j].Temp;
                if (mintemp > Lake[j].Temp) mintemp = Lake[j].Temp;
            }

            if (Lake[surfLayer].Height >= (fac_range_lower-Base)) {
                if (MIXwithdraw && i != NumOut && (Outflows[i+1].Draw * Outflows[i+1].Factor) > 0) {
                    AED_REAL MASSdepvarwith = 0, MASSbotout = 0, TEMPbotout = 0;
                    AED_REAL TEMPdepvarwith = 0, TARGETtemp = -1.;

                    if (Outflows[i].Type == 5 && maxtemp < MINlaketemp)
                        DrawHeight = NewDrawHeightTempdep(jday,i,idx_dep,(fac_range_upper-Base),lWithdrawalTemp, maxtemp, mintemp,0,999,1,0, sc_files[i]);
                    else {
                        MASSdepvarwith = (Outflows[i].Draw * Outflows[i].Factor); //mass of water for depth-variable withdrawal
                        MASSbotout = (Outflows[i+1].Draw * Outflows[i+1].Factor); //mass of water for bottom outlet

                        idx_lay = find_layer_by_height(Outflows[i+1].OLev);

                        TEMPbotout = Lake[idx_lay].Temp; //temperature of bottom outlet layer
                        //# calculate temperature for depth-variable withdrawal from mixing temperature
                        TEMPdepvarwith = (lWithdrawalTemp*(MASSdepvarwith+MASSbotout)-(MASSbotout*TEMPbotout)) / MASSdepvarwith;
                        TARGETtemp = TEMPdepvarwith;
                        if (TARGETtemp > maxtemp) // check if lake temperature is lower than target temperature
                            DrawHeight = NewDrawHeightTempdep(jday,i,idx_dep,(fac_range_upper-Base),lWithdrawalTemp, maxtemp, mintemp,0,999,1,0, sc_files[i]);
                        else if (TARGETtemp < mintemp) // check if lake temperature is higher than target temperature
                            DrawHeight = NewDrawHeightTempdep(jday,i,idx_dep,(fac_range_lower-Base),lWithdrawalTemp, maxtemp, mintemp,0,999,0,1, sc_files[i]);
                        else {
                            mindifftemp = 100;
                            for (j = 0; j < NumLayers; j++) {
                                laketemp = fabs(Lake[j].Temp-TARGETtemp);
                                if (mindifftemp > laketemp)
                                    mindifftemp = laketemp;
                                if (laketemp == mindifftemp)
                                    targetlyr = j;
                            }
                            if (Lake[targetlyr].Height >= (fac_range_upper-Base))
                                DrawHeight = NewDrawHeightTempdep(jday,i,idx_dep,(fac_range_upper-Base),lWithdrawalTemp, maxtemp, mintemp,1,0,1,0, sc_files[i]);
                            else if (Lake[targetlyr].Height <= (fac_range_lower-Base))
                                DrawHeight = NewDrawHeightTempdep(jday,i,idx_dep,(fac_range_lower-Base),lWithdrawalTemp, maxtemp, mintemp,1,0,0,1, sc_files[i]);
                            else {
                                DrawHeight = Lake[targetlyr].MeanHeight; // combine layer index with mean height
                                if (DrawHeight <= (fac_range_lower-Base))
                                    DrawHeight = NewDrawHeightTempdep(jday,i,idx_dep,(fac_range_lower-Base),lWithdrawalTemp, maxtemp, mintemp,1,0,0,1, sc_files[i]);
                                else
                                    NewDrawHeightTempdep(jday,i,idx_dep,DrawHeight,lWithdrawalTemp, maxtemp, mintemp,1,1,0,0, sc_files[i]);
                            }
                        }
                    }
                } else {
                    if (Outflows[i].Type == 5 && maxtemp < MINlaketemp)
                        DrawHeight = NewDrawHeightTempdep(jday,i,idx_dep,(fac_range_upper-Base),lWithdrawalTemp, maxtemp, mintemp,0,999,1,0, sc_files[i]);
                    else {
                        if (lWithdrawalTemp > maxtemp) // check if lake temperature is lower than target temperature
                            DrawHeight = NewDrawHeight(jday,i,idx_dep,(fac_range_upper-Base),lWithdrawalTemp, maxtemp, mintemp,0,999,1,0, sc_files[i]);
                        else if (lWithdrawalTemp < mintemp) // check if lake temperature is higher than target temperature
                            DrawHeight = NewDrawHeight(jday,i,idx_dep,(fac_range_lower-Base),lWithdrawalTemp, maxtemp, mintemp,0,999,0,1, sc_files[i]);
                        else {
                            mindifftemp = 100;
                            for (j = 0; j < NumLayers; j++) {
                                laketemp = fabs(Lake[j].Temp-lWithdrawalTemp);
                                if (mindifftemp > laketemp)
                                    mindifftemp = laketemp;
                                if (laketemp == mindifftemp)
                                    targetlyr = j;
                            }
                            if (Lake[targetlyr].Height >= (fac_range_upper-Base))
                                DrawHeight = NewDrawHeight(jday,i,idx_dep,(fac_range_upper-Base),lWithdrawalTemp, maxtemp, mintemp,1,0,1,0, sc_files[i]);
                            else if (Lake[targetlyr].Height <= (fac_range_lower-Base))
                                DrawHeight = NewDrawHeight(jday,i,idx_dep,(fac_range_lower-Base),lWithdrawalTemp, maxtemp, mintemp,1,0,0,1, sc_files[i]);
                            else {
                                DrawHeight = Lake[targetlyr].MeanHeight; // combine layer index with mean height
                                if (DrawHeight <= (fac_range_lower-Base))
                                    DrawHeight = NewDrawHeight(jday,i,idx_dep,(fac_range_lower-Base),lWithdrawalTemp, maxtemp, mintemp,1,0,0,1, sc_files[i]);
                                else
                                    NewDrawHeight(jday,i,idx_dep,DrawHeight,lWithdrawalTemp, maxtemp, mintemp,1,1,0,0, sc_files[i]);
                            }
                        }
                    }
                }
            } else {
                DrawHeight = Outflows[i].Hcrit-Base;      //# withdraw at another fixed height (e.g. bottom outlet)
                if (DrawHeight > Lake[surfLayer].Height)  // get TargetLayerTemp
                    DrawHeight = Lake[surfLayer].Height;
            }
        }
    }
    return DrawHeight;
}



/******************************************************************************
 * Function to calculate the proportion of fluid withdrawn from any layer,    *
 * given the depth of its top and bottom, using a curve which fits the region *
 * of water drawn in a given time to decide which set of withdrawal curves to *
 * use. If large withdrawal use first set, otherwise the 2nd.                 *
 ******************************************************************************/
AED_REAL delta_volume(AED_REAL z1, AED_REAL z2, AED_REAL da, AED_REAL avdel,
                      AED_REAL hh, AED_REAL DeltaTop, AED_REAL DeltaBot)
{
    AED_REAL a, da4, da7, s1, s2, s3,
             tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, z3;
    AED_REAL dV = 0.0;

/*----------------------------------------------------------------------------*/
    if (da >= 0.9*avdel) {
        // Curves for large withdrawal
        s1 = (z1 - hh) / avdel;
        s2 = (z2 - hh) / avdel;
        a = da / avdel;

        // If top and bottom of layer fall within lower curve
        if (z1 <= 0.25*a*avdel+hh) {
            tmp1 = s1 + 2.0 / a / (1.0 + 0.5*a*(DeltaBot/avdel + s1));
            tmp2 = s2 + 2.0 / a / (1.0 + 0.5*a*(DeltaBot/avdel + s2));
            dV  = tmp1 - tmp2;
        } else if (z2 >= hh+0.25*avdel*a) {
            //  If layer falls within upper curve
            tmp1 = exp(-1 * sqrt(2.0) * (DeltaTop/avdel + 0.75*a));
            tmp2 = 1.0 + 0.5 * a * DeltaBot / avdel + sqr(a)/8.0;
            tmp3 = (1 - 1.0 / sqr(tmp2)) * sqr((1 + tmp1) / (1 - tmp1));
            tmp4 = sqrt(2.0) * (s1 - DeltaTop / avdel - a);
            tmp5 = sqrt(2.0) * (s2 - DeltaTop / avdel - a);
            tmp6 = 4.0 / (1 + exp(tmp4)) + tmp4;
            tmp7 = 4.0 / (1 + exp(tmp5)) + tmp5;
            dV = (tmp6 - tmp7) / sqrt(2.0) * tmp3;
        } else {
            //  If join of curves lies within the layer
            s3 = 0.25*a;
            tmp2 = s2 + 2.0 / a / (1.0 + 0.5*a*(DeltaBot/avdel + s2));
            tmp1 = s3 + 2.0 / a / (1.0 + 0.5*a*(DeltaBot/avdel + s3));
            dV = tmp1 - tmp2;
            tmp1 = exp(-1 * sqrt(2.0) * (DeltaTop / avdel + 0.75*a));
            tmp2 = 1.0 + 0.5 * a * DeltaBot / avdel + sqr(a)/8.0;
            tmp3 = (1 - 1.0 / sqr(tmp2)) * sqr((1+tmp1) / (1-tmp1));
            tmp4 = sqrt(2.0) * (s1 - DeltaTop/avdel - a);
            tmp5 = sqrt(2.0) * (s3 - DeltaTop/avdel - a);
            tmp6 = 4.0 / (1 + exp(tmp4)) + tmp4;
            tmp7 = 4.0 / (1 + exp(tmp5)) + tmp5;
            dV = dV + (tmp6 - tmp7) / sqrt(2.0) * tmp3;
        }
    } else {
        //  Curves for small withdrawal
        da4 = da / 4.0;
        da7 = da * 0.75;
        if (z1 <= da4 + hh) {
            //  If top and bottom fall within lower curve
            tmp1 = z1 + (DeltaBot + da4) / Pi * sin(Pi * (z1 - hh - da4) / (DeltaBot + da4));
            tmp2 = z2 + (DeltaBot + da4) / Pi * sin(Pi * (z2 - hh - da4) / (DeltaBot + da4));
            dV = tmp1 - tmp2;
        } else if (z2 >= hh+da4) {
            //  If layer falls within upper curve
            tmp1 = z1 + (DeltaTop + da7) / Pi * sin(Pi * (z1 - hh - da4) / (DeltaTop + da7));
            tmp2 = z2 + (DeltaTop + da7) / Pi * sin(Pi * (z2 - hh - da4) / (DeltaTop + da7));
            dV = tmp1 - tmp2;
        } else {
            // If join of curves falls within the layer
            z3 = 0.25 * da + hh;
            tmp1 = z3 + (DeltaBot + da4) / Pi * sin(Pi * (z3 - hh - da4) / (DeltaBot + da4));
            tmp2 = z2 + (DeltaBot + da4) / Pi * sin(Pi * (z2 - hh - da4) / (DeltaBot + da4));
            dV = tmp1 - tmp2;
            tmp3 = z3 + (DeltaTop + da7) / Pi * sin(Pi * (z3 - hh - da4) / (DeltaTop + da7));
            tmp4 = z1 + (DeltaTop + da7) / Pi * sin(Pi * (z1 - hh - da4) / (DeltaTop + da7));
            dV = tmp4 - tmp3 + dV;
        }
    }

    return dV;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
