/******************************************************************************
 *                                                                            *
 * glm_bubbler.c                                                              *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 2018 -  The University of Western Australia                      *
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
#include <math.h>

#include "glm.h"
#include "glm_globals.h"
#include "glm_util.h"
#include "glm_layers.h"
#include "glm_mixu.h"
#include "glm_deep.h"
#include "glm_input.h"
#include "glm_const.h"
#include "glm_bubbler.h"

#include "aed_time.h"

#define SQRT_Pi sqrt(Pi_)
#define SQRT_Two sqrt(2.0)

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

#ifdef _VISUAL_C_
#define ALLOCATE1D_I(dst, cnt) dst=malloc(sizeof(int) * (cnt))
#define ALLOCATE1D_R(dst, cnt) dst=malloc(sizeof(AED_REAL) * (cnt))
#define ALLOCATE2D(dst, cnt1, cnt2) dst=malloc(sizeof(AED_REAL) * (cnt1) * (cnt2))
#endif


/******************************************************************************/

LOGICAL     bubbler_on        = FALSE;
char       *bubbler_data_file = NULL;
AED_REAL    bubbler_aflow     = 0.0;
int         bubbler_nports    = 0;
AED_REAL    bubbler_bublen    = 0.0;
AED_REAL    bubbler_bdepth    = 0.0;
LOGICAL     bubbler_opopt     = FALSE;
AED_REAL    bubbler_ton       = 0.0;
AED_REAL    bubbler_toff      = 0.0;
LOGICAL     bubbler_intopt    = FALSE;
char       *bubbler_start     = NULL;
char       *bubbler_stop      = NULL;
LOGICAL     bubbler_eff       = FALSE;
char       *bubbler_eff_file  = NULL;



//# Bubbler efficiency variables
static AED_REAL tMean, tMeanS;
static AED_REAL tEffic1, tEffic2, tEffic3, tPn, tTDiff, tHTRise, tFlowT;
static AED_REAL aEffic1, aEffic2,          aPn, aTDiff, aHTRise, aFlowT;

static AED_REAL timeOn;
static int  bSecs, nTimeSteps;
static FILE *ueff = NULL, *uair = NULL;

static AED_REAL PeI, PeF, xMeI, xMeF, zCI, zCF;
static AED_REAL cee, WISO, effic1, effic2, effic3;
static AED_REAL xPerI, xPerF, xVolR;
static AED_REAL htRise;

static AED_REAL lm2, lambda, pWidth, alp2, ub2d, ub2;

static const AED_REAL one = 1.0, ten = 10.0;

/*############################################################################*/

static AED_REAL aFlow;
static AED_REAL bDepth;
static AED_REAL bubLen;
static AED_REAL tOn;
static AED_REAL tOff;
static AED_REAL cdifBeg;

static int nPorts;
static int iDayStart;
static int iDayFinish;
static int bLevel;

static int ibub;
static LOGICAL opopt, eff, intopt, bubOn, bubprint;

static char *fair = NULL;
static char *feff = NULL;

#define nStepY   10

/*============================================================================*/
static void pe3(AED_REAL *_xp, AED_REAL *_zc, AED_REAL *_xmas, AED_REAL BDE,
                                       AED_REAL *_xper, AED_REAL *_xvolr);

static void inBubble(LOGICAL intopt, int MaxNSy);

static void InActRK4(AED_REAL *Y, AED_REAL *DYDX, int N, AED_REAL H,
              AED_REAL *YOUT, AED_REAL B0, AED_REAL SQLAM, AED_REAL DRH,
              AED_REAL ALP, AED_REAL W0, AED_REAL FLOW, AED_REAL PA,
              AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0, int DIM);

static void InDepRK4(AED_REAL *Y, AED_REAL *DYDX, int N, AED_REAL H,
              AED_REAL *YOUT, AED_REAL B0, AED_REAL SQLAM, AED_REAL DRH,
              AED_REAL ALP, AED_REAL W0, AED_REAL FLOW, AED_REAL PA,
              AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0);

static void InActDrv(AED_REAL X, AED_REAL *Y, AED_REAL *DYDX, AED_REAL B0,
              AED_REAL SQLAM, AED_REAL DRH, AED_REAL ALP, AED_REAL W0,
              AED_REAL FLOW, AED_REAL PA, AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0, int DIM);

static void InDepDrv(AED_REAL X, AED_REAL *Y, AED_REAL *DYDX, AED_REAL B0,
              AED_REAL SQLAM, AED_REAL DRH, AED_REAL ALP, AED_REAL W0,
              AED_REAL FLOW, AED_REAL PA, AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0);


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void init_bubbler()
{
    int secs1, secs2;
    bSecs = 0;

    //# Variables for changes to the efficiency files
    nTimeSteps = 0;
    bubprint = FALSE;

    ibub = -1;

    //# Initialise bubbler details
    intopt = FALSE;
    timeOn = 0.0;
    aFlow = 0.0001;
    bDepth = 0.0;
    bubLen = 1.0;
    nPorts = 1;

    if (bubbler_on) {
        printf("Bubbler is ON\n");
        bubprint = TRUE;

        bubOn  = bubbler_on;
        fair   = bubbler_data_file;
        aFlow  = bubbler_aflow;
        nPorts = bubbler_nports;
        bubLen = bubbler_bublen;
        bDepth = bubbler_bdepth;
        opopt  = bubbler_opopt;
        tOn    = bubbler_ton;
        tOff   = bubbler_toff;
        intopt = bubbler_intopt;
        read_time_string(bubbler_start, &iDayStart, &secs1);
        read_time_string(bubbler_stop, &iDayFinish, &secs2);
        eff    = bubbler_eff;
        feff   = bubbler_eff_file;

        ibub = 0;

        //# Bubbler efficiency file setup
        if (eff) {
            ueff = fopen(feff, "w");

            fprintf(ueff, "JDAY, EFFECT, T DIFF, C, HT RISE, aFlow, EFFIC\n");
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void do_bubbler(int jday, int iclock)
{
    AED_REAL tDiff;
    int i;

    if ( ! bubOn ) return;

    //# Bubbler calculations
    if (ibub >= 0) {
        /**********************************************************************
         * Install bubble, get pe etc                                         *
         * Test for operation of bubbler                                      *
         * Test that JDAY falls within period of bubbler operation            *
         * Calculate the time step                                            *
         **********************************************************************/
        for (i = surfLayer; i > botmLayer; i--) {
            if (Lake[i].Height >= bDepth) bLevel = i;
        }

        if (bLevel == surfLayer)
            tDiff = 0.0;
        else
            tDiff = Lake[surfLayer].Temp - Lake[bLevel+1].Temp;

        /**********************************************************************
         * Bubbler optimisation                                               *
         * if tDiff GT Delta T (BUBBLER ON) { TURN BUBBLER ON                 *
         * if tDiff LT Delta T (BUBBLER OFF) { TURN BUBBLER OFF               *
         * } else { LEAVE AS FOR PREVIOUS DAY                                 *
         * Note: pe's from s/r pe3 are based on area calculations             *
         *       from TABLE.  Units of area (in TABLE) are 10**6 sq.m         *
         *       Therefore this multiplier needs to be used on                *
         *       everything that comes out of pe3                             *
         **********************************************************************/
        if (jday >= iDayStart && jday <= iDayFinish) {
            if (opopt) {
                if (! bubOn && tDiff > tOn) bubOn = TRUE;
                if (tDiff < tOff) bubOn = FALSE;
            } else {
                bubOn = TRUE;
            }
            if (aFlow <= 0.0) bubOn = FALSE;
        } else {
            bubOn = FALSE;
        }

        if (iclock == 0) {
            bSecs = 0;
            tEffic1 = 0.0;
            tEffic2 = 0.0;
            tEffic3 = 0.0;
            tPn = 0.0;
            tTDiff = 0.0;
            tHTRise = 0.0;
            tFlowT = 0.0;
        }

        if (bubOn) {
            effic1 = 0.;
            effic2 = 0.;
            effic3 = 0.;
            ibub = ibub+1;

            pe3(&PeI, &zCI, &xMeI, bDepth, &xPerI, &xVolR);

            if (ibub == 1)
                cdifBeg = PeI;

            inBubble(intopt, nStepY*MaxLayers + 1);   //# +1: sub-layers are 1-based, top index = (surfLayer+1)*nStepY

            check_layer_stability();
            check_layer_thickness();
            pe3(&PeF, &zCF, &xMeF, bDepth, &xPerF, &xVolR);

            //# Change in PE as we proceed
            //# Mechanical efficiency
            //# Isothermal compression work of compressor
            WISO = 2.303 * (1.0) * aFlow * noSecs *
                              log10((Lake[surfLayer].Height+10.2)/10.2) * 1.E5;
            effic2 = (xPerF-xPerI) * 9.81 * 1.E6 / WISO;

            //# change in pe based on initial stratification
            //# effic3 = (xPerF*1.e6-peibeg)/(pembeg-peibeg)
            //# try using difference between c.o.m and c.o.v
            effic3 = PeF/cdifBeg;

            bSecs = bSecs + noSecs;
            tEffic1 = tEffic1 + effic1 * (noSecs);
            tEffic2 = tEffic2 + effic2 * (noSecs);
            tEffic3 = tEffic3 + effic3 * (noSecs);
            tPn = tPn + cee * (noSecs);
            tTDiff = tTDiff + tDiff * (noSecs);
            tHTRise = tHTRise + htRise * (noSecs);
            tFlowT = tFlowT + aFlow * (noSecs);

            //# Calculate mean and variance of tDiff
            timeOn = timeOn + (noSecs) / 86400.0;
            nTimeSteps = nTimeSteps + 1;
            tMean = tMean + tDiff;
            tMeanS = tMeanS + tDiff * tDiff;
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void read_bubbler(int jday)
{
    //# Input from bubbler file if required
    if ( bubOn && jday >= iDayStart && jday <= iDayFinish) {
        read_bubble_data(jday, &aFlow, &nPorts, &bDepth, &bubLen);
        if (aFlow <= 0.0) bubOn = FALSE;
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void write_bubbler(int jday)
{
    if ( ! bubOn ) return;

    //# Bubbler efficiency file writing
    //# Write daily bubbler stuff if required
    if (eff && ueff!=NULL) {
        aEffic1 = tEffic1/(bSecs);
        aEffic2 = tEffic2/(bSecs);
        aPn = tPn/(bSecs);
        aTDiff = tTDiff/(bSecs);
        aHTRise = tHTRise/(bSecs);
        aFlowT = tFlowT;

        fprintf(ueff, "%d, %12.3f, %12.3f, %12.3f, %12.3f, %12.3f, %12.3f\n",
                          jday, aEffic1, aTDiff, aPn, aHTRise, aFlowT, aEffic2);
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void print_bubbler()
{
    if ( ! bubOn ) return;

    if (bubprint && eff) {
        //# Check whether an efficiency file exists, and is opened
        //# If so, then need to write out all the data read in
        if ( ueff != NULL ) {
            //# Check whether writing bubbler data from keyboard or file
            if (uair != NULL )
                fprintf(ueff,
                       "Bubbler Information\nInput File is \"%s\" Time On = %8.2f days\n",
                                  fair, timeOn);
            else
                fprintf(ueff, " Air Flow = %10.4f No. Ports = %5d Time On = %8.2f days\n",
                                  aFlow, nPorts, timeOn);
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void done_bubbler()
{
    if ( ! bubOn ) return;

    if (uair != NULL) fclose(uair);
    if (ueff != NULL) fclose(ueff);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *  This gets the PE stored above the level of BDEP (Includes the layer       *
 *  containing BDEP)                                                          *
 *                                                                            *
 * NOTE: All totals are relative to BDEP                                      *
 *                                                                            *
 * Variable Dictionary                                                        *
 *                                                                            *
 *     XP=(Centre of mass-Centre of vol)*g                                    *
 *     ZC=Centre of volume                                                    *
 *     XMAS=Total mass above bubbler                                          *
 *     XVOL=Total vol above bubbler                                           *
 *     XPE=Potential energy/g                                                 *
 *                                                                            *
 ******************************************************************************/
void pe3(AED_REAL *_xp, AED_REAL *_zc, AED_REAL *_xmas, AED_REAL BDE,
                                       AED_REAL *_xper, AED_REAL *_xvolr)
{
    AED_REAL  XP = *_xp, ZC = *_zc, XMAS = *_xmas;
    AED_REAL  XPER = *_xper, xVolR = *_xvolr;

    DOUBLETYPE X, Y, ZCVP;
    DOUBLETYPE AB, HB, DA, HT, DENS;
    DOUBLETYPE YB, XVOL, XVOLP, XMASSP;
    DOUBLETYPE XPE, ZCV, XMASS, BDEP;

    int i, ij, il, ix, L, IBDEP, ik;

/*----------------------------------------------------------------------------*/

    IBDEP = 0;
    DENS = 0.;

    BDEP = BDE;

    //# Get the layer containing BDEP
    for (i = surfLayer; i > botmLayer; i--)
        if (Lake[i].Height > BDEP) IBDEP = i;

    /**************************************************************************
     * Get centre of volume and potential energy, total mass above bdep etc   *
     * locate surface in relation to storage table                            *
     * ij is the last storage entry before the surface, y is the distance     *
     * to the surface from there                                              *
     **************************************************************************/
    X = Lake[surfLayer].Height * ten;
    Y = AMOD(X, one);
    ij = (X - Y);
    if (ij >= X) ij--;
    if (ij > Nmorph) {
        Y += (ij - Nmorph);
        ij = Nmorph;
    }
    ij--; // offset for 0 based index
    Y /= ten;

    //# Locate BDEP in relation to storage table
    //# ik is the last storage entry before BDEP, YB is the distance to
    //# BDEP from there; NOTE ik could = 0
    X = BDEP * 10.0;
    YB = AMOD(X, one);
    ik = (X - YB);
    if (ik > Nmorph) {
        YB += (ik - Nmorph);
        ik = Nmorph;
    }
    ik--; // offset for 0 based index
    YB /= 10.0;

    //# Initialise the loop, il is the current layer number,
    //# DA the area gradient
    //# starting at IBDEP
    //#
    ZCV = 0.0;
    XPE = 0.0;
    XMASS = 0.0;
    XVOL = 0.0;
    il = IBDEP;
    ix = 0;
    if (ik < 0) {
        DA = (MphLevelArea[0]) * 10.0;
        AB = DA * YB;
    } else {
        DA = (dMphLevelArea[ik]) * 10.0;
        AB = (MphLevelArea[ik]) + DA * YB;
    }
    HB = 0.0;

    //# Loop
    if (ik != Nmorph) {
        ik++;
        for (i = ik; i < ij; i++) {
            ix = 1;
            while (ix == 1) {
                DENS = (Lake[il].Density);
                HT = (i)/10.0 - BDEP;
                if (i != 0) DA = (dMphLevelArea[i-1])*10.0;
                ix = 0;
                if (HT>=((Lake[il].Height) - BDEP)) {
                    HT = (Lake[il].Height) - BDEP;
                    il = il+1;
                    ix = 1;
                }
                ZCVP = AB*(SQR(HT)-SQR(HB))/2.0+DA*(CUBE(HT)/3.0-HB*SQR(HT)/2.0+CUBE(HB)/6.0);
                XVOLP = (HT-HB)*(AB+DA*(HT-HB)/2.0);
                XMASSP = XVOLP*(1000.0+DENS);
                ZCV = ZCV+ZCVP;
                XPE = XPE+(1000.0+DENS)*ZCVP;
                XMASS = XMASS+XMASSP;
                XVOL = XVOL+XVOLP;
                AB = DA*(HT-HB)+AB;
                HB = HT;
            }
        }
    }

    //# This takes us to the last storage table entry
    //# Now go to the surface
    if (il == surfLayer) {
        HT = (Lake[surfLayer].Height) - BDEP;
        ZCVP = AB*(SQR(HT)-SQR(HB))/2.0+DA*(CUBE(HT)/3.0-HB*SQR(HT)/2.0+CUBE(HB)/6.0);
        XVOLP = (HT-HB)*(AB+DA*(HT-HB)/2.0);
        XMASSP = XVOLP*(1000.0+DENS);
        ZCV = ZCV+ZCVP;
        XPE = XPE+(1000.0+DENS)*ZCVP;
        XMASS = XMASS+XMASSP;
        XVOL = XVOL+XVOLP;
    } else {
        for (L = il; L <= surfLayer; L++) {
           DENS = (Lake[L].Density);
           HT = (Lake[L].Height) - BDEP;
           ZCVP = AB*(SQR(HT)-SQR(HB))/2.0+DA*(CUBE(HT)/3.0-HB*SQR(HT)/2.0+CUBE(HB)/6.0);
           XVOLP = (HT-HB)*(AB+DA*(HT-HB)/2.0);
           XMASSP = XVOLP*(1000.0+DENS);
           ZCV = ZCV+ZCVP;
           XPE = XPE+(1000.0+DENS)*ZCVP;
           XMASS = XMASS+XMASSP;
           XVOL = XVOL+XVOLP;
           AB = (Lake[L].LayerArea);
           HB = HT;
       }
    }

    //# Note that the units of XP are millions of kg M**2/sec**2,
    //# units of XMAS are millions of kg
    ZC = (ZCV/XVOL);
    XP = ((XPE/XMASS-ZCV/XVOL)*9.810);
    XMAS = (XMASS);
    XPER = (XPE);
    xVolR = (XVOL);

    // set returned values
    *_xp = XP; *_zc = ZC; *_xmas = XMAS;
    *_xper = XPER; *_xvolr = xVolR;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
static int inner_loop(int iLevel, int level, int *pCountI, int *retI, int NSY, int *iTop,
                                AED_REAL *yDepth, AED_REAL *yDepthM, AED_REAL *TH,
                                AED_REAL *W, AED_REAL *B, AED_REAL *ENT, AED_REAL *yDen)
{
    const AED_REAL SQLAM=0.09, US=0.3E0, HPRES=10.2E0;
    const AED_REAL R2=1000., ALP=0.083;
    const int threeI=3;

    int DIM0, DIM;

    AED_REAL GRADMO;
    AED_REAL DYDX[3], Y[3], YOUT[3];
    AED_REAL BIGH, BIGM, Z, DZ, SMALLZ;
    AED_REAL B0, W0, TH0;
    AED_REAL LPP;
    AED_REAL Q0, BU0, XM0, DQ, DRH, DW;
    AED_REAL B2, W2, TH2, XM2, XM1, BU01, BU02;
    AED_REAL DQDZ, DXMDZ, DBUDZ, ZZ, DW0DZ;
    AED_REAL PA, UB, HLAY;
    AED_REAL FLOW;
    AED_REAL RHO0;

    int i, pCount = *pCountI;

    *retI = level;

    LPP = bubLen/(2.0*(nPorts));
    RHO0 = Lake[iLevel].Density;
    PA = RHO0 * g * HPRES;
    FLOW = aFlow/(nPorts);
    UB = US * (1.E0 + SQLAM);

//  HLAY = Lake[surfLayer].Height - BD;

    /********************%%%%%%%%%%%%%%%%%%%%%%%%%%%%**********************/
    //# Start of outer loop
    //# Initilize parameters for plume
    do  {
        bDepth = yDepth[level];

        DIM = 1;

        //# Origin of bubbles is level; locate depth and set initial conditions
        //# Based on series solutions McDougall
        HLAY = yDepth[NSY] - bDepth;
        BIGH = HLAY + HPRES;
        //# McDougall source-strength parameter M (dimensionless);
        //# Patterson & Imberger 1989 eq (11) - slip velocity is CUBED
        BIGM = FLOW * PA * (1.E0 + SQLAM) /
                       (4.E0 * Pi * SQR(ALP) * RHO0 * SQR(BIGH) * CUBE(UB));

        //#  Change initial height to 1 layer depth
        if (level != NSY)
            Z = yDepth[level+1]-yDepth[level];
        else
            Z = yDepth[level]-yDepth[level-1];

        ZZ = Z/BIGH;
        B0 = 2.E0 * ALP * Z *
               (0.6E0 + 0.01719E0 * pow(ZZ/BIGM, 0.333E0) -
                  0.002527E0 * pow(ZZ/BIGM, 0.667E0) +
                  ZZ * (0.04609E0 - 0.000031E0/BIGM));  //# P&I89 eq (10): +x(0.04609 - 0.000031/M)
        W0 = UB * pow(BIGM/ZZ, 0.333E0) * (1.609E0 - 0.3195E0 * pow(ZZ/BIGM, 0.333E0) +
              0.06693E0 * pow(ZZ/BIGM, 0.667E0) + ZZ * (0.4536E0 - 0.0105E0 / BIGM));

        if (intopt)  //# NOTE TH0 no longer has g in it
            TH0 = FLOW * PA * (1.E0+SQLAM)/(SQLAM*(1.-ZZ)*BIGH*Pi*RHO0*(W0+UB)*B0*B0);
        else
            TH0 = FLOW * PA * (1.E0+SQLAM)/(SQLAM*(1.-ZZ)*BIGH*Pi*RHO0*(W0+UB)*B0*B0*g);

        //# Use 0 for initial gradient
        DW0DZ = 0.0;

        if (intopt) {
            //# Initial check for plume interaction
            if (B0>=LPP)
                DIM = 2;
            else
                DIM = 1;

            //#  **  Initial conditions for line source conditions
            //#       however presently it is assumed to be similar to point source
            //#       but if aerator is set up with single ports this should be ok
            //#       because it will not go through these calculations anyhow.
            //#       should be modified to be consistent.
            if (DIM==2) {
                BIGM = FLOW * PA * (1.0+lm2) / (4.0 * Pi * SQR(alp2) * RHO0 * SQR(BIGH) * CUBE(UB));
                B0 = 2.0 * alp2 * Z * (0.6E0 + 0.01719E0 * pow(ZZ/BIGM, 0.333E0) -
                       0.002527E0 * pow(ZZ / BIGM, 0.667E0) + ZZ * (0.04609E0 - 0.000031E0 / BIGM));
                W0 = ub2 * pow(BIGM/ZZ, 0.333E0) * (1.609E0 - 0.3195E0 *
                       pow(ZZ/BIGM, 0.333E0) + 0.06693E0 * pow(ZZ/BIGM, 0.667E0) +
                                  ZZ * (0.4536E0 - 0.0105E0 / BIGM));

                //# NOTE TH0 no longer has g in it
                TH0 = FLOW*PA*(1.0+lm2)/(lm2*(1.-ZZ)*BIGH*Pi* RHO0*(W0+UB)*B0*B0);

                Q0 = B0*W0*pWidth* sqrt(Pi);
                BU0 = RHO0*lambda*B0*TH0* sqrt(Pi)*W0*(pWidth)/sqrt(1.E0+lm2);
                XM0 = RHO0*W0*W0*sqrt(Pi)*B0*(pWidth)/sqrt(2.E0);
            } else {
                //# Initial fluxes
                Q0 = B0 * B0 * W0 * Pi;
                BU0 = RHO0 * SQLAM * B0 * B0 * TH0 * Pi * W0 / (1.E0 + SQLAM);
                XM0 = RHO0 * W0 * W0 * Pi * B0 * B0 / 2.E0;
            }
        } else {
            //# Initial fluxes
            Q0 = B0*B0*W0;
            BU0 = Q0*TH0;
            XM0 = Q0*W0;
        }

        //# Start inner loop
        iTop[pCount] = level;

        for (i = level; i <= NSY; i++) {
            if (W0 <= 0.0) break;   //# a no-rise plume (i == level) is handled below

            if (i == 1) {
                //# index 0 is the virtual bottom node (all zeros) - use values
                //# one sub-layer up rather than differencing against it
                DZ = yDepth[1] - yDepth[0];
                DRH = (yDen[2] - yDen[1]) / (yDepthM[2] - yDepthM[1]);
            } else {
                DZ = yDepth[i] - yDepth[i-1];

                if (i == NSY) {
                    DRH = (yDen[NSY]-yDen[NSY-1]) / (yDepthM[NSY]-yDepthM[NSY-1]);
                } else {
                    DRH = (yDen[i+1]-yDen[i-1]) / (yDepthM[i+1]-yDepthM[i-1]);
                }
            }

            SMALLZ = BIGH - (yDepth[i] - bDepth);
            DRH = DRH / RHO0;

            //#  Calculate initial derivatives
            if (intopt) {
                if (DIM == 2) {
                    DQDZ = 2.0 * alp2 * pWidth * W0;
                    DXMDZ = sqrt(Pi) * RHO0 * TH0 * lambda * B0 * pWidth;
                    DBUDZ = RHO0 * sqrt(Pi) * B0 * pWidth * W0 * g * DRH + FLOW*PA *
                              (DW0DZ*(1.E0-(W0 / (W0+ub2d)))+(W0/SMALLZ))/
                              (SMALLZ*(W0+ub2d));
                } else {
                    DQDZ = 2.0 * ALP * B0 * W0 * Pi;
                    DXMDZ = Pi * RHO0 * TH0 * SQLAM * B0 * B0;
                    DBUDZ = RHO0 * Pi * B0 * B0 * W0 * g * DRH + FLOW * PA *
                              (DW0DZ*(1.-W0/(W0+UB))+ W0/SMALLZ)/(SMALLZ*(W0+UB));
                }
            } else {
                DQDZ = 2.0 * ALP * B0 * W0;
                DXMDZ = 2.0 * g * SQLAM * B0 * B0 * TH0;
                DBUDZ = (1.0 + SQLAM) * DRH * Q0 / SQLAM +
                          ((1.+SQLAM)*FLOW*PA/(g*Pi*SQLAM*RHO0*SMALLZ*(W0+UB))) *
                                             (W0/SMALLZ + (DW0DZ)*(1.-W0/(W0+UB)));
            }

            DYDX[0] = DQDZ;
            DYDX[1] = DXMDZ;
            DYDX[2] = DBUDZ;
            Y[0] = Q0;
            Y[1] = XM0;
            Y[2] = BU0;

            if (intopt) {
                if (i != level){
                    if (DIM==2) {
                        W0 = sqrt(2.0)*Y[1]/(RHO0*Y[0]);
                        TH0 = Y[2]*(sqrt(1.0+lm2))/(Y[0]*lambda*R2);
                    } else {
                        W0 = 2.0*Y[1]/(RHO0*Y[0]);
                        TH0 = Y[2]*(1.0+SQLAM)/(SQLAM*Y[0]*R2);
                    }
                }
            } else {
                if (i != level){
                    W0 = Y[1]/Y[0];
                    TH0 = Y[2]/Y[0];
                }
            }

            if (intopt)
                InActRK4(Y, DYDX, threeI, DZ, YOUT, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0, DIM);
            else
                InDepRK4(Y, DYDX, threeI, DZ, YOUT, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0);

            if (YOUT[1] <= 0.0) break;

            Q0 = YOUT[0];
            XM0 = YOUT[1];
            BU0 = YOUT[2];

            DQ = YOUT[0]-Y[0];

            //# Use test for change in sign of momentum flux and also
            //# rate of change of momentume flux.  the latter is necessary as
            //# erroneous results will be produced if the momentum just falls
            //# short of going negative
            GRADMO = Y[1]/YOUT[1];
            GRADMO = fabs(GRADMO);

            if (XM0 > 0.0 && Q0 > 0.0) {
                if (intopt) {
                    if (DIM==2)
                        W[i] = sqrt(2.0)*XM0/(RHO0*Q0);
                    else
                        W[i] = 2.0*XM0/((RHO0)*Q0);
                } else
                    W[i] = XM0 / Q0;
            } else
                break;

            if (GRADMO > 1.5)
                break;

            if (W[i] > 0 && Q0 > 0.) {
                if (intopt) {
                    if (DIM==2)
                        B[i] = Q0/(W[i]*pWidth*sqrt(Pi));
                    else
                        B[i] = sqrt(Q0/(W[i]*Pi));
                } else
                    B[i] = sqrt(Q0 / W[i]);
            } else //#  Change to avoid divison by zero and sqrt(-ve)
                B[i] = B0;

            DW = W[i] - W0;

            if (intopt) {
                if (Q0 != 0.0) {
                    if (DIM==2)
                        TH[i] = BU0*(sqrt(1.0+lm2))/(Q0*lambda*R2);
                    else
                        TH[i] = BU0*(1.0+SQLAM)/(SQLAM*Q0*R2);
                } else {
                    if (DIM==2)
                        TH[i] = BU0*(sqrt(1.0+lm2))/(Y[0]*lambda*R2);
                    else
                        TH[i] = BU0*(1.0+SQLAM)/(SQLAM*Y[0]*R2);
                }

                if (i==level)
                    ENT[i] = Q0;
                else
                    ENT[i] = DQ;

                DIM0 = DIM;

                if (B[i]>=LPP) DIM = 2;

                if (DIM0==1 && DIM==2) {
                    B2 = (sqrt(2.0*Pi))*(SQR(B[i]))/pWidth;
                    W2 = W[i]/sqrt(2.0);
                    TH2 = TH[i] * SQLAM * (sqrt(1.0 + lambda * lambda))/(lambda * (SQLAM + 1.0));
                    XM2 = W2 * W2 * sqrt(Pi) * B2 * pWidth / sqrt(2.0);
                    XM1 = W[i] * W[i] * Pi * B[i] * B[i] / 2.0;
                    XM0 = XM0 * XM2 / XM1;
                    BU01 = RHO0 * SQLAM * B[i] * B[i] * TH[i] * Pi * W[i] / (1.0 + SQLAM);
                    BU02 = RHO0 * lambda * B2 * TH2 * sqrt(Pi) * W2 * pWidth / (sqrt(1.0+lm2));
                    BU0 = BU02 / BU01 * BU0;
                    B[i] = B2;
                    W[i] = W2;
                    TH[i] = TH2;
                }
            } else {
                if (Q0 != 0.0)
                    TH[i] = BU0 / Q0;
                else
                    TH[i] = BU0 / Y[i];

                if (i==level)
                    ENT[i] = Pi * Q0;
                else
                    ENT[i] = Pi * DQ;
            }

            B0 = B[i];
            W0 = W[i];
            TH0 = TH[i];
            DW0DZ = DW/DZ;

        }

        if ( i < NSY ) {
            //# Here when plume has reached top
            //# set iTop(pCount) to level before last one processed
            //# this is because solution is inaccurate as W --> 0
            //# height of rise of initial plume
            iTop[pCount] = i-1;

            if (pCount == 1)
                htRise = yDepth[i-1];

        } else {
            //# Top of water column met
            iTop[pCount] = NSY;
            if (pCount == 1) htRise = yDepth[NSY];
        }

        //# plume has stopped at layer iTop(pCount)
        //# total discharge is Pi*Q0, with ENT(j) from
        //# each layer between level and iTop(pCount).
        //# Q0 includes entrainment from iTop
        //# NOTE that ENT(j) comes from each port; need nPorts*ENT(j) from layer j

        //# Check to see if iTop(pCount) equals surfLayer (IE. top of water column)
        //# If not set level to iTop(pCount) + 1
        if (iTop[pCount] < NSY) {
            //# The restart level must strictly increase: a plume that failed to
            //# rise even one sub-layer leaves iTop[pCount] = level-1, and would
            //# otherwise freeze level while pCount grows unbounded (stack overflow)
            if (iTop[pCount] < level) {
                pCount--;               //# discard the failed plume
                *pCountI = pCount;
                if (pCount == 0) return FALSE;   //# no plume formed at all
                break;
            }
            level = iTop[pCount] + 1;
            pCount++;
            *pCountI = pCount;
        } else
            break;
    } while (1);
    /********************%%%%%%%%%%%%%%%%%%%%%%%%%%%%**********************/
    *retI = i;
    return TRUE;
}


/******************************************************************************
 *                                                                            *
 *                                                                            *
 ******************************************************************************/
void inBubble(LOGICAL intopt, int MaxNSy)
{
//CONSTANTS
    const AED_REAL US=0.3E0, HPRES=10.2E0;

    int level, iLevel, i, ii, j, jj, iRiver, ll, IUP;
    int NINT, NINT0, pCount, iSkip;
    int NSY, indexY;

    AED_REAL BSL;
    AED_REAL FLOW;
    AED_REAL RHO0;
    AED_REAL QVOLT, QVOL;
    AED_REAL botH;

#ifndef _VISUAL_C_
    int iTop[MaxNSy];
    AED_REAL T[MaxNSy], S[MaxNSy], QF[MaxNSy];
    AED_REAL PLMDEN[MaxNSy], BFSQ[MaxNSy];
    AED_REAL yDepth[MaxNSy], yDepthM[MaxNSy], YVOL[MaxNSy];
    AED_REAL ENT[MaxNSy], W[MaxNSy], B[MaxNSy], TH[MaxNSy];
    AED_REAL yTemp[MaxNSy], ySal[MaxNSy], yDen[MaxNSy];
    AED_REAL ywqual[Num_WQ_Vars * MaxNSy];
    AED_REAL xw[Num_WQ_Vars * MaxNSy];
    AED_REAL wqx[Num_WQ_Vars];
#else
    int *iTop;
    AED_REAL *T, *S, *QF;
    AED_REAL *PLMDEN, *BFSQ;
    AED_REAL *yDepth, *yDepthM, *YVOL;
    AED_REAL *ENT, *W, *B, *TH;
    AED_REAL *yTemp, *ySal, *yDen;
    AED_REAL *xw, *wqx, *ywqual;
#endif

#define _ywqual(i, j) ywqual[_IDX_2d(Num_WQ_Vars, MaxNSy, i, j)]
#define _xw(i, j)         xw[_IDX_2d(Num_WQ_Vars, MaxNSy, i, j)]

    AED_REAL BD, XMom0, HL, HU, DENTOP, DENBOT;

    AED_REAL tWidth;
    LOGICAL  done;

/*----------------------------------------------------------------------------*/

    //# Calculated parameters and initialization
    IUP = 0;
    HU = 0.;
    HL = 0.;

#ifdef _VISUAL_C_
    ALLOCATE1D_I(iTop, MaxNSy);

    ALLOCATE1D_R(T, MaxNSy);      ALLOCATE1D_R(S, MaxNSy);    ALLOCATE1D_R(QF, MaxNSy);
    ALLOCATE1D_R(PLMDEN, MaxNSy); ALLOCATE1D_R(BFSQ, MaxLayers);   //# BFSQ is layer-indexed
    ALLOCATE1D_R(yDepth, MaxNSy);
    ALLOCATE1D_R(yDepthM, MaxNSy); ALLOCATE1D_R(YVOL, MaxNSy); ALLOCATE1D_R(ENT, MaxNSy);
    ALLOCATE1D_R(W, MaxNSy);     ALLOCATE1D_R(B, MaxNSy);    ALLOCATE1D_R(TH, MaxNSy);
    ALLOCATE1D_R(yTemp, MaxNSy); ALLOCATE1D_R(ySal, MaxNSy); ALLOCATE1D_R(yDen, MaxNSy);
    ALLOCATE1D_R(wqx,  Num_WQ_Vars);

    ALLOCATE2D(ywqual, Num_WQ_Vars, MaxNSy);
    ALLOCATE2D(xw,     Num_WQ_Vars, MaxNSy);
#endif

    htRise = 0.;
    cee = 0.;
    DENTOP = 0.;
    DENBOT = 0.;
    lambda = 0.3;
    alp2 = 0.102;
    iSkip = 0;
    if (Lake[surfLayer].Temp == Lake[botmLayer].Temp) iSkip = 1;
    lm2 = lambda*lambda;
    pWidth = bubLen/(nPorts);

    BD = bDepth;
    FLOW = aFlow/(nPorts);

    //# Basic bubbler time step = 15 Min
    NINT0 = noSecs;
    NINT = 900;

    //# Initialize sublayer variables
    for (i = 0; i < MaxNSy; i++) {
        yDepth[i] = 0.0;
        yDepthM[i] = 0.0;
        YVOL[i] = 0.0;
        yTemp[i] = 0.0;
        ySal[i] = 0.0;
        yDen[i] = 0.0;
    }

    ub2 = US * (1.0+lm2);
    ub2d = US * sqrt(1.0+lm2);

    do  {
        //# Find center of buoyancy; XMoment1 is 1st moment, XMom0 is 0th moment
        //# use this to estimate two-layer equivalent structure for purposes
        //# of calculating cee
        if (iSkip != 1) {
            for (i = 0; i < MaxLayers; i++) BFSQ[i] = 0.;

            XMoment1 = 0.;
            XMom0 = 0.;

            for (i = 1; i <= surfLayer; i++) {
                BFSQ[i] = gprime(Lake[i].Density, Lake[i-1].Density)/(Lake[i].MeanHeight-Lake[i-1].MeanHeight);
                XMoment1 += Lake[i-1].Height * BFSQ[i] * (Lake[i].LayerVol + Lake[i-1].LayerVol) / 2.;
                XMom0 +=                       BFSQ[i] * (Lake[i].LayerVol + Lake[i-1].LayerVol) / 2.;
            }

            //# Define length scales, XMoment1 is the center of buoyancy ( M above bot
            //#                       HU is the thickness of the upper mixed layer
            if (XMom0 != 0.) {
                HL = XMoment1 / XMom0;
                HU = Lake[surfLayer].Height - HL;
            } else {
                HL = Lake[surfLayer].Height;
                HU = 0.;
            }

            //# Find layer containing the centre of buoyancy (layer i)
            for (i = botmLayer; i <= surfLayer; i++)
                if (Lake[i].Height > HL) break;

            DENBOT = 0.;
            DENTOP = 0.;

            IUP = i;
            if (IUP != surfLayer) {
                for (i = botmLayer; i < IUP; i++)
                    DENBOT += Lake[i].Density*Lake[i].LayerVol;
                DENBOT = DENBOT/Lake[IUP].Vol1;

                for (i = IUP; i <= surfLayer; i++)
                    DENTOP += Lake[i].Density*Lake[i].LayerVol;
                DENTOP = DENTOP/(Lake[surfLayer].Vol1-Lake[IUP].Vol1);
            } else
                cee = 0.0;
        }

        //# Set first value of level
        iLevel = surfLayer;

        for (i = surfLayer; i >= botmLayer; i--)
            if (Lake[i].Height >= bDepth) iLevel = i;

        //# Set conditions at starting level
        //# (GLM stores absolute density; the +1000 sigma-t offset was a DYRESM-ism)
        RHO0 = Lake[iLevel].Density;

        //# pCount is number of plumes
        pCount = 1;
        level = iLevel;

        //# plume number for 2-layer approximation - assumes bubbler on bottom
        if (IUP == surfLayer)
            cee = 0.;
        else {
            if (DENBOT > DENTOP) {
                cee = pow(((DENBOT-DENTOP)/RHO0)*HL, 1.5)*Lake[surfLayer].Height *
                                      pow(HU/Lake[surfLayer].Height, 3.)/FLOW;
                cee = cee * ((Lake[surfLayer].Height + HPRES) / HPRES) * RHO0 / 1000.;
            } else
                cee = 0.0;
        }

        //# Now get into bubbler internal substeps
        //# Sub-layer convention (1-based, from the DYRESM Fortran): layer i
        //# (botmLayer = 0) owns sub-layers [i*nStepY+1 .. (i+1)*nStepY], so the
        //# top of the column is sub-layer NSY.  Index 0 is a virtual node at
        //# the lake bottom (all arrays zeroed above).  Every sub-layer of every
        //# layer must be filled or volumes are lost in the write-back below.
        NSY = (surfLayer+1) * nStepY;

        for (i = surfLayer; i >= botmLayer; i--) {
            botH = (i == botmLayer) ? 0.0 : Lake[i-1].Height;
            for (j = 1; j <= nStepY; j++) {
                indexY = (i+1)*nStepY-j+1;
                yDepth[indexY] = Lake[i].Height - (Lake[i].Height - botH) *
                                    (j-1)/(nStepY);
                yDepthM[indexY] = Lake[i].Height - (Lake[i].Height - botH) *
                                    ((j-1) + 0.5)/(nStepY);
                YVOL[indexY] = Lake[i].LayerVol/(nStepY);
                yTemp[indexY] = Lake[i].Temp;
                ySal[indexY] = Lake[i].Salinity;
                yDen[indexY] = Lake[i].Density;

                for (ii = 0; ii < Num_WQ_Vars; ii++)
                    _ywqual(ii, indexY) = _WQ_Vars(ii, i);
            }
        }

        //# Plume origin: top sub-layer of the layer containing the diffuser
        level = (iLevel+1) * nStepY;

        for (i = 0; i < MaxNSy; i++) {
            ENT[i] = 0.0;
            W[i] = 0.0;
            B[i] = 0.0;
            TH[i] = 0.0;
        }

        done = FALSE;
        do  {
            if ( !done && inner_loop(iLevel, level, &pCount, &i, NSY, iTop, yDepth, yDepthM, TH, W, B, ENT, yDen) ) {
                //# Multiply ENT(j) by the number of ports
                for (j = 1; j <= NSY; j++)
                    ENT[j] = (nPorts) * ENT[j];

                //# Check if new time step required
                if (NINT0 < NINT) NINT = NINT0;

                //# Fix entrainments
                //# ENT comes back from inner_loop as m3/s per port; convert to a
                //# volume over the NINT substep.  (DYRESM scaled by 1/1000 here
                //# because its volumes were in ML; GLM is SI so no such factor.)
                for (i = 1; i <= NSY; i++) {
                    ENT[i] *= (NINT);
                    YVOL[i] -= ENT[i];

                    if (YVOL[i] <= 0.0) {
                        YVOL[i] += ENT[i];
                        if (i != NSY)   //# push the unmet demand up, back in m3/s
                            ENT[i+1] += (ENT[i] - 0.9*YVOL[i])/(NINT);

                        ENT[i] = 0.9 * YVOL[i];
                        YVOL[i] -= ENT[i];
                    }
                }

                //# Fix properties of plume
                //# -incoming
                i = 1;
                level = 1;
                done = TRUE;
            }

            if ( done ) {
                T[i] = 0.0;
                S[i] = 0.0;
                PLMDEN[i] = 0.0;

                for (ii = 0; ii < Num_WQ_Vars; ii++) _xw(ii, i) = 0.0;

                //# -entrainment
                QF[i] = 0.0;

                for (j = level; j <= iTop[i]; j++) {   //# ENT is set up to and including iTop
                    if (ENT[j] > 0.0) {
                        T[i] = combine(T[i], QF[i], PLMDEN[i], yTemp[j], ENT[j], yDen[j]);
                        S[i] = combine(S[i], QF[i], PLMDEN[i],  ySal[j], ENT[j], yDen[j]);
                        PLMDEN[i] = calculate_density(T[i], S[i]);

                        for (ii = 0; ii < Num_WQ_Vars; ii++)
                            _xw(ii, i) = combine_vol(_xw(ii, i), QF[i], _ywqual(ii, j), ENT[j]);

                        QF[i] += ENT[j];
                    }
                }
                done = FALSE;
            }

            //# Check for next plume
            //# Restart at layer iTop+1 and go to top of next plume
            if (i < pCount) {
                level = iTop[i] + 1;
                i++;
                done = TRUE;
            }
        } while (done);

        //# Now swap back to the regular layer size
        //# Volume is only physical property we need change
        for (ii = botmLayer; ii <= surfLayer; ii++) {
            Lake[ii].LayerVol = 0.;
            for (jj = 1; jj <= nStepY; jj++)
                Lake[ii].LayerVol += YVOL[ii*nStepY+jj];

            Lake[ii].Temp = yTemp[(ii+1)*nStepY];
            Lake[ii].Salinity = ySal[(ii+1)*nStepY];
            Lake[ii].Density = yDen[(ii+1)*nStepY];

            for (jj = 0; jj < Num_WQ_Vars; jj++)
                _WQ_Vars(jj, ii) = _ywqual(jj, (ii+1)*nStepY);
        }

        //# Find the level of intrusion of the water when it falls back
        //# and insert it using the routine from the DoInFlow routine.
        //# assume the basin length at the bubbler is determined by the slope of
        //# the major DoInFlow valley, I.E river 1.
        iRiver = 0;
        BSL = Inflows[iRiver].Phi * Pi / 180.0;

        //# Make insertion WIDTH 0, then insert will use reservoir width at discha
        QVOLT = 0.0;

        for (i = 1; i <= pCount; i++) {   //# plumes are indexed 1..pCount
            for (jj = 0; jj < Num_WQ_Vars; jj++) wqx[jj] = _xw(jj, i);

            QVOLT = QVOLT+QF[i];
            QVOL = QF[i];

            tWidth = 1.E-8;
            insert(QVOL, PLMDEN[i], BSL, T[i], S[i], wqx, NINT, &tWidth, &ll);
        }

        //# Fix depths with resint
        resize_internals(2, 1);

        //# Fix layers
        check_layer_thickness();

        //# Reinitialize bubbler depth to original (needed by pe3)
        bDepth = BD;

        //# Check if number of seconds left
        NINT0 = NINT0 - NINT;
        if (NINT0 > 0)
            iSkip = 1;
        else
            break;
    } while (1);//# Finished!

#ifdef _VISUAL_C_
    free(iTop);
    free(T);       free(S);      free(QF);
    free(PLMDEN);  free(BFSQ);   free(yDepth);
    free(yDepthM); free(YVOL);   free(ENT);
    free(W);       free(B);      free(TH);
    free(yTemp);   free(ySal);   free(yDen);
    free(wqx);     free(ywqual); free(xw);
#endif
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void InActRK4(AED_REAL *Y, AED_REAL *DYDX, int N, AED_REAL H,
              AED_REAL *YOUT, AED_REAL B0, AED_REAL SQLAM, AED_REAL DRH,
              AED_REAL ALP, AED_REAL W0, AED_REAL FLOW, AED_REAL PA,
              AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0, int DIM)
{
    AED_REAL YT[3], DYT[3], DYM[3], HH, H6;
    int i;

/*----------------------------------------------------------------------------*/

    HH = H*0.5;
    H6 = H/6;

    for (i = 0; i < N; i++)
        YT[i] = Y[i] + HH * DYDX[i];

    if (YT[1] <= 0.0) {
        YOUT[0] = YT[0];
        YOUT[1] = YT[1];
        YOUT[2] = YT[2];
        return;
    }

    InActDrv(HH, YT, DYT, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0, DIM);

    for (i = 0; i < N; i++)
        YT[i] = Y[i] + HH * DYT[i];

    if (YT[1] <= 0.0) {
        YOUT[0] = YT[0];
        YOUT[1] = YT[1];
        YOUT[2] = YT[2];
        return;
    }

    InActDrv(HH, YT, DYM, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0, DIM);

    for (i = 0; i < N; i++) {
        YT[i] = Y[i] + H * DYM[i];
        DYM[i] = DYT[i] + DYM[i];
    }
    if (YT[1] <= 0.0) {
        YOUT[0] = YT[0];
        YOUT[1] = YT[1];
        YOUT[2] = YT[2];
        return;
    }

    InActDrv(HH, YT, DYT, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0, DIM);

    for (i = 0; i < N; i++)
        YOUT[i] = Y[i] + H6*(DYDX[i] + DYT[i] + 2.*DYM[i]);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void InActDrv(AED_REAL X, AED_REAL *Y, AED_REAL *DYDX, AED_REAL B0,
              AED_REAL SQLAM, AED_REAL DRH, AED_REAL ALP, AED_REAL W0,
              AED_REAL FLOW, AED_REAL PA, AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0, int DIM)
{
//CONSTANTS
    const AED_REAL R2 = 1000.;

//LOCALS
    AED_REAL W, B, TH;

/*----------------------------------------------------------------------------*/

    if (fabs(Y[0]) <= 1.E-7) {
        W = W0;
        TH = TH0;
    } else {
        if (DIM==2) {
            W = (sqrt(2.0)) * Y[1] / (RHO0 * Y[0]);
            TH = Y[2] * (sqrt(1.0+lm2)) / (Y[0] * lambda * R2);
        } else {
            W = 2.0 * Y[1] / (RHO0*Y[0]);
            TH = Y[2] * (1.0 + SQLAM) / (SQLAM * Y[0] * R2);
        }
    }

    if (W > 0. && Y[0] > 0.) {
        if (DIM==2)
            B = Y[0] / (W * pWidth * sqrt(Pi));
        else
            B = sqrt(Y[0] / (W * Pi));
    } else
        B = B0;

    if (DIM==2) {
        DYDX[0] = 2.0 * alp2 * pWidth * W;
        DYDX[1] = sqrt(Pi) * RHO0 * TH * lambda * B * pWidth;
        DYDX[2] = sqrt(Pi) * B * pWidth * W * g * DRH * RHO0 +
                  FLOW * PA * (((W-W0)/X) * (1. - W/(W + ub2d)) + W/SMALLZ) /
                  (SMALLZ * (W+ub2d));
    } else {
        DYDX[0] = 2.0 * ALP * B * W * Pi;
        DYDX[1] = Pi * RHO0 * TH * SQLAM * B * B;
        DYDX[2] = Pi * B * B * W * g * DRH * RHO0 +
                  FLOW * PA * (((W-W0)/X) * (1. - W/(W + UB)) + W/SMALLZ) /
                  (SMALLZ * (W+UB));
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void InDepRK4(AED_REAL *Y, AED_REAL *DYDX, int N, AED_REAL H,
              AED_REAL *YOUT, AED_REAL B0, AED_REAL SQLAM, AED_REAL DRH,
              AED_REAL ALP, AED_REAL W0, AED_REAL FLOW, AED_REAL PA,
              AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0)
{
    AED_REAL YT[3], DYT[3], DYM[3], HH, H6;

    int i;

/*----------------------------------------------------------------------------*/

    HH = H/2;
    H6 = H/6;

    for (i = 0; i < N; i++)
        YT[i] = Y[i] + HH * DYDX[i];

    if (YT[1] <= 0.0) {
        YOUT[0] = YT[0];
        YOUT[1] = YT[1];
        YOUT[2] = YT[2];
        return;
    }

    InDepDrv(HH, YT, DYT, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0);

    for (i = 0; i < N; i++)
        YT[i] = Y[i] + HH*DYT[i];

    if (YT[1] <= 0.0) {
        YOUT[0] = YT[0];
        YOUT[1] = YT[1];
        YOUT[2] = YT[2];
        return;
    }

    InDepDrv(HH, YT, DYM, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0);

    for (i = 0; i < N; i++) {
        YT[i] = Y[i] + H*DYM[i];
        DYM[i] = DYT[i] + DYM[i];
    }
    if (YT[1] <= 0.0) {
        YOUT[0] = YT[0];
        YOUT[1] = YT[1];
        YOUT[2] = YT[2];
        return;
    }

    InDepDrv(H, YT, DYT, B0, SQLAM, DRH, ALP, W0, FLOW, PA, UB, SMALLZ, RHO0, TH0);

    for (i = 0; i < N; i++)
        YOUT[i] = Y[i] + H6*(DYDX[i] + DYT[i] + 2.*DYM[i]);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void InDepDrv(AED_REAL X, AED_REAL *Y, AED_REAL *DYDX, AED_REAL B0,
              AED_REAL SQLAM, AED_REAL DRH, AED_REAL ALP, AED_REAL W0,
              AED_REAL FLOW, AED_REAL PA, AED_REAL UB, AED_REAL SMALLZ,
              AED_REAL RHO0, AED_REAL TH0)
{
    AED_REAL W, B, TH;

/*----------------------------------------------------------------------------*/

    if (fabs(Y[0]) <= 1.E-7) {
        W = W0;
        TH = TH0;
    } else {
        W = Y[1]/Y[0];
        TH = Y[2]/Y[0];
    }

    if (W > 0. && Y[0] > 0.)
        B = sqrt(Y[0] / W);
    else
        B = B0;

    DYDX[0] = 2.*ALP*B*W;
    DYDX[1] = 2.*9.81*SQLAM*SQR(B)*TH;
    DYDX[2] = (1.+SQLAM)*DRH*Y[0]/SQLAM
             +((1.+SQLAM)*FLOW*PA/(g*Pi*SQLAM*RHO0*SMALLZ*(W+UB)))*(W/SMALLZ +
                                                      ((W-W0)/X)*(1.-W/(W+UB)));
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
