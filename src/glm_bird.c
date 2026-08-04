/******************************************************************************
 *                                                                            *
 * glm_bird.c                                                                 *
 *                                                                            *
 * Bird and Hulstrom's Solar Irradiance Model                                 *
 * Richard E. Bird and Roland L. Hulstrom                                     *
 * A Simplified Clear Sky Model for Direct and Diffuse Insolation on          *
 *                                                       Horizontal Surfaces. *
 * SERI/TR-642-761                                                            *
 * Solar Energy Research Institute                                            *
 * Golden, Colorado, USA, February 1981.                                      *
 *                                                                            *
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
#include <math.h>

#include "glm.h"
#include "glm_const.h"
#include "glm_util.h"
#include "glm_bird.h"
#include "glm_globals.h"

#include <aed_time.h>
#include <namelist.h>


static AED_REAL AP = 973;          //# Atmospheric Pressure in milibars
static AED_REAL Oz = 0.279;        //# Ozone concentration in atm-cm
static AED_REAL WatVap = 1.1;      //# Total Precipitable water vapor in atm-cm
static AED_REAL AOD500 = 0.033;    //# Dimensionless Aerosol Optical Depth at wavelength 500 nm
static AED_REAL AOD380 = 0.038;    //# Dimensionless Aerosol Optical Depth at wavelength 380 nm
static AED_REAL Albedo = 0.2;      //# Albedo value kept at default


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
int config_bird(int namlst)
{
    /*---------------------------------------------
     * bird_model
     *-------------------------------------------*/
    AED_REAL t_AP     = MISVAL;
    AED_REAL t_Oz     = MISVAL;
    AED_REAL t_WatVap = MISVAL;
    AED_REAL t_AOD500 = MISVAL;
    AED_REAL t_AOD380 = MISVAL;
    AED_REAL t_Albedo = MISVAL;
    /*-------------------------------------------*/

    NAMELIST bird_model[] = {
          { "bird_model",  TYPE_START,            NULL         },
          { "AP",          TYPE_DOUBLE,           &t_AP        },
          { "Oz",          TYPE_DOUBLE,           &t_Oz        },
          { "WatVap",      TYPE_DOUBLE,           &t_WatVap    },
          { "AOD500",      TYPE_DOUBLE,           &t_AOD500    },
          { "AOD380",      TYPE_DOUBLE,           &t_AOD380    },
          { "Albedo",      TYPE_DOUBLE,           &t_Albedo    },
          { NULL,          TYPE_END,              NULL         }
    };
    if ( glm_write_nml_mode ) {
        //# show the module's real defaults, not the MISVAL "unset" sentinels
        t_AP = AP; t_Oz = Oz; t_WatVap = WatVap;
        t_AOD500 = AOD500; t_AOD380 = AOD380; t_Albedo = Albedo;
        write_namelist(glm_example_fp, bird_model, glm_nml_describe);
        return 1;
    }

    if ( get_namelist(namlst, bird_model) )
        return 1;

    if ( t_AP     != MISVAL ) AP = t_AP;
    if ( t_Oz     != MISVAL ) Oz = t_Oz;
    if ( t_WatVap != MISVAL ) WatVap = t_WatVap;
    if ( t_AOD500 != MISVAL ) AOD500 = t_AOD500;
    if ( t_AOD380 != MISVAL ) AOD380 = t_AOD380;
    if ( t_Albedo != MISVAL ) Albedo = t_Albedo;
    return 0;
}


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
AED_REAL calc_bird(AED_REAL lon, AED_REAL lat, int jday, int iclock, AED_REAL TZ)
{
    AED_REAL phi_day;          // Day Angle :- Position of the earth in sun's orbit
    AED_REAL ETR;              // Extra Terrestrial Beam Intensity
                               // Correction of Earth Sun Distance based on elliptical path of the sun
    AED_REAL ZenithAngle;      // Zenith Angle
    AED_REAL AirMass = 0.;     // Air Mass
    AED_REAL AMp = 0.;
    AED_REAL Trayleigh = 0.;   // Rayleigh Scattering
    AED_REAL OzAM = 0.;        // Ozone Scattering
    AED_REAL Toz = 0.;
    AED_REAL Tm = 0.;          // Scattering due to mixed gases
    AED_REAL Wm = 0.;
    AED_REAL Twater = 0.;      // Scattering due to Water Vapor
    AED_REAL TauA = 0.;        // Scattering due to Aerosols
    AED_REAL Ta = 0.;
    AED_REAL Taa = 0.;
    AED_REAL Tas = 0.;
    AED_REAL rs = 0.;          // Scattered Radiation
    AED_REAL phi_db = 0.;      // Direct Beam Horizontal Radiation
    AED_REAL phi_as = 0.;
    AED_REAL GHI = 0.;         // Global Horizontal Irradiation

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    int day = day_of_year(jday);

    ZenithAngle = zenith_angle(lon, lat, day, iclock, TZ);

    // Day Angle :- Position of the earth in sun's orbit
    phi_day = (two_Pi*(day-1)/365);

    // Extra Terrestrial Beam Intensity
    // Correction of Earth Sun Distance based on elliptical path of the sun
    // where did this come from ?
    ETR = I_sc * (1.00011 + 0.034221*cos(   two_Pi*phi_day) +
                            0.00128 *sin(   two_Pi*phi_day) +
                            0.000719*cos(2*(two_Pi*phi_day)) +
                            0.000077*sin(2*(two_Pi*phi_day)) );

    // Air Mass
    if (ZenithAngle < 89)
        AirMass = 1. / (cos(ZenithAngle * deg2rad) + 0.15 / pow(93.885-ZenithAngle, 1.25) );

    if ( AirMass > 0) {
        AMp = (AirMass*AP) / 1013;

        // Rayleigh Scattering
        Trayleigh = exp(-0.0903 * pow(AMp,0.84) * (1 + AMp - pow(AMp,1.01))) ;

        // Ozone Scattering
        OzAM = Oz * AirMass;
        Toz = 1 - 0.1611 * OzAM * pow(1. + 139.48 * OzAM,-0.3035) -
                      0.002715 * OzAM / (1. + 0.044 * OzAM + 0.0003 * pow(OzAM,2));

        // Scattering due to mixed gases
        Tm = exp(-0.0127 * pow(AMp, 0.26)) ;

        // Scattering due to Water Vapor
        Wm = AirMass * WatVap;
        Twater = 1 - 2.4959 * Wm / ((pow(1. + 79.034 * Wm, 0.6828)) + 6.385 * Wm) ;

        // Scattering due to Aerosols
        TauA = 0.2758 * AOD380 + 0.35 * AOD500;
        Ta = exp((-pow(TauA,0.873)) * (1.+TauA-(pow(TauA,0.7088)))*pow(AirMass,0.9108));

        Taa = 1. - 0.1 * (1 - AirMass + pow(AirMass, 1.06)) * (1 - Ta);

        Tas = Ta / Taa;

        // Scattered Radiation
        rs = 0.0685 + (1 - 0.84) * (1 - Tas);

        // Direct Beam Radiation (Extra-Terrestrial)
        if (ZenithAngle < 90)
            phi_db = 0.9662 * ETR * Trayleigh * Toz * Tm * Twater * Ta * cos(ZenithAngle * deg2rad) ;

        phi_as = 0.79 * ETR * Toz * Tm * Twater * Taa * cos(ZenithAngle * deg2rad) *
                   (0.5 * (1-Trayleigh) + 0.84*(1.-Tas)) / (1-AirMass + pow(AirMass, 1.02));

        // Global Horizontal Irradiation
        GHI = (phi_db + phi_as)/(1 - Albedo * rs) ;
    }

    return GHI;
}


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
AED_REAL clouded_bird(AED_REAL GHI, AED_REAL cloud)
{
   // f_C = (0.66182*(cloud * cloud) - 1.5236 * cloud + 0.98475)
    return GHI * (0.66182*(cloud * cloud) - 1.5236 * cloud + 0.98475);
}


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
AED_REAL cloud_from_bird(AED_REAL GHI, AED_REAL Solar)
{
#define a   0.66182
#define b  (-1.5236)
#define c   0.98475

    AED_REAL f_C;          // Fraction of clear sky solar radiation
    AED_REAL cloud = 0.;   // Cloud cover as fraction of one
    AED_REAL tmp;

    if (Solar > 0. && GHI > 0.)
        f_C = Solar / GHI;
    else
        f_C = 1;  //Night assume CC = 1}

    tmp = (b * b) - 4 * a * (c - f_C);
    if  (tmp > 0)   // Assume fully clouded if f_C <<<< 1
        cloud = (-b - sqrt(tmp)) / 2 * a;
    else
        cloud = 1.;

    if (cloud < 0.) cloud = 0.;

    return cloud;
}
