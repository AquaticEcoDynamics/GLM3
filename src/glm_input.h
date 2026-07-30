/******************************************************************************
 *                                                                            *
 * glm_input.h                                                                *
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
#ifndef _GLM_INPUT_H_
#define _GLM_INPUT_H_

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void open_inflow_file(int inf_id, const char *fname, const char *timefmt);
void index_inflow_file(int inf_id, int nvars, const char *vars[]);
void read_daily_inflow(int julian, int NumInf, AED_REAL *flow,
                                               AED_REAL *temp, AED_REAL *salt,
                                                  AED_REAL *elev, AED_REAL *wq);
void open_gw_file(int inf_id, const char *fname,
                            int nvars, const char *vars[], const char *timefmt);
void read_daily_gw(int julian, int NumGW, AED_REAL *flow);

void open_outflow_file(int i, const char *fname, const char *timefmt);
void read_daily_outflow(int julian, int NumOut, AED_REAL *drw, AED_REAL *elev, AED_REAL *heat_flux, AED_REAL *wq);

void open_withdrtemp_file(const char *fname, const char *timefmt);
void open_evap_file(const char *fname, const char *timefmt);
void read_daily_evap(int julian, AED_REAL *evap);

void read_daily_withdraw_temp(int julian, AED_REAL *withdrTemp);

void open_met_file(const char *fname, int snow_sw, int rain_sw,
                                                           const char *timefmt);
void read_daily_met(int julian, MetDataType *met);
void open_kw_file(const char *fname, const char *timefmt);
void read_daily_kw(int julian, AED_REAL *kwout);

AED_REAL get_fetch(AED_REAL windDir);
void read_sub_daily_met(int julian,int iclock, MetDataType *met);

void close_evap_files(void);
void close_kw_files(void);
void close_met_files(void);
void close_inflow_files(void);
void close_outflow_files(void);
void close_withdrtemp_files(void);

void open_bubbler_file(const char *fname, const char *timefmt);
void read_bubble_data(int julian, AED_REAL *aFlow, int *nPorts,
                                           AED_REAL *bDepth, AED_REAL *bLength);

extern int lw_ind;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif
