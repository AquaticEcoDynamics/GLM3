/******************************************************************************
 *                                                                            *
 * glm_main.c                                                                 *
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
#include <time.h>

#include "glm.h"

#include "glm_types.h"
#include "glm_globals.h"

#ifdef PLOTS
   #include <libplot.h>
   extern FLOGICAL do_plots;
   extern CLOGICAL saveall;
   extern char *plots_nml_name;
#ifdef XPLOTS
   extern int xdisp;
#endif
#endif
#if API
   #include <aed_api.h>
#endif
#if AED
 #ifdef WITH_AED_PLUS
   #include <aed+.h>
 #else
   #include <aed.h>
 #endif
#endif

#include "glm_debug.h"

char *all_plots_name = NULL;

extern char glm_nml_file[];
extern void run_model(void);

#ifdef PLOTS
/*----------------------------------------------------------------------------*/
static char *plot_filename()
{
    static char nam[256];
    strcpy(nam, "GLM.png");
    return nam;
}
#endif

/******************************************************************************/
int main(int argc, char *argv[])
{
    char *nmlfile = NULL;
    int all_ok = 1, show_options = 0, show_vers = 0;

    // Set random number seed
    srand ((unsigned int)time(NULL));

#ifdef PLOTS
    init_plotter_main(argv[0]);

    saveall = 0;
#ifdef XPLOTS
    xdisp = 0;
    int no_gui = 0;
#endif
#endif

    argv++; argc--;
    while (argc > 0) {
        if ( strcmp(*argv, "--nml") == 0) {
            argv++; argc--;
            nmlfile = *argv;
        }
        else if (strcmp(*argv, "--debug") == 0) {
            _glm_dbg_on();
        }
        else if (strcmp(*argv, "--dbgmix") == 0) {
            _mix_dbg_on();
        }
#ifdef PLOTS
#ifdef XPLOTS
        else if (strcmp(*argv, "--xdisp") == 0) {
            if ( !no_gui ) xdisp = 1;
            if ( argc > 1 && strncmp(argv[1], "--", 2) != 0 ) {
                argv++; argc--;
                plots_nml_name = *argv;
            }
        }
        else if (strcmp(*argv, "--no-gui") == 0) {
            no_gui = 1; xdisp = 0;
        }
#endif
        else if (strcmp(*argv, "--saveall") == 0) {
            saveall = TRUE;
        }
        else if (strcmp(*argv, "--save-all-in-one") == 0) {
            saveall = TRUE;
            all_plots_name = plot_filename();
            if ( argc > 1 && strncmp(argv[1], "--", 2) != 0 ) {
                argv++; argc--;
                all_plots_name = *argv;
            }
        }
#endif
        else if (strcmp(*argv, "--write_nml") == 0) {
            glm_write_nml_mode = TRUE;
            glm_example_fname = "glm4.nml";
            if ( argc > 1 && strncmp(argv[1], "--", 2) != 0 ) {
                argv++; argc--;
                glm_example_fname = *argv;
            }
        }
        else if (strcmp(*argv, "--quiet") == 0) {
            quiet = 5;
            if ( argc > 1 && strncmp(argv[1], "--", 2) != 0 ) {
                argv++; argc--;
                quiet = atoi(*argv);
            }
        }
        else if ( (*argv)[0] != '-' ) {
            // assume its the run file name
            nmlfile = *argv;
        }
        else {
            if (strcmp(*argv, "--version") == 0) {
                show_vers = 1;
            } else if (strcmp(*argv, "--help") == 0) {
                show_options = 1;
            } else {
                fprintf(stderr, "Unknown flag %s\n", *argv);
            }
            all_ok = 0;
        }
        argc--; argv++;
    }

#ifdef PLOTS
# ifdef XPLOTS
    do_plots = xdisp || saveall;
# else
    do_plots = saveall;
# endif
#endif

    if (quiet < 10 || show_vers || show_options ) {
        printf("\n");
        printf("    -------------------------------------------------------\n");
        printf("    |  General Lake Model (GLM)   Version %-12s    |\n", GLM_VERSION);
        printf("    -------------------------------------------------------\n");
        printf("\n");

#ifdef __GNUC__
        printf("     glm built using gcc version %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#elif defined(_MSC_VER) && !defined(__INTEL_COMPILER)
        printf("     glm built using MSC version %ld\n", _MSC_VER);
#endif
#ifndef _WIN32
        printf("     build date %s\n", BUILDDATE);
#endif
    }
    if ( show_vers ) {
        printf("\n");
#if API
        printf("        Includes API version %s\n", AED_API_VERSION);
#endif
#if AED
        printf("        Includes AED version %s\n", AED_VERSION);
#ifdef AED_PLUS_VERSION
        printf("        Includes AED+ version %s\n", AED_PLUS_VERSION);
#endif
#endif
    } else if ( show_options ) {
        printf("\n");
        printf("     --help  : show this blurb\n");
        printf("     --version  : report version numbers\n");
        printf("\n");
        printf("     --nml <nmlfile> : get parameters from nmlfile\n");
        printf("\n");
        printf("     --write_nml <nmlfile> : write a baseline/example config and exit (default 'glm4.nml')\n");
        printf("\n");
#ifdef PLOTS
#ifdef XPLOTS
        printf("     --xdisp : display temp/salt and selected others in x-window\n");
        printf("     --xdisp <plotsfile> : like --xdisp, but use <plotsfile> instead of plots.nml\n");
#endif
        printf("\n");
        printf("     --saveall : save plots to png files\n");
        printf("     --save-all-in-one : save all plots to png file\n");
        printf("     --save-all-in-one <destfile> : save all plots to png file <destfile>\n");
        printf("\n");
        printf("     --quiet   : less messages\n");
        printf("     --quiet <level> : set quiet level (1-10)\n");
#endif
    } else if ( all_ok ) {
        if ( nmlfile != NULL ) strncpy(glm_nml_file, nmlfile, 256);

        run_model();

        if (quiet < 10) {
            printf("    Model Run Complete\n");
            printf("    -------------------------------------------------------\n\n");
        }
    }

    exit(0);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
