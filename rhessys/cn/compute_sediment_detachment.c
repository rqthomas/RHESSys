/*--------------------------------------------------------------*/
/*                                                              */
/*              compute_sediment_detachment                     */
/*                                                              */
/*  NAME                                                        */
/*  compute_sediment_detachment - computes daily sediment       */
/*  mobilization per patch from rainsplash detachment.          */
/*  Adds the detached mass to patch.surface_sediment for        */
/*  subsequent routing by the overland-flow routing routines.   */
/*                                                              */
/*  SYNOPSIS                                                    */
/*  void compute_sediment_detachment(                           */
/*              struct patch_object *,                          */
/*              struct zone_object  *,                          */
/*              struct command_line_object *)                   */
/*                                                              */
/*  DESCRIPTION                                                 */
/*  Rainsplash detachment:                                      */
/*    KE (MJ/m2/day) = rain_mm * 0.029  (simplified USLE KE)  */
/*    cover_frac = 1 - exp(-0.5 * LAI)                         */
/*    C_factor: from soil defaults if in (0,1], else 1-cover   */
/*    rainsplash [kg/m2] = KE * soil_erodibility_K * C_factor  */
/*                                                              */
/*  patch.surface_sediment is the accumulating mobilised        */
/*  sediment pool.  Overland-flow routing (update_drainage_*)   */
/*  carries it downslope each timestep; the residual stays in   */
/*  the pool for future routing.                                */
/*--------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include "rhessys.h"

void compute_sediment_detachment(
    struct patch_object        *patch,
    struct zone_object         *zone,
    struct command_line_object *command_line)
{
    double rain_m;          /* daily rainfall (m/day) */
    double rain_mm;         /* daily rainfall (mm/day) */
    double KE;              /* rainfall kinetic energy (MJ/m2/day) */
    double cover_frac;      /* canopy cover fraction (0-1, dim) */
    double C_factor;        /* USLE C-factor (0-1) */
    double rainsplash;      /* rainsplash detachment (kg/m2) */

    rain_m = zone[0].rain;
    if (rain_m <= 0.0)
        return;

    rain_mm = rain_m * 1000.0;

    /* Simplified Wischmeier & Smith kinetic energy proxy    */
    KE = rain_mm * 0.029; /* MJ/m2/day */

    /* Canopy cover fraction attenuates KE to soil surface   */
    cover_frac = 1.0 - exp(-0.5 * patch[0].lai);

    /* C-factor: user-supplied or derived from canopy cover  */
    C_factor = patch[0].soil_defaults[0][0].cover_and_management_C;
    if (C_factor <= 0.0 || C_factor > 1.0)
        C_factor = 1.0 - cover_frac; /* low cover -> more erosion */

    /* Detachment adds to the surface sediment pool [kg/m2]  */
    rainsplash = KE
        * patch[0].soil_defaults[0][0].soil_erodibility_K
        * C_factor;

    /*fprintf(stderr, "rain_mm: %lf mm/day, KE: %lf MJ/m2/day,C_factor: %lf, soil_erodibility_K: %lf, rainsplash: %lf kg/m2/day\n",
        rain_mm, KE, C_factor, patch[0].soil_defaults[0][0].soil_erodibility_K, rainsplash);
    */

    patch[0].surface_sediment += rainsplash;

} /* end compute_sediment_detachment */

