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
/*    cover_frac = 1 - exp(-k * LAI), k = cover_canopy_k       */
/*      (soil def, default 0.5; larger k gives forest-like C  */
/*      under full canopy, standing in for litter cover)       */
/*    C_factor: from soil defaults if in (0,1], else           */
/*      (1-cover) * exp(-cover_litter_a * litter C)            */
/*    rainsplash [kg/m2] = KE * soil_erodibility_K * C_factor  */
/*                                                              */
/*  patch.surface_sediment is the mobilised sediment pool.      */
/*  Overland-flow routing (update_drainage_*) carries it        */
/*  downslope; whatever is left when the next day's detachment  */
/*  runs settles back (pool set to 0, its organic C/N returned  */
/*  to the soil pools).                                         */
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

    /* redeposition: detached material that overland flow did not carry away
       since yesterday's detachment settles back on the patch - the mineral
       sediment leaves the surface pool and its organic C/N returns to the soil
       pools (split by their current C), so soil organic matter is lost only
       with sediment that actually leaves the patch */
    if ((patch[0].surface_sedC > 0.0) || (patch[0].surface_sedN > 0.0)) {
        double *sc[4] = { &(patch[0].soil_cs.soil1c), &(patch[0].soil_cs.soil2c),
                          &(patch[0].soil_cs.soil3c), &(patch[0].soil_cs.soil4c) };
        double *sn[4] = { &(patch[0].soil_ns.soil1n), &(patch[0].soil_ns.soil2n),
                          &(patch[0].soil_ns.soil3n), &(patch[0].soil_ns.soil4n) };
        double totC = 0.0, w, addC[4], addN[4], sumN = 0.0, r;
        int k;
        for (k = 0; k < 4; k++) if (*sc[k] > 0.0) totC += *sc[k];
        /* C split by the pools' C; N at each pool's own N:C so pool C:N is kept
           (soil pools decompose at fixed C:N); any remaining N goes to mineral NH4,
           or, if the sediment carries less N than that, the additions are scaled down */
        for (k = 0; k < 4; k++) {
            w = (totC > 0.0) ? ((*sc[k] > 0.0) ? *sc[k] / totC : 0.0) : ((k == 3) ? 1.0 : 0.0);
            addC[k] = w * patch[0].surface_sedC;
            addN[k] = ((*sc[k] > 0.0) && (*sn[k] > 0.0)) ? addC[k] * (*sn[k]) / (*sc[k]) : 0.0;
            sumN += addN[k];
        }
        r = (sumN > patch[0].surface_sedN) ? ((sumN > 0.0) ? patch[0].surface_sedN / sumN : 0.0) : 1.0;
        for (k = 0; k < 4; k++) {
            *sc[k] += addC[k];
            *sn[k] += r * addN[k];
        }
        patch[0].soil_ns.sminn += patch[0].surface_sedN - r * sumN;   /* >= 0 */
        patch[0].surface_sedC = 0.0;
        patch[0].surface_sedN = 0.0;
    }
    patch[0].surface_sediment = 0.0;

    /* -spinmode: no erosion during the spin-up period (bare cold-start soil
       would otherwise lose its organic matter before vegetation establishes) */
    if (command_line[0].spin_active == 1)
        return;

    rain_m = zone[0].rain;
    if (rain_m <= 0.0)
        return;

    rain_mm = rain_m * 1000.0;

    /* Simplified Wischmeier & Smith kinetic energy proxy    */
    KE = rain_mm * 0.029; /* MJ/m2/day */

    /* Canopy cover fraction attenuates KE to soil surface   */
    cover_frac = 1.0 - exp(-patch[0].soil_defaults[0][0].cover_canopy_k * patch[0].lai);

    /* C-factor: user-supplied or derived from canopy cover  */
    C_factor = patch[0].soil_defaults[0][0].cover_and_management_C;
    if (C_factor <= 0.0 || C_factor > 1.0) {
        C_factor = 1.0 - cover_frac; /* low cover -> more erosion */
        /* forest-floor litter as an independent ground cover (residue-cover
           relation, Gregory 1982): x exp(-a * litter C), a = cover_litter_a
           (m2/kgC, soil def, default 0 = canopy only) */
        if (patch[0].soil_defaults[0][0].cover_litter_a > 0.0) {
            double litrc = patch[0].litter_cs.litr1c + patch[0].litter_cs.litr2c
                         + patch[0].litter_cs.litr3c + patch[0].litter_cs.litr4c;
            if (litrc > 0.0)
                C_factor *= exp(-patch[0].soil_defaults[0][0].cover_litter_a * litrc);
        }
    }

    /* Detachment adds to the surface sediment pool [kg/m2]  */
    rainsplash = KE
        * patch[0].soil_defaults[0][0].soil_erodibility_K
        * C_factor;

    /* splash under surface water: when the soil is saturated to the surface
       (same test as the ndays_sat output) rain falls on a sheet of saturation
       runoff, which absorbs drop impact (detachment ~ exp(-h/h0), negligible
       beyond ~3 drop diameters; Torri et al. 1987) */
    if ((patch[0].soil_defaults[0][0].splash_saturated_factor < 1.0) &&
        (patch[0].sat_deficit - patch[0].unsat_storage - patch[0].rz_storage <= 0.0))
        rainsplash *= patch[0].soil_defaults[0][0].splash_saturated_factor;

    /*fprintf(stderr, "rain_mm: %lf mm/day, KE: %lf MJ/m2/day,C_factor: %lf, soil_erodibility_K: %lf, rainsplash: %lf kg/m2/day\n",
        rain_mm, KE, C_factor, patch[0].soil_defaults[0][0].soil_erodibility_K, rainsplash);
    */

    patch[0].surface_sediment += rainsplash;

    /* soil organic matter detached with the sediment (sediment_OC_frac kgC/kg), taken
       from the soil pools in proportion to their C, with each pool's own N; it travels
       with surface_sediment (surface_sedC/N) and enters the stream as POC/PON */
    /* With sediment_OC_mixing_mass M (kg soil/m2, the topsoil mass the eroded
       sediment comes from) the sediment carries the fraction rainsplash/M of each
       soil pool, so its organic content follows the soil's (C/M = ~4% for 5 kg C/m2
       and M = 130) and falls toward zero on depleted soil; otherwise the fixed
       sediment_OC_frac (kgC/kg sediment) is used. */
    if (((patch[0].soil_defaults[0][0].sediment_OC_mixing_mass > 0.0) ||
         (patch[0].soil_defaults[0][0].sediment_OC_frac > 0.0)) && (rainsplash > 0.0)) {
        double *sc[4] = { &(patch[0].soil_cs.soil1c), &(patch[0].soil_cs.soil2c),
                          &(patch[0].soil_cs.soil3c), &(patch[0].soil_cs.soil4c) };
        double *sn[4] = { &(patch[0].soil_ns.soil1n), &(patch[0].soil_ns.soil2n),
                          &(patch[0].soil_ns.soil3n), &(patch[0].soil_ns.soil4n) };
        double totC = 0.0, oc, f;
        int k;
        for (k = 0; k < 4; k++) if (*sc[k] > 0.0) totC += *sc[k];
        if (patch[0].soil_defaults[0][0].sediment_OC_mixing_mass > 0.0)
            oc = totC * rainsplash / patch[0].soil_defaults[0][0].sediment_OC_mixing_mass;
        else
            oc = rainsplash * patch[0].soil_defaults[0][0].sediment_OC_frac;
        if (oc > 0.5 * totC) oc = 0.5 * totC;
        if ((totC > 0.0) && (oc > 0.0)) {
            f = oc / totC;
            for (k = 0; k < 4; k++) {
                double dc = (*sc[k] > 0.0) ? f * (*sc[k]) : 0.0;
                double dn = (*sn[k] > 0.0) ? f * (*sn[k]) : 0.0;
                *sc[k] -= dc;  *sn[k] -= dn;
                patch[0].surface_sedC += dc;
                patch[0].surface_sedN += dn;
            }
        }
    }

} /* end compute_sediment_detachment */


/*--------------------------------------------------------------*/
/*  sediment_export_scale                                       */
/*  Slope-dependent overland-flow transport capacity (soil def  */
/*  sediment_transport_slope_exp = gamma > 0):                  */
/*    Tc = sediment_transport_capacity_c * Qout^1.5             */
/*         * sin(slope)^gamma * exp(-b * surface litter C)      */
/*                                  (kg/m2 per routing step)    */
/*  (stream-power form Tc = k q^b S^g, Prosser & Rustomji 2000; */
/*  litter/mulch slows overland flow, b =                       */
/*  sediment_transport_litter_b, so fire that consumes litter   */
/*  raises Tc without fire-specific parameters).                */
/*  Returns the factor (0-1) by which the requested export      */
/*  flow_frac * surface_sediment must be scaled to stay within  */
/*  Tc; what stays behind redeposits at the next detachment.    */
/*  gamma <= 0 returns 1 (legacy behaviour).                    */
/*--------------------------------------------------------------*/
double sediment_export_scale(
    struct patch_object *patch,
    double Qout,          /* overland flow leaving this step (m) */
    double flow_frac)     /* requested fraction of surface_sediment */
{
    double gam = patch[0].soil_defaults[0][0].sediment_transport_slope_exp;
    double want, Tc;
    if ((gam <= 0.0) || (patch[0].surface_sediment <= 0.0) || (flow_frac <= 0.0))
        return(1.0);
    want = flow_frac * patch[0].surface_sediment;
    Tc = patch[0].soil_defaults[0][0].sediment_transport_capacity_c
        * pow(max(Qout, 0.0), 1.5) * pow(max(sin(patch[0].slope), 0.0), gam);
    if (patch[0].soil_defaults[0][0].sediment_transport_litter_b > 0.0) {
        double litrc = patch[0].litter_cs.litr1c + patch[0].litter_cs.litr2c
                     + patch[0].litter_cs.litr3c + patch[0].litter_cs.litr4c;
        Tc *= exp(-patch[0].soil_defaults[0][0].sediment_transport_litter_b * max(litrc, 0.0));
    }
    if (want <= Tc) return(1.0);
    return(Tc / want);
}
