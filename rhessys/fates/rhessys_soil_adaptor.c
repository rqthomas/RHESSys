/*--------------------------------------------------------------------------*/
/* rhessys_soil_adaptor.c -- see include/rhessys_soil_adaptor.h             */
/* (ED/FATES design doc 9.4d, Phase 1). Diagnostic only: reads the RHESSys  */
/* stores, never changes them.                                              */
/*--------------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>

#include "rhessys.h"
#include "rhessys_soil_adaptor.h"

double	compute_delta_water(int, double, double, double, double, double);
double	compute_z_final(int, double, double, double, double, double);
double	compute_soil_water_potential(int, int, double, double, double, double,
		double, double, double, double, double, double);

/* driest suction handed to FATES (MPa, RHESSys units); CLM's smpmin is -1e8 mm */
#define SA_PSI_FLOOR	-1000.0
/* smallest saturation used for suction, as CLM (s_node >= 0.01) */
#define SA_S_MIN	0.01
/* thinnest layer the root-zone split may produce (m); else split evenly */
#define SA_MIN_DZ	0.05

/* pore volume (m water) between depths top < bot */
static double pore_between(const struct soil_adaptor_params *par, double top, double bot)
{
	if (bot <= top)
		return(0.0);
	return(compute_delta_water(0, par[0].porosity_0, par[0].porosity_decay,
		par[0].soil_depth, bot, top));
}

int soil_adaptor_geometry(int nlev, double rz_depth0, double soil_depth, double *zi)
{
	int	j, nrest;
	double	top;

	if ((nlev < 1) || (nlev > SA_NLEV_MAX) || (soil_depth <= 0.0))
		return(1);
	zi[0] = 0.0;
	/* root-zone split only if no layer comes out thinner than
	   max(SA_MIN_DZ, 10 % of an even layer) -- e.g. not when the root zone
	   fills the soil (CCR: rootzone.depth ~ soil_depth) */
	top = max(SA_MIN_DZ, 0.1 * soil_depth / nlev);
	if ((nlev >= 2) && (rz_depth0 >= top) && ((soil_depth - rz_depth0) / (nlev - 1) >= top)) {
		zi[1] = rz_depth0;
		top = rz_depth0;
		nrest = nlev - 1;
		for (j = 1; j <= nrest; j++)
			zi[1 + j] = top + (soil_depth - top) * (double) j / (double) nrest;
	}
	else {
		for (j = 1; j <= nlev; j++)
			zi[j] = soil_depth * (double) j / (double) nlev;
	}
	zi[nlev] = soil_depth;
	return(0);
}

void soil_adaptor_profile(const struct soil_adaptor_params *par,
			  const struct soil_adaptor_state *st,
			  int nlev, const double *zi,
			  struct soil_adaptor_out *out)
{
	int	j, k, nzone;
	double	depth = par[0].soil_depth;
	double	zwt, zrz, pz, ztop[3], zbot[3], Sz[3];
	double	dz, s_node, psi, pA, pB, wA, wB, move;

	/*--------------------------------------------------------------*/
	/* water table from sat_deficit, as patch_daily_F does		*/
	/*--------------------------------------------------------------*/
	if (st[0].sat_deficit <= 0.0)
		zwt = 0.0;
	else
		zwt = compute_z_final(0, par[0].porosity_0, par[0].porosity_decay,
			depth, 0.0, -1.0 * st[0].sat_deficit);
	zwt = min(max(zwt, 0.0), depth);
	zrz = min(max(st[0].rz_depth, 0.0), depth);
	out[0].zwt = zwt;
	out[0].S_rz = -1.0;
	out[0].S_unsat = -1.0;

	/*--------------------------------------------------------------*/
	/* saturation zones (RHESSys definitions)			*/
	/*--------------------------------------------------------------*/
	nzone = 0;
	if (zwt > zrz) {
		/* root zone, then unsaturated zone down to the water table */
		pA = pore_between(par, 0.0, zrz);
		pB = pore_between(par, zrz, zwt);
		wA = st[0].rz_storage;
		wB = st[0].unsat_storage;
		/* a zone cannot hold more than its pore space: spill the excess into
		   the other zone above the water table (column water unchanged) */
		if (wB > pB) { wA += wB - pB; wB = pB; }
		if (wA > pA) { move = min(wA - pA, pB - wB); wA -= move; wB += move; }
		if (zrz > 0.0) {
			ztop[nzone] = 0.0; zbot[nzone] = zrz;
			Sz[nzone] = (pA > 0.0) ? wA / pA : 0.0;
			out[0].S_rz = Sz[nzone];
			nzone++;
		}
		else
			wB += wA;		/* no root zone: all of it is unsaturated zone */
		ztop[nzone] = zrz; zbot[nzone] = zwt;
		Sz[nzone] = (pB > 0.0) ? wB / pB : 0.0;
		out[0].S_unsat = Sz[nzone];
		nzone++;
	}
	else if (zwt > 0.0) {
		/* water table inside the root zone: RHESSys moves unsat into rz */
		pA = pore_between(par, 0.0, zwt);
		ztop[nzone] = 0.0; zbot[nzone] = zwt;
		Sz[nzone] = (pA > 0.0) ? (st[0].rz_storage + st[0].unsat_storage) / pA : 0.0;
		out[0].S_rz = Sz[nzone];
		nzone++;
	}
	/* saturated below the water table */
	ztop[nzone] = zwt; zbot[nzone] = depth; Sz[nzone] = 1.0;
	nzone++;

	/*--------------------------------------------------------------*/
	/* per-layer water, porosity, saturation, suction		*/
	/*--------------------------------------------------------------*/
	out[0].col_water = 0.0;
	for (j = 0; j < nlev; j++) {
		dz = zi[j+1] - zi[j];
		out[0].pore[j] = pore_between(par, zi[j], zi[j+1]);
		out[0].water[j] = 0.0;
		for (k = 0; k < nzone; k++)
			out[0].water[j] += Sz[k] * pore_between(par,
				max(zi[j], ztop[k]), min(zi[j+1], zbot[k]));
		out[0].col_water += out[0].water[j];
		out[0].theta[j] = (dz > 0.0) ? out[0].water[j] / dz : 0.0;
		out[0].watsat[j] = (dz > 0.0) ? out[0].pore[j] / dz : 0.0;
		out[0].S[j] = (out[0].pore[j] > 0.0) ? out[0].water[j] / out[0].pore[j] : 0.0;

		s_node = max(out[0].S[j], SA_S_MIN);
		/* soil curve only: wet cap = value at saturation, dry floor = SA_PSI_FLOOR */
		psi = compute_soil_water_potential(0, par[0].theta_psi_curve, st[0].Tsoil,
			-0.01 * par[0].psi_air_entry, SA_PSI_FLOOR,
			par[0].psi_air_entry, par[0].pore_size_index, par[0].p3, par[0].p4,
			par[0].porosity_0, par[0].porosity_decay, s_node);
		out[0].smp_mm[j] = psi * SA_MM_PER_MPA;
		/* the same call patch_daily_I makes for patch.psi (clipped), for comparison */
		out[0].psi_patch[j] = compute_soil_water_potential(0, par[0].theta_psi_curve,
			st[0].Tsoil, -1.0 * par[0].psi_max, -10.0,
			par[0].psi_air_entry, par[0].pore_size_index, par[0].p3, par[0].p4,
			par[0].porosity_0, par[0].porosity_decay, out[0].S[j]);
	}

	/*--------------------------------------------------------------*/
	/* patch.psi from the layers: patch_daily_I sets patch.S =	*/
	/* (rz + unsat) / sat_deficit (mean saturation above the water	*/
	/* table, 1 if ponded) and evaluates the curve there. Rebuild	*/
	/* that S from the layer water (minus the saturated part) so it	*/
	/* checks the layering against RHESSys's own psi.		*/
	/*--------------------------------------------------------------*/
	if (st[0].sat_deficit < ZERO)
		out[0].S_patch_equiv = 1.0;
	else
		out[0].S_patch_equiv = (out[0].col_water - pore_between(par, zwt, depth))
			/ st[0].sat_deficit;
	out[0].psi_patch_equiv = compute_soil_water_potential(0, par[0].theta_psi_curve,
		st[0].Tsoil, -1.0 * par[0].psi_max, -10.0,
		par[0].psi_air_entry, par[0].pore_size_index, par[0].p3, par[0].p4,
		par[0].porosity_0, par[0].porosity_decay, out[0].S_patch_equiv);

	/*--------------------------------------------------------------*/
	/* RHESSys column water: water above the water table plus the	*/
	/* saturated pore volume below it. Ponded water (sat_deficit < 0)	*/
	/* is not soil water; RHESSys lets sat_deficit exceed		*/
	/* soil_water_cap (water table below the soil bottom), and then	*/
	/* there is no saturated soil, only rz + unsat.			*/
	/*--------------------------------------------------------------*/
	out[0].ref_water = max(par[0].soil_water_cap - max(st[0].sat_deficit, 0.0), 0.0)
		+ st[0].rz_storage + st[0].unsat_storage;
}
