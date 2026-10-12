/*--------------------------------------------------------------------------*/
/* rhessys_soil_adaptor.h -- diagnostic RHESSys -> FATES soil-layer adaptor */
/*                                                                          */
/* ED/FATES design doc 9.4d, Phase 1. RHESSys keeps its bucket hydrology    */
/* (rz_storage / unsat_storage / sat_deficit); each day this adaptor        */
/* presents it to FATES as a fixed set of soil layers:                      */
/*   - layer pore volumes from RHESSys's own porosity profile               */
/*     (compute_delta_water: porosity_0 * exp(-z / porosity_decay)),        */
/*   - three saturation zones from the RHESSys stores, using RHESSys's own  */
/*     definitions: root zone S = rz / pore(0..rz_depth), unsaturated zone  */
/*     S = unsat / pore(rz_depth..water table), saturated below the water   */
/*     table (or, with the water table inside the root zone,                */
/*     S = (rz + unsat) / sat_deficit above it); water beyond a zone's pore */
/*     space spills into the other zone above the water table,             */
/*   - suction from the same Clapp-Hornberger curve as patch.psi            */
/*     (compute_soil_water_potential).                                      */
/* Layer water sums to the RHESSys column water by construction; the        */
/* soil_water_cap check (computed independently at soil-default            */
/* construction) verifies it.                                               */
/*--------------------------------------------------------------------------*/
#ifndef _RHESSYS_SOIL_ADAPTOR_H_
#define _RHESSYS_SOIL_ADAPTOR_H_

#define SA_NLEV_MAX 10

struct soil_adaptor_params {
	double	porosity_0;		/* (DIM) surface porosity */
	double	porosity_decay;		/* (m) e-folding depth of porosity (>= 999 = uniform) */
	double	soil_depth;		/* (m) */
	double	soil_water_cap;		/* (m water) pore volume of the column, from soil defaults */
	int	theta_psi_curve;	/* 1 Clapp-Hornberger, 2 van Genuchten, 3 power law */
	double	psi_air_entry;		/* (m) */
	double	pore_size_index;	/* (DIM) Clapp-Hornberger 1/b */
	double	p3, p4;			/* curve 3 coefficients */
	double	psi_max;		/* (m) wet-end cap used for patch.psi (MPa = -psi_max in that call) */
};

struct soil_adaptor_state {
	double	rz_depth;		/* (m) current root-zone depth */
	double	rz_storage;		/* (m water) */
	double	unsat_storage;		/* (m water) */
	double	sat_deficit;		/* (m water); < 0 = ponded */
	double	Tsoil;			/* (deg C) */
};

struct soil_adaptor_out {
	double	zwt;			/* (m) water-table depth, recomputed from sat_deficit */
	double	S_rz, S_unsat;		/* zone saturations used, after spilling (-1 = zone absent) */
	double	S_patch_equiv;		/* patch_daily_I's patch.S = (rz+unsat)/sat_deficit, rebuilt from the layers */
	double	psi_patch_equiv;	/* (MPa) patch.psi's curve + clipping at S_patch_equiv */
	double	pore[SA_NLEV_MAX];	/* (m water) layer pore volume */
	double	water[SA_NLEV_MAX];	/* (m water) layer water */
	double	theta[SA_NLEV_MAX];	/* (m3/m3) layer water / dz (unclamped) */
	double	watsat[SA_NLEV_MAX];	/* (m3/m3) layer-mean porosity */
	double	S[SA_NLEV_MAX];		/* (DIM) layer saturation water/pore (unclamped) */
	double	smp_mm[SA_NLEV_MAX];	/* (mm, < 0) suction for FATES: soil curve only */
	double	psi_patch[SA_NLEV_MAX];	/* (MPa) same curve with patch.psi's clipping (diagnostic) */
	double	col_water;		/* (m water) sum of layer water */
	double	ref_water;		/* (m water) RHESSys column water: max(cap - max(sat_def,0), 0) + rz + unsat */
};

/* Fixed layer interfaces zi[0..nlev] (zi[0] = 0, zi[nlev] = soil_depth): the
   first layer is the root zone at construction (when 0 < rz_depth0 < depth),
   the rest split the remainder evenly. Returns 0 on success. */
int soil_adaptor_geometry(int nlev, double rz_depth0, double soil_depth, double *zi);

/* Daily layered state from the RHESSys stores. */
void soil_adaptor_profile(const struct soil_adaptor_params *par,
			  const struct soil_adaptor_state *st,
			  int nlev, const double *zi,
			  struct soil_adaptor_out *out);

/* RHESSys convention used by compute_soil_water_potential: 1 m H2O = 0.01 MPa,
   so suction in mm = psi(MPa) * 1e5. */
#define SA_MM_PER_MPA 1.0e5

#endif /* _RHESSYS_SOIL_ADAPTOR_H_ */
