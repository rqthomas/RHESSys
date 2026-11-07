/******************************************************************************
 * rhessys_veg.h -- pluggable vegetation-backend interface for RHESSys        *
 *                                                                            *
 * RHESSys is the HOST; the vegetation module is a backend selected at        *
 * runtime, mirroring the GLM<->AED coupling (GLM3/src/glm_wqual.h):          *
 *                                                                            *
 *   backend "legacy" : the built-in big-leaf code (default; runs inline in   *
 *                       the daily cycle) -- selected when demog_flag == 0.    *
 *   backend "fates"  : the coupled FATES Fortran library (demog_flag > 0).   *
 *                                                                            *
 * Only the BRANCH (which backend runs) is gated on demog_flag; shared        *
 * helpers may be reused/refactored by both paths (design doc 6.2 opt i,9.2). *
 * Off-mode (legacy) must stay bit-identical (design doc 6.5).                 *
 *                                                                            *
 * Layers: the backend entry points (p_veg_*, Fortran BIND(C) for FATES) only *
 * see flat interoperable structs (struct veg_bc / struct veg_time). The C    *
 * host driver (fates/rhessys_veg_host.c, veg_host_*) walks the RHESSys world,*
 * runs the soil adaptor and fills those structs. This header is dual-use:    *
 * included from C and (for the BIND(C) names) from rhessys_fates.F90; the    *
 * Fortran mirror of the structs lives there and must match field order.      *
 ******************************************************************************/
#ifndef _RHESSYS_VEG_H_
#define _RHESSYS_VEG_H_

/*--------------------------------------------------------------------------*/
/* Shared C <-> Fortran BIND(C) entry-point names (used by BOTH sides).     */
/*--------------------------------------------------------------------------*/
#define _VEG_INIT      "rhessys_veg_init"
#define _VEG_INIT_SITE "rhessys_veg_init_site"
#define _VEG_COLDSTART "rhessys_veg_coldstart"
#define _VEG_SET_TIME  "rhessys_veg_set_time"
#define _VEG_SET_BC    "rhessys_veg_set_bc"
#define _VEG_DO        "rhessys_veg_do"
#define _VEG_DO_ALL    "rhessys_veg_do_all"
#define _VEG_CLEAN     "rhessys_veg_clean"
#define _VEG_DUMP      "rhessys_veg_dump_cohorts"
#define _VEG_RESTART_WRITE "rhessys_veg_restart_write"
#define _VEG_RESTART_READ  "rhessys_veg_restart_read"

/* struct sizes shared with the Fortran mirror */
#define VEG_NLEV_MAX 10      /* max soil layers handed to the backend */
#define VEG_NSWB     2       /* shortwave bands: 1 = visible, 2 = near-IR */
#define VEG_NSUB_MAX 12      /* max daylight flux sub-steps (-demog_diurnal) */

#ifndef _FORTRAN_SOURCE_
/*==========================================================================*/
/* C side                                                                   */
/*==========================================================================*/
#include <stddef.h>

/*--------------------------------------------------------------------------*/
/* Daily boundary conditions for one site (= one RHESSys patch), in the     */
/* backend's (FATES/CLM) units. Ints first, then doubles: interoperable     */
/* with the BIND(C) type veg_bc_t in rhessys_fates.F90 -- keep the order.   */
/*--------------------------------------------------------------------------*/
struct veg_bc {
	int	nlevsoil;			/* soil layers used (<= VEG_NLEV_MAX) */
	int	max_rooting_depth_index;	/* deepest layer roots may use */
	/* soil column (fixed at site init) */
	double	zi[VEG_NLEV_MAX + 1];		/* (m) interface depths, zi[0] = 0 */
	double	dz[VEG_NLEV_MAX];		/* (m) layer thickness */
	double	z[VEG_NLEV_MAX];		/* (m) layer node (mid) depth */
	double	watsat[VEG_NLEV_MAX];		/* (m3/m3) layer-mean porosity */
	/* soil state (soil adaptor, daily) */
	double	h2o_liqvol[VEG_NLEV_MAX];	/* (m3/m3) liquid water, <= watsat */
	double	eff_porosity[VEG_NLEV_MAX];	/* (m3/m3) porosity minus ice */
	double	smp[VEG_NLEV_MAX];		/* (mm, < 0) soil matric potential */
	double	tempk[VEG_NLEV_MAX];		/* (K) soil temperature */
	/* atmosphere (daylight means) */
	double	forc_pbot;			/* (Pa) */
	double	tgcm;				/* (K) air temperature */
	double	t_veg;				/* (K) vegetation temperature */
	double	eair;				/* (Pa) canopy-air vapour pressure */
	double	esat_tv;			/* (Pa) saturation vapour pressure at t_veg */
	double	oair;				/* (Pa) O2 partial pressure */
	double	cair;				/* (Pa) CO2 partial pressure */
	double	rb;				/* (s/m) leaf boundary-layer resistance */
	double	dayl_factor;			/* (0-1) daylength scaling */
	double	coszen;				/* (0-1) daylight-mean cosine solar zenith */
	double	solad[VEG_NSWB];		/* (W/m2) direct beam, vis / nir */
	double	solai[VEG_NSWB];		/* (W/m2) diffuse, vis / nir */
	double	albgr_dir[VEG_NSWB];		/* (0-1) ground albedo, direct */
	double	albgr_dif[VEG_NSWB];		/* (0-1) ground albedo, diffuse */
	double	snow_depth;			/* (m) */
	double	frac_sno_eff;			/* (0-1) ground snow cover */
	double	fcansno;			/* (0-1) canopy snow cover */
	double	precip24;			/* (mm/s) */
	double	relhumid24;			/* (%) */
	double	wind24;				/* (m/s) */
	/* day/night sub-steps of the daily flux step (Phase 2) */
	double	dayl;				/* (s) daylight length */
	double	t_night;			/* (K) night air/vegetation temperature */
	double	esat_night;			/* (Pa) saturation vapour pressure at t_night */
	double	eair_night;			/* (Pa) night canopy-air vapour pressure (<= esat_night) */
	double	t_avg;				/* (K) daily mean temperature (running means) */
	/* decomposition scalars for FATES litter fragmentation (Phase 3) */
	double	w_scalar;			/* (0-1) RHESSys decomposition moisture scalar */
	double	t_scalar;			/* (DIM) RHESSys decomposition temperature scalar (1 at 25 C) */
	/* mode 3 (Phase 4): plant N uptake granted by RHESSys's soil competition
	   (resolve_sminn_competition) on the previous day, handed to FATES before
	   its dynamics (as CLM integrates uptake over the day) */
	double	n_uptake_nh4;			/* (gN/m2/day) */
	double	n_uptake_no3;			/* (gN/m2/day) */
	/* Phase 6: today's disturbance from TEC events (mode 3).
	   dist_kind 0 = none, 1 = fire (uniform severity), 2 = harvest (clear-cut),
	   3 = partial removal -> proportional thinning (thin_frac) */
	double	dist_kind;
	double	fire_mort[2];			/* (0-1) cohort mortality, [0] evergreen, [1] deciduous PFTs */
	double	crown_burn[2];			/* (0-1) share of crowns burned (consumed) */
	double	fuel_burn_leaf;			/* (0-1) FATES dead-leaf litter consumed */
	double	fuel_burn_twig;			/* (0-1) FATES twig + small-branch CWD consumed */
	double	thin_frac[2];			/* (0-1) stems removed, [0] evergreen, [1] deciduous PFTs */
	/* dist_kind 4 = WMFire (execute_firespread_event, -firespread): burned fraction of the
	   patch, mean burn intensity of its burned cells, RHESSys + FATES litter C consumed per
	   m2 burned (kgC/m2), and WMFire's severity parameters (stratum / soil defaults) */
	double	wm_frac, wm_burn, wm_litter_consumed;
	double	wm_under_mort, wm_consumption, wm_over_k1, wm_over_k2, wm_under_h, wm_over_h;
	/* -demog_nmine: soil organic N the host granted yesterday to C-paid acquisition (gN/m2/day) */
	double	n_uptake_mine;
	/* -demog_nforage: share of yesterday's granted NH4 + NO3 won by AM foraging (to forager cohorts) */
	double	n_uptake_am_frac;
	/* -demog_nfun: yesterday's root-zone mineral N available to plants (gN/m2), the N term of the
	   FUN-style uptake cost */
	double	n_avail;
	/* -demog_diurnal: daylight flux sub-steps (<= 1: one daylight-mean step) and the
	   diurnal shape inputs: cos(zenith) = cz_a + cz_b cos(hour angle), sunset hour angle
	   hhalf (rad); daily tmax, tmin (K) for the temperature cycle */
	double	n_diurnal;
	double	cz_a, cz_b, hhalf;
	double	tmax, tmin;
};

/*--------------------------------------------------------------------------*/
/* Daily flux-step results for one site (backend -> host). Area-weighted    */
/* over the site's FATES patches. Same interop rules as struct veg_bc.      */
/*--------------------------------------------------------------------------*/
struct veg_out {
	int	npatches;			/* FATES patches in the site */
	int	ncohorts;			/* FATES cohorts in the site */
	double	tlai, elai, tsai;		/* (m2/m2) total / exposed leaf, stem area index */
	double	htop;				/* (m) canopy top height (area-weighted) */
	double	laisun, laisha;			/* (m2/m2) daylight sunlit / shaded leaf area */
	double	rssun, rssha;			/* (s/m) daylight stomatal resistance per leaf area */
	double	gc_canopy;			/* (m/s) daylight canopy conductance:	*/
						/*   laisun/(rb+rssun) + laisha/(rb+rssha) */
	double	swabs_canopy;			/* (W/m2) daylight-mean shortwave absorbed by the canopy */
	double	btran;				/* (0-1) soil water stress factor */
	double	gpp_day, gpp_night;		/* (kgC/m2) GPP in each sub-step */
	double	resp_m;				/* (kgC/m2/day) maintenance respiration, day + night */
	double	biomass_c;			/* (kgC/m2) live plant carbon (after today's dynamics) */
	double	litter_c;			/* (kgC/m2) FATES litter + CWD carbon */
	double	lai_evg_frac;			/* (0-1) evergreen share of LAI (leaf C x SLA by leaf habit) */
	/* Phase 4 (mode 3) */
	double	ncap;				/* (gN/m2/day) RD uptake capacity: sum fine-root C x (Vmax NH4 + NO3) */
	double	n_uptake_recv;			/* (gN/m2/day) N uptake FATES received today (cohorts) */
	double	lit_c_lab, lit_c_cel, lit_c_lig;	/* (gC/m2/day) FATES litter to host: labile, cellulose, lignin */
	double	lit_n_lab, lit_n_cel, lit_n_lig;	/* (gN/m2/day) */
	double	c_stock, n_stock;		/* (kg/m2) FATES total C, N (plants + FATES litter + seeds) */
	double	mb_n_uptake, mb_n_frag, mb_n_err;	/* (kgN/m2/day) FATES's own N mass balance record:
						   net root uptake, fragmentation out, error */
	/* litter handoff split by where it enters RHESSys (g/m2/day), mirroring
	   FATES FluxIntoLitterPools: surface = leaf fines + above-ground CWD + seed
	   decay, below ground = fine roots + below-ground CWD. lit_* above = sums. */
	double	ls_c_lab, ls_c_cel, ls_c_lig, ls_n_lab, ls_n_cel, ls_n_lig;	/* surface */
	double	lb_c_lab, lb_c_cel, lb_c_lig, lb_n_lab, lb_n_cel, lb_n_lig;	/* below ground */
	double	lit_split_err;			/* (g/m2/day) |split sum - FATES litt_flux|, C+N */
	/* FATES litter still unfragmented, seen by RHESSys as forest floor (kgC/m2) */
	double	leaf_litter_c;			/* above-ground dead leaves (leaf fines) */
	double	ag_cwd_c;			/* above-ground coarse woody debris */
	/* root efflux (g/m2/day): FATES books it out of the plants and puts it in
	   root_fines_frag, which CWDOut then overwrites before FluxIntoLitterPools
	   (sci.1.93.2), so it never reaches the host flux -> delivered by the host */
	double	efflux_c, efflux_n;
	/* Phase 6: FATES disturbance losses today (g/m2/day), from its mass balance */
	double	fire_c, fire_n;			/* burned to the atmosphere */
	double	harvest_c, harvest_n;		/* exported as wood products */
	/* (kgN/m2/day) FATES's other N terms in its own balance: seed in - out,
	   generic in - out, patch_resize_err (mass FATES gains/loses when it
	   resizes patches, e.g. at logging), - herbivory, - land-use products */
	double	mb_n_other, mb_n_resize;
	/* (gN/m2/day) symbiotic N fixation by FATES cohorts today (fates_cnp_nfix1 > 0): an N input
	   to FATES that is not taken from RHESSys soil mineral N */
	double	n_fix;
	/* -demog_nmine (gN/m2/day): today's demand = excess C of mining PFTs / cost; mined N the
	   cohorts received today */
	double	n_mine_demand, n_mine_recv;
	/* -demog_nforage (gN/m2/day): extra mineral-N claim = excess C of foraging PFTs / cost; N the
	   forager cohorts received today through it */
	double	n_am_demand, n_am_recv;
	/* -demog_diurnal: number of daylight sub-steps run today (0 or 1: the single
	   daylight-mean step) and per sub-step canopy conductance (m/s), absorbed
	   shortwave (W/m2), air temperature (C) and VPD (Pa) for the host transpiration */
	double	n_sub;
	double	sub_gc[VEG_NSUB_MAX], sub_swabs[VEG_NSUB_MAX], sub_t[VEG_NSUB_MAX], sub_vpd[VEG_NSUB_MAX];
};

/* Model time (FATES SetFatesTime); day_of_year is 1-based (1 = 1 Jan) as in CLM. */
struct veg_time {
	int	year, month, day;
	int	day_of_year;
	int	days_per_year;
	int	current_date;			/* yyyymmdd */
	int	reference_date;			/* yyyymmdd (run start) */
	double	model_day;			/* days since run start */
};

/* Backend entry-point signatures (host -> vegetation backend). Scalars by  */
/* pointer (GLM convention) so the Fortran dummies need no VALUE attribute. */
typedef void (*veg_init_t)(const char *paramfile, size_t *len, const char *invfile,
			   size_t *invlen, const char *logfile, size_t *loglen, int *nlevsoil,
			   int *nlevdecomp, int *demog_flag, int *nsites, int *is_restart);
/* FATES restart vectors <-> binary file (Phase 5) */
typedef void (*veg_restart_write_t)(const char *file, size_t *flen);
typedef void (*veg_restart_read_t)(const char *file, size_t *flen, struct veg_out *out_all);
/* forest management (rhessys_fates_mgmt.F90): driver file + RHESSys patch/zone/hill
   ID of each site; empty file = TEC thinning only */
/* per-patch inventories: inventory index (0-based, as listed in the generated FATES
   control file) of every site, before the cold start */
typedef void (*veg_inv_assign_t)(int *nsites, int *idx);
/* plant N demand mode (bit mask, see rhessys_fates.F90 ndemand_mode) and the -demog_ncap need multiplier */
typedef void (*veg_set_ndemand_t)(int *mode, double *ncap_mult);
/* FATES stomatal conductance model, set before veg_init: 1 = Ball-Berry, 2 = Medlyn */
typedef void (*veg_set_stomatal_t)(int *model);
/* FATES regeneration model, set before veg_init: 1 = default, 2 = TRS with seedling dynamics */
typedef void (*veg_set_regen_t)(int *model);
/* FATES cohort age tracking, set before veg_init: 0 = off, 1 = on */
typedef void (*veg_set_coage_t)(int *on);
/* FATES day-length factor on photosynthetic capacity, set before veg_init: 0 = off, 1 = on */
typedef void (*veg_set_dayl_t)(int *on);
/* FATES seed dispersal between sites: on/off before veg_init; site geometry (centroid x, y in m, area
   in m2, by FATES site index) after it */
typedef void (*veg_set_disperse_t)(int *on);
typedef void (*veg_set_disperse_geom_t)(int *nsites, double *x, double *y, double *area);
/* upper bound on FATES's dynamic leaf-to-fineroot ratio (<= 0: none) */
typedef void (*veg_set_l2fr_max_t)(double *l2fr_max);
typedef void (*veg_set_min_n_t)(double *min_n);
typedef void (*veg_set_nvcmax_t)(double *floor);
typedef void (*veg_set_nmine_t)(double *cost, int *mask);
typedef void (*veg_set_nforage_t)(double *cost, int *mask);
/* -demog_nfun: k[4] = kN_AM, kC_AM, kN_ECM, kC_ECM; charge 0/1; ECM PFT bit mask */
typedef void (*veg_set_nfun_t)(double *k, int *charge, int *mask);
/* -demog_nfun: today's site terms (g/m2/day): N to AM cohorts, N to ECM cohorts, sum N x cost AM,
   sum N x cost ECM, C charged to NPP */
typedef void (*veg_nfun_diag_t)(int *site, double *v);
/* -demog_mort_disturb: per-site fates_mort_disturb_frac override (FATES host patch; < 0 = parameter) */
typedef void (*veg_set_mort_disturb_t)(int *site, double *frac);
/* -demog_browse: rate (1/yr), height (m), PFT bit mask (FATES host patch) */
typedef void (*veg_set_browse_t)(double *rate, double *height, int *mask);
/* -demog_suppress: critical dbh growth (cm/yr), max dbh (cm), n PFTs, rate by PFT (1/yr), min canopy layer (FATES host patch) */
typedef void (*veg_set_suppress_t)(double *gcrit, double *dbhmax, int *npft, double *rate, int *minlayer);
/* site potential growth N need (kgN/m2/day) */
typedef void (*veg_site_nneed_t)(int *site, double *need);
/* today's FATES litter by source (g/m2/day): C leaf, seed, AG CWD, fine root, BG CWD; then N likewise */
typedef void (*veg_litter_src_t)(int *site, double *v);
/* stream-side sites: share of fresh leaf litter that falls into the channel (set daily) */
typedef void (*veg_set_leaf_export_t)(int *site, double *frac);
/* today's leaf litter taken out of FATES for the stream (g/m2/day): C, N, labile C, labile N */
typedef void (*veg_leaf_export_t)(int *site, double *v);
/* -demog_litter_direct on/off (FATES releases all fine litter on the day it falls) */
typedef void (*veg_set_litter_direct_t)(int *flag);
/* -demog_litter_direct: today's litter (g/m2/day) v[class + 3 source + 15 element],
   class 0 labile 1 cellulose 2 lignin; source 0 leaf fines 1 seeds 2 AG CWD 3 fine roots 4 BG CWD;
   element 0 C 1 N */
typedef void (*veg_litter_dcmp_t)(int *site, double *v);
/* FATES CWD fragmentation split: v[0] cellulose, v[1] lignin fraction */
typedef void (*veg_cwd_frac_t)(double *v);
/* FATES live leaf C of a site (kgC/m2), WMFire fuel */
typedef void (*veg_leaf_c_t)(int *site, double *leafc);
typedef void (*veg_mgmt_init_t)(const char *file, size_t *flen, int *nsites, int *pid,
				int *zid, int *hid);
typedef void (*veg_init_site_t)(int *site, int *patch_id, double *lat, struct veg_bc *bc);
typedef void (*veg_coldstart_t)(struct veg_out *out_all);	/* fills out_all[0..nsites-1] */
typedef void (*veg_set_time_t)(struct veg_time *t);
typedef void (*veg_set_bc_t)(int *site, struct veg_bc *bc, int *n_bad);
typedef void (*veg_do_t)(int *site, struct veg_bc *bc, struct veg_out *out);	/* daily flux step */
/* all sites at once (OpenMP over sites inside the backend) */
typedef void (*veg_do_all_t)(int *nsites, struct veg_bc *bcs, struct veg_out *outs);
typedef void (*veg_clean_t)(void);
/* append one site's cohorts (date-stamped) to a CSV; header when *first */
typedef void (*veg_dump_t)(int *site, int *patch_id, int *ymd, const char *file,
			   size_t *flen, int *first);

extern veg_init_t      p_veg_init;
extern veg_init_site_t p_veg_init_site;
extern veg_coldstart_t p_veg_coldstart;
extern veg_set_time_t  p_veg_set_time;
extern veg_set_bc_t    p_veg_set_bc;
extern veg_do_t        p_veg_do;
extern veg_do_all_t    p_veg_do_all;
extern veg_clean_t     p_veg_clean;
extern veg_dump_t      p_veg_dump;
extern veg_restart_write_t p_veg_restart_write;
extern veg_restart_read_t  p_veg_restart_read;
extern veg_mgmt_init_t     p_veg_mgmt_init;
extern veg_leaf_c_t        p_veg_leaf_c;
extern veg_inv_assign_t    p_veg_inv_assign;
extern veg_set_ndemand_t   p_veg_set_ndemand;
extern veg_set_stomatal_t  p_veg_set_stomatal;
extern veg_set_regen_t     p_veg_set_regen;
extern veg_set_coage_t     p_veg_set_coage;
extern veg_set_dayl_t      p_veg_set_dayl;
extern veg_set_disperse_t  p_veg_set_disperse;
extern veg_set_disperse_geom_t p_veg_set_disperse_geom;
extern veg_set_l2fr_max_t  p_veg_set_l2fr_max;
extern veg_set_min_n_t     p_veg_set_min_n;
extern veg_set_nvcmax_t    p_veg_set_nvcmax;
extern veg_set_nmine_t     p_veg_set_nmine;
extern veg_set_nforage_t   p_veg_set_nforage;
extern veg_set_nfun_t      p_veg_set_nfun;
extern veg_nfun_diag_t     p_veg_nfun_diag;
extern veg_set_mort_disturb_t p_veg_set_mort_disturb;
extern veg_set_browse_t   p_veg_set_browse;
extern veg_set_suppress_t p_veg_set_suppress;
extern veg_site_nneed_t    p_veg_site_nneed;
extern veg_litter_src_t    p_veg_litter_src;
extern veg_set_leaf_export_t p_veg_set_leaf_export;
extern veg_leaf_export_t   p_veg_leaf_export;
extern veg_set_litter_direct_t p_veg_set_litter_direct;
extern veg_litter_dcmp_t   p_veg_litter_dcmp;
extern veg_cwd_frac_t      p_veg_cwd_frac;

#define veg_init       (*p_veg_init)
#define veg_init_site  (*p_veg_init_site)
#define veg_coldstart  (*p_veg_coldstart)
#define veg_set_time   (*p_veg_set_time)
#define veg_set_bc     (*p_veg_set_bc)
#define veg_do         (*p_veg_do)
#define veg_do_all     (*p_veg_do_all)
#define veg_clean      (*p_veg_clean)
#define veg_dump       (*p_veg_dump)
#define veg_restart_write (*p_veg_restart_write)
#define veg_restart_read  (*p_veg_restart_read)
#define veg_mgmt_init     (*p_veg_mgmt_init)
#define veg_leaf_c        (*p_veg_leaf_c)
#define veg_inv_assign    (*p_veg_inv_assign)
#define veg_set_ndemand   (*p_veg_set_ndemand)
#define veg_set_stomatal  (*p_veg_set_stomatal)
#define veg_set_regen     (*p_veg_set_regen)
#define veg_set_coage     (*p_veg_set_coage)
#define veg_set_dayl      (*p_veg_set_dayl)
#define veg_set_disperse  (*p_veg_set_disperse)
#define veg_set_disperse_geom (*p_veg_set_disperse_geom)
#define veg_set_l2fr_max  (*p_veg_set_l2fr_max)
#define veg_set_min_n     (*p_veg_set_min_n)
#define veg_set_nvcmax    (*p_veg_set_nvcmax)
#define veg_set_nmine     (*p_veg_set_nmine)
#define veg_set_nforage   (*p_veg_set_nforage)
#define veg_set_nfun      (*p_veg_set_nfun)
#define veg_nfun_diag     (*p_veg_nfun_diag)
#define veg_set_mort_disturb (*p_veg_set_mort_disturb)
#define veg_set_browse    (*p_veg_set_browse)
#define veg_set_suppress  (*p_veg_set_suppress)
#define veg_site_nneed    (*p_veg_site_nneed)
#define veg_litter_src    (*p_veg_litter_src)
#define veg_set_leaf_export (*p_veg_set_leaf_export)
#define veg_leaf_export   (*p_veg_leaf_export)
#define veg_set_litter_direct (*p_veg_set_litter_direct)
#define veg_litter_dcmp   (*p_veg_litter_dcmp)
#define veg_cwd_frac      (*p_veg_cwd_frac)

/* Select the backend: "legacy" (default) or "fates". 0 on success. */
int prime_veg(const char *which);

/*--------------------------------------------------------------------------*/
/* C host driver (fates/rhessys_veg_host.c). Called only when demog_flag.  */
/*--------------------------------------------------------------------------*/
struct world_object;
struct basin_object;
struct zone_object;
struct patch_object;
struct command_line_object;
struct date;

void veg_host_init(struct world_object *, struct command_line_object *);
void veg_host_daily(struct world_object *, struct basin_object *,
		    struct zone_object *, struct patch_object *,
		    struct command_line_object *, struct date);
void veg_host_transpiration(struct zone_object *, struct patch_object *,
		    struct command_line_object *);
void veg_host_daily_end(struct zone_object *, struct patch_object *,
		    struct command_line_object *, struct date);
void veg_host_ndemand(struct patch_object *, struct command_line_object *);
/* FATES unfragmented above-ground leaf litter on this patch (kgC/m2): RHESSys's
   forest floor = its litter pools + this (sediment cover/transport, litter
   interception). 0 unless -demog mode 3 and the patch is a FATES site. */
double veg_host_surface_litter_c(struct patch_object *);
void veg_host_litter_interception(struct patch_object *);
/* Phase 6: a TEC thin event (redefine_world_thin_*) on a stratum of a FATES
   site in mode 3 is recorded for FATES instead of applied to the legacy
   stratum. Returns 1 if recorded (caller skips the legacy update). */
struct canopy_strata_object;
struct mortality_struct;
int veg_host_strata_event(struct patch_object *, struct canopy_strata_object *, int thintyp,
			  struct mortality_struct *mort);
/* write <statefile>.fates + <statefile>.fates_host next to a RHESSys state file */
void veg_host_write_restart(const char *statefile);
void veg_host_nuptake(struct patch_object *, struct command_line_object *);
void veg_host_finish(struct command_line_object *);
/* world_daily_F two-pass switch (0 = legacy single pass, 1 = zone weather
   only, 2 = zone patches only) and the batched FATES daily step */
extern int veg_zone_phase;
void veg_host_daily_batch(struct command_line_object *, struct date);

/* WMFire coupling (-demog 3 + -firespread): FATES fuels of a patch (kgC/m2; returns 0
   if the patch is not a FATES site) and one fire event per burned patch family */
struct patch_family_object;
int  veg_host_fire_fuels(struct patch_object *patch, double *veg, double *litter, double *cwd);
void veg_host_wmfire_family(struct patch_family_object *fam);
/* -demog_nneed: plant-priority soil N competition (returns 0 = use resolve_sminn_competition) */
int  veg_host_sminn_competition(struct patch_object *patch, struct command_line_object *command_line);
/* diagnostics: FATES litter C to this patch today by source (gC/m2/day); 0 = not a FATES site */
int  veg_host_litter_src(struct patch_object *patch, double *v);

#endif /* _FORTRAN_SOURCE_ */
#endif /* _RHESSYS_VEG_H_ */
