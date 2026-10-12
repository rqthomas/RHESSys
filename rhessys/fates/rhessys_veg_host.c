/*--------------------------------------------------------------------------*/
/* rhessys_veg_host.c -- C host driver for the FATES vegetation backend     */
/*                                                                          */
/* ED/FATES design doc 9.4a-d, 12 (Phase 1). Plays the role of CLM's        */
/* clmfates_interfaceMod on the host side:                                  */
/*   veg_host_init   : every RHESSys patch becomes one FATES site; fixed    */
/*                     soil layers per site (soil adaptor geometry); the    */
/*                     backend allocates sites + boundary conditions.       */
/*   veg_host_daily  : called at the start of patch_daily_F; packs today's  */
/*                     soil state (soil adaptor) and daylight-mean          */
/*                     atmosphere/radiation into struct veg_bc, hands it to */
/*                     the backend and runs the FATES daily flux step       */
/*                     (Phase 2); checks the soil-adaptor mass balance.     */
/*   veg_host_transpiration : after the legacy canopy, before the patch     */
/*                     sums the soil-water demand: FATES transpiration =    */
/*                     RHESSys Penman-Monteith with the FATES canopy        */
/*                     conductance (as CLM computes transpiration from      */
/*                     FATES stomatal resistance). demog_flag 2: it         */
/*                     replaces the legacy strata transpiration demand.     */
/*   veg_host_daily_end : after the patch water balance; legacy vs FATES    */
/*                     comparison and diagnostics.                          */
/*   veg_host_finish : summary + backend cleanup.                           */
/* Sites: every patch with at least one vegetated (non NON_VEG) stratum     */
/* (water/reservoir patches are not FATES sites, as CLM uses only soil      */
/* columns). demog_flag 1 leaves RHESSys state untouched.                   */
/*                                                                          */
/* Diagnostics: one patch's daily boundary conditions are written to        */
/* <output prefix>_fates_bc.csv (patch from env RHESSYS_FATES_DIAG_PATCH,   */
/* default the first site).                                                 */
/*--------------------------------------------------------------------------*/
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "rhessys.h"
#include "rhessys_veg.h"
#include "rhessys_soil_adaptor.h"

long	julday(struct date);
double	compute_saturation_vapor_pressure(double);
double	penman_monteith(int, double, double, double, double, double, double, int);

#define VH_NLEVSOIL	3	/* soil-adaptor layers (design doc 9.4d: start coarse) */
#define VH_NLEVDECOMP	1	/* RHESSys owns decomposition */
#define VH_EPAR		4.55	/* (umol/J) as EPAR in phys_constants.h (RHESSys PAR from Kdown) */
#define VH_TFRZ		273.15
#define VH_O2_FRAC	0.209	/* O2 mole fraction (CLM) */
#define VH_SNOW_ALB_VIS	0.85	/* ground albedo under snow, provisional (CLM-like) */
#define VH_SNOW_ALB_NIR	0.65
#define VH_MASS_TOL	1.0e-9	/* (m water) soil-adaptor mass-check tolerance */

struct site_map {
	struct patch_object	*patch;
	int			site;	/* 1-based FATES site index */
};

static struct site_map	*smap = NULL;
static int		nsites = 0;
static double		(*site_zi)[VEG_NLEV_MAX + 1] = NULL;
static struct veg_out	*site_out = NULL;	/* today's flux-step results per site */
static double		*site_tf = NULL;	/* today's FATES transpiration (m/day) */
static double		*site_tleg = NULL;	/* today's legacy transpiration demand (m/day) */
/* mode 3: N uptake granted by RHESSys (gN/m2/day) for the next FATES step, */
/* FATES N stock yesterday (kg/m2), legacy strata C+N fingerprint at start  */
static double		*site_nh4 = NULL, *site_no3 = NULL, *site_nprev = NULL, *site_legacy0 = NULL;
static double		tot_effc = 0.0, tot_effn = 0.0;
/* Phase 6: TEC disturbance events pending for the next FATES step, per site
   and leaf habit ([0] evergreen, [1] deciduous): removal (thin_harvest) and
   mortality (thin_remain / _fire_remain / _snags) fractions of above-ground C */
static double		(*ev_rem)[2] = NULL, (*ev_kill)[2] = NULL;
static int		(*ev_has)[2] = NULL;
static long		n_fire_events = 0, n_harvest_events = 0;
static double		resize_max = 0.0;	/* max |FATES patch_resize_err N| (kg/m2/day) */
static double		err_max = 0.0, res_err_max = 0.0;	/* FATES's own N balance error */
static double		tot_firec = 0.0, tot_firen = 0.0, tot_harvc = 0.0, tot_harvn = 0.0;
/* FATES fuel consumed in a fire: RHESSys's litter consumption fractions
   (litr1/2 100 %, litr4 71 %; compute_fire_effects) weighted by FATES leaf-
   litter chemistry (0.25 / 0.50 / 0.25 labile / cellulose / lignin) */
#define VH_FIRE_FUEL_BURN	(0.25 * 1.0 + 0.50 * 1.0 + 0.25 * 0.71)
static char		cohort_file[FILEPATH_LEN + 32] = "rhessys_fates_cohorts.csv";
static int		*site_pid = NULL;	/* RHESSys patch ID of each site (restart check) */
/* batch (world_daily_F pass 1): per site its patch/zone/basin, today's bc */
int			veg_zone_phase = 0;
static double		t_pack = 0.0, t_fates = 0.0;	/* wall seconds in the batch */
static double		sum_coh = 0.0, sum_pat = 0.0, n_batch_sd = 0.0;	/* cohorts, patches per site-day */
static double wall_now(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return(ts.tv_sec + 1.0e-9 * ts.tv_nsec);
}
static struct patch_object **site_patch = NULL;
static struct zone_object **site_zone = NULL;
static struct basin_object **site_basin = NULL;
static struct veg_bc	*site_bc = NULL;
static int		*site_nbad = NULL;
#define VH_HOST_MAGIC	"RHFHOST2"
/* FATES model day continues across restarts / spin-up segments */
static double		md_offset = 0.0, md_last = -1.0;
static double		nb_max = 0.0, nb_sum = 0.0, leg_max = 0.0, tot_up = 0.0, tot_recv = 0.0,
			tot_litc = 0.0, tot_litn = 0.0, tot_area_days = 0.0;
static long		nb_n = 0, nb_big = 0;
static int		host_mode = 0;		/* demog_flag of this run */
static double		split_err_max = 0.0;
static double		nb_max_vs_fates = 0.0, nb_max_d = 0, nb_max_recv = 0, nb_max_lit = 0,
			nb_max_mbup = 0, nb_max_mbfrag = 0;
static int		nb_max_patch = -1;
static long		last_julday = -1;
static long		start_julday = 0;
static int		reference_date = 0;

/* diagnostics */
static FILE		*diag = NULL;
static struct patch_object *diag_patch = NULL;
static double		max_resid = 0.0, sum_resid = 0.0;
static long		n_site_days = 0, n_bad_total = 0, n_resid_fail = 0;
static int		max_resid_patch = -1;
/* water balance (patch_daily_F's patch.water_balance), all sites */
static double		wb_max = 0, wb_sum = 0;
static long		wb_n = 0, wb_big = 0;
/* basin-mean daily comparison (area-weighted over sites) */
static FILE		*fbasin = NULL, *fflux = NULL;
static long		bday = -1;
static struct date	bdate;
static double		b_area, b_gpp_f, b_gpp_l, b_tf, b_tleg, b_tact, b_lai_f, b_lai_l,
			b_btran, b_gc, b_sw, b_wbmax, b_bio, b_lit, b_fleaf, b_fcwd, b_rlit, b_rlitbg;
static long		b_n;
/* psi tracking: patch.psi vs the same call on the adaptor's rebuilt patch.S */
static double		trk_n = 0, trk_x = 0, trk_y = 0, trk_xx = 0, trk_yy = 0, trk_xy = 0;
static double		trk_absdiff = 0, trk_maxdiff = 0, trk_close = 0;
static struct date	max_resid_date;

static int find_site(struct patch_object *patch);

static int cmp_map(const void *a, const void *b)
{
	const struct patch_object *pa = ((const struct site_map *) a)->patch;
	const struct patch_object *pb = ((const struct site_map *) b)->patch;
	return((pa > pb) - (pa < pb));
}

/* a FATES site: at least one vegetated stratum (not NON_VEG) */
static int is_site(struct patch_object *patch)
{
	int i;
	for (i = 0; i < patch[0].num_canopy_strata; i++)
		if (patch[0].canopy_strata[i][0].defaults[0][0].epc.veg_type != NON_VEG)
			return(1);
	return(0);
}

/*--------------------------------------------------------------------------*/
/* Mode 2 (Phase 3b): the vegetated legacy strata take FATES's canopy        */
/* structure for today's RHESSys physics (radiation, interception,           */
/* aerodynamics, evaporation), after update_phenology (patch_daily_I) set    */
/* the legacy values. FATES LAI is split by leaf habit: evergreen strata get  */
/* the evergreen share, deciduous strata the rest (all to the other habit if */
/* a patch lacks one); within a habit LAI is per stratum cover. Sunlit share */
/* = FATES daylight laisun/(laisun+laisha). Heights = FATES canopy top, then */
/* the layers are re-sorted (canopy_stratum_daily_F finds the top layer by    */
/* height). Legacy carbon pools run on as an unused shadow (Phase 4).        */
/*--------------------------------------------------------------------------*/
void	sort_patch_layers(struct patch_object *, int *);

static void apply_structure(struct patch_object *patch, struct veg_out *o)
{
	int	i, rec = 0, habit;
	double	cov[2] = {0.0, 0.0}, lai_h[2], sai_h[2], fe, fsun, lai, sai;
	struct canopy_strata_object *st;

	for (i = 0; i < patch[0].num_canopy_strata; i++) {
		st = patch[0].canopy_strata[i];
		if (st[0].defaults[0][0].epc.veg_type == NON_VEG)
			continue;
		habit = (st[0].defaults[0][0].epc.phenology_type == EVERGREEN) ? 0 : 1;
		cov[habit] += st[0].cover_fraction;
	}
	fe = o[0].lai_evg_frac;
	if (cov[0] <= ZERO) fe = 0.0;
	if (cov[1] <= ZERO) fe = 1.0;
	lai_h[0] = o[0].tlai * fe;  lai_h[1] = o[0].tlai * (1.0 - fe);
	sai_h[0] = o[0].tsai * fe;  sai_h[1] = o[0].tsai * (1.0 - fe);
	fsun = ((o[0].laisun + o[0].laisha) > ZERO) ? o[0].laisun / (o[0].laisun + o[0].laisha) : 0.0;

	for (i = 0; i < patch[0].num_canopy_strata; i++) {
		st = patch[0].canopy_strata[i];
		if (st[0].defaults[0][0].epc.veg_type == NON_VEG)
			continue;
		habit = (st[0].defaults[0][0].epc.phenology_type == EVERGREEN) ? 0 : 1;
		if (cov[habit] <= ZERO)
			continue;
		lai = lai_h[habit] / cov[habit];
		sai = sai_h[habit] / cov[habit];
		st[0].epv.proj_lai = lai;
		st[0].epv.all_lai = lai * st[0].defaults[0][0].epc.lai_ratio;
		st[0].epv.proj_lai_sunlit = fsun * lai;
		st[0].epv.proj_lai_shade = lai - st[0].epv.proj_lai_sunlit;
		if (st[0].defaults[0][0].epc.veg_type == TREE) {
			st[0].epv.proj_pai = lai + sai;
			st[0].epv.all_pai = st[0].epv.all_lai + sai;
		}
		else {
			st[0].epv.proj_pai = lai;
			st[0].epv.all_pai = st[0].epv.all_lai;
		}
		st[0].epv.height = o[0].htop;
	}
	sort_patch_layers(patch, &rec);
}

/*--------------------------------------------------------------------------*/
/* Mode 3 (Phase 4): FATES litter (labile / cellulose / lignin, already      */
/* fragmented from FATES CWD and leaf/root litter) goes to RHESSys litter    */
/* pools litr1 / litr2 / litr4 (g -> kg); RHESSys decomposition takes it     */
/* from there. Budget: FATES total N change must equal N uptake received    */
/* minus litter N sent to RHESSys (FATES has no other N exits with fire and  */
/* land use off); the interface itself is checked by comparing the uptake    */
/* handed to FATES with the uptake it distributed to cohorts.                */
/*--------------------------------------------------------------------------*/
static void fates_litter_and_budget(struct patch_object *patch, int site, double given)
{
	double	litn, res, res2;
	struct veg_out *o = &site_out[site - 1];

	/* surface fragments (leaves, above-ground CWD, seeds) -> forest-floor
	   litter; fine roots + below-ground CWD -> RHESSys below-ground litter
	   (decomposed by update_decomp_root) */
	patch[0].litter_cs.litr1c += o[0].ls_c_lab / 1000.0;
	patch[0].litter_cs.litr2c += o[0].ls_c_cel / 1000.0;
	patch[0].litter_cs.litr4c += o[0].ls_c_lig / 1000.0;
	patch[0].litter_ns.litr1n += o[0].ls_n_lab / 1000.0;
	patch[0].litter_ns.litr2n += o[0].ls_n_cel / 1000.0;
	patch[0].litter_ns.litr4n += o[0].ls_n_lig / 1000.0;
	patch[0].litter_cs.litr1c_bg += o[0].lb_c_lab / 1000.0;
	patch[0].litter_cs.litr2c_bg += o[0].lb_c_cel / 1000.0;
	patch[0].litter_cs.litr4c_bg += o[0].lb_c_lig / 1000.0;
	patch[0].litter_ns.litr1n_bg += o[0].lb_n_lab / 1000.0;
	patch[0].litter_ns.litr2n_bg += o[0].lb_n_cel / 1000.0;
	patch[0].litter_ns.litr4n_bg += o[0].lb_n_lig / 1000.0;
	split_err_max = max(split_err_max, o[0].lit_split_err);
	/* root efflux -> below-ground labile litter, where FATES meant to put it */
	patch[0].litter_cs.litr1c_bg += o[0].efflux_c / 1000.0;
	patch[0].litter_ns.litr1n_bg += o[0].efflux_n / 1000.0;
	tot_effc += patch[0].area * o[0].efflux_c;
	tot_effn += patch[0].area * o[0].efflux_n;
	litn = o[0].lit_n_lab + o[0].lit_n_cel + o[0].lit_n_lig + o[0].efflux_n;
	tot_firec += patch[0].area * o[0].fire_c;
	tot_firen += patch[0].area * o[0].fire_n;
	tot_harvc += patch[0].area * o[0].harvest_c;
	tot_harvn += patch[0].area * o[0].harvest_n;
	tot_litc += patch[0].area * (o[0].lit_c_lab + o[0].lit_c_cel + o[0].lit_c_lig);
	tot_litn += patch[0].area * litn;
	tot_up += patch[0].area * given;
	tot_recv += patch[0].area * o[0].n_uptake_recv;
	tot_area_days += patch[0].area;
	if (site_nprev[site - 1] >= 0.0) {
		/* FATES N change = uptake received - (litter + efflux) N sent
		   - N burned - N exported in harvested wood */
		res = (o[0].n_stock - site_nprev[site - 1])
			- (o[0].n_uptake_recv - litn - o[0].fire_n - o[0].harvest_n) / 1000.0
			- o[0].mb_n_other;	/* FATES-internal terms (seeds, patch resizing, ...) */
		resize_max = max(resize_max, fabs(o[0].mb_n_resize));
		err_max = max(err_max, fabs(o[0].mb_n_err));
		res_err_max = max(res_err_max, fabs(res + o[0].mb_n_err));
		/* same change against FATES's own daily record (root uptake, fragmentation) */
		res2 = (o[0].n_stock - site_nprev[site - 1]) - (o[0].mb_n_uptake - o[0].mb_n_frag);
		nb_sum += fabs(res);
		nb_n++;
		if (fabs(res) > 1.0e-9)
			nb_big++;
		nb_max_vs_fates = max(nb_max_vs_fates, fabs(res2));
		if (fabs(res) > nb_max) {
			nb_max = fabs(res);
			nb_max_patch = patch[0].ID;
			nb_max_d = (o[0].n_stock - site_nprev[site - 1]) * 1000.0;
			nb_max_recv = o[0].n_uptake_recv;
			nb_max_lit = litn;
			nb_max_mbup = o[0].mb_n_uptake * 1000.0;
			nb_max_mbfrag = o[0].mb_n_frag * 1000.0;
		}
	}
	site_nprev[site - 1] = o[0].n_stock;
}

/*--------------------------------------------------------------------------*/
/* What RHESSys sees of the FATES forest floor. FATES keeps dead leaves in   */
/* its own pool until they fragment (fates_frag_maxdecomp, ~1/yr for leaves) */
/* so in mode 3 RHESSys's litr1-4 alone would miss them. Above-ground CWD is */
/* not added: legacy RHESSys also kept dead wood out of litr1-4 (stratum     */
/* cwdc) until it fragmented, so cover_litter_a / sediment_transport_litter_b */
/* keep their calibrated meaning.                                           */
/*--------------------------------------------------------------------------*/
double veg_host_surface_litter_c(struct patch_object *patch)
{
	int site;
	if (host_mode != 3)
		return(0.0);
	site = find_site(patch);
	return((site > 0) ? max(site_out[site - 1].leaf_litter_c, 0.0) : 0.0);
}

/* mode 3: litter interception/depth from RHESSys litter + FATES dead leaves,
   same formula as update_litter_interception_capacity */
void veg_host_litter_interception(struct patch_object *patch)
{
	double extra, total;
	extra = veg_host_surface_litter_c(patch);
	if (extra <= 0.0)
		return;
	total = patch[0].litter_cs.litr1c + patch[0].litter_cs.litr2c
		+ patch[0].litter_cs.litr3c + patch[0].litter_cs.litr4c + extra;
	patch[0].litter.rain_capacity = max(total * patch[0].litter.moist_coef
		* patch[0].litter.cover_fraction, 0.0);
	patch[0].litter.proj_pai = 1.0 * patch[0].litter.cover_fraction;
	patch[0].litter.depth = total / patch[0].litter.density;
}

/* legacy plant C + N of a patch's vegetated strata (must stay constant in mode 3) */
static double legacy_plant_cn(struct patch_object *patch)
{
	int	i;
	double	t = 0.0;
	struct canopy_strata_object *st;
	for (i = 0; i < patch[0].num_canopy_strata; i++) {
		st = patch[0].canopy_strata[i];
		if (st[0].defaults[0][0].epc.veg_type == NON_VEG)
			continue;
		t += st[0].cs.leafc + st[0].cs.frootc + st[0].cs.live_stemc + st[0].cs.dead_stemc
			+ st[0].cs.live_crootc + st[0].cs.dead_crootc + st[0].cs.cpool
			+ st[0].ns.leafn + st[0].ns.frootn + st[0].ns.live_stemn + st[0].ns.dead_stemn
			+ st[0].ns.live_crootn + st[0].ns.dead_crootn + st[0].ns.npool
			+ st[0].cs.mr_deficit + st[0].cs.availc;
	}
	return(t);
}

/*--------------------------------------------------------------------------*/
/* Mode 3, before resolve_sminn_competition: the patch's plant N demand is   */
/* the FATES RD capacity (replaces the legacy strata sum).                    */
/*--------------------------------------------------------------------------*/
void veg_host_ndemand(struct patch_object *patch, struct command_line_object *command_line)
{
	int site;
	if (command_line[0].demog_flag != 3)
		return;
	site = find_site(patch);
	if (site == 0)
		return;
	patch[0].ndf.plant_potential_ndemand = site_out[site - 1].ncap / 1000.0;	/* kgN/m2/day */
}

/*--------------------------------------------------------------------------*/
/* Mode 3, after resolve_sminn_competition (legacy strata growth skipped):  */
/* the plant uptake RHESSys allows is debited from the soil by update_decomp */
/* (via ndf.sminn_to_npool, nitrate first) and handed to FATES tomorrow.     */
/*--------------------------------------------------------------------------*/
void veg_host_nuptake(struct patch_object *patch, struct command_line_object *command_line)
{
	int site;
	double u, no3;
	if (command_line[0].demog_flag != 3)
		return;
	site = find_site(patch);
	if (site == 0)
		return;
	u = max(patch[0].ndf.plant_avail_uptake, 0.0);			/* kgN/m2/day */
	patch[0].ndf.sminn_to_npool = u;
	no3 = min(max(patch[0].soil_ns.nitrate, 0.0), u);
	site_no3[site - 1] = 1000.0 * no3;
	site_nh4[site - 1] = 1000.0 * (u - no3);
}

/* the tallest vegetated stratum (aerodynamics + longwave of the canopy) */
static struct canopy_strata_object *overstory(struct patch_object *patch)
{
	int i;
	struct canopy_strata_object *top = NULL, *st;
	for (i = 0; i < patch[0].num_canopy_strata; i++) {
		st = patch[0].canopy_strata[i];
		if (st[0].defaults[0][0].epc.veg_type == NON_VEG)
			continue;
		if ((top == NULL) || (st[0].epv.height > top[0].epv.height))
			top = st;
	}
	return(top);
}

static int find_site(struct patch_object *patch)
{
	struct site_map key, *hit;
	if (smap == NULL)
		return(0);
	key.patch = patch;
	hit = bsearch(&key, smap, nsites, sizeof(struct site_map), cmp_map);
	return(hit ? hit->site : 0);
}

static void adaptor_inputs(struct patch_object *patch,
			   struct soil_adaptor_params *par,
			   struct soil_adaptor_state *st)
{
	struct soil_default *sd = patch[0].soil_defaults[0];
	par[0].porosity_0 = sd[0].porosity_0;
	par[0].porosity_decay = sd[0].porosity_decay;
	par[0].soil_depth = sd[0].soil_depth;
	par[0].soil_water_cap = sd[0].soil_water_cap;
	par[0].theta_psi_curve = sd[0].theta_psi_curve;
	par[0].psi_air_entry = sd[0].psi_air_entry;
	par[0].pore_size_index = sd[0].pore_size_index;
	par[0].p3 = sd[0].p3;
	par[0].p4 = sd[0].p4;
	par[0].psi_max = sd[0].psi_max;
	st[0].rz_depth = patch[0].rootzone.depth;
	st[0].rz_storage = patch[0].rz_storage;
	st[0].unsat_storage = patch[0].unsat_storage;
	st[0].sat_deficit = patch[0].sat_deficit;
	st[0].Tsoil = patch[0].Tsoil;
}

/* soil part of the boundary conditions from the adaptor */
static void soil_bc(struct patch_object *patch, int site, struct veg_bc *bc,
		    struct soil_adaptor_out *out, double tsoil)
{
	int j, n = VH_NLEVSOIL;
	struct soil_adaptor_params par;
	struct soil_adaptor_state st;
	const double *zi = site_zi[site - 1];

	adaptor_inputs(patch, &par, &st);
	st.Tsoil = tsoil;
	soil_adaptor_profile(&par, &st, n, zi, out);
	bc[0].nlevsoil = n;
	bc[0].max_rooting_depth_index = n;	/* soil_depth is the bedrock */
	bc[0].w_scalar = patch[0].cdf.decomp_w_scalar;	/* from patch_daily_I (today) */
	bc[0].t_scalar = patch[0].cdf.decomp_t_scalar;
	bc[0].zi[0] = zi[0];
	for (j = 0; j < n; j++) {
		bc[0].zi[j+1] = zi[j+1];
		bc[0].dz[j] = zi[j+1] - zi[j];
		bc[0].z[j] = 0.5 * (zi[j] + zi[j+1]);
		bc[0].watsat[j] = out[0].watsat[j];
		bc[0].h2o_liqvol[j] = min(max(out[0].theta[j], 0.0), out[0].watsat[j]);
		bc[0].eff_porosity[j] = max(out[0].watsat[j], 0.01);	/* no soil ice in RHESSys */
		bc[0].smp[j] = out[0].smp_mm[j];
		bc[0].tempk[j] = tsoil + VH_TFRZ;	/* zone soil temperature, all layers */
	}
}

/* daylength (s) for latitude phi and declination delta (radians) */
static double daylength(double phi, double delta, double *H)
{
	double c = -tan(phi) * tan(delta);
	c = min(max(c, -1.0), 1.0);
	*H = acos(c);
	return(86400.0 * (*H) / PI);
}

/* atmosphere/radiation part (daylight means) */
static void atm_bc(struct basin_object *basin, struct zone_object *zone,
		   struct patch_object *patch, struct date current_date, struct veg_bc *bc)
{
	int	i;
	long	doy0;
	double	phi, delta, H, Hmax, dayl, dayl_max, esat, kdir, kdif, pvis, pdif, swe;
	struct date jan1;

	jan1 = current_date; jan1.month = 1; jan1.day = 1;
	doy0 = julday(current_date) - julday(jan1);		/* 0-based */
	phi = basin[0].latitude * DtoR;
	delta = 23.44 * DtoR * sin(2.0 * PI * (284.0 + doy0 + 1) / 365.0);
	daylength(phi, delta, &H);
	dayl_max = daylength(fabs(phi), 23.44 * DtoR, &Hmax);
	dayl = zone[0].metv.dayl;

	/* daylight-mean cosine of the solar zenith angle */
	if (H > 0.0)
		bc[0].coszen = sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(H) / H;
	else
		bc[0].coszen = 0.0;
	bc[0].coszen = min(max(bc[0].coszen, 0.0), 1.0);
	/* CLM's daylength factor */
	bc[0].dayl_factor = (dayl_max > 0.0) ?
		min(1.0, max(0.01, (dayl * dayl) / (dayl_max * dayl_max))) : 1.0;

	/* shortwave: kJ/m2/day and umol/m2/day -> daylight-mean W/m2, vis = PAR */
	if (dayl > 1.0) {
		kdir = 1000.0 * zone[0].Kdown_direct / dayl;
		kdif = 1000.0 * zone[0].Kdown_diffuse / dayl;
		pvis = (zone[0].PAR_direct >= 0.0) ? zone[0].PAR_direct / VH_EPAR / dayl : 0.5 * kdir;
		pdif = (zone[0].PAR_diffuse >= 0.0) ? zone[0].PAR_diffuse / VH_EPAR / dayl : 0.5 * kdif;
	}
	else {
		kdir = kdif = pvis = pdif = 0.0;
	}
	bc[0].solad[0] = max(pvis, 0.0);
	bc[0].solad[1] = max(kdir - pvis, 0.0);
	bc[0].solai[0] = max(pdif, 0.0);
	bc[0].solai[1] = max(kdif - pdif, 0.0);

	/* air: daylight temperature, vapour pressure from the dew point */
	bc[0].forc_pbot = zone[0].metv.pa;
	bc[0].tgcm = zone[0].metv.tday + VH_TFRZ;
	bc[0].t_veg = bc[0].tgcm;
	esat = compute_saturation_vapor_pressure(zone[0].metv.tday);
	bc[0].esat_tv = esat;
	bc[0].eair = min(max(zone[0].e_dewpoint, 1.0), esat);
	bc[0].oair = VH_O2_FRAC * bc[0].forc_pbot;
	bc[0].cair = zone[0].CO2 * 1.0e-6 * bc[0].forc_pbot;
	/* CLM leaf boundary layer: rb = 1 / (0.01 sqrt(u / 0.04 m)); provisional */
	bc[0].rb = 1.0 / (0.01 * sqrt(max(zone[0].wind, 0.1) / 0.04));

	/* ground albedo and snow */
	swe = patch[0].snowpack.water_equivalent_depth;
	for (i = 0; i < VEG_NSWB; i++) {
		bc[0].albgr_dir[i] = patch[0].soil_defaults[0][0].albedo;
		bc[0].albgr_dif[i] = patch[0].soil_defaults[0][0].albedo;
	}
	if (swe > 0.0) {
		bc[0].albgr_dir[0] = bc[0].albgr_dif[0] = VH_SNOW_ALB_VIS;
		bc[0].albgr_dir[1] = bc[0].albgr_dif[1] = VH_SNOW_ALB_NIR;
	}
	bc[0].snow_depth = max(patch[0].snowpack.height, 0.0);
	bc[0].frac_sno_eff = (swe > 0.0) ? 1.0 : 0.0;
	bc[0].fcansno = 0.0;
	for (i = 0; i < patch[0].num_canopy_strata; i++)
		if (patch[0].canopy_strata[i][0].snow_stored > 0.0)
			bc[0].fcansno += patch[0].canopy_strata[i][0].cover_fraction;
	bc[0].fcansno = min(max(bc[0].fcansno, 0.0), 1.0);

	/* "24 h" means = today's daily values (daily host step) */
	bc[0].precip24 = max(zone[0].rain + zone[0].snow, 0.0) * 1000.0 / 86400.0;
	bc[0].relhumid24 = 100.0 * zone[0].relative_humidity;
	bc[0].wind24 = zone[0].wind;

	/* night sub-step + daily mean (Phase 2) */
	bc[0].dayl = min(max(dayl, 0.0), 86400.0);
	bc[0].t_night = zone[0].metv.tnight + VH_TFRZ;
	bc[0].esat_night = compute_saturation_vapor_pressure(zone[0].metv.tnight);
	bc[0].eair_night = min(max(zone[0].e_dewpoint, 1.0), bc[0].esat_night);
	bc[0].t_avg = zone[0].metv.tavg + VH_TFRZ;
}

static void write_diag_header(void)
{
	int j;
	fprintf(diag, "year,month,day,patchID,zwt,rz_depth,rz_storage,unsat_storage,sat_deficit,"
		"ref_water,col_water,resid,patch_S,rootzone_S,patch_psi,S_rz,S_unsat,"
		"S_patch_equiv,psi_patch_equiv");
	for (j = 1; j <= VH_NLEVSOIL; j++)
		fprintf(diag, ",zi_%d,theta_%d,watsat_%d,S_%d,smp_mm_%d,psi_patch_%d", j, j, j, j, j, j);
	fprintf(diag, ",tgcm,t_veg,eair,esat_tv,forc_pbot,oair,cair,rb,coszen,dayl_factor,"
		"solad_vis,solad_nir,solai_vis,solai_nir,albgr_vis,albgr_nir,snow_depth,frac_sno,"
		"fcansno,precip24,relhumid24,wind24,tempk_1,n_bad,w_scalar,t_scalar\n");
}

static void write_diag_row(struct patch_object *patch, struct date d, struct veg_bc *bc,
			   struct soil_adaptor_out *out, int n_bad)
{
	int j;
	fprintf(diag, "%ld,%ld,%ld,%d,%.6f,%.6f,%.9f,%.9f,%.9f,%.12f,%.12f,%.3e,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f",
		d.year, d.month, d.day, patch[0].ID, out[0].zwt, patch[0].rootzone.depth,
		patch[0].rz_storage, patch[0].unsat_storage, patch[0].sat_deficit,
		out[0].ref_water, out[0].col_water, out[0].col_water - out[0].ref_water,
		patch[0].S, patch[0].rootzone.S, patch[0].psi, out[0].S_rz, out[0].S_unsat,
		out[0].S_patch_equiv, out[0].psi_patch_equiv);
	for (j = 0; j < VH_NLEVSOIL; j++)
		fprintf(diag, ",%.4f,%.6f,%.6f,%.6f,%.6g,%.6f", bc[0].zi[j+1], out[0].theta[j],
			out[0].watsat[j], out[0].S[j], bc[0].smp[j], out[0].psi_patch[j]);
	fprintf(diag, ",%.3f,%.3f,%.2f,%.2f,%.1f,%.1f,%.3f,%.3f,%.5f,%.5f,%.3f,%.3f,%.3f,%.3f,"
		"%.3f,%.3f,%.4f,%.1f,%.3f,%.4e,%.2f,%.3f,%.3f,%d,%.4f,%.4f\n",
		bc[0].tgcm, bc[0].t_veg, bc[0].eair, bc[0].esat_tv, bc[0].forc_pbot, bc[0].oair,
		bc[0].cair, bc[0].rb, bc[0].coszen, bc[0].dayl_factor, bc[0].solad[0], bc[0].solad[1],
		bc[0].solai[0], bc[0].solai[1], bc[0].albgr_dir[0], bc[0].albgr_dir[1],
		bc[0].snow_depth, bc[0].frac_sno_eff, bc[0].fcansno, bc[0].precip24,
		bc[0].relhumid24, bc[0].wind24, bc[0].tempk[0], n_bad, bc[0].w_scalar, bc[0].t_scalar);
}

/* FATES model time for a RHESSys date (day_of_year 0-based, real calendar) */
static void fill_time(struct date d, struct veg_time *t)
{
	long jd;
	struct date jan1, jan1n;
	jd = julday(d);
	jan1 = d; jan1.month = 1; jan1.day = 1;
	jan1n = jan1; jan1n.year += 1;
	t[0].year = (int) d.year;
	t[0].month = (int) d.month;
	t[0].day = (int) d.day;
	t[0].day_of_year = (int) (jd - julday(jan1));
	t[0].days_per_year = (int) (julday(jan1n) - julday(jan1));
	t[0].current_date = t[0].year * 10000 + t[0].month * 100 + t[0].day;
	t[0].reference_date = reference_date;
	t[0].model_day = md_offset + (double) (jd - start_julday);
	md_last = t[0].model_day;
}

/*--------------------------------------------------------------------------*/
/* Restart (Phase 5): <statefile>.fates = FATES restart vectors (Fortran),   */
/* <statefile>.fates_host = this driver's carry-over state per site: the    */
/* last flux-step outputs (struct veg_out; patch_daily_I reads the FATES     */
/* structure and dead-leaf litter before the next flux step) and the N      */
/* uptake RHESSys granted but FATES has not yet received. Written at        */
/* output_current_state, so a run restarted from that state file continues  */
/* as the uninterrupted run would.                                          */
/*--------------------------------------------------------------------------*/
void veg_host_write_restart(const char *statefile)
{
	int	s, n = nsites, sz = (int) sizeof(struct veg_out);
	size_t	len;
	char	file[FILEPATH_LEN + 160];
	FILE	*f;

	if (smap == NULL)
		return;
	sprintf(file, "%s.fates", statefile);
	len = strlen(file);
	veg_restart_write(file, &len);
	sprintf(file, "%s.fates_host", statefile);
	if ((f = fopen(file, "wb")) == NULL) {
		fprintf(stderr, "FATAL ERROR: cannot write %s\n", file);
		exit(EXIT_FAILURE);
	}
	fwrite(VH_HOST_MAGIC, 1, 8, f);
	fwrite(&n, sizeof(int), 1, f);
	fwrite(&sz, sizeof(int), 1, f);
	fwrite(&md_last, sizeof(double), 1, f);		/* model day of the last day run */
	fwrite(&reference_date, sizeof(int), 1, f);
	for (s = 0; s < nsites; s++) {
		fwrite(&site_pid[s], sizeof(int), 1, f);
		fwrite(&site_out[s], sizeof(struct veg_out), 1, f);
		fwrite(&site_nh4[s], sizeof(double), 1, f);
		fwrite(&site_no3[s], sizeof(double), 1, f);
	}
	fclose(f);
	fprintf(stderr, "FATES host: restart written: %s.fates(+_host)\n", statefile);
}

/* model-day count + reference date only (needed before the FATES restart read) */
static void read_host_header(const char *file)
{
	int	n, sz, refd;
	double	mdl;
	char	magic[8];
	FILE	*f;
	if (((f = fopen(file, "rb")) == NULL) || (fread(magic, 1, 8, f) != 8)
	    || (memcmp(magic, VH_HOST_MAGIC, 8) != 0) || (fread(&n, sizeof(int), 1, f) != 1)
	    || (fread(&sz, sizeof(int), 1, f) != 1) || (fread(&mdl, sizeof(double), 1, f) != 1)
	    || (fread(&refd, sizeof(int), 1, f) != 1)) {
		fprintf(stderr, "FATAL ERROR: cannot read FATES host restart %s\n", file);
		exit(EXIT_FAILURE);
	}
	fclose(f);
	md_offset = mdl + 1.0;		/* the restarted run starts the next day */
	reference_date = refd;
}

static void read_host_restart(const char *file)
{
	int	s, n, sz, pid, refd;
	double	mdl;
	char	magic[8];
	FILE	*f;

	if ((f = fopen(file, "rb")) == NULL) {
		fprintf(stderr, "FATAL ERROR: cannot open FATES host restart %s\n", file);
		exit(EXIT_FAILURE);
	}
	if ((fread(magic, 1, 8, f) != 8) || (memcmp(magic, VH_HOST_MAGIC, 8) != 0)
	    || (fread(&n, sizeof(int), 1, f) != 1) || (fread(&sz, sizeof(int), 1, f) != 1)
	    || (n != nsites)
	    || ((sz != (int) sizeof(struct veg_out)) && (sz != (int) offsetof(struct veg_out, fire_c)))
	    || (fread(&mdl, sizeof(double), 1, f) != 1) || (fread(&refd, sizeof(int), 1, f) != 1)) {
		fprintf(stderr, "FATAL ERROR: FATES host restart %s does not match this run\n", file);
		exit(EXIT_FAILURE);
	}
	md_offset = mdl + 1.0;		/* the restarted run starts the next day */
	reference_date = refd;
	for (s = 0; s < nsites; s++) {
		if ((fread(&pid, sizeof(int), 1, f) != 1) || (pid != site_pid[s])
		    || (memset(&site_out[s], 0, sizeof(struct veg_out)) == NULL)
		    /* pre-Phase-6 sidecars stop before fire_c (the per-day disturbance
		       fluxes, zero-filled here) */
		    || (fread(&site_out[s], (size_t) sz, 1, f) != 1)
		    || (fread(&site_nh4[s], sizeof(double), 1, f) != 1)
		    || (fread(&site_no3[s], sizeof(double), 1, f) != 1)) {
			fprintf(stderr, "FATAL ERROR: FATES host restart %s: site %d (patch %d) mismatch\n",
				file, s + 1, site_pid[s]);
			exit(EXIT_FAILURE);
		}
	}
	fclose(f);
}

/*--------------------------------------------------------------------------*/
void veg_host_init(struct world_object *world, struct command_line_object *command_line)
{
	int	b, h, z, p, s, nlevsoil = VH_NLEVSOIL, nlevdecomp = VH_NLEVDECOMP, diag_id = -1;
	int	n_nonveg;
	size_t	pf_len, inv_len, log_len;
	double	*lat;
	char	*env, path[FILEPATH_LEN + 32], path2[FILEPATH_LEN + 32];
	struct basin_object *basin;
	struct zone_object *zone;
	struct patch_object *patch;
	struct veg_bc bc;
	struct soil_adaptor_out out;
	struct date d0;

	/*--- every vegetated patch is a site ---*/
	nsites = 0;
	n_nonveg = 0;
	for (b = 0; b < world[0].num_basin_files; b++)
		for (h = 0; h < world[0].basins[b][0].num_hillslopes; h++)
			for (z = 0; z < world[0].basins[b][0].hillslopes[h][0].num_zones; z++) {
				zone = world[0].basins[b][0].hillslopes[h][0].zones[z];
				for (p = 0; p < zone[0].num_patches; p++) {
					if (is_site(zone[0].patches[p]))
						nsites++;
					else
						n_nonveg++;
				}
			}
	if (nsites == 0) {
		fprintf(stderr, "FATAL ERROR: veg_host_init: no vegetated patches for FATES\n");
		exit(EXIT_FAILURE);
	}
	smap = (struct site_map *) calloc(nsites, sizeof(struct site_map));
	site_zi = calloc(nsites, sizeof(*site_zi));
	lat = (double *) calloc(nsites, sizeof(double));
	site_pid = (int *) calloc(nsites, sizeof(int));
	site_patch = calloc(nsites, sizeof(*site_patch));
	site_zone = calloc(nsites, sizeof(*site_zone));
	site_basin = calloc(nsites, sizeof(*site_basin));
	site_bc = (struct veg_bc *) calloc(nsites, sizeof(struct veg_bc));
	site_nbad = (int *) calloc(nsites, sizeof(int));
	ev_rem = calloc(nsites, sizeof(*ev_rem));
	ev_kill = calloc(nsites, sizeof(*ev_kill));
	ev_has = calloc(nsites, sizeof(*ev_has));
	if ((site_patch == NULL) || (site_zone == NULL) || (site_basin == NULL) || (site_bc == NULL)
	    || (site_nbad == NULL) || (ev_rem == NULL) || (ev_kill == NULL) || (ev_has == NULL)) {
		fprintf(stderr, "FATAL ERROR: veg_host_init: out of memory\n");
		exit(EXIT_FAILURE);
	}
	if ((smap == NULL) || (site_zi == NULL) || (lat == NULL)) {
		fprintf(stderr, "FATAL ERROR: veg_host_init: out of memory (%d sites)\n", nsites);
		exit(EXIT_FAILURE);
	}
	s = 0;
	for (b = 0; b < world[0].num_basin_files; b++) {
		basin = world[0].basins[b];
		for (h = 0; h < basin[0].num_hillslopes; h++)
			for (z = 0; z < basin[0].hillslopes[h][0].num_zones; z++) {
				zone = basin[0].hillslopes[h][0].zones[z];
				for (p = 0; p < zone[0].num_patches; p++) {
					patch = zone[0].patches[p];
					if (!is_site(patch))
						continue;
					smap[s].patch = patch;
					smap[s].site = s + 1;
					site_pid[s] = patch[0].ID;
					site_patch[s] = patch;
					site_zone[s] = zone;
					site_basin[s] = basin;
					lat[s] = basin[0].latitude;
					if (soil_adaptor_geometry(nlevsoil, patch[0].rootzone.depth,
							patch[0].soil_defaults[0][0].soil_depth, site_zi[s]) != 0) {
						fprintf(stderr, "FATAL ERROR: veg_host_init: bad soil column "
							"for patch %d (soil_depth %g)\n", patch[0].ID,
							patch[0].soil_defaults[0][0].soil_depth);
						exit(EXIT_FAILURE);
					}
					s++;
				}
			}
	}

	/*--- backend globals + sites ---*/
	pf_len = strlen(command_line[0].fates_paramfile);
	inv_len = strlen(command_line[0].fates_invfile);
	sprintf(path, "%s_fates.log",
		(command_line[0].output_prefix != NULL) ? command_line[0].output_prefix : "rhessys");
	log_len = strlen(path);
	veg_init(command_line[0].fates_paramfile, &pf_len, command_line[0].fates_invfile,
		&inv_len, path, &log_len, &nlevsoil, &nlevdecomp, &command_line[0].demog_flag,
		&nsites, &command_line[0].demog_restart);
	for (s = 0; s < nsites; s++) {
		memset(&bc, 0, sizeof(bc));
		soil_bc(smap[s].patch, smap[s].site, &bc, &out, smap[s].patch[0].Tsoil);
		veg_init_site(&smap[s].site, &smap[s].patch[0].ID, &lat[s], &bc);
	}
	/*--- cold start: patches + cohorts (inventory or seedlings) ---*/
	site_out = (struct veg_out *) calloc(nsites, sizeof(struct veg_out));
	site_tf = (double *) calloc(nsites, sizeof(double));
	site_tleg = (double *) calloc(nsites, sizeof(double));
	site_nh4 = (double *) calloc(nsites, sizeof(double));
	site_no3 = (double *) calloc(nsites, sizeof(double));
	site_nprev = (double *) calloc(nsites, sizeof(double));
	site_legacy0 = (double *) calloc(nsites, sizeof(double));
	for (s = 0; s < nsites; s++) {
		site_nprev[s] = -1.0;
		site_legacy0[s] = -1.0;
	}
	if ((site_out == NULL) || (site_tf == NULL) || (site_tleg == NULL) || (site_nh4 == NULL)
	    || (site_no3 == NULL) || (site_nprev == NULL) || (site_legacy0 == NULL)) {
		fprintf(stderr, "FATAL ERROR: veg_host_init: out of memory\n");
		exit(EXIT_FAILURE);
	}
	d0 = world[0].start_date;
	start_julday = julday(d0);
	reference_date = (int) (d0.year * 10000 + d0.month * 100 + d0.day);
	if (command_line[0].demog_restart) {
		/* restart: FATES state + this driver's carry-over from the state file's companions */
		struct veg_time t0;
		size_t rl;
		char rfile[FILEPATH_LEN + 160];
		sprintf(rfile, "%s.fates_host", command_line[0].world_filename);
		read_host_header(rfile);		/* model-day count, reference date */
		fill_time(d0, &t0);
		veg_set_time(&t0);
		sprintf(rfile, "%s.fates", command_line[0].world_filename);
		rl = strlen(rfile);
		veg_restart_read(rfile, &rl, site_out);
		sprintf(rfile, "%s.fates_host", command_line[0].world_filename);
		read_host_restart(rfile);		/* host carry-over (overrides site_out) */
		if (command_line[0].demog_flag >= 2)
			for (s = 0; s < nsites; s++)
				apply_structure(smap[s].patch, &site_out[s]);
	}
	else
		veg_coldstart(site_out);

	/*--- diagnostics: pick the patch before sorting the map ---*/
	env = getenv("RHESSYS_FATES_DIAG_PATCH");
	if (env != NULL)
		diag_id = atoi(env);
	diag_patch = smap[0].patch;
	if (diag_id >= 0) {
		diag_patch = NULL;
		for (s = 0; s < nsites; s++)
			if (smap[s].patch[0].ID == diag_id) { diag_patch = smap[s].patch; break; }
		if (diag_patch == NULL) {
			fprintf(stderr, "WARNING: RHESSYS_FATES_DIAG_PATCH %d not found; using patch %d\n",
				diag_id, smap[0].patch[0].ID);
			diag_patch = smap[0].patch;
		}
	}
	qsort(smap, nsites, sizeof(struct site_map), cmp_map);
	free(lat);

	sprintf(path, "%s_fates_bc.csv",
		(command_line[0].output_prefix != NULL) ? command_line[0].output_prefix : "rhessys");
	diag = fopen(path, "w");
	if (diag == NULL)
		fprintf(stderr, "WARNING: cannot open FATES diagnostics file %s\n", path);
	else
		write_diag_header();
	sprintf(path2, "%s_fates_basin.csv",
		(command_line[0].output_prefix != NULL) ? command_line[0].output_prefix : "rhessys");
	fbasin = fopen(path2, "w");
	if (fbasin != NULL)
		fprintf(fbasin, "year,month,day,n_sites,gpp_fates,gpp_legacy,transp_fates,"
			"transp_legacy_demand,transp_actual,lai_fates,lai_legacy,btran,gc_canopy,"
			"swabs_canopy,wb_absmax,biomass_c,litter_c,fates_leaf_litter_c,fates_ag_cwd_c,"
			"rhessys_litter_c,rhessys_litter_bg_c\n");
	sprintf(cohort_file, "%s_fates_cohorts.csv",
		(command_line[0].output_prefix != NULL) ? command_line[0].output_prefix : "rhessys");
	sprintf(path2, "%s_fates_flux.csv",
		(command_line[0].output_prefix != NULL) ? command_line[0].output_prefix : "rhessys");
	fflux = fopen(path2, "w");
	if (fflux != NULL)
		fprintf(fflux, "year,month,day,patchID,npatches,ncohorts,tlai,elai,htop,laisun,laisha,"
			"rssun,rssha,gc_canopy,swabs_canopy,btran,gpp_day,gpp_night,resp_m,"
			"transp_fates,transp_legacy_demand,transp_actual,gpp_legacy,lai_legacy,"
			"water_balance,biomass_c,litter_c,n_stock,n_uptake_recv,lit_n,mb_n_uptake,mb_n_frag,mb_n_err,ncap\n");

	host_mode = command_line[0].demog_flag;
	last_julday = -1;
	fprintf(stderr, "FATES host: %d sites (vegetated patches; %d non-vegetated patches "
		"skipped), %d soil layers, mode %d (%s); diagnostics for patch %d -> %s\n",
		nsites, n_nonveg, nlevsoil, command_line[0].demog_flag,
		(command_line[0].demog_flag == 3) ? "FATES water, structure, N uptake, litter" :
		(command_line[0].demog_flag == 2) ? "FATES transpiration + structure drive water" : "read-only",
		diag_patch[0].ID, path);
	fprintf(stderr, "FATES host: %s (%s): site 1 has %d patches, %d cohorts, "
		"LAI %.2f, canopy height %.1f m\n",
		command_line[0].demog_restart ? "restart" : "cold start",
		command_line[0].demog_restart ? command_line[0].world_filename :
		(command_line[0].fates_invfile[0] != '\0') ? command_line[0].fates_invfile : "seedlings",
		site_out[0].npatches, site_out[0].ncohorts, site_out[0].tlai, site_out[0].htop);
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/* Phase 6: TEC thin events on FATES sites (mode 3). The legacy CCR events  */
/* are pairs on one day: thin_harvest (removal = consumed by fire, or a     */
/* clear-cut) and thin_remain / thin_fire_remain (killed, left on site).    */
/* Resolved at the next FATES step (same day: events run at H1/H2, before   */
/* world_daily_F):                                                          */
/*   removal + kill  -> fire: mortality c + k(1-c), crowns burned c         */
/*   removal only    -> harvest (clear-cut) if >= 0.99, else fire (c, c)    */
/*   kill only       -> fire-type disturbance without combustion            */
/* Above-ground fractions are used (FATES kills whole trees; the legacy     */
/* separate root fractions have no equivalent).                             */
/*--------------------------------------------------------------------------*/
int veg_host_strata_event(struct patch_object *patch, struct canopy_strata_object *st,
			  int thintyp, struct mortality_struct *mort)
{
	int site, h;
	double f;
	if (host_mode != 3)
		return(0);
	site = find_site(patch);
	if (site == 0)
		return(0);
	h = (st[0].defaults[0][0].epc.phenology_type == EVERGREEN) ? 0 : 1;
	f = min(max(max(mort[0].mort_leafc, max(mort[0].mort_livestemc, mort[0].mort_deadstemc)), 0.0), 1.0);
	if (thintyp == 2)
		ev_rem[site - 1][h] = max(ev_rem[site - 1][h], f);
	else
		ev_kill[site - 1][h] = max(ev_kill[site - 1][h], f);
	ev_has[site - 1][h] = 1;
	return(1);
}

static void resolve_disturbance(int s, struct veg_bc *bc)
{
	int h, any_rem = 0, any_kill = 0;
	double c, k;
	bc[0].dist_kind = 0.0;
	if (!ev_has[s][0] && !ev_has[s][1])
		return;
	for (h = 0; h < 2; h++)			/* a patch without one habit's stratum: same event for it */
		if (!ev_has[s][h]) {
			ev_rem[s][h] = ev_rem[s][1 - h];
			ev_kill[s][h] = ev_kill[s][1 - h];
		}
	for (h = 0; h < 2; h++) {
		if (ev_rem[s][h] > 0.0) any_rem = 1;
		if (ev_kill[s][h] > 0.0) any_kill = 1;
	}
	if (any_rem && !any_kill && (max(ev_rem[s][0], ev_rem[s][1]) >= 0.99)) {
		bc[0].dist_kind = 2.0;			/* clear-cut harvest */
		n_harvest_events++;
	}
	else if (any_rem || any_kill) {
		bc[0].dist_kind = 1.0;			/* fire, uniform severity */
		for (h = 0; h < 2; h++) {
			c = ev_rem[s][h];
			k = ev_kill[s][h];
			bc[0].fire_mort[h] = c + k * (1.0 - c);
			bc[0].crown_burn[h] = c;
		}
		bc[0].fuel_burn_leaf = any_rem ? VH_FIRE_FUEL_BURN : 0.0;
		bc[0].fuel_burn_twig = any_rem ? VH_FIRE_FUEL_BURN : 0.0;
		n_fire_events++;
	}
	ev_rem[s][0] = ev_rem[s][1] = ev_kill[s][0] = ev_kill[s][1] = 0.0;
	ev_has[s][0] = ev_has[s][1] = 0;
}

/*--------------------------------------------------------------------------*/
/* Batched FATES daily step (world_daily_F pass 1, after every zone's daily  */
/* weather): pack each site's boundary conditions (soil adaptor + daylight-  */
/* mean forcing; patch state is what patch_daily_F starts from -- nothing   */
/* else touches a patch before its own patch_daily_F) and run all sites'    */
/* flux sub-steps + dynamics in one call (OpenMP over sites in the backend). */
/*--------------------------------------------------------------------------*/
void veg_host_daily_batch(struct command_line_object *command_line, struct date current_date)
{
	int	s, site, n_bad;
	double	resid, tsoil;
	struct veg_time t;
	struct soil_adaptor_out out;
	struct patch_object *patch;
	struct zone_object *zone;

	double	w0, w1;
	if (smap == NULL)
		return;
	w0 = wall_now();
	fill_time(current_date, &t);
	veg_set_time(&t);
	last_julday = julday(current_date);

	for (s = 0; s < nsites; s++) {
		site = s + 1;
		patch = site_patch[s];
		zone = site_zone[s];
		/* patch_daily_F sets patch Tsoil from the zone at its start */
		tsoil = (command_line[0].surface_energy_flag == 0) ? zone[0].metv.tsoil : patch[0].Tsoil;
		memset(&site_bc[s], 0, sizeof(struct veg_bc));
		soil_bc(patch, site, &site_bc[s], &out, tsoil);
		atm_bc(site_basin[s], zone, patch, current_date, &site_bc[s]);

		/*--- soil-adaptor mass check ---*/
		resid = fabs(out.col_water - out.ref_water);
		sum_resid += resid;
		n_site_days++;
		if (resid > VH_MASS_TOL) {
			n_resid_fail++;
			if (n_resid_fail <= 10)
				fprintf(stderr, "FATES host: mass check %ld-%02ld-%02ld patch %d: adaptor %.9f "
					"RHESSys %.9f m\n", current_date.year, current_date.month,
					current_date.day, patch[0].ID, out.col_water, out.ref_water);
		}
		if (resid > max_resid) {
			max_resid = resid;
			max_resid_patch = patch[0].ID;
			max_resid_date = current_date;
		}
		/*--- psi tracking (patch.psi from patch_daily_I vs the adaptor) ---*/
		trk_n += 1.0;
		trk_x += patch[0].psi; trk_y += out.psi_patch_equiv;
		trk_xx += patch[0].psi * patch[0].psi; trk_yy += out.psi_patch_equiv * out.psi_patch_equiv;
		trk_xy += patch[0].psi * out.psi_patch_equiv;
		trk_absdiff += fabs(patch[0].psi - out.psi_patch_equiv);
		trk_maxdiff = max(trk_maxdiff, fabs(patch[0].psi - out.psi_patch_equiv));
		if (fabs(patch[0].psi - out.psi_patch_equiv) <= 0.01)
			trk_close += 1.0;

		n_bad = 0;
		veg_set_bc(&site, &site_bc[s], &n_bad);
		site_nbad[s] = n_bad;
		n_bad_total += n_bad;
		/* mode 3: N uptake RHESSys granted yesterday */
		site_bc[s].n_uptake_nh4 = site_nh4[s];
		site_bc[s].n_uptake_no3 = site_no3[s];
		/* Phase 6: today's TEC disturbance, if any */
		resolve_disturbance(s, &site_bc[s]);
	}
	/*--- all FATES sites: day + night flux sub-steps + daily dynamics ---*/
	w1 = wall_now();
	veg_do_all(&nsites, site_bc, site_out);
	t_pack += w1 - w0;
	t_fates += wall_now() - w1;
	for (s = 0; s < nsites; s++) {
		sum_coh += site_out[s].ncohorts;
		sum_pat += site_out[s].npatches;
	}
	n_batch_sd += nsites;
	/* yearly timing line: FATES wall time and mean cohorts per site over the year */
	{
		static double t_yr0 = 0.0, coh_yr0 = 0.0, sd_yr0 = 0.0;
		if ((current_date.month == 12) && (current_date.day == 31)) {
			fprintf(stderr, "FATES host year %ld: FATES step %.1f s wall, mean %.2f cohorts per site\n",
				current_date.year, t_fates - t_yr0,
				(n_batch_sd > sd_yr0) ? (sum_coh - coh_yr0) / (n_batch_sd - sd_yr0) : 0.0);
			t_yr0 = t_fates; coh_yr0 = sum_coh; sd_yr0 = n_batch_sd;
		}
	}
}

/*--------------------------------------------------------------------------*/
/* Start of patch_daily_F: per-patch use of the batch results.              */
/*--------------------------------------------------------------------------*/
void veg_host_daily(struct world_object *world, struct basin_object *basin,
		    struct zone_object *zone, struct patch_object *patch,
		    struct command_line_object *command_line, struct date current_date)
{
	int	site;
	struct veg_bc bc;
	struct soil_adaptor_out out;

	(void) world; (void) basin; (void) zone;
	site = find_site(patch);
	if (site == 0)
		return;					/* non-vegetated patch */

	/*--- diagnostics: cohort structure of the diagnostic patch, 1 Jan + 1 Jul ---*/
	if ((patch == diag_patch) && (current_date.day == 1)
	    && ((current_date.month == 1) || (current_date.month == 7))) {
		static int first = 1;
		int ymd = (int) (current_date.year * 10000 + current_date.month * 100 + current_date.day);
		size_t fl = strlen(cohort_file);
		veg_dump(&site, &patch[0].ID, &ymd, cohort_file, &fl, &first);
		first = 0;
	}

	/*--- modes 2-3: FATES canopy structure drives RHESSys canopy physics ---*/
	if (command_line[0].demog_flag >= 2)
		apply_structure(patch, &site_out[site - 1]);

	/*--- mode 3: FATES litter -> RHESSys litter pools; FATES N budget ---*/
	if (command_line[0].demog_flag == 3)
		fates_litter_and_budget(patch, site,
			site_bc[site - 1].n_uptake_nh4 + site_bc[site - 1].n_uptake_no3);
	site_tf[site - 1] = 0.0;
	site_tleg[site - 1] = 0.0;

	if ((diag != NULL) && (patch == diag_patch)) {
		/* adaptor details for the row (state unchanged since the batch) */
		soil_bc(patch, site, &bc, &out, site_bc[site - 1].tempk[0] - VH_TFRZ);
		write_diag_row(patch, current_date, &site_bc[site - 1], &out, site_nbad[site - 1]);
	}
}

/*--------------------------------------------------------------------------*/
/* FATES transpiration: RHESSys penman_monteith (as the legacy canopy uses)  */
/* with the FATES daylight canopy conductance (incl. leaf boundary layer)    */
/* and FATES canopy-absorbed shortwave, plus the overstory stratum's net     */
/* longwave and aerodynamic conductance from today's legacy canopy step.     */
/* Daylight only, minus the wet-canopy (rain) time, as canopy_stratum_daily_F*/
/* does. demog_flag 2: rescale every stratum's transpiration demand so the   */
/* patch total equals the FATES value (RHESSys's soil-water supply limits,   */
/* sat/unsat partition and accounting then apply unchanged).                 */
/*--------------------------------------------------------------------------*/
void veg_host_transpiration(struct zone_object *zone, struct patch_object *patch,
			    struct command_line_object *command_line)
{
	int	i, site;
	double	dayl, rain_day, rnet, rate, tf, tleg, f;
	struct veg_out *o;
	struct canopy_strata_object *st, *top;

	site = find_site(patch);
	if (site == 0)
		return;
	o = &site_out[site - 1];
	dayl = zone[0].metv.dayl;
	top = overstory(patch);

	/* legacy demand (cover-weighted, as patch_daily_F sums it) */
	tleg = 0.0;
	for (i = 0; i < patch[0].num_canopy_strata; i++) {
		st = patch[0].canopy_strata[i];
		tleg += st[0].cover_fraction * (st[0].transpiration_unsat_zone + st[0].transpiration_sat_zone);
	}

	tf = 0.0;
	if ((dayl > 1.0) && (top != NULL) && (o[0].gc_canopy > 0.0) && (top[0].ga > 0.0)) {
		rnet = o[0].swabs_canopy + 1000.0 * top[0].Lstar_day / dayl;
		if (rnet > 0.0) {
			rate = penman_monteith(command_line[0].verbose_flag, zone[0].metv.tday,
				zone[0].metv.pa, zone[0].metv.vpd_day, rnet,
				1.0 / o[0].gc_canopy, 1.0 / top[0].ga, 2);		/* m/s */
			rain_day = zone[0].rain_duration * dayl / 86400.0;
			tf = max(rate, 0.0) * max(dayl - rain_day, 0.0);
		}
	}
	site_tf[site - 1] = tf;
	site_tleg[site - 1] = tleg;

	if (command_line[0].demog_flag < 2)
		return;
	if (tleg > ZERO) {
		f = tf / tleg;
		for (i = 0; i < patch[0].num_canopy_strata; i++) {
			st = patch[0].canopy_strata[i];
			st[0].transpiration_unsat_zone *= f;
			st[0].transpiration_sat_zone *= f;
		}
	}
	else if ((tf > 0.0) && (top != NULL) && (top[0].cover_fraction > ZERO)) {
		/* no legacy demand to rescale: the overstory draws it, from the
		   saturated zone if the water table is in the root zone */
		if (patch[0].sat_deficit_z > patch[0].rootzone.depth)
			top[0].transpiration_unsat_zone = tf / top[0].cover_fraction;
		else
			top[0].transpiration_sat_zone = tf / top[0].cover_fraction;
	}
}

static void flush_basin_day(void)
{
	if ((fbasin != NULL) && (b_n > 0) && (b_area > 0.0))
		fprintf(fbasin, "%ld,%ld,%ld,%ld,%.6f,%.6f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.4f,%.3f,%.3e,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
			bdate.year, bdate.month, bdate.day, b_n,
			1000.0 * b_gpp_f / b_area, 1000.0 * b_gpp_l / b_area,		/* gC/m2/day */
			1000.0 * b_tf / b_area, 1000.0 * b_tleg / b_area,		/* mm/day */
			1000.0 * b_tact / b_area, b_lai_f / b_area, b_lai_l / b_area,
			b_btran / b_area, 1000.0 * b_gc / b_area, b_sw / b_area, b_wbmax,
			b_bio / b_area, b_lit / b_area, b_fleaf / b_area, b_fcwd / b_area,
			b_rlit / b_area, b_rlitbg / b_area);
	b_area = b_gpp_f = b_gpp_l = b_tf = b_tleg = b_tact = b_lai_f = b_lai_l = 0.0;
	b_btran = b_gc = b_sw = b_wbmax = b_bio = b_lit = 0.0;
	b_fleaf = b_fcwd = b_rlit = b_rlitbg = 0.0;
	b_n = 0;
}

/*--------------------------------------------------------------------------*/
/* End of patch_daily_F (after the water balance): legacy vs FATES.          */
/*--------------------------------------------------------------------------*/
void veg_host_daily_end(struct zone_object *zone, struct patch_object *patch,
			struct command_line_object *command_line, struct date current_date)
{
	int	i, site;
	long	jd;
	double	gpp_l, lai_l, tact, wb, a;
	struct veg_out *o;
	struct canopy_strata_object *st;

	(void) zone;
	site = find_site(patch);
	if (site == 0)
		return;
	o = &site_out[site - 1];

	gpp_l = lai_l = 0.0;
	for (i = 0; i < patch[0].num_canopy_strata; i++) {
		st = patch[0].canopy_strata[i];
		gpp_l += st[0].cover_fraction * st[0].cdf.psn_to_cpool;
		lai_l += st[0].cover_fraction * st[0].epv.proj_lai;
	}
	tact = patch[0].transpiration_sat_zone + patch[0].transpiration_unsat_zone;
	if (command_line[0].demog_flag == 3) {
		a = legacy_plant_cn(patch);
		if (site_legacy0[site - 1] < 0.0)
			site_legacy0[site - 1] = a;
		leg_max = max(leg_max, fabs(a - site_legacy0[site - 1]));
	}
	wb = fabs(patch[0].water_balance);
	wb_max = max(wb_max, wb);
	wb_sum += wb;
	wb_n++;
	if (wb > 1.0e-6)
		wb_big++;

	/* basin-mean daily row (area-weighted), written when the day changes */
	jd = julday(current_date);
	if (jd != bday) {
		flush_basin_day();
		bday = jd;
		bdate = current_date;
	}
	a = patch[0].area;
	b_area += a;
	b_n++;
	b_gpp_f += a * (o[0].gpp_day + o[0].gpp_night);
	b_gpp_l += a * gpp_l;
	b_tf += a * site_tf[site - 1];
	b_tleg += a * site_tleg[site - 1];
	b_tact += a * tact;
	b_lai_f += a * o[0].tlai;
	b_lai_l += a * lai_l;
	b_btran += a * o[0].btran;
	b_gc += a * o[0].gc_canopy;
	b_sw += a * o[0].swabs_canopy;
	b_bio += a * o[0].biomass_c;
	b_lit += a * o[0].litter_c;
	b_fleaf += a * o[0].leaf_litter_c;
	b_fcwd += a * o[0].ag_cwd_c;
	b_rlit += a * (patch[0].litter_cs.litr1c + patch[0].litter_cs.litr2c
		+ patch[0].litter_cs.litr3c + patch[0].litter_cs.litr4c);
	b_rlitbg += a * (patch[0].litter_cs.litr1c_bg + patch[0].litter_cs.litr2c_bg
		+ patch[0].litter_cs.litr3c_bg + patch[0].litter_cs.litr4c_bg);
	b_wbmax = max(b_wbmax, wb);

	if ((fflux != NULL) && (patch == diag_patch))
		fprintf(fflux, "%ld,%ld,%ld,%d,%d,%d,%.4f,%.4f,%.3f,%.4f,%.4f,%.2f,%.2f,%.6f,%.3f,%.4f,"
			"%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.4f,%.3e,%.5f,%.5f,%.8e,%.6e,%.6e,%.6e,%.6e,%.3e,%.5e\n",
			current_date.year, current_date.month, current_date.day, patch[0].ID,
			o[0].npatches, o[0].ncohorts, o[0].tlai, o[0].elai, o[0].htop, o[0].laisun,
			o[0].laisha, o[0].rssun, o[0].rssha, o[0].gc_canopy, o[0].swabs_canopy,
			o[0].btran, o[0].gpp_day, o[0].gpp_night, o[0].resp_m, site_tf[site - 1],
			site_tleg[site - 1], tact, gpp_l, lai_l, patch[0].water_balance,
			o[0].biomass_c, o[0].litter_c, o[0].n_stock, o[0].n_uptake_recv,
			o[0].lit_n_lab + o[0].lit_n_cel + o[0].lit_n_lig, o[0].mb_n_uptake * 1000.0,
			o[0].mb_n_frag * 1000.0, o[0].mb_n_err, o[0].ncap);
}

/*--------------------------------------------------------------------------*/
void veg_host_finish(struct command_line_object *command_line)
{
	(void) command_line;
	if (diag != NULL) {
		fclose(diag);
		diag = NULL;
	}
	flush_basin_day();
	if (fbasin != NULL) { fclose(fbasin); fbasin = NULL; }
	if (fflux != NULL) { fclose(fflux); fflux = NULL; }
	if (command_line[0].demog_flag == 3)
		fprintf(stderr, "FATES host N/litter (mode 3, area-weighted site means per day): uptake "
			"granted %.4f gN/m2, received by FATES %.4f; litter to RHESSys %.4f gC/m2, %.5f gN/m2; "
			"FATES N budget |dN - (uptake - litter N - efflux N)| max %.3e mean %.3e kgN/m2 (%ld site-days); "
			"legacy strata C+N max change %.3e kg/m2\n",
			tot_area_days ? tot_up / tot_area_days : 0.0, tot_area_days ? tot_recv / tot_area_days : 0.0,
			tot_area_days ? tot_litc / tot_area_days : 0.0, tot_area_days ? tot_litn / tot_area_days : 0.0,
			nb_max, nb_n ? nb_sum / nb_n : 0.0, nb_n, leg_max);
	if (command_line[0].demog_flag == 3)
		fprintf(stderr, "FATES host disturbance (TEC events on FATES sites): %ld site fire events, "
			"%ld site harvests; burned %.2f gC/m2 %.4f gN/m2, harvested wood %.2f gC/m2 %.4f gN/m2 "
			"(site-area-weighted totals over the run); max FATES patch-resize N error %.3e kg/m2; "
			"FATES's own N balance error max %.3e kg/m2, host residual after adding it max %.3e\n",
			n_fire_events, n_harvest_events,
			tot_area_days ? tot_firec / (tot_area_days / max(nb_n / (double) max(nsites, 1), 1.0)) : 0.0,
			tot_area_days ? tot_firen / (tot_area_days / max(nb_n / (double) max(nsites, 1), 1.0)) : 0.0,
			tot_area_days ? tot_harvc / (tot_area_days / max(nb_n / (double) max(nsites, 1), 1.0)) : 0.0,
			tot_area_days ? tot_harvn / (tot_area_days / max(nb_n / (double) max(nsites, 1), 1.0)) : 0.0,
			resize_max, err_max, res_err_max);
	if (command_line[0].demog_flag == 3)
		fprintf(stderr, "FATES host litter split (surface vs below ground) vs FATES merged flux: "
			"max error %.3e g/m2/day; root efflux delivered to below-ground labile litter "
			"%.5f gC/m2/day, %.6f gN/m2/day (site-area mean)\n", split_err_max,
			tot_area_days ? tot_effc / tot_area_days : 0.0, tot_area_days ? tot_effn / tot_area_days : 0.0);
	if (command_line[0].demog_flag == 3)
		fprintf(stderr, "FATES host N budget detail: %ld site-days > 1e-9 kgN/m2; worst patch %d: "
			"dN %.6f recv %.6f litterN %.6f | FATES record uptake %.6f frag %.6f gN/m2; "
			"max |dN - FATES record| %.3e kgN/m2\n", nb_big, nb_max_patch, nb_max_d, nb_max_recv,
			nb_max_lit, nb_max_mbup, nb_max_mbfrag, nb_max_vs_fates);
	fprintf(stderr, "FATES host timing: FATES step (all sites) %.1f s, boundary-condition packing %.1f s; "
		"mean %.2f FATES patches, %.2f cohorts per site\n",
		t_fates, t_pack, (n_batch_sd > 0) ? sum_pat / n_batch_sd : 0.0,
		(n_batch_sd > 0) ? sum_coh / n_batch_sd : 0.0);
	fprintf(stderr, "FATES host water balance (patch_daily_F patch.water_balance, sites): "
		"max |wb| = %.3e m, mean %.3e m, %ld of %ld site-days > 1e-6 m\n",
		wb_max, wb_n ? wb_sum / wb_n : 0.0, wb_big, wb_n);
	fprintf(stderr, "FATES host summary: %d sites, %ld site-days; soil-adaptor mass check "
		"max |resid| = %.3e m", nsites, n_site_days, max_resid);
	if (max_resid_patch >= 0)
		fprintf(stderr, " (patch %d, %ld-%02ld-%02ld)", max_resid_patch,
			max_resid_date.year, max_resid_date.month, max_resid_date.day);
	fprintf(stderr, ", mean %.3e m, %ld site-days > %.0e m; bc_in validation failures: %ld\n",
		n_site_days ? sum_resid / n_site_days : 0.0, n_resid_fail, VH_MASS_TOL, n_bad_total);
	if (trk_n > 1.0) {
		double vx = trk_xx / trk_n - pow(trk_x / trk_n, 2.0);
		double vy = trk_yy / trk_n - pow(trk_y / trk_n, 2.0);
		double cxy = trk_xy / trk_n - (trk_x / trk_n) * (trk_y / trk_n);
		fprintf(stderr, "FATES host psi tracking: patch.psi vs adaptor (same curve, rebuilt "
			"patch.S): r = %.4f, mean |diff| = %.4f MPa, max %.4f MPa, %.1f%% within 0.01 MPa\n",
			(vx > 0.0 && vy > 0.0) ? cxy / sqrt(vx * vy) : 1.0, trk_absdiff / trk_n,
			trk_maxdiff, 100.0 * trk_close / trk_n);
	}
	veg_clean();
	free(smap); smap = NULL;
	free(site_zi); site_zi = NULL;
	free(site_out); site_out = NULL;
	free(site_pid); site_pid = NULL;
	free(site_patch); free(site_zone); free(site_basin); free(site_bc); free(site_nbad);
	free(ev_rem); free(ev_kill); free(ev_has); ev_rem = NULL; ev_kill = NULL; ev_has = NULL;
	site_patch = NULL; site_zone = NULL; site_basin = NULL; site_bc = NULL; site_nbad = NULL;
	free(site_tf); site_tf = NULL;
	free(site_tleg); site_tleg = NULL;
	free(site_nh4); free(site_no3); free(site_nprev); free(site_legacy0);
	site_nh4 = site_no3 = site_nprev = site_legacy0 = NULL;
}
