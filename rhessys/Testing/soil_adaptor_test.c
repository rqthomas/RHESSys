/*--------------------------------------------------------------------------*/
/* soil_adaptor_test.c -- standalone check of the RHESSys -> FATES soil     */
/* adaptor (design doc 9.4d / 12 pre-work spike). No FATES, no world file.  */
/*                                                                          */
/* build + run (from rhessys/):                                             */
/*   cc -I include Testing/soil_adaptor_test.c fates/rhessys_soil_adaptor.c */
/*      hydro/compute_delta_water.c hydro/compute_z_final.c                 */
/*      hydro/compute_soil_water_potential.c -lm -o /tmp/sa_test && /tmp/sa_test */
/*--------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rhessys.h"
#include "rhessys_soil_adaptor.h"

double	compute_delta_water(int, double, double, double, double, double);
double	compute_soil_water_potential(int, int, double, double, double, double,
		double, double, double, double, double, double);

static int n_fail = 0, n_check = 0;

static void check(int ok, const char *what, double got, double want)
{
	n_check++;
	if (!ok) {
		n_fail++;
		printf("  FAIL %-48s got %.12g want %.12g\n", what, got, want);
	}
}
#define NEAR(a, b, tol) (fabs((a) - (b)) <= (tol) * (1.0 + fabs(b)))

static struct soil_adaptor_params make_par(double p0, double decay, double depth)
{
	struct soil_adaptor_params par;
	par.porosity_0 = p0;
	par.porosity_decay = decay;
	par.soil_depth = depth;
	par.soil_water_cap = compute_delta_water(0, p0, decay, depth, depth, 0.0);
	par.theta_psi_curve = 1;
	par.psi_air_entry = 0.218;
	par.pore_size_index = 0.204;
	par.p3 = 0.0; par.p4 = 0.0;
	par.psi_max = 0.01;
	return(par);
}

/* patch.psi exactly as patch_daily_I computes it */
static double patch_psi(const struct soil_adaptor_params *p, double Tsoil, double S)
{
	return(compute_soil_water_potential(0, p[0].theta_psi_curve, Tsoil,
		-1.0 * p[0].psi_max, -10.0, p[0].psi_air_entry, p[0].pore_size_index,
		p[0].p3, p[0].p4, p[0].porosity_0, p[0].porosity_decay, S));
}

int main(void)
{
	int	j, i, nlev = 3;
	double	zi[SA_NLEV_MAX + 1], pot_sat, b, want, maxres;
	struct soil_adaptor_params par;
	struct soil_adaptor_state st;
	struct soil_adaptor_out out;

	/*--- 1. geometry ---*/
	printf("1. geometry\n");
	check(soil_adaptor_geometry(3, 0.5, 2.0, zi) == 0, "geometry ok", 0, 0);
	check(NEAR(zi[1], 0.5, 1e-15) && NEAR(zi[2], 1.25, 1e-15) && NEAR(zi[3], 2.0, 1e-15),
		"zi = 0,0.5,1.25,2", zi[2], 1.25);
	soil_adaptor_geometry(3, 0.0, 1.5, zi);
	check(NEAR(zi[1], 0.5, 1e-15) && NEAR(zi[3], 1.5, 1e-15), "no roots -> even split", zi[1], 0.5);
	check(soil_adaptor_geometry(SA_NLEV_MAX + 1, 0.5, 2.0, zi) != 0, "nlev too big rejected", 0, 0);
	soil_adaptor_geometry(3, 1.4999, 1.5, zi);	/* CCR: root zone fills the soil */
	check(NEAR(zi[1], 0.5, 1e-15) && NEAR(zi[2], 1.0, 1e-15), "no slivers -> even split", zi[1], 0.5);

	/*--- 2. uniform porosity, water table below root zone ---*/
	printf("2. uniform porosity, water table below root zone\n");
	par = make_par(0.45, 4000.0, 2.0);
	soil_adaptor_geometry(nlev, 0.5, par.soil_depth, zi);
	st.rz_depth = 0.5; st.sat_deficit = 0.6; st.rz_storage = 0.15;
	st.unsat_storage = 0.2; st.Tsoil = 10.0;
	soil_adaptor_profile(&par, &st, nlev, zi, &out);
	check(NEAR(out.zwt, 0.6 / 0.45, 1e-12), "zwt", out.zwt, 0.6 / 0.45);
	check(NEAR(out.S[0], 0.15 / (0.45 * 0.5), 1e-12), "S layer 1 = rootzone.S", out.S[0], 0.15 / 0.225);
	pot_sat = 0.45 * 0.5;
	want = 0.2 / (0.6 - pot_sat);			/* RHESSys patch.S */
	check(NEAR(out.S[1], want, 1e-12), "S layer 2 = patch.S", out.S[1], want);
	want = (want * 0.45 * (0.6 / 0.45 - 1.25) + 0.45 * (2.0 - 0.6 / 0.45)) / (0.45 * 0.75);
	check(NEAR(out.S[2], want, 1e-12), "S layer 3 (unsat + saturated)", out.S[2], want);
	check(NEAR(out.col_water, out.ref_water, 1e-12), "mass: sum(theta dz) = RHESSys water", out.col_water, out.ref_water);
	check(NEAR(out.watsat[0], 0.45, 1e-12), "watsat uniform", out.watsat[0], 0.45);
	/* layer 2 lies wholly in the unsaturated zone: same psi as patch.psi */
	/* patch_daily_I: patch.S = (rz + unsat) / sat_deficit, then patch.psi */
	check(NEAR(out.S_patch_equiv, 0.35 / 0.6, 1e-12), "S_patch_equiv = (rz+unsat)/sat_deficit", out.S_patch_equiv, 0.35 / 0.6);
	check(NEAR(out.psi_patch_equiv, patch_psi(&par, st.Tsoil, 0.35 / 0.6), 1e-12),
		"psi_patch_equiv = patch.psi", out.psi_patch_equiv, patch_psi(&par, st.Tsoil, 0.35 / 0.6));
	check(NEAR(out.psi_patch[1], patch_psi(&par, st.Tsoil, 0.2 / (0.6 - pot_sat)), 1e-12),
		"layer-2 psi == patch.psi", out.psi_patch[1], patch_psi(&par, st.Tsoil, 0.2 / (0.6 - pot_sat)));

	/*--- 3. suction curve ---*/
	printf("3. suction (Clapp-Hornberger, RHESSys units)\n");
	b = 1.0 / par.pore_size_index;
	want = -1000.0 * par.psi_air_entry * pow(out.S[1], -b);
	check(NEAR(out.smp_mm[1], want, 1e-10), "smp = -psi_ae S^-b (mm)", out.smp_mm[1], want);
	st.sat_deficit = 0.0; st.rz_storage = 0.0; st.unsat_storage = 0.0;
	soil_adaptor_profile(&par, &st, nlev, zi, &out);
	check(NEAR(out.smp_mm[0], -1000.0 * par.psi_air_entry, 1e-12), "saturated -> air entry", out.smp_mm[0], -218.0);
	st.sat_deficit = 0.9; st.rz_storage = 0.0; st.unsat_storage = 0.0;	/* bone dry above WT */
	soil_adaptor_profile(&par, &st, nlev, zi, &out);
	/* S floored at 0.01; the curve there (-1.4e12 mm) is below the -1000 MPa
	   floor, i.e. -1e8 mm = CLM's smpmin */
	want = max(-1000.0 * par.psi_air_entry * pow(0.01, -b), -1.0e8);
	check(NEAR(out.smp_mm[0], want, 1e-10), "dry: S>=0.01, smp>=-1e8 mm (CLM smpmin)", out.smp_mm[0], want);
	st.Tsoil = -1.0; st.sat_deficit = 0.6; st.rz_storage = 0.15; st.unsat_storage = 0.2;
	soil_adaptor_profile(&par, &st, nlev, zi, &out);
	want = -2000.0 * par.psi_air_entry * pow(out.S[1], -b);
	check(NEAR(out.smp_mm[1], want, 1e-10), "frozen soil doubles suction (RHESSys)", out.smp_mm[1], want);

	/*--- 4. water table inside the root zone ---*/
	printf("4. water table inside the root zone\n");
	st.Tsoil = 10.0; st.rz_depth = 1.0; st.sat_deficit = 0.2;
	st.rz_storage = 0.08; st.unsat_storage = 0.0;
	soil_adaptor_geometry(nlev, 1.0, par.soil_depth, zi);
	soil_adaptor_profile(&par, &st, nlev, zi, &out);
	check(NEAR(out.S_rz, 0.08 / 0.2, 1e-12), "S above WT = rz / sat_deficit", out.S_rz, 0.4);
	check(out.S_unsat < 0.0, "no unsat zone", out.S_unsat, -1);
	check(NEAR(out.col_water, out.ref_water, 1e-12), "mass", out.col_water, out.ref_water);

	/*--- 4a. over-full unsaturated zone spills into the root zone ---*/
	printf("4a. over-full zone spills, mass kept\n");
	par = make_par(0.45, 4000.0, 2.0);
	soil_adaptor_geometry(nlev, 0.5, par.soil_depth, zi);
	st.Tsoil = 10.0; st.rz_depth = 0.5; st.sat_deficit = 0.27;	/* zwt 0.6: unsat pore 0.045 */
	st.rz_storage = 0.10; st.unsat_storage = 0.06;
	soil_adaptor_profile(&par, &st, nlev, zi, &out);
	check(NEAR(out.S_unsat, 1.0, 1e-12), "unsat zone capped at S = 1", out.S_unsat, 1.0);
	check(NEAR(out.S_rz, (0.10 + 0.015) / 0.225, 1e-12), "excess moved to root zone", out.S_rz, 0.115 / 0.225);
	check(NEAR(out.S_patch_equiv, 0.16 / 0.27, 1e-12), "S_patch_equiv unchanged by spill", out.S_patch_equiv, 0.16 / 0.27);
	check(NEAR(out.col_water, out.ref_water, 1e-12), "mass", out.col_water, out.ref_water);
	for (j = 0; j < nlev; j++)
		check(out.theta[j] <= out.watsat[j] + 1e-12, "theta <= watsat after spill", out.theta[j], out.watsat[j]);

	/*--- 4b. water table below the soil bottom (sat_deficit > soil_water_cap) ---*/
	printf("4b. water table below the soil bottom\n");
	par = make_par(0.485, 4000.0, 1.5);		/* CCR-like: cap 0.7275 */
	soil_adaptor_geometry(nlev, 1.4999, par.soil_depth, zi);
	st.rz_depth = 1.4999; st.sat_deficit = 0.732183602; st.rz_storage = 0.485453399;
	st.unsat_storage = 0.0; st.Tsoil = 5.0;
	soil_adaptor_profile(&par, &st, nlev, zi, &out);
	check(NEAR(out.zwt, par.soil_depth, 1e-12), "zwt clamped to soil depth", out.zwt, par.soil_depth);
	check(NEAR(out.ref_water, 0.485453399, 1e-12), "soil water = rz + unsat", out.ref_water, 0.485453399);
	check(NEAR(out.col_water, out.ref_water, 1e-12), "mass", out.col_water, out.ref_water);

	/*--- 4c. ponded (sat_deficit < 0) ---*/
	printf("4c. ponded\n");
	st.sat_deficit = -0.01; st.rz_storage = 0.0; st.unsat_storage = 0.0;
	soil_adaptor_profile(&par, &st, nlev, zi, &out);
	check(NEAR(out.col_water, par.soil_water_cap, 1e-12), "ponded: column saturated", out.col_water, par.soil_water_cap);
	check(NEAR(out.col_water, out.ref_water, 1e-12), "mass (ponded water excluded)", out.col_water, out.ref_water);

	/*--- 5. exponential porosity, random sweep ---*/
	printf("5. exponential porosity, 20000 random states x 3 profiles\n");
	srand(12345);
	maxres = 0.0;
	for (i = 0; i < 3; i++) {
		double decays[3] = {0.5, 1.2, 4000.0};
		par = make_par(0.5, decays[i], 1.5);
		soil_adaptor_geometry(nlev, 0.4, par.soil_depth, zi);
		for (j = 0; j < 20000; j++) {
			double u = rand() / (double) RAND_MAX, pa, pb;
			st.rz_depth = 0.05 + 1.3 * (rand() / (double) RAND_MAX);
			st.sat_deficit = (1.05 * u - 0.02) * par.soil_water_cap;	/* incl. ponded and below-bottom */
			st.Tsoil = 5.0;
			/* storages within the pore space they occupy */
			pa = compute_delta_water(0, par.porosity_0, par.porosity_decay, par.soil_depth,
				min(st.rz_depth, 1.5), 0.0);
			pb = min(st.sat_deficit, par.soil_water_cap) - pa;	/* pore space below the root zone */
			if (st.sat_deficit <= 0.0) {		/* ponded: nothing above the WT */
				st.rz_storage = 0.0;
				st.unsat_storage = 0.0;
			} else if (pb > 0) {
				st.rz_storage = pa * (rand() / (double) RAND_MAX);
				st.unsat_storage = pb * (rand() / (double) RAND_MAX);
			} else {
				st.rz_storage = st.sat_deficit * (rand() / (double) RAND_MAX);
				st.unsat_storage = 0.0;
			}
			soil_adaptor_profile(&par, &st, nlev, zi, &out);
			maxres = max(maxres, fabs(out.col_water - out.ref_water));
			check(out.S[0] >= 0.0 && out.S[0] <= 1.0 + 1e-12 && out.smp_mm[0] < 0.0,
				"0<=S<=1, smp<0", out.S[0], 0.5);
		}
	}
	check(maxres < 1e-12, "max |mass residual| (m)", maxres, 0.0);
	printf("   max |sum(theta dz) - RHESSys water| = %.3e m\n", maxres);

	printf("%s: %d checks, %d failures\n", n_fail ? "FAIL" : "PASS", n_check, n_fail);
	return(n_fail ? 1 : 0);
}
