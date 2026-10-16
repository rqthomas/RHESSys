/*--------------------------------------------------------------*/
/*								*/
/*	compute_stream_temperature				*/
/*								*/
/*	Daily water temperature of one stream reach (-strtemp).	*/
/*								*/
/*	Each day the reach is at the steady state of mixing	*/
/*	plus surface heat exchange (CCR reaches flush within a	*/
/*	day, so reach storage is not carried):			*/
/*								*/
/*	  0 = Q_up T_up + Q_gw T_gw + Q_sub T_sub - Q_in T	*/
/*	      + A_s/rc H(T)					*/
/*								*/
/*	Q_in = Q_up + Q_gw + Q_sub (m3/day), rc = 4186 kJ/m3/K.	*/
/*	H(T) (kJ/m2/day) is linearized once about the mixing	*/
/*	temperature T_mix of the inflows: K = -dH/dT,		*/
/*	T_e = T_mix + H(T_mix)/K, so				*/
/*	  T = (Q_in T_mix + A_s K T_e/rc) / (Q_in + A_s K/rc)	*/
/*	T is kept at >= 0 C (no ice). Lake reaches only mix.	*/
/*								*/
/*	Source temperatures (running means of max(Tair, 0),	*/
/*	M += (T_d - M)(1 - e^(-1/tau))):			*/
/*	  hillslope groundwater   T_gw = M(tau_gw) + dT_gw	*/
/*	  patch streamflow (subsurface and return flow)		*/
/*	                          T_sub = M(tau_sub)		*/
/*								*/
/*	Surface flux, kJ/m2/day:				*/
/*	  H(T) = (1 - a_w)(1 - S) Kdown				*/
/*	       + e_w [(1 - S) Ldown + S e_v s (Ta+273.15)^4]	*/
/*	       - e_w s (T+273.15)^4				*/
/*	       - f(U) (e_s(T) - e_a)		evaporation	*/
/*	       - 0.61 (P/1013.25) f(U) (T - Ta)	sensible	*/
/*	f(U) = mult (596.3 + 29.81 U^2) kJ/m2/day/hPa		*/
/*	(Brady, Graves & Geyer 1969, as in QUAL2K).		*/
/*								*/
/*	Canopy shade from the stream-side (lateral) patches:	*/
/*	  S = min(1, c_oh H_can / w_s) (1 - exp(-k PAI))	*/
/*	PAI = area-weighted strata cover x proj_pai, k = strata	*/
/*	ext_coef, H_can = strata height, w_s = water surface	*/
/*	width.							*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

#define RHO_C_WATER 4186.0     /* kJ/m3/K */
#define ALBEDO_W    0.08       /* water shortwave albedo */
#define EMISS_W     0.97       /* water longwave emissivity */
#define EMISS_VEG   0.97       /* canopy longwave emissivity (canopy at air temperature) */

int stream_temp_is_lake(int reach_ID, struct stream_temp_defaults *p)
{
	int i;
	if ((p == NULL) || (p->num_lake_reaches <= 0)) return(0);
	for (i = 0; i < p->num_lake_reaches; i++)
		if (p->lake_reaches[i] == reach_ID) return(1);
	return(0);
}

static double sat_vp(double T)       /* hPa, Magnus */
{
	return(6.1078 * exp(17.27 * T / (T + 237.3)));
}

/* net surface heat flux into the water (kJ/m2/day) at water temperature T;
   sets dHdT and the components sw, lw (net), evap and sens (+ = into the water) */
static double surface_flux(double T, struct stream_temp_input *in, double S,
			   struct stream_temp_defaults *p, double *dHdT,
			   double *sw, double *lw, double *evap, double *sens)
{
	double sig, TK, TaK, fU, es, ea, Td, P, bow, Kd, Ld;

	sig = SBC * SECONDS_PER_DAY / 1000.0;          /* kJ/m2/day/K4 */
	TK = T + 273.15;
	TaK = in->tavg + 273.15;
	Kd = (in->Kdown > 0.0) ? in->Kdown : 0.0;
	Ld = (in->Ldown > 0.0) ? in->Ldown : 0.80 * sig * pow(TaK, 4.0);   /* missing: clear-ish sky */
	Td = (in->tdew > -100.0) ? ((in->tdew < in->tavg) ? in->tdew : in->tavg) : in->tavg - 2.0;
	fU = p->wind_fn_mult * (596.3 + 29.81 * ((in->wind > 0.0) ? in->wind * in->wind : 0.0));
	es = sat_vp(T);
	ea = sat_vp(Td);
	P = 1013.25 * pow(1.0 - 2.25577e-5 * ((in->z > 0.0) ? in->z : 0.0), 5.25588);
	bow = 0.61 * P / 1013.25;

	*sw = (1.0 - ALBEDO_W) * (1.0 - S) * Kd;
	*lw = EMISS_W * ((1.0 - S) * Ld + S * EMISS_VEG * sig * pow(TaK, 4.0))
		- EMISS_W * sig * pow(TK, 4.0);
	*evap = -fU * (es - ea);
	*sens = -bow * fU * (T - in->tavg);
	*dHdT = -4.0 * EMISS_W * sig * pow(TK, 3.0)
		- fU * es * 17.27 * 237.3 / ((T + 237.3) * (T + 237.3))
		- bow * fU;
	return(*sw + *lw + *evap + *sens);
}

void compute_stream_temperature(struct stream_network_object *reach,
				struct stream_temp_input *in,
				struct stream_temp_defaults *p)
{
	double drive, T_gw, T_sub, Q_lat, heat_lat, Q_in, T_mix, T;
	double bw, w_s, d, A_s, S, f_oh, H, dHdT, Kex, Teq;
	double sw, lw, ev, se;

	/*--------------------------------------------------------------*/
	/*	source temperatures from the running means		*/
	/*--------------------------------------------------------------*/
	drive = (in->tavg > 0.0) ? in->tavg : 0.0;
	reach->T_sub_mem += (drive - reach->T_sub_mem) * (1.0 - exp(-1.0 / p->tau_sub));
	reach->T_gw_mem += (drive - reach->T_gw_mem) * (1.0 - exp(-1.0 / p->tau_gw));
	T_sub = reach->T_sub_mem;
	T_gw = reach->T_gw_mem + p->dT_gw;
	if (T_gw < 0.0) T_gw = 0.0;

	/*--------------------------------------------------------------*/
	/*	mixing of the inflows (m3/day, m3/day x C)		*/
	/*--------------------------------------------------------------*/
	Q_lat = in->Q_gw + in->Q_sub;
	heat_lat = in->Q_gw * T_gw + in->Q_sub * T_sub;
	Q_in = in->Q_up + Q_lat;
	T_mix = (Q_in > 0.0) ? (in->heat_up + heat_lat) / Q_in : in->tavg;

	/*--------------------------------------------------------------*/
	/*	shade: water surface width of the trapezoid holding	*/
	/*	the cross-section A_x					*/
	/*--------------------------------------------------------------*/
	bw = (reach->bottom_width > 0.0) ? reach->bottom_width : 1.0;
	if ((in->tangent > 0.0) && (in->A_x > 0.0)) {
		d = (-bw + sqrt(bw * bw + 4.0 * in->tangent * in->A_x)) / (2.0 * in->tangent);
		w_s = bw + 2.0 * in->tangent * d;
	}
	else w_s = bw;
	A_s = w_s * reach->length;
	f_oh = p->c_overhang * in->h_can / w_s;
	if (f_oh > 1.0) f_oh = 1.0;
	S = f_oh * (1.0 - exp(-in->k_ext * in->pai));

	/*--------------------------------------------------------------*/
	/*	steady state with surface exchange, linearized once	*/
	/*--------------------------------------------------------------*/
	reach->H_sw = reach->H_lw = reach->H_evap = reach->H_sens = 0.0;
	if (reach->temp_is_lake || (Q_in <= 0.0)) T = T_mix;
	else {
		H = surface_flux(T_mix, in, S, p, &dHdT, &sw, &lw, &ev, &se);
		Kex = (-dHdT > 1.0) ? -dHdT : 1.0;
		Teq = T_mix + H / Kex;
		T = (Q_in * T_mix + A_s * Kex * Teq / RHO_C_WATER) / (Q_in + A_s * Kex / RHO_C_WATER);
		surface_flux(T, in, S, p, &dHdT, &reach->H_sw, &reach->H_lw, &reach->H_evap, &reach->H_sens);
	}
	if (T < 0.0) T = 0.0;     /* no ice */

	/*--------------------------------------------------------------*/
	/*	store state and diagnostics				*/
	/*--------------------------------------------------------------*/
	reach->T_mean = T;
	reach->T_gw = T_gw;
	reach->T_sub = T_sub;
	reach->T_lat = (Q_lat > 0.0) ? heat_lat / Q_lat : -9999.0;
	reach->T_up = (in->Q_up > 0.0) ? in->heat_up / in->Q_up : -9999.0;
	reach->Q_lat_gw = in->Q_gw;
	reach->Q_lat_sub = in->Q_sub;
	reach->shade = S;
	reach->pai_rip = in->pai;
	return;
} /*end compute_stream_temperature*/
