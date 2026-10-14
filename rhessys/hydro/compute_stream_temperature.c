/*--------------------------------------------------------------*/
/*								*/
/*	compute_stream_temperature				*/
/*								*/
/*	Daily water temperature of one stream reach (-strtemp).	*/
/*								*/
/*	The reach is one well-mixed water column of volume V	*/
/*	(alfa Q^0.6 L, the same V as -strbgc):			*/
/*								*/
/*	V dT/dt = Q_up T_up + sum Q_k T_k - Q_in T		*/
/*	        + A_s/rc H(T) + A_b/rc k_bed (T_gw - T)		*/
/*								*/
/*	Q_in = Q_up + sum Q_k (m3/day), rc = 4186 kJ/m3/K, so	*/
/*	with no exchange T is the flow-weighted mean of the	*/
/*	inputs. H(T) (kJ/m2/day) is linearized about T_l,	*/
/*	K = -dH/dT, T_e = T_l + H(T_l)/K, and the linear ODE	*/
/*	is solved exactly over the day (stable for any V):	*/
/*	  T* = a/b, b = [Q_in + (A_s K + A_b k_bed)/rc] / V	*/
/*	  T_end = T* + (T0 - T*) e^-b				*/
/*	  T_mean = T* + (T0 - T*)(1 - e^-b)/b			*/
/*	twice, T_l = T0 then T_l = T_mean. Below V_min or for	*/
/*	b > 50 the reach is at its steady state T*. T is kept	*/
/*	at >= 0 C (no ice); the heat removed is H_ice.		*/
/*								*/
/*	Lateral source temperatures:				*/
/*	  hillslope groundwater  T_gw = running mean (tau_gw)	*/
/*	                         + dT_gw, or T_gw_const	*/
/*	  patch subsurface flow  T_sub = running mean (tau_sub)	*/
/*	  return flow            rain/snowmelt of the day at T_P	*/
/*	                         (from the caller), the rest	*/
/*	                         exfiltrated soil water at T_sub */
/*	Running means M += (T_d - M)(1 - e^(-1/tau)) are driven	*/
/*	by T_d = max(Tair, 0), and by 0 C on the snow-covered	*/
/*	share (stream_T_snow_insulate).				*/
/*								*/
/*	Surface flux H(T), daily, kJ/m2/day:			*/
/*	  (1 - a_w)(1 - S) Kdown				*/
/*	+ e_w [(1 - S_v) Ldown + S_v e_v s (Ta+273.15)^4]	*/
/*	- e_w s (T+273.15)^4					*/
/*	- f(U) (e_s(T) - e_a)		evaporation		*/
/*	- 0.61 (P/1013.25) f(U) (T - Ta)	sensible (Bowen)	*/
/*	f(U) = mult (19.0 + 0.95 U^2) cal/cm2/d/mmHg		*/
/*	     = mult (596.3 + 29.81 U^2) kJ/m2/day/hPa		*/
/*	(Brady, Graves & Geyer 1969, as in QUAL2K), U sheltered	*/
/*	by the canopy: U (1 - c_ws S_v).			*/
/*								*/
/*	Shading from the riparian (lateral-input) patches:	*/
/*	  S_v = f_oh (1 - exp(-k PAI)),				*/
/*	  f_oh = min(1, c_oh H_can / w_s),			*/
/*	  PAI = area-weighted sum of strata cover x proj_pai	*/
/*	        (+ stream_SAI), k = strata ext_coef or		*/
/*	        stream_k_shade, H_can = strata height,		*/
/*	  S = 1 - (1 - S_topo)(1 - S_v)				*/
/*	(shade_override >= 0 from the reach file replaces S_v).	*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

#define RHO_C_WATER 4186.0     /* kJ/m3/K */

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
static double surface_flux(double T, struct stream_temp_input *in, double S, double S_v,
			   struct stream_temp_defaults *p, double *dHdT,
			   double *sw, double *lw, double *evap, double *sens)
{
	double sig, TK, TaK, U, fU, es, ea, Td, P, bow, Kd, Ld;

	sig = SBC * SECONDS_PER_DAY / 1000.0;          /* kJ/m2/day/K4 */
	TK = T + 273.15;
	TaK = in->tavg + 273.15;
	Kd = (in->Kdown > 0.0) ? in->Kdown : 0.0;
	Ld = (in->Ldown > 0.0) ? in->Ldown : 0.80 * sig * pow(TaK, 4.0);   /* missing: clear-ish sky */
	Td = (in->tdew > -100.0) ? ((in->tdew < in->tavg) ? in->tdew : in->tavg) : in->tavg - 2.0;
	U = (in->wind > 0.0) ? in->wind * (1.0 - p->c_wind_shelter * S_v) : 0.0;
	fU = p->wind_fn_mult * (596.3 + 29.81 * U * U);
	es = sat_vp(T);
	ea = sat_vp(Td);
	P = 1013.25 * pow(1.0 - 2.25577e-5 * ((in->z > 0.0) ? in->z : 0.0), 5.25588);
	bow = 0.61 * P / 1013.25;

	*sw = (1.0 - p->albedo_w) * (1.0 - S) * Kd;
	*lw = p->emiss_w * ((1.0 - S_v) * Ld + S_v * p->emiss_veg * sig * pow(TaK, 4.0))
		- p->emiss_w * sig * pow(TK, 4.0);
	*evap = -fU * (es - ea);
	*sens = -bow * fU * (T - in->tavg);
	*dHdT = -4.0 * p->emiss_w * sig * pow(TK, 3.0)
		- fU * es * 17.27 * 237.3 / ((T + 237.3) * (T + 237.3))
		- bow * fU;
	return(*sw + *lw + *evap + *sens);
}

void compute_stream_temperature(struct stream_network_object *reach,
				struct stream_temp_input *in,
				struct stream_temp_defaults *p)
{
	double drive, T_gw, T_sub, Q_lat, heat_lat, Q_in, heat_adv;
	double bw, w_s, d, A_s, A_b, S_v, S, k, pai, f_oh, kb;
	double T0, Tl, T_mean, T_end, Tstar, num, den, b, e, H, dHdT, Kex, Teq;
	double sw, lw, ev, se;
	int pass, exch;

	/*--------------------------------------------------------------*/
	/*	source temperatures from the running means		*/
	/*--------------------------------------------------------------*/
	drive = (in->tavg > 0.0) ? in->tavg : 0.0;
	if (p->snow_insulate) drive *= (1.0 - in->snow_frac);
	reach->T_sub_mem += (drive - reach->T_sub_mem) * (1.0 - exp(-1.0 / p->tau_sub));
	reach->T_gw_mem += (drive - reach->T_gw_mem) * (1.0 - exp(-1.0 / p->tau_gw));
	T_sub = (reach->T_sub_mem > 0.0) ? reach->T_sub_mem : 0.0;
	T_gw = (p->T_gw_const > -998.0) ? p->T_gw_const : reach->T_gw_mem + p->dT_gw;
	if (T_gw < 0.0) T_gw = 0.0;

	/*--------------------------------------------------------------*/
	/*	advected heat (m3/day x C)				*/
	/*--------------------------------------------------------------*/
	Q_lat = in->Q_gw + in->Q_sub + in->Q_ret_event + in->Q_ret_soil;
	heat_lat = in->Q_gw * T_gw + (in->Q_sub + in->Q_ret_soil) * T_sub + in->heat_ret_event;
	Q_in = in->Q_up + Q_lat;
	heat_adv = in->heat_up + heat_lat;

	/*--------------------------------------------------------------*/
	/*	geometry: water surface width of the trapezoid holding	*/
	/*	the cross-section A_x					*/
	/*--------------------------------------------------------------*/
	bw = (reach->bottom_width > 0.0) ? reach->bottom_width : 1.0;
	if ((in->tangent > 0.0) && (in->A_x > 0.0)) {
		d = (-bw + sqrt(bw * bw + 4.0 * in->tangent * in->A_x)) / (2.0 * in->tangent);
		w_s = bw + 2.0 * in->tangent * d;
	}
	else w_s = bw;
	A_s = w_s * reach->length;
	A_b = bw * reach->length;

	/*--------------------------------------------------------------*/
	/*	shading							*/
	/*--------------------------------------------------------------*/
	pai = in->pai + p->SAI;
	if (reach->shade_override >= 0.0) S_v = reach->shade_override;
	else {
		k = (p->k_shade >= 0.0) ? p->k_shade : in->k_ext;
		f_oh = reach->c_oh * in->h_can / w_s;
		if (f_oh > 1.0) f_oh = 1.0;
		S_v = f_oh * (1.0 - exp(-k * pai));
	}
	S = 1.0 - (1.0 - reach->S_topo) * (1.0 - S_v);

	/*--------------------------------------------------------------*/
	/*	exact daily step, linearized twice			*/
	/*--------------------------------------------------------------*/
	exch = (p->exchange != 0) && (reach->temp_is_lake == 0);
	kb = exch ? p->k_bed : 0.0;
	T0 = reach->T_water;
	Tl = T0;
	T_mean = T_end = T0;
	for (pass = 0; pass < 2; pass++) {
		Kex = 0.0; Teq = 0.0;
		if (exch) {
			H = surface_flux(Tl, in, S, S_v, p, &dHdT, &sw, &lw, &ev, &se);
			Kex = (-dHdT > 1.0) ? -dHdT : 1.0;
			Teq = Tl + H / Kex;
		}
		num = heat_adv + (A_s * Kex * Teq + A_b * kb * T_gw) / RHO_C_WATER;   /* m3/day x C */
		den = Q_in + (A_s * Kex + A_b * kb) / RHO_C_WATER;                   /* m3/day */
		if (den <= 0.0) break;            /* no inflow and no exchange: T unchanged */
		Tstar = num / den;
		if ((in->V < p->V_min) || (den / in->V > 50.0)) {
			T_end = T_mean = Tstar;
		}
		else {
			b = den / in->V;
			e = exp(-b);
			T_end = Tstar + (T0 - Tstar) * e;
			T_mean = Tstar + (T0 - Tstar) * (1.0 - e) / b;
		}
		if (!exch) break;                 /* linear already */
		Tl = T_mean;
	}
#ifdef STREAM_HEAT_BUDGET
	/* storage change = advection + exchange over the day (both in m3 x C) */
	if ((den > 0.0) && (in->V >= p->V_min) && (den / in->V <= 50.0)) {
		double resid = in->V * (T_end - T0) - (num - den * T_mean);
		if (fabs(resid) > 1e-6 * (fabs(num) + den * fabs(T_mean) + 1.0))
			fprintf(stderr, "STREAM_HEAT_BUDGET reach %d residual %e\n", reach->reach_ID, resid);
	}
#endif

	/*--------------------------------------------------------------*/
	/*	no ice: heat taken to keep T >= 0, per m2 of surface	*/
	/*--------------------------------------------------------------*/
	reach->H_ice = 0.0;
	if ((T_mean < 0.0) || (T_end < 0.0)) {
		reach->H_ice = RHO_C_WATER * (Q_in * ((T_mean < 0.0) ? -T_mean : 0.0)
			+ in->V * ((T_end < 0.0) ? -T_end : 0.0)) / A_s;
		if (T_mean < 0.0) T_mean = 0.0;
		if (T_end < 0.0) T_end = 0.0;
	}

	/*--------------------------------------------------------------*/
	/*	store state and diagnostics				*/
	/*--------------------------------------------------------------*/
	reach->T_water = T_end;
	reach->T_mean = T_mean;
	reach->T_gw = T_gw;
	reach->T_sub = T_sub;
	reach->T_lat = (Q_lat > 0.0) ? heat_lat / Q_lat : -9999.0;
	reach->T_up = (in->Q_up > 0.0) ? in->heat_up / in->Q_up : -9999.0;
	reach->Q_lat_gw = in->Q_gw;
	reach->Q_lat_sub = in->Q_sub + in->Q_ret_soil;
	reach->Q_lat_ret = in->Q_ret_event;
	reach->shade = S;
	reach->pai_rip = pai;
	if (exch) {
		surface_flux(T_mean, in, S, S_v, p, &dHdT, &sw, &lw, &ev, &se);
		reach->H_sw = sw; reach->H_lw = lw; reach->H_evap = ev; reach->H_sens = se;
		reach->H_bed = kb * (T_gw - T_mean);
	}
	else {
		reach->H_sw = reach->H_lw = reach->H_evap = reach->H_sens = reach->H_bed = 0.0;
	}
	return;
} /*end compute_stream_temperature*/
