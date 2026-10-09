/*--------------------------------------------------------------*/
/*								*/
/*	compute_stream_bgc					*/
/*								*/
/*	In-stream storage and processing for one reach and day	*/
/*	(only with -strbgc). Loads in kg/day, stores in kg.	*/
/*								*/
/*	Water column: each reach is a well-mixed store of water	*/
/*	volume V = A x length (A = alfa Q^0.6, the kinematic-wave	*/
/*	cross-section) holding NO3, NH4, DON, DOC and suspended	*/
/*	POC/PON (state). Each day a store with start mass M0	*/
/*	receives its inputs I as a steady inflow and loses mass to	*/
/*	outflow (rate r = Q/V) and a first-order process (rate k)	*/
/*	- exact solution of dM/dt = I - (r + k) M over the day:	*/
/*	  M1 = M0 e^-(r+k) + I (1 - e^-(r+k)) / (r+k)		*/
/*	  loss M0 + I - M1, split r : k into outflow : process.	*/
/*	High flow (V/Q << 1 day) passes everything through the	*/
/*	same day; slow reaches hold and process mass for days.	*/
/*	Process rates k = v fT A_bed / V (v uptake/settling		*/
/*	velocity, m/day; A_bed = bottom width x length):		*/
/*	  suspended POC/PON settle (v_dep) to the streambed store;	*/
/*	  DOC/DON decompose (vf_DOM); NO3 denitrifies (vf_denit).	*/
/*	Streambed store (benthic POC/PON): leaf litter in the	*/
/*	channel; breakdown (k_frag) - a share 1 - f_mic to DOC/DON	*/
/*	(fragmentation/leaching) and f_mic decomposed in place by	*/
/*	attached microbes (below) - and entrainment			*/
/*	e_max (Q/Q_bf)^p to suspended POC are inputs to the water	*/
/*	column; settled POC and microbial biomass return to it.	*/
/*	Attached microbes (leaf packs): keep CUE x C at CN_mic on	*/
/*	the bed, respire the rest; N surplus -> NH4, N deficit	*/
/*	immobilized from NH4 then NO3 in the reach water out of	*/
/*	what the DOM microbes leave of the same uptake share; if	*/
/*	short, only that part is processed (Bed_N_limit) and the	*/
/*	rest stays on the bed. N-poor autumn leaves (C:N >		*/
/*	CN_mic/CUE) draw N from the water; N-rich bed releases it.	*/
/*	Decomposed DOM: microbes keep CUE x C as biomass at C:N	*/
/*	CN_mic (-> streambed), respire the rest; N surplus -> NH4	*/
/*	store; N deficit immobilized from NH4 then NO3 in the reach	*/
/*	water (at most the uptake share kN/(r+kN) of what passed	*/
/*	through, kN = vf_N fT A_bed/V), otherwise			*/
/*	the decomposition rate is scaled (N_limit) and the DOC/DON	*/
/*	stores recomputed, so unused DOM flows on.			*/
/*	fT = theta^(Tstream - 20), Tstream =			*/
/*	max(2.99 + 0.72 Tair, 0).					*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"

#define STREAM_TINY 1.0e-9

/* one water-column store over a day: start mass *Ms, input Iin (kg/day), outflow
   rate r and process rate k (1/day); r < 0 means no water left to hold mass
   (everything flows out). Returns outflow; *proc gets the processed mass. */
static double store_step(double *Ms, double Iin, double r, double k, double *proc)
{
	double lam, e, M1, loss;
	*proc = 0.0;
	if (r < 0.0) {			/* flush: no storage */
		double out = *Ms + Iin;
		*Ms = 0.0;
		return(out);
	}
	lam = r + k;
	if (lam <= STREAM_TINY) {	/* no outflow, no processing: hold */
		*Ms += Iin;
		return(0.0);
	}
	e = exp(-lam);
	M1 = (*Ms) * e + Iin * (1.0 - e) / lam;
	if (M1 < 0.0) M1 = 0.0;
	loss = *Ms + Iin - M1;
	if (loss < 0.0) loss = 0.0;
	*Ms = M1;
	*proc = loss * k / lam;
	return(loss * r / lam);
}

/* immobilize amount (kg N) from the reach water: NH4 first, then NO3, out of the
   outflow and the store in proportion (uptake lowers both what leaves and what stays) */
static void take_N(struct stream_network_object *reach, double *outNH4, double *outNO3, double amount)
{
	double tNH4 = *outNH4 + reach->wc_NH4, tNO3 = *outNO3 + reach->wc_NO3, f, take;
	if (amount <= 0.0) return;
	take = (amount < tNH4) ? amount : tNH4;
	if ((tNH4 > 0.0) && (take > 0.0)) { f = 1.0 - take / tNH4; *outNH4 *= f; reach->wc_NH4 *= f; }
	take = amount - take;
	if (take > tNO3) take = tNO3;
	if ((tNO3 > 0.0) && (take > 0.0)) { f = 1.0 - take / tNO3; *outNO3 *= f; reach->wc_NO3 *= f; }
}

/* 1 if the reach is listed as inside a lake/reservoir (no in-stream processing) */
int stream_bgc_is_lake(int reach_ID, struct stream_bgc_defaults *p)
{
	int i;
	if ((p == NULL) || (p->num_lake_reaches <= 0)) return(0);
	for (i = 0; i < p->num_lake_reaches; i++)
		if (p->lake_reaches[i] == reach_ID) return(1);
	return(0);
}

void compute_stream_bgc(struct stream_network_object *reach,
			double NO3, double NH4, double DON, double DOC,   /* water-column inputs (upstream + lateral), kg/day */
			double POC_in, double PON_in,                    /* suspended POC/PON inputs (upstream + eroded sediment) */
			double lat_POC, double lat_PON,                  /* leaf litter into the channel (to the streambed) */
			double Q, double V, double T_air,                /* outflow m3/day, water volume m3, air temperature C */
			struct stream_bgc_defaults *p)
{
	double Tstream, fT, A_bed, r, kdep, kdec, kden, kN, fN;
	double fragC, fragN, fexp, entC, entN, rr;
	double outNO3, outNH4, outDON, outDOC, outPOC, outPON;
	double depC, depN, decC, decN, den, dummy;
	double netN, need, supply, nlim, immob, bioC, bioN, M0DOC, M0DON;
	double brkC, brkN, micC, micN, supply_tot, netNb, needb, availb, fb;

	Tstream = T_air * 0.72 + 2.99;
	if (Tstream < 0.0) Tstream = 0.0;
	fT = pow(p->theta, Tstream - 20.0);
	A_bed = reach->bottom_width * reach->length;
	reach->V_water = V;

	/* water-column turnover and process rates (1/day) */
	if (Q <= STREAM_TINY) r = 0.0;                /* no flow: hold */
	else if (V <= STREAM_TINY) r = -1.0;          /* flow but no stored water: pass through */
	else r = Q / V;
	if (V > STREAM_TINY) {
		kdep = p->v_dep * fT * A_bed / V;
		kdec = p->vf_DOM * fT * A_bed / V;
		kden = p->vf_denit * fT * A_bed / V;
		kN = p->vf_N * fT * A_bed / V;
	} else kdep = kdec = kden = kN = 0.0;

	/* 1. streambed store: leaf litter in, fragmentation and entrainment out */
	reach->benthic_POC += lat_POC;
	reach->benthic_PON += lat_PON;
	brkC = p->k_frag * fT * reach->benthic_POC;
	if (brkC > reach->benthic_POC) brkC = reach->benthic_POC;
	brkN = (reach->benthic_POC > 0.0) ? brkC * reach->benthic_PON / reach->benthic_POC : 0.0;
	if (brkN > reach->benthic_PON) brkN = reach->benthic_PON;
	if (brkN < 0.0) brkN = 0.0;
	/* share f_mic is set aside for the attached microbes (step 4b), the rest fragments to DOM */
	micC = p->f_mic * brkC;  micN = p->f_mic * brkN;
	fragC = brkC - micC;     fragN = brkN - micN;
	reach->benthic_POC -= brkC;  reach->benthic_PON -= brkN;
	rr = (reach->Q_bf > 0.0) ? ((Q < reach->Q_bf) ? Q : reach->Q_bf) / reach->Q_bf : 0.0;
	fexp = p->e_max * pow(rr, p->export_p);
	if (fexp > 1.0) fexp = 1.0;
	entC = fexp * reach->benthic_POC;  entN = fexp * reach->benthic_PON;
	reach->benthic_POC -= entC;  reach->benthic_PON -= entN;

	/* 2. water-column stores: inflow, outflow and first-order processes */
	outPOC = store_step(&(reach->wc_POC), POC_in + entC, r, kdep, &depC);
	outPON = store_step(&(reach->wc_PON), PON_in + entN, r, kdep, &depN);
	M0DOC = reach->wc_DOC;  M0DON = reach->wc_DON;
	outDOC = store_step(&(reach->wc_DOC), DOC + fragC, r, kdec, &decC);
	outDON = store_step(&(reach->wc_DON), DON + fragN, r, kdec, &decN);
	outNO3 = store_step(&(reach->wc_NO3), NO3, r, kden, &den);
	outNH4 = store_step(&(reach->wc_NH4), NH4, r, 0.0, &dummy);

	/* 3. settled POC/PON to the streambed */
	reach->benthic_POC += depC;  reach->benthic_PON += depN;

	/* N available to microbes (DOM and streambed) during the day: the uptake share of
	   all NH4 + NO3 in the reach water (what flowed out plus what is stored), uptake
	   (rate kN) competing with outflow like the other processes */
	if ((r < 0.0) || (V <= STREAM_TINY)) fN = 0.0;
	else if (r + kN > STREAM_TINY) fN = kN / (r + kN) * (1.0 - exp(-(r + kN)));
	else fN = 0.0;
	supply_tot = fN * (outNH4 + reach->wc_NH4 + outNO3 + reach->wc_NO3);

	/* 4. decomposed DOM: microbial stoichiometry, N limitation */
	netN = decN - p->CUE * decC / p->CN_mic;
	nlim = 1.0;  immob = 0.0;  reach->N_mineral = 0.0;
	if (netN >= 0.0) {
		reach->wc_NH4 += netN;              /* net mineralization */
		reach->N_mineral = netN;
	} else {
		need = -netN;
		supply = supply_tot;
		if (supply < need) {
			/* N-limited: slow decomposition (rate kdec x nlim) and redo the
			   DOC/DON stores, so DOM that microbes cannot use flows on
			   (a few passes: the slower rate changes the decomposed amount) */
			int it;
			double kd;
			nlim = (need > 0.0) ? supply / need : 1.0;
			for (it = 0; it < 4; it++) {
				kd = kdec * nlim;
				reach->wc_DOC = M0DOC;  reach->wc_DON = M0DON;
				outDOC = store_step(&(reach->wc_DOC), DOC + fragC, r, kd, &decC);
				outDON = store_step(&(reach->wc_DON), DON + fragN, r, kd, &decN);
				need = p->CUE * decC / p->CN_mic - decN;
				if ((need <= supply) || (need <= 0.0)) break;
				nlim *= supply / need;
			}
			if (need > supply) {   /* remaining small excess: return it to the store */
				double f = (need > 0.0) ? supply / need : 1.0;
				reach->wc_DOC += (1.0 - f) * decC;  reach->wc_DON += (1.0 - f) * decN;
				decC *= f;  decN *= f;  need = supply;
				nlim *= f;
			}
			if (need < 0.0) need = 0.0;
		}
		immob = need;
		take_N(reach, &outNH4, &outNO3, immob);
	}
	bioC = p->CUE * decC;
	bioN = bioC / p->CN_mic;               /* = decN + immob - mineral */
	reach->benthic_POC += bioC;  reach->benthic_PON += bioN;

	/* 4b. streambed breakdown by attached microbes (share f_mic of step 1) */
	reach->Bed_N_mineral = 0.0;  reach->Bed_N_immob = 0.0;  reach->Bed_N_limit = 1.0;
	netNb = micN - p->CUE * micC / p->CN_mic;
	if (netNb >= 0.0) {
		reach->wc_NH4 += netNb;            /* N-rich bed material: net mineralization */
		reach->Bed_N_mineral = netNb;
	} else {
		needb = -netNb;
		availb = supply_tot - immob;       /* what the DOM microbes left of the uptake share */
		if (availb < 0.0) availb = 0.0;
		fb = (needb > availb) ? availb / needb : 1.0;
		if (fb < 1.0) {                    /* N-limited: the unprocessed part stays on the bed */
			reach->benthic_POC += (1.0 - fb) * micC;
			reach->benthic_PON += (1.0 - fb) * micN;
			micC *= fb;  micN *= fb;  needb *= fb;
		}
		take_N(reach, &outNH4, &outNO3, needb);
		reach->Bed_N_immob = needb;
		reach->Bed_N_limit = fb;
	}
	reach->benthic_POC += p->CUE * micC;   /* attached microbial biomass */
	reach->benthic_PON += p->CUE * micC / p->CN_mic;   /* = micN + immob - mineral */
	reach->Bed_mic_C = micC;

	/* outlet loads and diagnostics */
	reach->NO3_out = outNO3;  reach->NH4_out = outNH4;
	reach->DON_out = outDON;  reach->DOC_out = outDOC;
	reach->POC_out = outPOC;  reach->PON_out = outPON;
	reach->POC_deposit = depC;   reach->POC_entrain = entC;
	reach->Frag_C = fragC;       reach->Frag_N = fragN;
	reach->DOM_dec_C = decC;     reach->N_immob = immob;
	reach->N_limit = nlim;       reach->stream_CO2 = decC - bioC + (1.0 - p->CUE) * micC;
	reach->stream_denitrif = den;
	return;
} /*end compute_stream_bgc*/
