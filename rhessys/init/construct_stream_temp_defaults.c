/*--------------------------------------------------------------*/
/*								*/
/*	construct_stream_temp_defaults				*/
/*								*/
/*	Reads the parameter file given with -strtemp <file> for	*/
/*	stream temperature (compute_stream_temperature.c).	*/
/*	Same key/value format as the other def files; missing	*/
/*	keys take the defaults below.				*/
/*								*/
/*	Optional files named in it:				*/
/*	stream_temp_lake_reaches: reach IDs inside a lake or	*/
/*	  reservoir (whitespace/comma separated, '#' comment	*/
/*	  lines); these only mix their inputs.			*/
/*	stream_temp_reach_file: one line per reach,		*/
/*	  reach_ID S_topo shade_override c_overhang		*/
/*	  (shade_override < 0 = LAI shading model; c_overhang	*/
/*	  < 0 = stream_c_overhang).				*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "rhessys.h"
#include "params.h"
#include <string.h>

struct stream_temp_defaults *construct_stream_temp_defaults(char *filename)
{
	void *alloc(size_t, char *, char *);
	param *paramPtr = NULL;
	int paramCnt = 0;
	struct stream_temp_defaults *p;

	p = (struct stream_temp_defaults *) alloc(sizeof(struct stream_temp_defaults),
		"stream_temp_defaults", "construct_stream_temp_defaults");
	printf("Reading %s\n", filename);
	paramPtr = readParamFile(&paramCnt, filename);

	p->albedo_w       = getDoubleParam(&paramCnt, &paramPtr, "stream_albedo_w", "%lf", 0.08, 1);       /* (DIM) */
	p->emiss_w        = getDoubleParam(&paramCnt, &paramPtr, "stream_emiss_w", "%lf", 0.97, 1);        /* (DIM) */
	p->emiss_veg      = getDoubleParam(&paramCnt, &paramPtr, "stream_emiss_veg", "%lf", 0.97, 1);      /* (DIM) */
	p->k_shade        = getDoubleParam(&paramCnt, &paramPtr, "stream_k_shade", "%lf", -1.0, 1);        /* < 0 = strata ext_coef */
	p->SAI            = getDoubleParam(&paramCnt, &paramPtr, "stream_SAI", "%lf", 0.0, 1);             /* m2/m2 added to proj_pai */
	p->c_overhang     = getDoubleParam(&paramCnt, &paramPtr, "stream_c_overhang", "%lf", 1.0, 1);      /* (DIM) */
	p->c_wind_shelter = getDoubleParam(&paramCnt, &paramPtr, "stream_c_wind_shelter", "%lf", 0.5, 1);  /* (DIM) */
	p->wind_fn_mult   = getDoubleParam(&paramCnt, &paramPtr, "stream_wind_fn_mult", "%lf", 1.0, 1);    /* (DIM) */
	p->k_bed          = getDoubleParam(&paramCnt, &paramPtr, "stream_k_bed", "%lf", 0.0, 1);           /* kJ/m2/day/K */
	p->tau_gw         = getDoubleParam(&paramCnt, &paramPtr, "stream_tau_gw", "%lf", 365.0, 1);        /* days */
	p->dT_gw          = getDoubleParam(&paramCnt, &paramPtr, "stream_dT_gw", "%lf", 1.0, 1);           /* C */
	p->T_gw_const     = getDoubleParam(&paramCnt, &paramPtr, "stream_T_gw_const", "%lf", -999.0, 1);   /* C; -999 = running mean */
	p->tau_sub        = getDoubleParam(&paramCnt, &paramPtr, "stream_tau_sub", "%lf", 60.0, 1);        /* days */
	p->T_init         = getDoubleParam(&paramCnt, &paramPtr, "stream_T_init", "%lf", 10.0, 1);         /* C */
	p->T_gw_init      = getDoubleParam(&paramCnt, &paramPtr, "stream_T_gw_init", "%lf", 12.5, 1);      /* C */
	p->V_min          = getDoubleParam(&paramCnt, &paramPtr, "stream_V_min", "%lf", 1.0, 1);           /* m3 */
	p->snow_insulate  = getIntParam(&paramCnt, &paramPtr, "stream_T_snow_insulate", "%d", 1, 1);
	p->exchange       = getIntParam(&paramCnt, &paramPtr, "stream_temp_exchange", "%d", 1, 1);

	/* reaches inside a lake/reservoir */
	p->num_lake_reaches = 0;
	p->lake_reaches = NULL;
	{
		char *lf = getStrParam(&paramCnt, &paramPtr, "stream_temp_lake_reaches", "%s", "none", 1);
		if (strcmp(lf, "none") != 0) {
			FILE *f = fopen(lf, "r");
			char line[4096], *tk;
			int n = 0, cap = 64;
			if (f == NULL) {
				fprintf(stderr, "FATAL ERROR: cannot open stream_temp_lake_reaches file %s\n", lf);
				exit(EXIT_FAILURE);
			}
			p->lake_reaches = (int *) malloc(cap * sizeof(int));
			while (fgets(line, sizeof(line), f) != NULL) {
				if (line[0] == '#') continue;
				for (tk = strtok(line, " \t\r\n,"); tk != NULL; tk = strtok(NULL, " \t\r\n,")) {
					if (n == cap) { cap *= 2; p->lake_reaches = (int *) realloc(p->lake_reaches, cap * sizeof(int)); }
					p->lake_reaches[n++] = atoi(tk);
				}
			}
			fclose(f);
			p->num_lake_reaches = n;
			printf("Stream temperature: %d lake/reservoir reaches only mix their inputs (%s)\n", n, lf);
		}
	}

	/* per-reach shading table */
	p->num_reach_params = 0;
	p->rp_ID = NULL;
	p->rp_S_topo = p->rp_shade = p->rp_c_oh = NULL;
	{
		char *rf = getStrParam(&paramCnt, &paramPtr, "stream_temp_reach_file", "%s", "none", 1);
		if (strcmp(rf, "none") != 0) {
			FILE *f = fopen(rf, "r");
			char line[4096];
			int n = 0, cap = 64, id;
			double st, sh, co;
			if (f == NULL) {
				fprintf(stderr, "FATAL ERROR: cannot open stream_temp_reach_file %s\n", rf);
				exit(EXIT_FAILURE);
			}
			p->rp_ID = (int *) malloc(cap * sizeof(int));
			p->rp_S_topo = (double *) malloc(cap * sizeof(double));
			p->rp_shade = (double *) malloc(cap * sizeof(double));
			p->rp_c_oh = (double *) malloc(cap * sizeof(double));
			while (fgets(line, sizeof(line), f) != NULL) {
				if ((line[0] == '#') || (sscanf(line, "%d %lf %lf %lf", &id, &st, &sh, &co) != 4)) continue;
				if ((st < 0.0) || (st > 1.0) || (sh > 1.0)) {
					fprintf(stderr, "FATAL ERROR: in %s reach %d needs 0 <= S_topo <= 1 and shade_override <= 1\n", rf, id);
					exit(EXIT_FAILURE);
				}
				if (n == cap) {
					cap *= 2;
					p->rp_ID = (int *) realloc(p->rp_ID, cap * sizeof(int));
					p->rp_S_topo = (double *) realloc(p->rp_S_topo, cap * sizeof(double));
					p->rp_shade = (double *) realloc(p->rp_shade, cap * sizeof(double));
					p->rp_c_oh = (double *) realloc(p->rp_c_oh, cap * sizeof(double));
				}
				p->rp_ID[n] = id; p->rp_S_topo[n] = st; p->rp_shade[n] = sh; p->rp_c_oh[n] = co;
				n++;
			}
			fclose(f);
			p->num_reach_params = n;
			printf("Stream temperature: shading parameters for %d reaches (%s)\n", n, rf);
		}
	}

	if ((p->albedo_w < 0.0) || (p->albedo_w > 1.0) || (p->emiss_w <= 0.0) || (p->emiss_w > 1.0)
		|| (p->emiss_veg < 0.0) || (p->emiss_veg > 1.0) || (p->c_wind_shelter < 0.0) || (p->c_wind_shelter > 1.0)) {
		fprintf(stderr, "FATAL ERROR: in %s stream_albedo_w, stream_emiss_w, stream_emiss_veg and stream_c_wind_shelter must be 0-1\n", filename);
		exit(EXIT_FAILURE);
	}
	if ((p->tau_gw < 1.0) || (p->tau_sub < 1.0) || (p->k_bed < 0.0) || (p->wind_fn_mult < 0.0)
		|| (p->c_overhang < 0.0) || (p->SAI < 0.0) || (p->V_min < 0.0)) {
		fprintf(stderr, "FATAL ERROR: in %s stream_tau_gw and stream_tau_sub must be >= 1 day; stream_k_bed, stream_wind_fn_mult, stream_c_overhang, stream_SAI, stream_V_min >= 0\n", filename);
		exit(EXIT_FAILURE);
	}
	if (paramPtr != NULL) free(paramPtr);
	return(p);
} /*end construct_stream_temp_defaults*/
