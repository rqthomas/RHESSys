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

	p->tau_gw       = getDoubleParam(&paramCnt, &paramPtr, "stream_tau_gw", "%lf", 365.0, 1);        /* days */
	p->dT_gw        = getDoubleParam(&paramCnt, &paramPtr, "stream_dT_gw", "%lf", 1.0, 1);           /* C */
	p->tau_sub      = getDoubleParam(&paramCnt, &paramPtr, "stream_tau_sub", "%lf", 60.0, 1);        /* days */
	p->c_overhang   = getDoubleParam(&paramCnt, &paramPtr, "stream_c_overhang", "%lf", 1.0, 1);      /* (DIM) */
	p->wind_fn_mult = getDoubleParam(&paramCnt, &paramPtr, "stream_wind_fn_mult", "%lf", 1.0, 1);    /* (DIM) */
	p->T_init       = getDoubleParam(&paramCnt, &paramPtr, "stream_T_init", "%lf", 10.0, 1);         /* C */
	p->T_gw_init    = getDoubleParam(&paramCnt, &paramPtr, "stream_T_gw_init", "%lf", 12.5, 1);      /* C */

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

	if ((p->tau_gw < 1.0) || (p->tau_sub < 1.0) || (p->wind_fn_mult < 0.0) || (p->c_overhang < 0.0)) {
		fprintf(stderr, "FATAL ERROR: in %s stream_tau_gw and stream_tau_sub must be >= 1 day; stream_wind_fn_mult and stream_c_overhang >= 0\n", filename);
		exit(EXIT_FAILURE);
	}
	if (paramPtr != NULL) free(paramPtr);
	return(p);
} /*end construct_stream_temp_defaults*/
