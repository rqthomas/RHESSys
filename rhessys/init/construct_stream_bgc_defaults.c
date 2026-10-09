/*--------------------------------------------------------------*/
/*								*/
/*	construct_stream_bgc_defaults				*/
/*								*/
/*	Reads the parameter file given with -strbgc <file> for	*/
/*	in-stream storage and processing (compute_stream_bgc.c).	*/
/*	Same key/value format as the other def files; missing	*/
/*	keys take the defaults below.				*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "rhessys.h"
#include "params.h"
#include <string.h>

struct stream_bgc_defaults *construct_stream_bgc_defaults(char *filename)
{
	void *alloc(size_t, char *, char *);
	param *paramPtr = NULL;
	int paramCnt = 0;
	struct stream_bgc_defaults *p;

	p = (struct stream_bgc_defaults *) alloc(sizeof(struct stream_bgc_defaults),
		"stream_bgc_defaults", "construct_stream_bgc_defaults");
	printf("Reading %s\n", filename);
	paramPtr = readParamFile(&paramCnt, filename);

	p->theta    = getDoubleParam(&paramCnt, &paramPtr, "stream_theta", "%lf", 1.047, 1);     /* temperature coefficient */
	p->v_dep    = getDoubleParam(&paramCnt, &paramPtr, "stream_v_dep", "%lf", 10.0, 1);      /* m/day, POC/PON deposition */
	p->k_frag   = getDoubleParam(&paramCnt, &paramPtr, "stream_k_frag", "%lf", 0.01, 1);     /* 1/day at 20 C */
	p->e_max    = getDoubleParam(&paramCnt, &paramPtr, "stream_e_max", "%lf", 0.5, 1);       /* 1/day at bankfull */
	p->export_p = getDoubleParam(&paramCnt, &paramPtr, "stream_export_p", "%lf", 2.0, 1);    /* (DIM) */
	p->vf_DOM   = getDoubleParam(&paramCnt, &paramPtr, "stream_vf_DOM", "%lf", 0.1, 1);      /* m/day */
	p->CUE      = getDoubleParam(&paramCnt, &paramPtr, "stream_CUE", "%lf", 0.3, 1);         /* (DIM) */
	p->CN_mic   = getDoubleParam(&paramCnt, &paramPtr, "stream_CN_mic", "%lf", 10.0, 1);     /* kgC/kgN */
	p->vf_N     = getDoubleParam(&paramCnt, &paramPtr, "stream_vf_N", "%lf", 2.0, 1);        /* m/day */
	p->vf_denit = getDoubleParam(&paramCnt, &paramPtr, "stream_vf_denit", "%lf", 0.1, 1);    /* m/day */
	p->f_mic    = getDoubleParam(&paramCnt, &paramPtr, "stream_f_mic", "%lf", 0.0, 1);       /* (DIM) microbial share of bed breakdown */

	/* optional list of reaches inside a lake/reservoir (whitespace-separated reach IDs,
	   '#' starts a comment line): no in-stream processing there, loads pass through */
	p->num_lake_reaches = 0;
	p->lake_reaches = NULL;
	{
		char *lf = getStrParam(&paramCnt, &paramPtr, "stream_bgc_lake_reaches", "%s", "none", 1);
		if (strcmp(lf, "none") != 0) {
			FILE *f = fopen(lf, "r");
			char line[4096], *tk;
			int n = 0, cap = 64;
			if (f == NULL) {
				fprintf(stderr, "FATAL ERROR: cannot open stream_bgc_lake_reaches file %s\n", lf);
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
			printf("Stream bgc: %d lake/reservoir reaches pass loads through without processing (%s)\n", n, lf);
		}
	}

	if ((p->CUE < 0.0) || (p->CUE > 1.0) || (p->CN_mic <= 0.0) || (p->theta <= 0.0) || (p->f_mic < 0.0) || (p->f_mic > 1.0)) {
		fprintf(stderr, "FATAL ERROR: in %s stream_CUE and stream_f_mic must be 0-1, stream_CN_mic and stream_theta > 0\n", filename);
		exit(EXIT_FAILURE);
	}
	if (paramPtr != NULL) free(paramPtr);
	return(p);
} /*end construct_stream_bgc_defaults*/
