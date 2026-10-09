/*--------------------------------------------------------------*/
/* 																*/
/*					output_stream_routing						*/
/*																*/
/*	output_stream_routing - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_stream_routing - outputs current contents of streamflow routing.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_stream_routing( int routing_flag,										*/	
/*					struct	basin_object	*basin,				*/
/*					struct	date	date,  						*/
/*					FILE 	*outfile)							*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	outputs streamflow routing results according to commandline			*/
/*	specifications to specific files							*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	We only permit one fileset per spatial modelling level.     */
/*	Each fileset has one file for each timestep.  				*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"

void	output_stream_routing(			
					 struct	stream_network_object  *stream_network,
					 struct command_line_object *command_line,
					 struct	date	date,
					 FILE *outfile)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*      Output streamflow routing results.			*/
	/*	Columns: day month year reach_ID Qout(m3/day)		*/
	/*	         lateral_input(m3/day) Qin(m3/day)		*/
	/*	         water_depth(m) reservoir_storage(ha-m)		*/
	/*	When grow_flag > 0, four additional columns are appended:*/
	/*	         NO3_out(kgN/day) NH4_out(kgN/day)		*/
	/*	         DON_out(kgN/day) DOC_out(kgC/day)		*/
	/*	Then (always): sediment_out(kg/day)                      */
	/*	When grow_flag > 0, lateral input loads follow:		*/
	/*	         lateral_NO3(kgN/day) lateral_NH4(kgN/day)	*/
	/*	         lateral_DON(kgN/day) lateral_DOC(kgC/day)	*/
	/*	Then (always): lateral_sediment(kg/day), sediment from	*/
	/*	         lateral input patches only			*/
	/*	Last (always): POC_out PON_out lateral_POC lateral_PON	*/
	/*	         lateral_POC_labile lateral_PON_labile (kg/day):	*/
	/*	         leaf litter from stream-side patches		*/
	/*--------------------------------------------------------------*/

	if (command_line[0].grow_flag > 0) {
		fprintf(outfile, "%d %d %d %d %lf %lf %lf %lf %lf %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e",
			date.day,
			date.month,
			date.year,
			stream_network[0].reach_ID,
			stream_network[0].Qout * 86400,
			stream_network[0].previous_lateral_input * stream_network[0].length * 86400,
			stream_network[0].previous_Qin * 86400,
			stream_network[0].water_depth,
			stream_network[0].reservoir.initial_storage/10000,
			stream_network[0].NO3_out,
			stream_network[0].NH4_out,
			stream_network[0].DON_out,
			stream_network[0].DOC_out,
			stream_network[0].sediment_out,
			stream_network[0].lateral_NO3,
			stream_network[0].lateral_NH4,
			stream_network[0].lateral_DON,
			stream_network[0].lateral_DOC,
			stream_network[0].lateral_sediment,
			stream_network[0].POC_out, stream_network[0].PON_out,
			stream_network[0].lateral_POC, stream_network[0].lateral_PON,
			stream_network[0].lateral_POC_labile, stream_network[0].lateral_PON_labile,
			stream_network[0].lateral_POC_sed, stream_network[0].lateral_PON_sed);
	} else {
		fprintf(outfile, "%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %e %e %e %e %e %e %e %e", 
			date.day,
			date.month,
			date.year,
			stream_network[0].reach_ID,
			stream_network[0].Qout * 86400,
			stream_network[0].previous_lateral_input * stream_network[0].length * 86400,
			stream_network[0].previous_Qin * 86400,
			stream_network[0].water_depth,
			stream_network[0].reservoir.initial_storage/10000,
			stream_network[0].sediment_out,
			stream_network[0].lateral_sediment,
			stream_network[0].POC_out, stream_network[0].PON_out,
			stream_network[0].lateral_POC, stream_network[0].lateral_PON,
			stream_network[0].lateral_POC_labile, stream_network[0].lateral_PON_labile,
			stream_network[0].lateral_POC_sed, stream_network[0].lateral_PON_sed);
	}

	if (command_line[0].stream_bgc_flag == 1)
		fprintf(outfile, " %e %e %e %e %e %e %e %e %e %f %e %e %e %e %e %e %e %e %e %e %e %e %f",
			stream_network[0].benthic_POC, stream_network[0].benthic_PON,
			stream_network[0].POC_deposit, stream_network[0].POC_entrain,
			stream_network[0].Frag_C, stream_network[0].Frag_N,
			stream_network[0].DOM_dec_C, stream_network[0].N_mineral,
			stream_network[0].N_immob, stream_network[0].N_limit,
			stream_network[0].stream_CO2, stream_network[0].stream_denitrif,
			stream_network[0].V_water, stream_network[0].wc_DOC, stream_network[0].wc_DON,
			stream_network[0].wc_NO3, stream_network[0].wc_NH4, stream_network[0].wc_POC, stream_network[0].wc_PON,
			stream_network[0].Bed_mic_C, stream_network[0].Bed_N_mineral,
			stream_network[0].Bed_N_immob, stream_network[0].Bed_N_limit);
	fprintf(outfile, "\n");
	return;
} /*end output_stream_routing*/
