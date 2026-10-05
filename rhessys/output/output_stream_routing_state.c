/*--------------------------------------------------------------*/
/* 																*/
/*					output_stream_routing_state					*/
/*																*/
/*	output_stream_routing_state - appends the stream reach		*/
/*	routing state to a world state file.						*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_stream_routing_state(						*/
/*					struct	world_object	*world,				*/
/*					FILE	*outfile)							*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	Writes the values compute_stream_routing carries from one	*/
/*	day to the next, so a restart from the state file matches	*/
/*	a continuous run.  The block follows all basins:			*/
/*																*/
/*	stream_routing_state										*/
/*	<n>  num_stream_state_basins								*/
/*	  <ID> basin_ID												*/
/*	  <n>  num_stream_reaches									*/
/*	  reach_ID initial_flow previous_lateral_input water_depth	*/
/*	           reservoir_storage Qout previous_Qin   (one line per reach) */
/*																*/
/*	Read back by read_stream_routing_state.						*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"

void	output_stream_routing_state(
									struct	world_object	*world,
									FILE	*outfile)
{
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	int b, i, num_basins;
	struct stream_list_object *stream_list;
	struct stream_network_object *reach;

	num_basins = 0;
	for (b=0; b < world[0].num_basin_files; b++)
		if (world[0].basins[b][0].stream_list.stream_network != NULL) num_basins++;

	fprintf(outfile, "\nstream_routing_state");
	fprintf(outfile, "\n%-30d %s", num_basins, "num_stream_state_basins");
	for (b=0; b < world[0].num_basin_files; b++) {
		stream_list = &(world[0].basins[b][0].stream_list);
		if (stream_list->stream_network == NULL) continue;
		fprintf(outfile, "\n   %-30d %s", world[0].basins[b][0].ID, "basin_ID");
		fprintf(outfile, "\n   %-30d %s", stream_list->num_reaches, "num_stream_reaches");
		for (i=0; i < stream_list->num_reaches; i++) {
			reach = &(stream_list->stream_network[i]);
			fprintf(outfile, "\n   %d %.17g %.17g %.17g %.17g %.17g %.17g",
				reach->reach_ID,
				reach->initial_flow,
				reach->previous_lateral_input,
				reach->water_depth,
				(reach->reservoir_ID != 0) ? reach->reservoir.initial_storage : 0.0,
				reach->Qout,
				reach->previous_Qin);
		}
	}
	fprintf(outfile, "\n");
	return;
} /*end output_stream_routing_state*/
