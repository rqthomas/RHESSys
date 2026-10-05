/*--------------------------------------------------------------*/
/* 																*/
/*					read_stream_routing_state					*/
/*																*/
/*	read_stream_routing_state - restores stream reach routing	*/
/*	state from the end of a world state file.					*/
/*																*/
/*	SYNOPSIS													*/
/*	void	read_stream_routing_state(							*/
/*					struct	world_object	*world,				*/
/*					FILE	*world_file)						*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	Called after all basins are constructed.  Looks for the		*/
/*	stream_routing_state block written by						*/
/*	output_stream_routing_state; if the file has none (an older	*/
/*	worldfile) the reaches keep their cold-start values.		*/
/*	Reaches are matched by basin_ID and reach_ID.				*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include "rhessys.h"

void	read_stream_routing_state(
								  struct	world_object	*world,
								  FILE	*world_file)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/
	int	read_record(FILE *, char *);
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	int b, i, k, s, num_basins, basin_ID, num_reaches, reach_ID, num_restored, num_missing;
	double initial_flow, previous_lateral_input, water_depth, reservoir_storage, Qout, previous_Qin;
	char	token[MAXSTR];
	char	record[MAXSTR];
	struct stream_list_object *stream_list;
	struct stream_network_object *reach;

	/*--------------------------------------------------------------*/
	/*	find the block; an older worldfile has none					*/
	/*--------------------------------------------------------------*/
	while (fscanf(world_file, "%s", token) == 1) {
		if (strcmp(token, "stream_routing_state") == 0) break;
	}
	if (feof(world_file) || strcmp(token, "stream_routing_state") != 0) {
		printf("\nNo stream_routing_state in worldfile; stream reaches start empty\n");
		return;
	}

	if (fscanf(world_file, "%d", &num_basins) != 1) {
		fprintf(stderr, "FATAL ERROR: in read_stream_routing_state, cannot read num_stream_state_basins\n");
		exit(EXIT_FAILURE);
	}
	read_record(world_file, record);

	for (s=0; s < num_basins; s++) {
		if (fscanf(world_file, "%d", &basin_ID) != 1) {
			fprintf(stderr, "FATAL ERROR: in read_stream_routing_state, cannot read basin_ID\n");
			exit(EXIT_FAILURE);
		}
		read_record(world_file, record);
		if (fscanf(world_file, "%d", &num_reaches) != 1) {
			fprintf(stderr, "FATAL ERROR: in read_stream_routing_state, cannot read num_stream_reaches\n");
			exit(EXIT_FAILURE);
		}
		read_record(world_file, record);

		stream_list = NULL;
		for (b=0; b < world[0].num_basin_files; b++)
			if (world[0].basins[b][0].ID == basin_ID) stream_list = &(world[0].basins[b][0].stream_list);
		if ((stream_list != NULL) && (stream_list->stream_network == NULL)) stream_list = NULL;
		if (stream_list == NULL)
			fprintf(stderr, "WARNING: stream_routing_state for basin %d ignored (no stream network)\n", basin_ID);

		num_restored = 0;
		num_missing = 0;
		for (i=0; i < num_reaches; i++) {
			if (fscanf(world_file, "%d %lf %lf %lf %lf %lf %lf", &reach_ID,
				&initial_flow, &previous_lateral_input, &water_depth,
				&reservoir_storage, &Qout, &previous_Qin) != 7) {
				fprintf(stderr, "FATAL ERROR: in read_stream_routing_state, bad reach record %d of basin %d\n",
					i, basin_ID);
				exit(EXIT_FAILURE);
			}
			if (stream_list == NULL) continue;
			reach = NULL;
			for (k=0; k < stream_list->num_reaches; k++)
				if (stream_list->stream_network[k].reach_ID == reach_ID) {
					reach = &(stream_list->stream_network[k]);
					break;
				}
			if (reach == NULL) {
				num_missing++;
				continue;
			}
			reach->initial_flow = initial_flow;
			reach->previous_lateral_input = previous_lateral_input;
			reach->water_depth = water_depth;
			if (reach->reservoir_ID != 0) reach->reservoir.initial_storage = reservoir_storage;
			reach->Qout = Qout;
			reach->previous_Qin = previous_Qin;
			num_restored++;
		}
		if (stream_list != NULL) {
			printf("\nStream routing state restored for %d of %d reaches in basin %d\n",
				num_restored, stream_list->num_reaches, basin_ID);
			if ((num_missing > 0) || (num_restored != stream_list->num_reaches))
				fprintf(stderr, "WARNING: stream_routing_state does not match the stream table "
					"(%d saved reaches not found, %d reaches start empty)\n",
					num_missing, stream_list->num_reaches - num_restored);
		}
	}
	return;
} /*end read_stream_routing_state*/
