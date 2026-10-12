/*--------------------------------------------------------------*/
/* 																*/
/*							world_daily_F							*/
/*																*/
/*	NAME														*/
/*	world_daily_F													*/
/*				 - performs cycling and output of a world		*/
/*					with parallel hillslope simulations			*/
/*																*/
/*																*/
/*	SYNOPSIS													*/
/*	void		world_daily_F(									*/
/*					long   	day,								*/
/*					struct world_object *world,					*/
/*					struct command_line_object *command_line,	*/
/*					struct tec_entry *event,					*/
/*					struct date current_date)					*/
/* 																*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	Performs daily simulation of the world with hillslope and	*/
/*	lower level executed in parallel.							*/
/*																*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	At the moment this routine only handles setting Io and 		*/
/*	declination for the world.  								*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "rhessys_veg.h"
#include <openmp.h>

void world_daily_F(
				   long   	day,
				   struct world_object *world,
				   struct command_line_object *command_line,
				   struct tec_entry *event,
				   struct date current_date)
{
	/*--------------------------------------------------------------*/
	/*  Local Function Declarations.                                */
	/*--------------------------------------------------------------*/
	void	basin_daily_F(
		long	,
		struct world_object *,
		struct basin_object *,
		struct command_line_object *,
		struct tec_entry *,
		struct date);
	/*--------------------------------------------------------------*/
	/*  Local variable definition.                                  */
	/*--------------------------------------------------------------*/
	void	zone_daily_F(long, struct world_object *, struct basin_object *,
		struct hillslope_object *, struct zone_object *,
		struct command_line_object *, struct tec_entry *, struct date);
	int	basin, h, z;
	struct basin_object *bas;
	/*--------------------------------------------------------------*/
	/*	Simulate over all of the basins.							*/
	/*--------------------------------------------------------------*/
  world[0].target_status = 1;
	/*--------------------------------------------------------------*/
	/*	FATES backend (-demog): every zone's daily weather first,	*/
	/*	then all FATES sites in one (parallel) batch, then the	*/
	/*	usual pass with zone_daily_F doing only its patches.	*/
	/*--------------------------------------------------------------*/
	if (command_line[0].demog_flag) {
		veg_zone_phase = 1;
		for (basin = 0; basin < world[0].num_basin_files; basin++) {
			bas = world[0].basins[basin];
			for (h = 0; h < bas[0].num_hillslopes; h++)
				for (z = 0; z < bas[0].hillslopes[h][0].num_zones; z++)
					zone_daily_F(day, world, bas, bas[0].hillslopes[h],
						bas[0].hillslopes[h][0].zones[z], command_line, event, current_date);
		}
		veg_host_daily_batch(command_line, current_date);
		veg_zone_phase = 2;
	}
	for ( basin = 0; basin < world[0].num_basin_files; basin++ ){
		basin_daily_F(	day,
			world,
			world[0].basins[basin],
			command_line,
			event,
			current_date);
	}
	veg_zone_phase = 0;
	return;
} /*end world_daily_F.c*/
