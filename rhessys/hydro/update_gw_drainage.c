/*--------------------------------------------------------------*/
/* 								*/
/*		update_gw_drainage					*/
/*								*/
/*								*/
/*	NAME							*/
/*	update_gw_drainage -  					*/
/* 		drainage shallow subsurface saturation zone	*/
/*		from each patch to a regional (hillslope scale)	*/
/*		groundwater store				*/
/*		nitrogen is also drained using assumption 	*/
/*		of an exponential decay of N with depth		*/
/*								*/
/*	SYNOPSIS						*/
/*	int update_gw_drainage(					*/
/*			struct patch_object *			*/
/*			struct hillslope_object *		*/
/*			struct command_line_object *		*/
/*			struct date,				*/
/*			)					*/
/*								*/
/*	returns:						*/
/*								*/
/*	OPTIONS							*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*	preset code just uses a user assigned loading rate	*/
/*	and all of it is nitrate				*/
/*								*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

int update_gw_drainage(
				  struct  patch_object   *patch,
				  struct  hillslope_object *hillslope,
				  struct  zone_object *zone,
				  struct  command_line_object *command_line,
				  struct	date	current_date)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.			*/
	/*------------------------------------------------------*/
	
	double  compute_z_final(
		int,
		double,
		double,
		double,
		double,
		double);


	/*------------------------------------------------------*/
	/*	Local Variable Definition. 			*/
	/*------------------------------------------------------*/
	int ok = 1;
	double drainage,sat_store,N_loss;
	double preday_sat_deficit_z, add_field_capacity;
	double sat_to_gw_coeff;
	/*------------------------------------------------------*/
	/*		assume percent of incoming precip	*/
	/*------------------------------------------------------*/
	if (zone[0].hourly_rain_flag==1){
	  sat_to_gw_coeff = patch[0].soil_defaults[0][0].sat_to_gw_coeff / 24;
	}
	else sat_to_gw_coeff = patch[0].soil_defaults[0][0].sat_to_gw_coeff; 

	/*------------------------------------------------------*/
	/* multiply by Ksat vertical which is surface Ksat 	*/
	/* so that impervious areas don't infiltrate to deep gw */
	/*------------------------------------------------------*/
	sat_to_gw_coeff = sat_to_gw_coeff * patch[0].Ksat_vertical;

	drainage = sat_to_gw_coeff * patch[0].detention_store;
	patch[0].detention_store -= drainage;
	patch[0].gw_drainage = drainage;
	hillslope[0].gw.storage += (drainage * patch[0].area / hillslope[0].area);
	/*------------------------------------------------------*/
	/*	soil DOM carried by the recharge water: a fraction	*/
	/*	(gw_DOM_recharge_frac) of the soil-water DOC/DON	*/
	/*	concentration, debited from the soil pools. Without	*/
	/*	it groundwater only receives surface DOM, so		*/
	/*	baseflow carries almost no DOC between storms.		*/
	/*------------------------------------------------------*/
	if ((patch[0].soil_defaults[0][0].gw_DOM_recharge_frac > ZERO) && (drainage > ZERO)) {
		double soil_water = patch[0].soil_defaults[0][0].soil_water_cap - patch[0].sat_deficit
			+ patch[0].rz_storage + patch[0].unsat_storage;
		if (soil_water > ZERO) {
			double f = min(1.0, patch[0].soil_defaults[0][0].gw_DOM_recharge_frac * drainage / soil_water);
			double doc = (patch[0].soil_cs.DOC > ZERO) ? f * patch[0].soil_cs.DOC : 0.0;
			double don = (patch[0].soil_ns.DON > ZERO) ? f * patch[0].soil_ns.DON : 0.0;
			/* cap the recharge DOC concentration (mg C/L = g/m3; drainage in m) and move
			   DON in the same proportion; what is not moved stays in the soil DOM pools */
			if ((patch[0].soil_defaults[0][0].gw_DOC_recharge_max > 0.0) && (doc > ZERO)) {
				double doc_max = patch[0].soil_defaults[0][0].gw_DOC_recharge_max * drainage / 1000.0; /* kg C/m2 */
				if (doc > doc_max) {
					don *= doc_max / doc;
					doc = doc_max;
				}
			}
			if (doc > 0.0) {
				hillslope[0].gw.DOC += (doc * patch[0].area / hillslope[0].area);
				patch[0].cdf.DOC_to_gw += doc;
				patch[0].soil_cs.DOC -= doc;
			}
			if (don > 0.0) {
				hillslope[0].gw.DON += (don * patch[0].area / hillslope[0].area);
				patch[0].ndf.DON_to_gw += don;
				patch[0].soil_ns.DON -= don;
			}
		}
	}

	/*------------------------------------------------------*/
	/*	determine associated N leached			*/
	/*------------------------------------------------------*/
	if (patch[0].surface_DON > ZERO) {
		N_loss = sat_to_gw_coeff * patch[0].surface_DON;
		hillslope[0].gw.DON += (N_loss * patch[0].area / hillslope[0].area);
		patch[0].ndf.DON_to_gw += N_loss;
		patch[0].surface_DON -= N_loss;
		}
	if (patch[0].surface_DOC > ZERO) {
		N_loss = sat_to_gw_coeff * patch[0].surface_DOC;
		hillslope[0].gw.DOC += (N_loss * patch[0].area / hillslope[0].area);
		patch[0].cdf.DOC_to_gw += N_loss;
		patch[0].surface_DOC -= N_loss;
		}
	
	
	if (patch[0].surface_NH4 > ZERO) {
		N_loss = sat_to_gw_coeff * patch[0].surface_NH4;
		hillslope[0].gw.NH4 += (N_loss * patch[0].area / hillslope[0].area);
		patch[0].ndf.N_to_gw += N_loss;
		patch[0].surface_NH4 -= N_loss;
		}
	
	if (patch[0].surface_NO3 > ZERO) {
		N_loss = sat_to_gw_coeff * patch[0].surface_NO3;
		hillslope[0].gw.NO3 += (N_loss * patch[0].area / hillslope[0].area);
		patch[0].ndf.N_to_gw += N_loss;
		patch[0].surface_NO3 -= N_loss;
		}

	return (!ok);
} /* end update_gw_drainage.c */
