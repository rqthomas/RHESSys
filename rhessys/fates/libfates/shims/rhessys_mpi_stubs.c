/* MPI entry points referenced by FATES DetermineGridCellNeighbors (seed
   dispersal). RHESSys never enables it (seeddisp_cadence = 0); abort loudly
   if that ever changes rather than silently returning garbage. */
#include <stdio.h>
#include <stdlib.h>
static void die(const char *f) {
	fprintf(stderr, "FATAL: %s called -- FATES seed dispersal/MPI is not "
		"supported in the RHESSys host build\n", f);
	abort();
}
void mpi_allgather_(void)  { die("MPI_Allgather");  }
void mpi_allgatherv_(void) { die("MPI_Allgatherv"); }
