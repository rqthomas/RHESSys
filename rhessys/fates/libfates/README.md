# libfates build recipe for the RHESSys FATES host

This folder holds everything needed to build `libfates.a` for RHESSys, apart from the FATES
source itself: `CMakeLists.txt` (standalone build), `shims/` (stand-ins for the CESM `shr_*` and
CLM seed-dispersal/MPI modules FATES expects from a climate model) and the host patch.

## Host patch

`rhessys_host_patches.patch` is applied to NGEET/fates at tag `sci.1.93.2_api.46.0.0`
(commit `ec37ae8b`) before libfates is built. It adds the hooks the RHESSys host
(`../rhessys_veg_host.c`, `../rhessys_fates.F90`, `../rhessys_fates_mgmt.F90`) uses:

- host-scheduled logging trigger, per-cohort harvest rates and whole above-ground harvest
  switch for the forest-management module (EDLoggingMortalityMod, EDPatchDynamicsMod,
  EDMortalityFunctionsMod)
- leaf export: the share of a site's fresh leaf litter the host takes out before it enters the
  FATES litter pool (e.g. stream-side leaves falling into the channel)
- canopy-structure fix: the scratch vectors are sized for the cohorts actually in the patch
  (EDCanopyStructureMod)
- phenology fix: growing degree days also accumulate in the `nevercold` state, so a site whose
  cold leaf-off comes more than 400 days after the previous one still leafs out (EDPhysiologyMod)

With the hooks unused, FATES behaves like the stock tag. Apply from the FATES source folder:

## Build

From the parent folder of the RHESSys repo:

    git clone --depth 1 --branch sci.1.93.2_api.46.0.0 https://github.com/NGEET/fates.git fates_lib/src
    (cd fates_lib/src && patch -p1 < ../../RHESSys/rhessys/fates/libfates/rhessys_host_patches.patch)
    cmake -S RHESSys/rhessys/fates/libfates -B fates_lib/build && cmake --build fates_lib/build -j8
    cd RHESSys/rhessys && make clobber PGM=__none__ && make fates=T PGM=rhessys7.5_fates_strtemp

`FATES_SRC` defaults to `<parent>/fates_lib/src`; pass `-DFATES_SRC=<path>` for another source
folder. The makefile's default `FATES_DIR` is `../../fates_lib/build`. The CCR kit's
`build.sh --fates` runs these steps.
