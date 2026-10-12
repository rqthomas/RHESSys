#define _FORTRAN_SOURCE_ 1
#include "rhessys_veg.h"
!===============================================================================
! rhessys_fates.F90 -- Fortran side of the RHESSys <-> FATES binding (Phase 0)
!
! Mirrors GLM3/src/glm_aed.F90: BIND(C) wrappers that the C host calls through
! function pointers (see rhessys_veg.c / rhessys_veg.h). The _VEG_* BIND(C)
! names come from rhessys_veg.h so both sides stay in sync.
!
! CONFIG TARGET (design doc 9.4b): coupled N (N-only, ECA), btran water stress,
! WMFire (SPITFIRE off), TEC management; FATES-Hydro/SP/LUH2/methane/P off.
!
! PINNED TARGET: NGEET/fates tag sci.1.93.2_api.46.0.0 (API v46, latest stable
! as of 2026-10-11). Entry-point names + ctrlparm tag strings below were verified
! against this tag. If you bump the tag, re-verify FatesInterfaceMod names, the
! set_fates_ctrlparms tag strings, and the required-parm list (design doc 11 #3).
!===============================================================================
module rhessys_fates

   use, intrinsic :: ISO_C_BINDING
   use FatesInterfaceMod,    only : FatesInterfaceInit,        &
                                    SetFatesGlobalElements1,    &
                                    SetFatesGlobalElements2,    &
                                    set_fates_ctrlparms,        &
                                    InitTimeAveragingGlobals
   use FatesConstantsMod,    only : ievergreen
   use FatesConstantsMod,    only : itrue, ifalse, default_regeneration, &
                                    lmrmodel_ryan_1991, cstarvation_model_lin
   use PRTGenericMod,        only : carbon_only
   use FatesRadiationMemMod, only : norman_solver
   use LeafBiophysicsMod,    only : FvCB1980
   use FatesHydraulicsMemMod,only : hydr_solver_1DTaylor
   use FatesInterfaceTypesMod,only: hlm_stepsize, numpft
   use FatesInterfaceMod,    only : fates_interface_type, allocate_bcin, allocate_bcout, &
                                    allocate_bcpconst, set_bcpconst, zero_bcs, SetFatesTime
   use FatesParameterDerivedMod, only : param_derived
   use EDInitMod,            only : init_site_vars, zero_site, set_site_properties, init_patches
   use EDMainMod,            only : ed_update_site
   use EDCanopyStructureMod, only : canopy_summarization, update_hlm_dynamics, &
                                    UpdateFatesAvgSnowDepth
   use FatesRadiationDriveMod, only : FatesNormalizedCanopyRadiation, FatesSunShadeFracs
   use EDBtranMod,           only : get_active_suction_layers, btran_ed
   use FatesPlantRespPhotosynthMod, only : FatesPlantRespPhotosynthDrive
   use EDAccumulateFluxesMod, only : AccumulateFluxes_ED
   use FatesInterfaceMod,    only : UpdateFatesRMeansTStep
   use EDMainMod,            only : ed_ecosystem_dynamics
   use FatesSoilBGCFluxMod,  only : UnPackNutrientAquisitionBCs
   use ChecksBalancesMod,    only : SiteMassStock
   use PRTParametersMod,     only : prt_params
   use PRTGenericMod,        only : leaf_organ, fnrt_organ, carbon12_element, &
                                    carbon_nitrogen_phosphorus
   use PRTGenericMod,        only : sapw_organ, struct_organ, store_organ, repro_organ
   use FatesInterfaceTypesMod,only: hlm_current_year
   use FatesConstantsMod,    only : i_term_mort_type_cstarv, i_term_mort_type_canlev, &
                                    i_term_mort_type_numdens
   use EDPftvarcon,          only : EDPftvarcon_inst
   use FatesLitterMod,       only : litter_type, ncwd, ilabile, icellulose, ilignin
   use FatesRestartInterfaceMod, only : fates_restart_interface_type
   use FatesIODimensionsMod, only : fates_bounds_type
   use FatesInterfaceTypesMod, only : fates_maxElementsPerSite
   use EDParamsMod,          only : ED_val_cwd_fcel, ED_val_cwd_flig
   use EDLoggingMortalityMod, only : host_logging_trigger    ! local FATES patch (fates/patches)
   use EDTypesMod,           only : area_inv
   use FatesPatchMod,        only : fates_patch_type
   use FatesCohortMod,       only : fates_cohort_type
   use shr_kind_mod,         only : r8 => shr_kind_r8
   implicit none
   private
   public :: rhessys_veg_do_all
   public :: rhessys_veg_init, rhessys_veg_init_site, rhessys_veg_coldstart, &
             rhessys_veg_set_time, rhessys_veg_set_bc, rhessys_veg_do, rhessys_veg_clean, &
             rhessys_veg_dump_cohorts, rhessys_veg_restart_write, rhessys_veg_restart_read

   ! --- Fortran mirrors of struct veg_bc / struct veg_time (rhessys_veg.h). ---
   ! Field order must match the C structs exactly.
   type, bind(C) :: veg_bc_t
      integer(c_int) :: nlevsoil
      integer(c_int) :: max_rooting_depth_index
      real(c_double) :: zi(VEG_NLEV_MAX+1)
      real(c_double) :: dz(VEG_NLEV_MAX)
      real(c_double) :: z(VEG_NLEV_MAX)
      real(c_double) :: watsat(VEG_NLEV_MAX)
      real(c_double) :: h2o_liqvol(VEG_NLEV_MAX)
      real(c_double) :: eff_porosity(VEG_NLEV_MAX)
      real(c_double) :: smp(VEG_NLEV_MAX)
      real(c_double) :: tempk(VEG_NLEV_MAX)
      real(c_double) :: forc_pbot, tgcm, t_veg, eair, esat_tv, oair, cair, rb
      real(c_double) :: dayl_factor, coszen
      real(c_double) :: solad(VEG_NSWB), solai(VEG_NSWB)
      real(c_double) :: albgr_dir(VEG_NSWB), albgr_dif(VEG_NSWB)
      real(c_double) :: snow_depth, frac_sno_eff, fcansno
      real(c_double) :: precip24, relhumid24, wind24
      real(c_double) :: dayl, t_night, esat_night, eair_night, t_avg
      real(c_double) :: w_scalar, t_scalar
      real(c_double) :: n_uptake_nh4, n_uptake_no3
      real(c_double) :: dist_kind, fire_mort(2), crown_burn(2), fuel_burn_leaf, fuel_burn_twig
   end type veg_bc_t

   type, bind(C) :: veg_out_t
      integer(c_int) :: npatches, ncohorts
      real(c_double) :: tlai, elai, tsai, htop
      real(c_double) :: laisun, laisha, rssun, rssha, gc_canopy, swabs_canopy
      real(c_double) :: btran, gpp_day, gpp_night, resp_m
      real(c_double) :: biomass_c, litter_c, lai_evg_frac
      real(c_double) :: ncap, n_uptake_recv
      real(c_double) :: lit_c_lab, lit_c_cel, lit_c_lig, lit_n_lab, lit_n_cel, lit_n_lig
      real(c_double) :: c_stock, n_stock
      real(c_double) :: mb_n_uptake, mb_n_frag, mb_n_err
      real(c_double) :: ls_c_lab, ls_c_cel, ls_c_lig, ls_n_lab, ls_n_cel, ls_n_lig
      real(c_double) :: lb_c_lab, lb_c_cel, lb_c_lig, lb_n_lab, lb_n_cel, lb_n_lig
      real(c_double) :: lit_split_err, leaf_litter_c, ag_cwd_c
      real(c_double) :: efflux_c, efflux_n
      real(c_double) :: fire_c, fire_n, harvest_c, harvest_n
      real(c_double) :: mb_n_other, mb_n_resize
   end type veg_out_t

   type, bind(C) :: veg_time_t
      integer(c_int) :: year, month, day, day_of_year, days_per_year
      integer(c_int) :: current_date, reference_date
      real(c_double) :: model_day
   end type veg_time_t

   ! All FATES sites (one per RHESSys patch) + their boundary conditions,
   ! the analogue of CLM's this%fates(nc) for a single clump.
   type(fates_interface_type), save, target :: fates
   integer, save :: nlevsoil_g = 0, nlevdecomp_g = 0
   integer, save :: n_setbc = 0, n_bad_total = 0
   logical, save :: cnp_coupled = .false.     ! mode 3: FATES CNP, N from RHESSys
   ! restart (Phase 5): FATES restart vectors written/read as a binary stream
   type(fates_restart_interface_type), save :: frest
   logical, save :: frest_ready = .false.
   logical, save :: restarting  = .false.     ! this run starts from a FATES restart
   character(len=8), parameter :: restart_magic = 'RHFATES1'
   integer, parameter :: max_bad_reports = 25

   logical, save :: fates_is_init = .false.
   integer, save :: fates_log_unit = -1        ! newunit= gives a NEGATIVE unit number
   logical, save :: fates_log_open = .false.

   ! Annual demography diagnostic (<output prefix>_fates_mort.csv): site-mean
   ! (unweighted over sites) stems, biomass, NPP, C lost per mortality mechanism,
   ! termination and impact mortality, by PFT x canopy layer (1 canopy, 2 understory).
   integer, parameter :: nmd = 21
   character(len=1024), save :: mort_file = ''
   logical, save :: mort_first = .true.
   integer, save :: mort_year = -1, mort_ndays = 0
   real(r8), allocatable, save :: mort_acc(:,:,:)      ! (nmd, pft, layer)
   real(r8), allocatable, save :: term_base(:,:,:)     ! (3 types + imort, pft, layer) cumulative at year start

   ! Codes FATES keeps private (LeafBiophysicsMod) -- values verified at the tag.
   integer, parameter :: ballberry_model = 1   ! stomatal_model
   integer, parameter :: net_assim_model = 1   ! stomatal_assim_model
   integer, parameter :: photosynth_acclim_model_none = 0
   ! SPITFIRE ignition-source codes: FATES only compares spitfire_mode against
   ! these, so they must be distinct; values follow CLM's fates_spitfire_mode
   ! namelist (0 no fire, 1 scalar lightning, 2 lightning data, 3 successful
   ! ignitions, 4 anthro ignitions). spitfire_mode = 0 keeps SPITFIRE off.
   integer, parameter :: sf_nofire = 0, sf_scalar_lightning = 1, &
                         sf_successful_ignitions = 3, sf_anthro_ignitions = 4

contains

!-------------------------------------------------------------------------------
! One-time init: FATES globals + parameter-file read + mode control parms.
! Phase 0 brings FATES up with ZERO sites/cohorts (no per-patch allocation yet).
!
! Order (constraints at the pinned tag, main/FatesInterfaceMod.F90):
!   FatesInterfaceInit -> flush_to_unset -> set ALL ctrlparms
!   -> SetFatesGlobalElements1 (reads the JSON paramfile; needs use_sp/nocomp/
!      luh2/masterproc/parteh_mode already set)
!   -> check_allset (reads EDPftvarcon_inst, so must follow the param read)
!   -> SetFatesGlobalElements2 -> InitTimeAveragingGlobals
!-------------------------------------------------------------------------------
subroutine rhessys_veg_init(c_paramfile, pf_len, c_invfile, inv_len, c_logfile, log_len, &
                            nlevsoil, nlevdecomp, demog_flag, nsites, is_restart) &
                            bind(C, name=_VEG_INIT)
   character(kind=c_char), intent(in) :: c_paramfile(*)
   integer(c_size_t),      intent(in) :: pf_len
   character(kind=c_char), intent(in) :: c_invfile(*)
   integer(c_size_t),      intent(in) :: inv_len
   character(kind=c_char), intent(in) :: c_logfile(*)
   integer(c_size_t),      intent(in) :: log_len
   integer(c_int),         intent(in) :: nlevsoil
   integer(c_int),         intent(in) :: nlevdecomp
   integer(c_int),         intent(in) :: demog_flag
   integer(c_int),         intent(in) :: nsites
   integer(c_int),         intent(in) :: is_restart   ! 1 = start from a FATES restart

   character(len=1024) :: paramfile, invfile, logfile
   integer :: i, s, numpft_surf, numcft

   ! --- C string -> Fortran string ---
   paramfile = ''
   do i = 1, min(int(pf_len), len(paramfile))
      paramfile(i:i) = c_paramfile(i)
   end do
   invfile = ''
   do i = 1, min(int(inv_len), len(invfile))
      invfile(i:i) = c_invfile(i)
   end do
   logfile = ''
   do i = 1, min(int(log_len), len(logfile))
      logfile(i:i) = c_logfile(i)
   end do

   ! FATES writes its (long) parameter dump + diagnostics to its own log file
   ! (<output prefix>_fates.log), keeping RHESSys stdout/stderr readable.
   open(newunit=fates_log_unit, file=trim(logfile), status='replace', action='write')
   fates_log_open = .true.
   i = index(logfile, '_fates.log', back=.true.)
   if (i > 0) then
      mort_file = logfile(1:i-1)//'_fates_mort.csv'
   else
      mort_file = trim(logfile)//'_mort.csv'
   end if
   call FatesInterfaceInit(fates_log_unit, .false.)

   ! --- control parameters: the full required set (flush -> set ALL) --------
   ! Mirrors CLM's CLMFatesGlobals1/2 (clmfates_interfaceMod.F90); every tag
   ! below is one set_fates_ctrlparms accepts at sci.1.93.2_api.46.0.0, and
   ! check_allset aborts if any required one is left unset.
   call set_fates_ctrlparms('flush_to_unset')

   ! process / run mode
   call set_fates_ctrlparms('masterproc', ival=itrue)            ! single process
   restarting = (is_restart == 1)
   if (restarting) then
      call set_fates_ctrlparms('is_restart', ival=itrue)
   else
      call set_fates_ctrlparms('is_restart', ival=ifalse)
   end if

   ! vegetation modes needed by SetFatesGlobalElements1
   call set_fates_ctrlparms('use_sp',            ival=ifalse)    ! satellite phenology OFF
   call set_fates_ctrlparms('use_nocomp',        ival=ifalse)    ! full competition
   call set_fates_ctrlparms('use_fixed_biogeog', ival=ifalse)
   call set_fates_ctrlparms('use_luh2',          ival=ifalse)    ! mgmt via RHESSys TEC (9.4b, 7)

   ! PARTEH allocation hypothesis. At this tag there is no N-only mode: the
   ! decided N-coupled config (9.4b) = carbon_nitrogen_phosphorus with P uptake
   ! prescribed in the paramfile (fates_cnp_prescribed_puptake = 1), as CLM
   ! does. Phase 0 runs carbon-only; Phase 4 switches.
   ! mode 3 (Phase 4): CNP with P uptake prescribed in the param file
   ! (fates_cnp_prescribed_puptake = 1) and N from RHESSys; modes 1-2 carbon-only.
   cnp_coupled = (demog_flag == 3)
   if (cnp_coupled) then
      call set_fates_ctrlparms('parteh_mode', ival=carbon_nitrogen_phosphorus)
   else
      call set_fates_ctrlparms('parteh_mode', ival=carbon_only)
   end if

   ! radiation bands (match the vis/NIR split RHESSys supplies in bc_in)
   call set_fates_ctrlparms('num_sw_bbands', ival=2)
   call set_fates_ctrlparms('vis_sw_index',  ival=1)
   call set_fates_ctrlparms('nir_sw_index',  ival=2)

   ! host identity + grid
   call set_fates_ctrlparms('hlm_name',       cval='CLM')        ! RHESSys masquerades as CLM host
   call set_fates_ctrlparms('num_lev_soil',   ival=nlevsoil)     ! soil-adaptor grid (9.4d)
   call set_fates_ctrlparms('hio_ignore_val', rval=-9999.0_c_double)
   call set_fates_ctrlparms('soilwater_ipedof', ival=0)          ! only used by FATES-Hydro (off)

   ! history (RHESSys owns output -> minimal FATES history, no seed dispersal)
   call set_fates_ctrlparms('seeddisp_cadence',    ival=0)
   call set_fates_ctrlparms('hist_hifrq_dimlevel', ival=0)
   call set_fates_ctrlparms('hist_dynam_dimlevel', ival=0)

   ! biogeochem: RHESSys owns decomposition -> FATES emits litter out (9.4c)
   call set_fates_ctrlparms('decomp_method', cval='NONE')
   call set_fates_ctrlparms('use_vertsoilc', ival=ifalse)
   call set_fates_ctrlparms('use_ch4',       ival=ifalse)        ! methane OFF
   ! N competition: relative demand (Phase 4 decision, 2026-10-12) -- the host splits
   ! mineral N between plants (capacity = fine-root C x Vmax) and microbes in
   ! RHESSys's resolve_sminn_competition. ECA possible later.
   call set_fates_ctrlparms('nu_com',        cval='RD')

   ! fire: WMFire drives fire -> SPITFIRE off (defs still required by check_allset)
   call set_fates_ctrlparms('spitfire_mode',               ival=sf_nofire)
   call set_fates_ctrlparms('sf_nofire_def',               ival=sf_nofire)
   call set_fates_ctrlparms('sf_scalar_lightning_def',     ival=sf_scalar_lightning)
   call set_fates_ctrlparms('sf_successful_ignitions_def', ival=sf_successful_ignitions)
   call set_fates_ctrlparms('sf_anthro_ignitions_def',     ival=sf_anthro_ignitions)
   call set_fates_ctrlparms('use_managed_fire',            ival=ifalse)

   ! stand structure / water stress
   call set_fates_ctrlparms('use_ed_st3',             ival=ifalse)  ! static-stand mode OFF
   call set_fates_ctrlparms('use_ed_prescribed_phys', ival=ifalse)
   call set_fates_ctrlparms('use_planthydro',         ival=ifalse)  ! FATES-Hydro OFF -> btran (9.4b)
   call set_fates_ctrlparms('hydr_solver',            ival=hydr_solver_1DTaylor) ! required even w/ hydro off
   call set_fates_ctrlparms('use_cohort_age_tracking',ival=ifalse)
   call set_fates_ctrlparms('use_tree_damage',        ival=ifalse)

   ! photosynthesis / physiology model selectors (FATES defaults)
   call set_fates_ctrlparms('photosynth_acclimation',      ival=photosynth_acclim_model_none)
   call set_fates_ctrlparms('use_daylength_factor_switch', ival=ifalse)  ! dayl_factor via bc_in instead
   call set_fates_ctrlparms('stomatal_model',          ival=ballberry_model)
   call set_fates_ctrlparms('stomatal_assim_model',    ival=net_assim_model)
   call set_fates_ctrlparms('maintresp_leaf_model',    ival=lmrmodel_ryan_1991)
   call set_fates_ctrlparms('mort_cstarvation_model',  ival=cstarvation_model_lin)
   call set_fates_ctrlparms('regeneration_model',      ival=default_regeneration)
   call set_fates_ctrlparms('radiation_model',         ival=norman_solver)
   call set_fates_ctrlparms('electron_transport_model',ival=FvCB1980)

   ! land use: management via RHESSys TEC, not LUH2 -> all OFF (9.4b, 7)
   call set_fates_ctrlparms('use_lu_harvest',            ival=ifalse)
   call set_fates_ctrlparms('num_lu_harvest_cats',       ival=0)
   ! logging on in mode 3 so host-scheduled harvests (TEC) can trigger it
   ! (fates_landuse_logging_event_code = 1 keeps FATES's own schedule off)
   if (demog_flag == 3) then
      call set_fates_ctrlparms('use_logging',            ival=itrue)
   else
      call set_fates_ctrlparms('use_logging',            ival=ifalse)
   end if
   call set_fates_ctrlparms('num_luh2_states',           ival=0)
   call set_fates_ctrlparms('num_luh2_transitions',      ival=0)
   call set_fates_ctrlparms('fates_lu_transition_logic', ival=1)   ! CLM default; must be 1-9 even with LUH off
   call set_fates_ctrlparms('use_fates_potentialveg',    ival=ifalse)

   ! initial stand: an inventory (-demog_inv, ED2-style pss/css files) or the
   ! FATES cold start (seedlings of every PFT at fates_recruit_init_density)
   if ((len_trim(invfile) > 0) .and. (.not. restarting)) then
      call set_fates_ctrlparms('use_inventory_init',  ival=itrue)
      call set_fates_ctrlparms('inventory_ctrl_file', cval=trim(invfile))
   else
      call set_fates_ctrlparms('use_inventory_init',  ival=ifalse)
      call set_fates_ctrlparms('inventory_ctrl_file', cval='none')
   end if
   call set_fates_ctrlparms('use_dbh_init',        ival=ifalse)

   ! --- parameter file (JSON at this tag) ---
   numpft_surf = 0 ; numcft = 0                        ! surface-dataset counts: SP mode only
   call SetFatesGlobalElements1(.true., numpft_surf, numcft, trim(paramfile))

   call set_fates_ctrlparms('check_allset')

   call SetFatesGlobalElements2(.true.)

   ! Host step size (s). Not a ctrlparm: CLM assigns the public global directly.
   ! RHESSys runs the FATES flux step once per day (design doc 9.4c option a);
   ! the running means (ema_24hr etc.) are defined from it.
   hlm_stepsize = 86400.0_c_double
   call InitTimeAveragingGlobals()

   ! --- sites + boundary conditions (CLM hlm_fates_interface_type%init) ---
   ! One site per RHESSys patch; RHESSys owns decomposition, so a single
   ! decomposition layer (use_vertsoilc off requires nlevdecomp = 1).
   call param_derived%Init(numpft)
   nlevsoil_g = nlevsoil
   nlevdecomp_g = nlevdecomp
   fates%nsites = nsites
   allocate(fates%sites(nsites))
   allocate(fates%bc_in(nsites))
   allocate(fates%bc_out(nsites))
   call allocate_bcpconst(fates%bc_pconst, nlevdecomp)
   call set_bcpconst(fates%bc_pconst, nlevdecomp)
   do s = 1, nsites
      ! no LU harvest / LUH2 / fixed biogeography: zero-size LU dims, surfpft bounds unused
      call allocate_bcin(fates%bc_in(s), nlevsoil, nlevdecomp, 0, 0, 0, 0, 0)
      call allocate_bcout(fates%bc_out(s), nlevsoil, nlevdecomp)
      call zero_bcs(fates, s)
   end do
   fates_is_init = .true.

   if (demog_flag < 1 .or. demog_flag > 3) then
      write(fates_log_unit,*) 'rhessys_fates: WARNING init called with demog_flag=', demog_flag
   end if
   write(fates_log_unit,*) 'rhessys_fates: FATES initialized: ', nsites, ' sites, ', nlevsoil, &
        ' soil layers; paramfile=', trim(paramfile), '; inventory=', trim(invfile)
   flush(fates_log_unit)
end subroutine rhessys_veg_init

!-------------------------------------------------------------------------------
! Static per-site setup (CLM init_soil_depths + init_coldstart's init_site_vars /
! zero_site). Soil geometry comes from the host's soil adaptor (design doc 9.4d).
! Patches/cohorts are not created yet (Phase 2-3).
!-------------------------------------------------------------------------------
subroutine rhessys_veg_init_site(site, patch_id, lat, bc) bind(C, name=_VEG_INIT_SITE)
   integer(c_int), intent(in) :: site
   integer(c_int), intent(in) :: patch_id
   real(c_double), intent(in) :: lat
   type(veg_bc_t), intent(in) :: bc
   integer :: n

   n = fates%bc_in(site)%nlevsoil
   if (bc%nlevsoil /= n) then
      write(fates_log_unit,*) 'rhessys_fates: site ', site, ' nlevsoil ', bc%nlevsoil, ' /= ', n
      error stop 'rhessys_veg_init_site: soil layer count mismatch'
   end if
   associate(bi => fates%bc_in(site))
      bi%zi_sisl(0:n) = bc%zi(1:n+1)
      bi%dz_sisl(1:n) = bc%dz(1:n)
      bi%z_sisl(1:n)  = bc%z(1:n)
      bi%decomp_id(1:n) = 1                  ! every soil layer -> the one decomposition layer
      bi%dz_decomp_sisl(1) = bc%zi(n+1)      ! that layer spans the column
      bi%max_rooting_depth_index_col = bc%max_rooting_depth_index
   end associate
   fates%sites(site)%h_gid = patch_id
   fates%sites(site)%lat = lat
   fates%sites(site)%lon = 0.0_r8           ! not used by FATES science
   ! on restart create_patchcohort_structure does this
   if (.not. restarting) then
      call init_site_vars(fates%sites(site), fates%bc_in(site), fates%bc_out(site))
      call zero_site(fates%sites(site))
   end if
end subroutine rhessys_veg_init_site

!-------------------------------------------------------------------------------
! Cold start of all sites (CLM init_coldstart): site phenology/disturbance
! state, patches + cohorts (inventory or seedlings), ed_update_site, canopy
! structure -> bc_out. out_all(1:nsites) gets each site's canopy structure.
!-------------------------------------------------------------------------------
subroutine rhessys_veg_coldstart(out_all) bind(C, name=_VEG_COLDSTART)
   type(veg_out_t), intent(out) :: out_all(*)
   integer :: s
   integer, allocatable :: fcolumn(:)

   do s = 1, fates%nsites
      out_all(s)%gpp_day = 0; out_all(s)%gpp_night = 0; out_all(s)%resp_m = 0
      out_all(s)%laisun = 0; out_all(s)%laisha = 0; out_all(s)%rssun = 0; out_all(s)%rssha = 0
      out_all(s)%gc_canopy = 0; out_all(s)%swabs_canopy = 0; out_all(s)%btran = 0
   end do

   call set_site_properties(fates%nsites, fates%sites, fates%bc_in)
   call init_patches(fates%nsites, fates%sites, fates%bc_in)
   do s = 1, fates%nsites
      fates%bc_in(s)%max_rooting_depth_index_col = fates%bc_in(s)%nlevsoil
      call ed_update_site(fates%sites(s), fates%bc_in(s), fates%bc_out(s), is_restarting=.false.)
   end do
   allocate(fcolumn(fates%nsites))
   fcolumn = [(s, s = 1, fates%nsites)]
   call UpdateFatesAvgSnowDepth(fates%sites, fates%bc_in)
   call canopy_summarization(fates%nsites, fates%sites, fates%bc_in)
   call update_hlm_dynamics(fates%nsites, fates%sites, fcolumn, fates%bc_out)
   deallocate(fcolumn)
   do s = 1, fates%nsites
      call sample_site(s, out_all(s), .false.)
   end do
   write(fates_log_unit,*) 'rhessys_fates: cold start done; site 1: ', out_all(1)%npatches, &
        ' patches, ', out_all(1)%ncohorts, ' cohorts, tlai ', out_all(1)%tlai, ' htop ', out_all(1)%htop
   flush(fates_log_unit)
end subroutine rhessys_veg_coldstart

!-------------------------------------------------------------------------------
! Area-weighted site summary over the FATES patches. with_fluxes: also the
! daylight radiation / conductance / btran values left by the last sub-step.
!-------------------------------------------------------------------------------
subroutine sample_site(s, out, with_fluxes)
   integer,         intent(in)    :: s
   type(veg_out_t), intent(inout) :: out
   logical,         intent(in)    :: with_fluxes
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort
   integer  :: ifp, ib
   real(r8) :: w, rb

   real(r8) :: tot, bio, lit, seed, wl, wl_evg
   wl = 0.0_r8; wl_evg = 0.0_r8
   out%npatches = 0; out%ncohorts = 0
   out%tlai = 0; out%elai = 0; out%tsai = 0; out%htop = 0
   if (with_fluxes) then
      out%laisun = 0; out%laisha = 0; out%rssun = 0; out%rssha = 0
      out%gc_canopy = 0; out%swabs_canopy = 0; out%btran = 0
   end if
   call SiteMassStock(fates%sites(s), 1, tot, bio, lit, seed)      ! element 1 = carbon
   out%biomass_c = bio * area_inv
   out%litter_c  = lit * area_inv
   associate(bi => fates%bc_in(s), bo => fates%bc_out(s))
      cpatch => fates%sites(s)%oldest_patch
      do while (associated(cpatch))
         ifp = cpatch%patchno
         w = cpatch%area * area_inv
         out%npatches = out%npatches + 1
         ccohort => cpatch%tallest
         do while (associated(ccohort))
            out%ncohorts = out%ncohorts + 1
            ! leaf area ~ leaf C x SLA (top), by leaf habit
            w = ccohort%prt%GetState(leaf_organ, carbon12_element) * ccohort%n * &
                prt_params%slatop(ccohort%pft)
            wl = wl + w
            if (prt_params%phen_leaf_habit(ccohort%pft) == ievergreen) wl_evg = wl_evg + w
            ccohort => ccohort%shorter
         end do
         w = cpatch%area * area_inv
         out%tlai = out%tlai + w * bo%tlai_pa(ifp)
         out%elai = out%elai + w * bo%elai_pa(ifp)
         out%tsai = out%tsai + w * bo%tsai_pa(ifp)
         out%htop = out%htop + w * bo%htop_pa(ifp)
         if (with_fluxes) then
            rb = bi%rb_pa(ifp)
            out%laisun = out%laisun + w * bo%laisun_pa(ifp)
            out%laisha = out%laisha + w * bo%laisha_pa(ifp)
            out%rssun  = out%rssun  + w * bo%rssun_pa(ifp)
            out%rssha  = out%rssha  + w * bo%rssha_pa(ifp)
            if (bo%laisun_pa(ifp) > 0.0_r8) &
               out%gc_canopy = out%gc_canopy + w * bo%laisun_pa(ifp) / (rb + bo%rssun_pa(ifp))
            if (bo%laisha_pa(ifp) > 0.0_r8) &
               out%gc_canopy = out%gc_canopy + w * bo%laisha_pa(ifp) / (rb + bo%rssha_pa(ifp))
            do ib = 1, VEG_NSWB
               out%swabs_canopy = out%swabs_canopy + w * &
                    (bi%solad_parb(ifp,ib) * bo%fabd_parb(ifp,ib) + &
                     bi%solai_parb(ifp,ib) * bo%fabi_parb(ifp,ib))
            end do
            out%btran = out%btran + w * bo%btran_pa(ifp)
         end if
         cpatch => cpatch%younger
      end do
   end associate
   out%lai_evg_frac = 0.0_r8
   if (wl > 0.0_r8) out%lai_evg_frac = wl_evg / wl
end subroutine sample_site

!-------------------------------------------------------------------------------
! Site GPP and maintenance respiration of the last sub-step (kgC/m2).
!-------------------------------------------------------------------------------
subroutine site_carbon(s, gpp, resp_m)
   integer,  intent(in)  :: s
   real(r8), intent(out) :: gpp, resp_m
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort
   gpp = 0.0_r8; resp_m = 0.0_r8
   cpatch => fates%sites(s)%oldest_patch
   do while (associated(cpatch))
      ccohort => cpatch%tallest
      do while (associated(ccohort))
         gpp    = gpp    + ccohort%gpp_tstep    * ccohort%n
         resp_m = resp_m + ccohort%resp_m_tstep * ccohort%n
         ccohort => ccohort%shorter
      end do
      cpatch => cpatch%younger
   end do
   gpp = gpp * area_inv
   resp_m = resp_m * area_inv
end subroutine site_carbon

!-------------------------------------------------------------------------------
! One flux sub-step of length dt (s) for site s, CLM order:
! canopy radiation -> sun/shade -> btran -> photosynthesis -> accumulate.
!-------------------------------------------------------------------------------
subroutine flux_substep(s, dt)
   integer,  intent(in) :: s
   real(r8), intent(in) :: dt
   call FatesNormalizedCanopyRadiation(fates%sites(s:s), fates%bc_in(s:s), fates%bc_out(s:s))
   call FatesSunShadeFracs(1, fates%sites(s:s), fates%bc_in(s:s), fates%bc_out(s:s))
   fates%bc_in(s)%filter_photo_pa(:) = 1
   call get_active_suction_layers(1, fates%sites(s:s), fates%bc_in(s:s), fates%bc_out(s:s))
   call btran_ed(1, fates%sites(s:s), fates%bc_in(s:s), fates%bc_out(s:s))
   fates%bc_in(s)%filter_photo_pa(:) = 2
   call FatesPlantRespPhotosynthDrive(1, fates%sites(s:s), fates%bc_in(s:s), fates%bc_out(s:s), dt)
   fates%bc_in(s)%filter_photo_pa(:) = 3
   call AccumulateFluxes_ED(1, fates%sites(s:s), fates%bc_in(s:s), fates%bc_out(s:s), dt)
end subroutine flux_substep

!-------------------------------------------------------------------------------
! Model time, once per day before the sites' boundary conditions.
!-------------------------------------------------------------------------------
subroutine rhessys_veg_set_time(t) bind(C, name=_VEG_SET_TIME)
   type(veg_time_t), intent(in) :: t
   call SetFatesTime(int(t%year), int(t%month), int(t%day), 0, &
                     int(t%current_date), int(t%reference_date), &
                     real(t%model_day, r8), int(t%day_of_year), &
                     int(t%days_per_year), 1.0_r8 / real(t%days_per_year, r8))
end subroutine rhessys_veg_set_time

!-------------------------------------------------------------------------------
! Daily boundary conditions for one site: copy struct veg_bc into bc_in, then
! validate every field the host owns (finite, physically plausible). n_bad =
! number of failed checks (design doc 12 Phase 1 exit: bc_in fully populated).
! Patch-indexed (_pa) arrays get the site value in every slot: one RHESSys
! patch = one forcing for all FATES patches of the site.
!-------------------------------------------------------------------------------
subroutine rhessys_veg_set_bc(site, bc, n_bad) bind(C, name=_VEG_SET_BC)
   integer(c_int), intent(in)  :: site
   type(veg_bc_t), intent(in)  :: bc
   integer(c_int), intent(out) :: n_bad
   integer :: n, j, ib

   n_bad = 0
   n_setbc = n_setbc + 1
   associate(bi => fates%bc_in(site))
      n = bi%nlevsoil
      call chk(bc%nlevsoil == n, 'nlevsoil', real(bc%nlevsoil, r8))

      ! --- soil (btran inputs; soil adaptor) ---
      bi%filter_btran = .true.
      bi%max_rooting_depth_index_col = bc%max_rooting_depth_index
      bi%smp_sl(1:n)          = bc%smp(1:n)
      bi%watsat_sl(1:n)       = bc%watsat(1:n)
      bi%eff_porosity_sl(1:n) = bc%eff_porosity(1:n)
      bi%h2o_liqvol_sl(1:n)   = bc%h2o_liqvol(1:n)
      bi%tempk_sl(1:n)        = bc%tempk(1:n)
      bi%t_soisno_sl(1:n)     = bc%tempk(1:n)
      ! decomposition scalars: FATES litter fragmentation always uses the host's
      ! (use_hlm_soil_scalar), so pass RHESSys's (patch cdf decomp_t/w_scalar)
      bi%w_scalar_sisl(1:n)   = bc%w_scalar
      bi%t_scalar_sisl(1:n)   = bc%t_scalar

      ! --- atmosphere / radiation (daylight means) ---
      bi%forc_pbot          = bc%forc_pbot
      bi%coszen             = bc%coszen
      bi%dayl_factor_pa(:)  = bc%dayl_factor
      bi%tgcm_pa(:)         = bc%tgcm
      bi%t_veg_pa(:)        = bc%t_veg
      bi%eair_pa(:)         = bc%eair
      bi%esat_tv_pa(:)      = bc%esat_tv
      bi%oair_pa(:)         = bc%oair
      bi%cair_pa(:)         = bc%cair
      bi%rb_pa(:)           = bc%rb
      do ib = 1, VEG_NSWB
         bi%solad_parb(:,ib) = bc%solad(ib)
         bi%solai_parb(:,ib) = bc%solai(ib)
         bi%albgr_dir_rb(ib) = bc%albgr_dir(ib)
         bi%albgr_dif_rb(ib) = bc%albgr_dif(ib)
      end do

      ! --- snow ---
      bi%snow_depth_si   = bc%snow_depth
      bi%frac_sno_eff_si = bc%frac_sno_eff
      bi%fcansno_pa(:)   = bc%fcansno

      ! --- 24 h means (fire BCs; SPITFIRE off, filled for completeness) ---
      bi%precip24_pa(:)   = bc%precip24
      bi%relhumid24_pa(:) = bc%relhumid24
      bi%wind24_pa(:)     = bc%wind24

      ! --- validation ---
      do j = 1, n
         call chk(bi%smp_sl(j) <= 0.0_r8 .and. bi%smp_sl(j) >= -1.0e8_r8*(1.0_r8+1.0e-12_r8), &
                  'smp_sl', bi%smp_sl(j))
         call chk(bi%watsat_sl(j) > 0.0_r8 .and. bi%watsat_sl(j) <= 1.0_r8, 'watsat_sl', bi%watsat_sl(j))
         call chk(bi%eff_porosity_sl(j) >= 0.01_r8 .and. bi%eff_porosity_sl(j) <= 1.0_r8, &
                  'eff_porosity_sl', bi%eff_porosity_sl(j))
         call chk(bi%h2o_liqvol_sl(j) >= 0.0_r8 .and. bi%h2o_liqvol_sl(j) <= bi%watsat_sl(j), &
                  'h2o_liqvol_sl', bi%h2o_liqvol_sl(j))
         call chk(bi%tempk_sl(j) > 200.0_r8 .and. bi%tempk_sl(j) < 340.0_r8, 'tempk_sl', bi%tempk_sl(j))
         call chk(bi%dz_sisl(j) > 0.0_r8, 'dz_sisl', bi%dz_sisl(j))
      end do
      call chk(bi%max_rooting_depth_index_col >= 1 .and. bi%max_rooting_depth_index_col <= n, &
               'max_rooting_depth_index_col', real(bi%max_rooting_depth_index_col, r8))
      call chk(bi%forc_pbot > 3.0e4_r8 .and. bi%forc_pbot < 1.1e5_r8, 'forc_pbot', bi%forc_pbot)
      call chk(bi%coszen >= 0.0_r8 .and. bi%coszen <= 1.0_r8, 'coszen', bi%coszen)
      call chk(bc%dayl_factor >= 0.0_r8 .and. bc%dayl_factor <= 1.0_r8, 'dayl_factor_pa', bc%dayl_factor)
      call chk(bc%tgcm > 200.0_r8 .and. bc%tgcm < 340.0_r8, 'tgcm_pa', bc%tgcm)
      call chk(bc%t_veg > 200.0_r8 .and. bc%t_veg < 340.0_r8, 't_veg_pa', bc%t_veg)
      call chk(bc%esat_tv > 0.0_r8, 'esat_tv_pa', bc%esat_tv)
      call chk(bc%eair > 0.0_r8 .and. bc%eair <= bc%esat_tv, 'eair_pa', bc%eair)
      call chk(bc%oair > 0.0_r8, 'oair_pa', bc%oair)
      call chk(bc%cair > 0.0_r8 .and. bc%cair < 0.01_r8*bc%forc_pbot, 'cair_pa', bc%cair)
      call chk(bc%rb > 0.0_r8 .and. bc%rb < 1.0e4_r8, 'rb_pa', bc%rb)
      do ib = 1, VEG_NSWB
         call chk(bc%solad(ib) >= 0.0_r8 .and. bc%solad(ib) < 1500.0_r8, 'solad_parb', bc%solad(ib))
         call chk(bc%solai(ib) >= 0.0_r8 .and. bc%solai(ib) < 1500.0_r8, 'solai_parb', bc%solai(ib))
         call chk(bc%albgr_dir(ib) >= 0.0_r8 .and. bc%albgr_dir(ib) <= 1.0_r8, 'albgr_dir_rb', bc%albgr_dir(ib))
         call chk(bc%albgr_dif(ib) >= 0.0_r8 .and. bc%albgr_dif(ib) <= 1.0_r8, 'albgr_dif_rb', bc%albgr_dif(ib))
      end do
      call chk(bi%snow_depth_si >= 0.0_r8, 'snow_depth_si', bi%snow_depth_si)
      call chk(bi%frac_sno_eff_si >= 0.0_r8 .and. bi%frac_sno_eff_si <= 1.0_r8, 'frac_sno_eff_si', bi%frac_sno_eff_si)
      call chk(bc%fcansno >= 0.0_r8 .and. bc%fcansno <= 1.0_r8, 'fcansno_pa', bc%fcansno)
      call chk(bc%precip24 >= 0.0_r8, 'precip24_pa', bc%precip24)
      call chk(bc%relhumid24 >= 0.0_r8 .and. bc%relhumid24 <= 100.0_r8 + 1.0e-9_r8, 'relhumid24_pa', bc%relhumid24)
      call chk(bc%wind24 >= 0.0_r8, 'wind24_pa', bc%wind24)
   end associate
   n_bad_total = n_bad_total + n_bad

contains

   ! a check fails if the condition is false or the value is not finite
   subroutine chk(ok, name, val)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      logical,          intent(in) :: ok
      character(len=*), intent(in) :: name
      real(r8),         intent(in) :: val
      if (ok .and. ieee_is_finite(val)) return
      n_bad = n_bad + 1
      if (n_bad_total + n_bad <= max_bad_reports) then
         write(fates_log_unit,'(a,i0,a,i0,3a,es14.6)') 'rhessys_fates: bc_in check failed: site ', site, &
              ' (patch ', fates%sites(site)%h_gid, ') ', trim(name), ' = ', val
      end if
   end subroutine chk

end subroutine rhessys_veg_set_bc

!-------------------------------------------------------------------------------
! Daily flux step for one site (design doc 9.4c option a, 12 Phase 2): RHESSys
! is a daily model, so the day is two FATES sub-steps -- daylight (daylight-
! mean forcing from rhessys_veg_set_bc, dt = daylength) and night (no
! shortwave, night temperature, dt = 86400 - daylength) -- then one update of
! the running means with the daily mean temperature. Demography (Phase 3) is
! not run: the stand and its leaf status stay as initialized.
!-------------------------------------------------------------------------------
subroutine rhessys_veg_do(site, bc, out) bind(C, name=_VEG_DO)
   integer(c_int),  intent(in)  :: site
   type(veg_bc_t),  intent(in)  :: bc
   type(veg_out_t), intent(out) :: out
   integer  :: s, ib
   real(r8) :: dt_day, dt_night, gpp, resp_m

   s = site
   if (s < 1 .or. s > fates%nsites) error stop 'rhessys_veg_do: bad site'
   dt_day   = min(max(bc%dayl, 0.0_r8), 86400.0_r8)
   dt_night = 86400.0_r8 - dt_day
   out%gpp_day = 0; out%gpp_night = 0; out%resp_m = 0
   call sample_site(s, out, .true.)               ! zero/structure before the steps

   ! --- daylight sub-step (bc_in holds the daylight forcing) ---
   if (dt_day > 0.0_r8) then
      call flux_substep(s, dt_day)
      call sample_site(s, out, .true.)            ! daylight radiation/conductance/btran
      call site_carbon(s, gpp, resp_m)
      out%gpp_day = gpp
      out%resp_m = out%resp_m + resp_m
   end if

   ! --- night sub-step ---
   if (dt_night > 0.0_r8) then
      associate(bi => fates%bc_in(s))
         bi%coszen          = 0.0_r8
         do ib = 1, VEG_NSWB
            bi%solad_parb(:,ib) = 0.0_r8
            bi%solai_parb(:,ib) = 0.0_r8
         end do
         bi%t_veg_pa(:)     = bc%t_night
         bi%tgcm_pa(:)      = bc%t_night
         bi%esat_tv_pa(:)   = bc%esat_night
         bi%eair_pa(:)      = bc%eair_night
      end associate
      call flux_substep(s, dt_night)
      call site_carbon(s, gpp, resp_m)
      out%gpp_night = gpp
      out%resp_m = out%resp_m + resp_m
   end if

   ! --- running means (hlm_stepsize = 1 day): daily mean temperature ---
   fates%bc_in(s)%t_veg_pa(:) = bc%t_avg
   call UpdateFatesRMeansTStep(fates%sites(s:s), fates%bc_in(s:s), fates%bc_out(s:s))

   ! --- daily demography (Phase 3; CLM dynamics_driv): phenology, growth,
   !     allocation, mortality, recruitment, disturbance, patch/cohort fusion.
   !     Carbon-only PARTEH: nutrients supplied (as CLM with carbon-only). ---
   if (cnp_coupled) then
      ! RHESSys-granted uptake -> cohorts (competitor = cohort), share by capacity
      call distribute_n_uptake(s, bc%n_uptake_nh4, bc%n_uptake_no3, out%n_uptake_recv)
      call UnPackNutrientAquisitionBCs(fates%sites(s:s), fates%bc_in(s:s), .false., .true.)
   else
      out%n_uptake_recv = 0.0_r8
      call UnPackNutrientAquisitionBCs(fates%sites(s:s), fates%bc_in(s:s), .true., .true.)
   end if

   ! Phase 6: host-imposed fire / harvest for today (always resets yesterday's)
   call apply_disturbance(s, bc)
   call ed_ecosystem_dynamics(fates%sites(s), fates%bc_in(s), fates%bc_out(s))
   call ed_update_site(fates%sites(s), fates%bc_in(s), fates%bc_out(s), is_restarting=.false.)
   call UpdateFatesAvgSnowDepth(fates%sites(s:s), fates%bc_in(s:s))
   call canopy_summarization(1, fates%sites(s:s), fates%bc_in(s:s))
   call update_hlm_dynamics(1, fates%sites(s:s), [s], fates%bc_out(s:s))
   call sample_site(s, out, .false.)              ! end-of-day canopy structure + stocks
   call nutrient_out(s, out)                      ! capacity, litter to host, total stocks
end subroutine rhessys_veg_do

!-------------------------------------------------------------------------------
! All sites' daily steps in one call, OpenMP over sites (sites are independent
! within a day; each thread works on its own site's FATES structures). Thread
! count from OMP_NUM_THREADS. Results do not depend on the thread count.
!-------------------------------------------------------------------------------
subroutine rhessys_veg_do_all(nsites, bcs, outs) bind(C, name=_VEG_DO_ALL)
   integer(c_int),  intent(in)    :: nsites
   type(veg_bc_t),  intent(in)    :: bcs(nsites)
   type(veg_out_t), intent(inout) :: outs(nsites)
   integer(c_int) :: s
   !$omp parallel do schedule(dynamic, 4) default(shared) private(s)
   do s = 1, nsites
      call rhessys_veg_do(s, bcs(s), outs(s))
   end do
   !$omp end parallel do
   call mort_accumulate(int(nsites))
end subroutine rhessys_veg_do_all

!-------------------------------------------------------------------------------
! Demography diagnostic: cumulative termination / impact-mortality C over all
! sites, kgC/m2 summed over sites (FATES zeroes these only in its history code,
! which RHESSys does not call, so they are cumulative since init).
!-------------------------------------------------------------------------------
subroutine mort_term_cum(nsites, cum)
   integer,  intent(in)  :: nsites
   real(r8), intent(out) :: cum(:,:,:)
   integer :: s, ft
   cum = 0.0_r8
   do s = 1, nsites
      do ft = 1, numpft
         cum(1,ft,1) = cum(1,ft,1) + fates%sites(s)%term_carbonflux_canopy(i_term_mort_type_cstarv,ft)  * area_inv
         cum(2,ft,1) = cum(2,ft,1) + fates%sites(s)%term_carbonflux_canopy(i_term_mort_type_canlev,ft)  * area_inv
         cum(3,ft,1) = cum(3,ft,1) + fates%sites(s)%term_carbonflux_canopy(i_term_mort_type_numdens,ft) * area_inv
         cum(1,ft,2) = cum(1,ft,2) + fates%sites(s)%term_carbonflux_ustory(i_term_mort_type_cstarv,ft)  * area_inv
         cum(2,ft,2) = cum(2,ft,2) + fates%sites(s)%term_carbonflux_ustory(i_term_mort_type_canlev,ft)  * area_inv
         cum(3,ft,2) = cum(3,ft,2) + fates%sites(s)%term_carbonflux_ustory(i_term_mort_type_numdens,ft) * area_inv
         ! imort_carbonflux: daily additions in gC/m2/s -> kgC/m2 per day (understory killed by treefall)
         cum(4,ft,2) = cum(4,ft,2) + fates%sites(s)%imort_carbonflux(ft) * 86400.0_r8 * 1.0e-3_r8
      end do
   end do
end subroutine mort_term_cum

subroutine mort_write(nsites)
   integer, intent(in) :: nsites
   real(r8) :: cum(4, numpft, 2), d, fd
   integer :: u, ft, il
   if (len_trim(mort_file) == 0 .or. mort_ndays == 0) return
   call mort_term_cum(nsites, cum)
   if (mort_first) then
      open(newunit=u, file=trim(mort_file), status='replace', action='write')
      write(u,'(a)') 'year,ndays,pft,layer,stems_per_ha,biomass_kgC_m2,npp_kgC_m2_yr,'// &
           'mort_bg,mort_cstarv,mort_hydr,mort_freeze,mort_senes,mort_agesen,mort_damage,'// &
           'term_cstarv,term_canlev,term_numdens,imort,'// &
           'store_frac,frac_C_lim,frac_N_lim,n_demand_gN_m2_yr,n_gain_gN_m2_yr,c_efflux_kgC_m2_yr,leafc_kgC_m2'
      mort_first = .false.
   else
      open(newunit=u, file=trim(mort_file), status='old', position='append', action='write')
   end if
   d  = real(mort_ndays, r8) * real(nsites, r8)    ! site-days: state means
   fd = 365.0_r8 / real(mort_ndays, r8) / real(nsites, r8)   ! flux sums -> per year, site mean
   do ft = 1, numpft
      do il = 1, 2
         associate(a => mort_acc(:, ft, il))
         write(u,'(i0,",",i0,",",i0,",",i0,22(",",es12.5))') mort_year, mort_ndays, ft, il, &
              a(1)/d*1.0e4_r8, a(2)/d, a(3)*fd, a(4:10)*fd, &
              (cum(1:4,ft,il) - term_base(1:4,ft,il)) * fd, &
              a(11)/max(a(2),1.0e-30_r8), a(12)/max(a(14),1.0e-30_r8), a(13)/max(a(14),1.0e-30_r8), &
              a(15)*fd*1.0e3_r8, a(16)*fd*1.0e3_r8, a(17)*fd, a(18)/d
         end associate
      end do
   end do
   close(u)
   term_base = cum
end subroutine mort_write

subroutine mort_accumulate(nsites)
   integer, intent(in) :: nsites
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort
   real(r8) :: nm2, ctot, leafc, fnrtc, storec
   integer  :: s, ft, il
   real(r8) :: cum(4, numpft, 2)

   if (len_trim(mort_file) == 0) return
   if (.not. allocated(mort_acc)) then
      allocate(mort_acc(nmd, numpft, 2), term_base(4, numpft, 2))
      mort_acc = 0.0_r8
      call mort_term_cum(nsites, cum)
      term_base = cum
      mort_year = hlm_current_year
   end if
   if (hlm_current_year /= mort_year) then
      call mort_write(nsites)
      mort_acc = 0.0_r8
      mort_ndays = 0
      mort_year = hlm_current_year
   end if
   mort_ndays = mort_ndays + 1
   do s = 1, nsites
      cpatch => fates%sites(s)%oldest_patch
      do while (associated(cpatch))
         ccohort => cpatch%tallest
         do while (associated(ccohort))
            ft = ccohort%pft
            il = min(max(ccohort%canopy_layer, 1), 2)
            nm2 = ccohort%n * area_inv                       ! plants per m2 of site
            leafc  = ccohort%prt%GetState(leaf_organ,  carbon12_element)
            fnrtc  = ccohort%prt%GetState(fnrt_organ,  carbon12_element)
            storec = ccohort%prt%GetState(store_organ, carbon12_element)
            ctot = leafc + fnrtc + storec + ccohort%prt%GetState(sapw_organ, carbon12_element) + &
                   ccohort%prt%GetState(struct_organ, carbon12_element) + &
                   ccohort%prt%GetState(repro_organ, carbon12_element)
            associate(a => mort_acc(:, ft, il))
            a(1) = a(1) + nm2
            a(2) = a(2) + nm2 * ctot
            a(3) = a(3) + nm2 * fin(ccohort%npp_acc_hold) / 365.0_r8
            a(4)  = a(4)  + nm2 * ctot * fin(ccohort%bmort)  / 365.0_r8
            a(5)  = a(5)  + nm2 * ctot * fin(ccohort%cmort)  / 365.0_r8
            a(6)  = a(6)  + nm2 * ctot * fin(ccohort%hmort)  / 365.0_r8
            a(7)  = a(7)  + nm2 * ctot * fin(ccohort%frmort) / 365.0_r8
            a(8)  = a(8)  + nm2 * ctot * fin(ccohort%smort)  / 365.0_r8
            a(9)  = a(9)  + nm2 * ctot * fin(ccohort%asmort) / 365.0_r8
            a(10) = a(10) + nm2 * ctot * fin(ccohort%dgmort) / 365.0_r8
            a(11) = a(11) + nm2 * storec
            a(14) = a(14) + nm2 * ctot
            a(18) = a(18) + nm2 * leafc
            if (cnp_coupled) then
               if (ccohort%cnp_limiter == 1) a(12) = a(12) + nm2 * ctot
               if (ccohort%cnp_limiter == 2) a(13) = a(13) + nm2 * ctot
               if (fin(ccohort%daily_n_demand) > 0.0_r8) a(15) = a(15) + nm2 * fin(ccohort%daily_n_demand)
               a(16) = a(16) + nm2 * fin(ccohort%daily_n_gain)
               a(17) = a(17) + nm2 * fin(ccohort%daily_c_efflux)
            end if
            end associate
            ccohort => ccohort%shorter
         end do
         cpatch => cpatch%younger
      end do
   end do
contains
   real(r8) function fin(x)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(r8), intent(in) :: x
      fin = 0.0_r8
      if (ieee_is_finite(x)) then
         if (abs(x) < 1.0e30_r8) fin = x
      end if
   end function fin
end subroutine mort_accumulate

!-------------------------------------------------------------------------------
! Split the site's daily N uptake (gN/m2) over the cohorts (FATES competitors in
! coupled N uptake mode: one per cohort, in oldest-patch / tallest-cohort order)
! in proportion to their RD capacity, fine-root C x (Vmax NH4 + Vmax NO3).
! recv = total actually placed (0 if no cohort has fine roots).
!-------------------------------------------------------------------------------
subroutine distribute_n_uptake(s, nh4, no3, recv)
   integer,  intent(in)  :: s
   real(r8), intent(in)  :: nh4, no3
   real(r8), intent(out) :: recv
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort
   real(r8) :: wsum, w
   integer  :: icomp, pass

   recv = 0.0_r8
   wsum = 0.0_r8
   associate(bi => fates%bc_in(s))
      bi%plant_nh4_uptake_flux(:,:) = 0.0_r8
      bi%plant_no3_uptake_flux(:,:) = 0.0_r8
      do pass = 1, 2
         icomp = 0
         cpatch => fates%sites(s)%oldest_patch
         do while (associated(cpatch))
            ccohort => cpatch%tallest
            do while (associated(ccohort))
               icomp = icomp + 1
               w = ccohort%prt%GetState(fnrt_organ, carbon12_element) * ccohort%n * &
                   (EDPftvarcon_inst%vmax_nh4(ccohort%pft) + EDPftvarcon_inst%vmax_no3(ccohort%pft))
               if (pass == 1) then
                  wsum = wsum + w
               else if (wsum > 0.0_r8 .and. icomp <= size(bi%plant_nh4_uptake_flux, 1)) then
                  bi%plant_nh4_uptake_flux(icomp,1) = nh4 * w / wsum
                  bi%plant_no3_uptake_flux(icomp,1) = no3 * w / wsum
                  recv = recv + (nh4 + no3) * w / wsum
               end if
               ccohort => ccohort%shorter
            end do
            cpatch => cpatch%younger
         end do
      end do
   end associate
end subroutine distribute_n_uptake

!-------------------------------------------------------------------------------
! End-of-day nutrient exchange terms: RD capacity for tomorrow's competition,
! today's FATES litter flux to the host (g/m3/s over the single decomposition
! layer -> g/m2/day), total C and N stocks (budget check).
!-------------------------------------------------------------------------------
subroutine nutrient_out(s, out)
   integer,         intent(in)    :: s
   type(veg_out_t), intent(inout) :: out
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort
   real(r8) :: f, tot, bio, lit, seed

   out%ncap = 0.0_r8
   cpatch => fates%sites(s)%oldest_patch
   do while (associated(cpatch))
      ccohort => cpatch%tallest
      do while (associated(ccohort))
         out%ncap = out%ncap + ccohort%prt%GetState(fnrt_organ, carbon12_element) * &
              ccohort%n * area_inv * 1000.0_r8 * 86400.0_r8 * &
              (EDPftvarcon_inst%vmax_nh4(ccohort%pft) + EDPftvarcon_inst%vmax_no3(ccohort%pft))
         ccohort => ccohort%shorter
      end do
      cpatch => cpatch%younger
   end do
   associate(bo => fates%bc_out(s), bi => fates%bc_in(s))
      f = bi%dz_decomp_sisl(1) * 86400.0_r8
      out%lit_c_lab = bo%litt_flux_lab_c_si(1) * f
      out%lit_c_cel = bo%litt_flux_cel_c_si(1) * f
      out%lit_c_lig = bo%litt_flux_lig_c_si(1) * f
      if (cnp_coupled) then
         out%lit_n_lab = bo%litt_flux_lab_n_si(1) * f
         out%lit_n_cel = bo%litt_flux_cel_n_si(1) * f
         out%lit_n_lig = bo%litt_flux_lig_n_si(1) * f
      else
         out%lit_n_lab = 0.0_r8; out%lit_n_cel = 0.0_r8; out%lit_n_lig = 0.0_r8
      end if
   end associate
   call litter_split(s, out)

   ! root efflux from FATES's mass-balance record: C net_root_uptake is only
   ! -efflux; N efflux = uptake received - net root uptake
   out%efflux_c = max(-fates%sites(s)%mass_balance(1)%net_root_uptake, 0.0_r8) * area_inv * 1000.0_r8
   out%efflux_n = 0.0_r8
   if (cnp_coupled) out%efflux_n = max(out%n_uptake_recv - &
        fates%sites(s)%mass_balance(2)%net_root_uptake * area_inv * 1000.0_r8, 0.0_r8)

   call SiteMassStock(fates%sites(s), 1, tot, bio, lit, seed)
   out%c_stock = tot * area_inv
   out%n_stock = 0.0_r8
   out%mb_n_uptake = 0.0_r8; out%mb_n_frag = 0.0_r8; out%mb_n_err = 0.0_r8
   if (cnp_coupled) then
      call SiteMassStock(fates%sites(s), 2, tot, bio, lit, seed)     ! element 2 = N in CNP
      out%n_stock = tot * area_inv
      associate(mb => fates%sites(s)%mass_balance(2))
         out%mb_n_uptake = mb%net_root_uptake * area_inv
         out%mb_n_frag   = mb%frag_out * area_inv
         out%mb_n_err    = mb%err_fates * area_inv
      end associate
   end if
   ! Phase 6: fire and harvest losses (FATES mass-balance record, kg/site/day)
   associate(mc => fates%sites(s)%mass_balance(1))
      out%fire_c    = sum(mc%burn_flux_to_atm(:)) * area_inv * 1000.0_r8
      out%harvest_c = sum(mc%wood_product_harvest(:)) * area_inv * 1000.0_r8
   end associate
   out%fire_n = 0.0_r8; out%harvest_n = 0.0_r8; out%mb_n_other = 0.0_r8; out%mb_n_resize = 0.0_r8
   if (cnp_coupled) then
      associate(mn => fates%sites(s)%mass_balance(2))
         out%fire_n    = sum(mn%burn_flux_to_atm(:)) * area_inv * 1000.0_r8
         out%harvest_n = sum(mn%wood_product_harvest(:)) * area_inv * 1000.0_r8
         out%mb_n_resize = mn%patch_resize_err * area_inv
         out%mb_n_other  = (mn%seed_in - mn%seed_out + mn%flux_generic_in - mn%flux_generic_out &
              + mn%patch_resize_err - mn%herbivory_flux_out - sum(mn%wood_product_landusechange(:))) &
              * area_inv
      end associate
   end if
end subroutine nutrient_out

!-------------------------------------------------------------------------------
! Phase 6: host-imposed disturbance for site s today, set just before
! ed_ecosystem_dynamics so FATES's own disturbance code (disturbance_rates ->
! spawn_patches -> TotalBalanceCheck) applies it. With SPITFIRE off nothing
! else sets or clears these fields, so every day they are first reset.
!  fire (uniform severity): every patch burns (frac_burnt 0.98, FATES's warning
!    limit), each cohort dies at fire_mort of its leaf habit, crowns of killed
!    and surviving trees burn at crown_burn, FATES dead leaves / twigs burn at
!    the fuel fractions; large CWD (classes 3-4) does not burn.
!  harvest: the host logging trigger (local FATES patch) makes today logging
!    time for this site; intensity from the logging parameters (clear-cut).
!-------------------------------------------------------------------------------
subroutine apply_disturbance(s, bc)
   integer,        intent(in) :: s
   type(veg_bc_t), intent(in) :: bc
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort
   integer :: kind, h
   kind = nint(bc%dist_kind)
   host_logging_trigger = (kind == 2)
   cpatch => fates%sites(s)%oldest_patch
   do while (associated(cpatch))
      cpatch%fire = 0
      cpatch%frac_burnt = 0.0_r8
      cpatch%fuel%frac_burnt(:) = 0.0_r8
      if (kind == 1) then
         cpatch%fire = 1
         cpatch%frac_burnt = 0.98_r8
         cpatch%fuel%frac_burnt(1:2) = bc%fuel_burn_twig
         cpatch%fuel%frac_burnt(5)   = bc%fuel_burn_leaf
         cpatch%fuel%frac_burnt(6)   = bc%fuel_burn_leaf
      end if
      ccohort => cpatch%tallest
      do while (associated(ccohort))
         ccohort%fire_mort = 0.0_r8
         ccohort%fraction_crown_burned = 0.0_r8
         if (kind == 1) then
            h = 2
            if (prt_params%phen_leaf_habit(ccohort%pft) == ievergreen) h = 1
            ccohort%fire_mort = min(max(bc%fire_mort(h), 0.0_r8), 1.0_r8)
            ccohort%fraction_crown_burned = min(max(bc%crown_burn(h), 0.0_r8), 1.0_r8)
         end if
         ccohort => ccohort%shorter
      end do
      cpatch => cpatch%younger
   end do
end subroutine apply_disturbance

!-------------------------------------------------------------------------------
! Split today's FATES litter handoff into surface and below-ground parts
! (g/m2/day), term by term as FatesSoilBGCFluxMod::FluxIntoLitterPools builds
! the merged litt_flux (one decomposition layer: surface profile = 1):
!   surface      = leaf fines frag (lab/cel/lig) + above-ground CWD frag
!                  (cel/lig by fates_frag_cwd_fcel/flig) + seed decay (leaf fractions)
!   below ground = fine-root fines frag (lab/cel/lig) + below-ground CWD frag
! lit_split_err checks the split against the merged flux. Also the standing
! unfragmented above-ground leaf litter and CWD stocks (kgC/m2).
!-------------------------------------------------------------------------------
subroutine litter_split(s, out)
   integer,         intent(in)    :: s
   type(veg_out_t), intent(inout) :: out
   type(fates_patch_type), pointer :: cpatch
   type(litter_type),      pointer :: litt
   real(r8) :: sur(3,2), blw(3,2), af, sd
   integer  :: el, ne, ic, j, ipft, nlev_eff

   sur = 0.0_r8; blw = 0.0_r8
   out%leaf_litter_c = 0.0_r8; out%ag_cwd_c = 0.0_r8
   ne = 1
   if (cnp_coupled) ne = 2
   nlev_eff = max(fates%bc_in(s)%max_rooting_depth_index_col, 1)
   do el = 1, ne
      cpatch => fates%sites(s)%oldest_patch
      do while (associated(cpatch))
         litt => cpatch%litter(el)
         af = cpatch%area * area_inv
         do ic = 1, ncwd
            sur(2,el) = sur(2,el) + litt%ag_cwd_frag(ic) * ED_val_cwd_fcel * af
            sur(3,el) = sur(3,el) + litt%ag_cwd_frag(ic) * ED_val_cwd_flig * af
            do j = 1, nlev_eff
               blw(2,el) = blw(2,el) + litt%bg_cwd_frag(ic,j) * ED_val_cwd_fcel * af
               blw(3,el) = blw(3,el) + litt%bg_cwd_frag(ic,j) * ED_val_cwd_flig * af
            end do
         end do
         sur(1,el) = sur(1,el) + litt%leaf_fines_frag(ilabile)    * af
         sur(2,el) = sur(2,el) + litt%leaf_fines_frag(icellulose) * af
         sur(3,el) = sur(3,el) + litt%leaf_fines_frag(ilignin)    * af
         do ipft = 1, numpft
            sd = (litt%seed_decay(ipft) + litt%seed_germ_decay(ipft)) * af
            sur(1,el) = sur(1,el) + sd * EDPftvarcon_inst%lf_flab(ipft)
            sur(2,el) = sur(2,el) + sd * EDPftvarcon_inst%lf_fcel(ipft)
            sur(3,el) = sur(3,el) + sd * EDPftvarcon_inst%lf_flig(ipft)
         end do
         do j = 1, nlev_eff
            blw(1,el) = blw(1,el) + litt%root_fines_frag(ilabile,j)    * af
            blw(2,el) = blw(2,el) + litt%root_fines_frag(icellulose,j) * af
            blw(3,el) = blw(3,el) + litt%root_fines_frag(ilignin,j)    * af
         end do
         if (el == 1) then
            out%leaf_litter_c = out%leaf_litter_c + sum(litt%leaf_fines(:)) * af
            out%ag_cwd_c      = out%ag_cwd_c + sum(litt%ag_cwd(:)) * af
         end if
         cpatch => cpatch%younger
      end do
   end do
   sur = sur * 1000.0_r8                     ! kg/m2/day -> g/m2/day
   blw = blw * 1000.0_r8
   out%ls_c_lab = sur(1,1); out%ls_c_cel = sur(2,1); out%ls_c_lig = sur(3,1)
   out%lb_c_lab = blw(1,1); out%lb_c_cel = blw(2,1); out%lb_c_lig = blw(3,1)
   out%ls_n_lab = sur(1,2); out%ls_n_cel = sur(2,2); out%ls_n_lig = sur(3,2)
   out%lb_n_lab = blw(1,2); out%lb_n_cel = blw(2,2); out%lb_n_lig = blw(3,2)
   out%lit_split_err = abs(out%ls_c_lab + out%lb_c_lab - out%lit_c_lab) + &
                       abs(out%ls_c_cel + out%lb_c_cel - out%lit_c_cel) + &
                       abs(out%ls_c_lig + out%lb_c_lig - out%lit_c_lig) + &
                       abs(out%ls_n_lab + out%lb_n_lab - out%lit_n_lab) + &
                       abs(out%ls_n_cel + out%lb_n_cel - out%lit_n_cel) + &
                       abs(out%ls_n_lig + out%lb_n_lig - out%lit_n_lig)
end subroutine litter_split

!-------------------------------------------------------------------------------
! Restart (Phase 5). FATES's restart vectors (FatesRestartInterfaceMod, one
! "column" per site, fates_maxElementsPerSite cohort slots per site) are written
! as a binary stream: magic, nsites, slots, nvars, then per variable its name,
! kind (0 = real, 1 = integer), length and data. CLM does the same through
! netCDF; RHESSys writes it next to its .state file at output_current_state.
!-------------------------------------------------------------------------------
subroutine setup_restart()
   type(fates_bounds_type) :: b
   integer :: s, n
   if (frest_ready) return
   n = fates%nsites
   b%cohort_begin = 1
   b%cohort_end   = n * fates_maxElementsPerSite
   b%column_begin = 1
   b%column_end   = n
   call frest%Init(1, b)
   call frest%SetThreadBoundsEach(1, b)
   allocate(frest%restart_map(1)%site_index(n))
   allocate(frest%restart_map(1)%cohort1_index(n))
   do s = 1, n
      frest%restart_map(1)%site_index(s)    = s
      frest%restart_map(1)%cohort1_index(s) = (s - 1) * fates_maxElementsPerSite + 1
   end do
   call frest%assemble_restart_output_types()
   call frest%initialize_restart_vars()
   frest_ready = .true.
end subroutine setup_restart

subroutine c2f_string(c, clen, f)
   character(kind=c_char), intent(in)  :: c(*)
   integer(c_size_t),      intent(in)  :: clen
   character(len=*),       intent(out) :: f
   integer :: i
   f = ''
   do i = 1, min(int(clen), len(f))
      f(i:i) = c(i)
   end do
end subroutine c2f_string

subroutine rhessys_veg_restart_write(c_file, flen) bind(C, name=_VEG_RESTART_WRITE)
   character(kind=c_char), intent(in) :: c_file(*)
   integer(c_size_t),      intent(in) :: flen
   character(len=1024) :: fname
   integer :: u, ivar, kind, n
   call c2f_string(c_file, flen, fname)
   call setup_restart()
   call frest%set_restart_vectors(1, fates%nsites, fates%sites, fates%bc_in)
   open(newunit=u, file=trim(fname), access='stream', form='unformatted', status='replace')
   write(u) restart_magic, int(fates%nsites, c_int), int(fates_maxElementsPerSite, c_int), &
            int(frest%num_restart_vars(), c_int)
   do ivar = 1, frest%num_restart_vars()
      associate(rv => frest%rvars(ivar))
         if (associated(rv%r81d)) then
            kind = 0; n = size(rv%r81d)
            write(u) rv%vname, int(kind, c_int), int(n, c_int)
            write(u) rv%r81d
         else
            kind = 1; n = size(rv%int1d)
            write(u) rv%vname, int(kind, c_int), int(n, c_int)
            write(u) rv%int1d
         end if
      end associate
   end do
   close(u)
   write(fates_log_unit,*) 'rhessys_fates: wrote restart ', trim(fname), ' (', frest%num_restart_vars(), ' variables)'
   flush(fates_log_unit)
end subroutine rhessys_veg_restart_write

subroutine rhessys_veg_restart_read(c_file, flen, out_all) bind(C, name=_VEG_RESTART_READ)
   character(kind=c_char), intent(in)  :: c_file(*)
   integer(c_size_t),      intent(in)  :: flen
   type(veg_out_t),        intent(out) :: out_all(*)
   character(len=1024) :: fname
   character(len=8)    :: magic
   character(len=32)   :: vname
   integer(c_int)      :: ns, nslot, nv, kind, n
   integer :: u, ivar, s, ios
   integer, allocatable :: fcolumn(:)

   call c2f_string(c_file, flen, fname)
   call setup_restart()
   open(newunit=u, file=trim(fname), access='stream', form='unformatted', status='old', &
        action='read', iostat=ios)
   if (ios /= 0) then
      write(*,*) 'FATAL: cannot open FATES restart ', trim(fname)
      error stop 'rhessys_veg_restart_read'
   end if
   read(u) magic, ns, nslot, nv
   if (magic /= restart_magic .or. ns /= fates%nsites .or. nslot /= fates_maxElementsPerSite &
       .or. nv /= frest%num_restart_vars()) then
      write(*,*) 'FATAL: FATES restart ', trim(fname), ' does not match this run: sites ', ns, &
           fates%nsites, ' slots ', nslot, fates_maxElementsPerSite, ' vars ', nv, frest%num_restart_vars()
      error stop 'rhessys_veg_restart_read'
   end if
   do ivar = 1, frest%num_restart_vars()
      read(u) vname, kind, n
      associate(rv => frest%rvars(ivar))
         if (trim(vname) /= trim(rv%vname)) then
            write(*,*) 'FATAL: FATES restart variable ', ivar, ' is ', trim(vname), ', expected ', trim(rv%vname)
            error stop 'rhessys_veg_restart_read'
         end if
         if (kind == 0) then
            read(u) rv%r81d(1:n)
         else
            read(u) rv%int1d(1:n)
         end if
      end associate
   end do
   close(u)

   ! CLM restart sequence (clmfates_interfaceMod%restart, flag 'read')
   call frest%create_patchcohort_structure(1, fates%nsites, fates%sites, fates%bc_in, fates%bc_out)
   call frest%get_restart_vectors(1, fates%nsites, fates%sites, fates%bc_in)
   do s = 1, fates%nsites
      fates%bc_in(s)%max_rooting_depth_index_col = fates%bc_in(s)%nlevsoil
      call ed_update_site(fates%sites(s), fates%bc_in(s), fates%bc_out(s), is_restarting=.true.)
   end do
   allocate(fcolumn(fates%nsites))
   fcolumn = [(s, s = 1, fates%nsites)]
   call canopy_summarization(fates%nsites, fates%sites, fates%bc_in)
   call update_hlm_dynamics(fates%nsites, fates%sites, fcolumn, fates%bc_out)
   deallocate(fcolumn)
   call frest%update_3dpatch_radiation(fates%nsites, fates%sites, fates%bc_out)
   do s = 1, fates%nsites
      out_all(s)%gpp_day = 0; out_all(s)%gpp_night = 0; out_all(s)%resp_m = 0
      out_all(s)%laisun = 0; out_all(s)%laisha = 0; out_all(s)%rssun = 0; out_all(s)%rssha = 0
      out_all(s)%gc_canopy = 0; out_all(s)%swabs_canopy = 0; out_all(s)%btran = 0
      call sample_site(s, out_all(s), .false.)
   end do
   write(fates_log_unit,*) 'rhessys_fates: restarted from ', trim(fname), '; site 1: ', &
        out_all(1)%npatches, ' patches, ', out_all(1)%ncohorts, ' cohorts'
   flush(fates_log_unit)
end subroutine rhessys_veg_restart_read

!-------------------------------------------------------------------------------
! Diagnostics: append every cohort of site s to a CSV (one row per cohort).
!-------------------------------------------------------------------------------
subroutine rhessys_veg_dump_cohorts(site, patch_id, ymd, c_file, flen, first) &
                                    bind(C, name=_VEG_DUMP)
   integer(c_int),         intent(in) :: site, patch_id, ymd, first
   character(kind=c_char), intent(in) :: c_file(*)
   integer(c_size_t),      intent(in) :: flen
   character(len=1024) :: fname
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort
   integer :: i, u, ipa

   fname = ''
   do i = 1, min(int(flen), len(fname))
      fname(i:i) = c_file(i)
   end do
   if (first /= 0) then
      open(newunit=u, file=trim(fname), status='replace', action='write')
      write(u,'(a)') 'date,patchID,fates_patch,patch_age_yr,patch_area_frac,pft,dbh_cm,height_m,'// &
           'stems_per_ha,canopy_layer,leafc_kg_per_plant,crown_area_m2'
   else
      open(newunit=u, file=trim(fname), status='old', position='append', action='write')
   end if
   ipa = 0
   cpatch => fates%sites(site)%oldest_patch
   do while (associated(cpatch))
      ipa = ipa + 1
      ccohort => cpatch%tallest
      do while (associated(ccohort))
         write(u,'(i8,",",i0,",",i0,",",f9.3,",",f8.5,",",i0,",",f8.3,",",f8.3,",",es12.5,",",i0,",",es12.5,",",f10.3)') &
              ymd, patch_id, ipa, cpatch%age, cpatch%area * area_inv, ccohort%pft, ccohort%dbh, &
              ccohort%height, ccohort%n * 1.0e4_r8 / cpatch%area * (cpatch%area * area_inv), &
              ccohort%canopy_layer, ccohort%prt%GetState(leaf_organ, carbon12_element), ccohort%c_area
         ccohort => ccohort%shorter
      end do
      cpatch => cpatch%younger
   end do
   close(u)
end subroutine rhessys_veg_dump_cohorts

!-------------------------------------------------------------------------------
subroutine rhessys_veg_clean() bind(C, name=_VEG_CLEAN)
   ! TODO(Phase 5): deallocate site structures (FATES has no fates_clean yet).
   if (allocated(mort_acc)) call mort_write(size(fates%sites))
   if (fates_log_open) then
      write(fates_log_unit,*) 'rhessys_fates: set_bc calls ', n_setbc, ', bc_in check failures ', n_bad_total
      close(fates_log_unit)
   end if
   fates_is_init = .false.
   fates_log_open = .false.
end subroutine rhessys_veg_clean

end module rhessys_fates
