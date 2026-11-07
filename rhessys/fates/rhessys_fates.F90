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
   use FatesConstantsMod,    only : ievergreen, pi_const, tfrz_const => t_water_freeze_k_1atm
   use FatesConstantsMod,    only : itrue, ifalse, default_regeneration, TRS_regeneration, &
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
   use PRTGenericMod,        only : sapw_organ, struct_organ, store_organ, repro_organ, nitrogen_element
   use PRTAllometricCNPMod,  only : stoich_growth_min, stoich_max
   use FatesInterfaceTypesMod,only: hlm_current_year, hlm_current_month
   use FatesInterfaceTypesMod,only: fates_dispersal_cadence_daily
   use FatesConstantsMod,    only : i_term_mort_type_cstarv, i_term_mort_type_canlev, &
                                    i_term_mort_type_numdens
   use EDPftvarcon,          only : EDPftvarcon_inst
   use FatesLitterMod,       only : litter_type, ncwd, ilabile, icellulose, ilignin, ndcmpy
   use FatesRestartInterfaceMod, only : fates_restart_interface_type
   use FatesIODimensionsMod, only : fates_bounds_type
   use FatesInterfaceTypesMod, only : fates_maxElementsPerSite
   use EDParamsMod,          only : ED_val_cwd_fcel, ED_val_cwd_flig
   use EDLoggingMortalityMod, only : host_logging_trigger    ! local FATES patch (fates/patches)
   use EDLoggingMortalityMod, only : host_leaf_export_frac, host_leaf_export   ! local FATES patch
   use EDLoggingMortalityMod, only : host_litter_passthrough                    ! local FATES patch
   use EDTypesMod,           only : area_inv
   use rhessys_fates_mgmt,   only : mgmt_init, mgmt_pre, mgmt_post, mgmt_write_day, mgmt_summary
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
      real(c_double) :: thin_frac(2)
      real(c_double) :: wm_frac, wm_burn, wm_litter_consumed
      real(c_double) :: wm_under_mort, wm_consumption, wm_over_k1, wm_over_k2, wm_under_h, wm_over_h
      real(c_double) :: n_uptake_mine
      real(c_double) :: n_uptake_am_frac
      real(c_double) :: n_avail
      real(c_double) :: n_diurnal
      real(c_double) :: cz_a, cz_b, hhalf
      real(c_double) :: tmax, tmin
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
      real(c_double) :: n_fix
      real(c_double) :: n_mine_demand, n_mine_recv
      real(c_double) :: n_am_demand, n_am_recv
      real(c_double) :: n_sub
      real(c_double) :: sub_gc(VEG_NSUB_MAX), sub_swabs(VEG_NSUB_MAX), sub_t(VEG_NSUB_MAX), sub_vpd(VEG_NSUB_MAX)
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

   ! Plant N demand handed to RHESSys's soil competition (mode 3):
   !   0 = RD uptake capacity only: fine-root C x (Vmax NH4 + Vmax NO3)
   ! A bit mask: bit 1 (value 1) = growth-based (-demog_nneed), bit 2 (value 2) = store-regulated
   ! (-demog_nstore_reg); both can be set.
   !   1 = growth-based (-demog_nneed): per cohort max(capacity, potential growth N need),
   !       need = potential C gain (GPP - maintenance resp, less growth resp) x tissue N:C,
   !       as legacy RHESSys sets plant demand from potential growth. The same weights
   !       share the granted N among cohorts.
   !   2 = store-regulated (-demog_nstore_reg): the claim times the cohort's room in its N store,
   !       max(0, 1 - N store / N store maximum), maximum = FATES's overflow level, target x
   !       (1 + fates_cnp_store_ovrflw_frac) (target = fates_cnp_nitr_store_ratio x full-leaf N), so cohorts
   !       with full N stores stop competing with soil microbes for mineral N. Not the target itself:
   !       FATES refills storage (allocation priority 2) before sapwood and structure, so a plant held at
   !       its store target gets no N for stature growth.
   !   4 = need-capped (-demog_ncap): the claim is at most the cohort's N need, potential growth N
   !       need (as bit 1) plus refilling its N store to target over ncap_refill_days:
   !       min(capacity, need + max(0, store target - store) / ncap_refill_days). Root capacity
   !       still limits sparse stands; a mature stand no longer claims N it cannot use, so
   !       decomposers are not starved of mineral N (litter build-up in mode 3).
   integer, save :: ndemand_mode = 0
   ! -demog_nsplit (mode bit 4, value 8): the granted plant N is shared among cohorts by max(claim, potential
   ! growth N need) instead of the claim alone; the site claim on soil N (out%ncap) is unchanged, so total plant
   ! uptake still comes from RHESSys's plant-microbe competition. Gives gap saplings, whose C gain per fine-root C
   ! is ~10x a canopy tree's, N in proportion to what they can use (2026-11-03).
   logical, save :: nsplit_need = .false.
   ! -demog_nsplit 2 (mode bit 5, value 16): share by min(claim, growth N need) instead - no cohort weighted above
   ! its root capacity; cohorts that cannot use N (C-limited) give up share (2026-11-03, after ns1's shrub takeover)
   logical, save :: nsplit_min = .false.
   integer, save :: stomatal_model_sel = 1   ! 1 = Ball-Berry (default), 2 = Medlyn (-demog_medlyn); set before init
   integer, save :: regen_model_sel = default_regeneration   ! or TRS_regeneration (-demog_trs); set before init
   integer, save :: coage_sel = 0            ! 1 = cohort age tracking (-demog_coage); set before init
   integer, save :: dayl_sel = 0             ! 1 = day-length factor on vcmax/jmax (-demog_dayl); set before init
   ! Seed dispersal between neighbouring sites (-demog_disperse; see disperse_seeds): FATES sends
   ! fates_seed_dispersal_fraction of each site's seed rain (the share crossing the patch boundary) to
   ! site%seed_out; the host splits it among the other sites whose edge gap g = centroid distance - r_i - r_j
   ! (r = sqrt(area/pi)) is <= fates_seed_dispersal_max_dist (m), in proportion to exp(-max(g,0)/a) r_j with
   ! a = fates_seed_dispersal_pdf_scale (m): touching patches take most, larger ones more. Shares sum to 1
   ! (seed falling on non-vegetated patches or across the domain edge is not lost); a site with no
   ! neighbour loses its dispersed seed. Received seed enters as site%seed_in the next day.
   integer, save :: disp_sel = 0
   integer, save :: disp_maxnb = 0
   integer,  allocatable, save :: disp_nnb(:), disp_nb(:,:)          ! (site), (neighbour, site)
   real(r8), allocatable, save :: disp_w(:,:,:)                      ! (pft, neighbour, source site)
   real(r8), allocatable, save :: disp_area(:)                       ! site area (m2)
   real(r8), allocatable, save :: disp_tot(:,:)                      ! (1 out, 2 in, 3 lost; pft) kgC this year
   integer, save :: disp_year = -1
   real(r8), parameter :: ncap_refill_days = 30.0_r8
   real(r8), save :: ncap_mult = 1.0_r8      ! -demog_ncap multiplier on the need
   ! -demog_nmine: C-paid N acquisition. Mining PFTs (bit k-1 of nmine_mask = PFT k) turn their
   ! excess respiration (C that N shortage kept from growth) into a demand for soil organic N at
   ! nmine_cost gC/gN; the host grants it from the soil (rate-limited) and the mined N is shared
   ! among the mining cohorts by their excess C. nmine_cost <= 0 = off.
   real(r8), save :: nmine_cost = -1.0_r8
   integer,  save :: nmine_mask = 0
   ! -demog_nforage: AM foraging. Foraging PFTs (bit mask) turn excess respiration into an extra
   ! claim on soil mineral N (excess C / nforage_cost); the share of the granted uptake that the
   ! extra claim won goes to the foraging cohorts by their excess C. nforage_cost <= 0 = off.
   real(r8), save :: nforage_cost = -1.0_r8
   integer,  save :: nforage_mask = 0
   ! -demog_nfun: FUN-style N acquisition cost (Fisher et al. 2010; CLM5 CNFUNMod active uptake):
   ! cost (gC/gN) = kN / root-zone mineral N (gN/m2) + kC / fine-root C per crown area (gC/m2),
   ! with AM (nfun_k(1:2)) or ECM (nfun_k(3:4), PFTs in nfun_mask) parameters. Each cohort's share of
   ! the plant mineral-N uptake RHESSys granted is weighted by 1 / cost (total unchanged); with
   ! nfun_charge = 1 the cohort pays N x cost from today's NPP (added to resp_m_acc, so it is
   ! booked as autotrophic respiration). nfun_k(1) <= 0 = off.
   real(r8), save :: nfun_k(4) = [-1.0_r8, 0.0_r8, 0.0_r8, 0.0_r8]
   integer,  save :: nfun_charge = 0
   integer,  save :: nfun_mask = 0
   ! today's site terms (g/m2/day): N to AM, N to ECM, N x cost AM, N x cost ECM, C charged
   real(r8), allocatable, save :: nfun_diag(:,:)       ! (5, site)

   ! litter C handed to the host today by FATES source (gC/m2/day), per site:
   ! 1 leaf fines, 2 seeds, 3 above-ground CWD, 4 fine roots, 5 below-ground CWD;
   ! 6-10 the same for N (gN/m2/day, mode 3)
   real(r8), allocatable, save :: lit_src(:,:)

   ! stream-side sites (local FATES patch, host_leaf_export): share of each site's fresh leaf
   ! litter that falls into the channel, set by the host; today's amount taken out of FATES,
   ! (4, site) in g/m2/day: C, N, labile C, labile N
   real(r8), allocatable, save :: leaf_exp_frac(:)
   real(r8), allocatable, save :: leaf_exp(:,:)

   ! -demog_litter_direct: FATES hands all fine litter to the host on the day it falls (local
   ! FATES patch, host_litter_passthrough), and the host files it into its litter pools with
   ! the legacy rules. lit_dc: today's litter to the host by decomposability class (1 labile,
   ! 2 cellulose, 3 lignin), source (1 leaf fines, 2 seeds, 3 above-ground CWD, 4 fine-root
   ! fines, 5 below-ground CWD) and element (1 C, 2 N), per site (g/m2/day)
   logical, save :: litter_direct = .false.
   real(r8), allocatable, save :: lit_dc(:,:,:,:)

   ! per-patch inventories: inventory index (0-based) of each site, cold start only
   integer, allocatable, save :: inv_idx(:)

   ! Annual demography diagnostic (<output prefix>_fates_mort.csv): site-mean
   ! (unweighted over sites) stems, biomass, NPP, C lost per mortality mechanism,
   ! termination and impact mortality, by PFT x canopy layer (1 canopy, 2 understory).
   integer, parameter :: nmd = 21
   integer, parameter :: n_aorg = 6
   integer, parameter :: aorg(n_aorg) = [leaf_organ, fnrt_organ, sapw_organ, struct_organ, store_organ, repro_organ]
   character(len=1024), save :: mort_file = ''
   ! Recruitment diagnostic (<output prefix>_fates_regen.csv, on with env RHESSYS_FATES_REGEN=1): per
   ! site (patch ID) x PFT and year, the seed -> seedling -> sapling pipeline, FATES-patch-area weighted:
   ! seed rain, germination and seedling mortality (kgC/m2/yr), recruits (stems/ha/yr, FATES's
   ! recruitment_rate, read and reset daily), mean seed bank and seedling pool (kgC/m2), and May-Sep means
   ! of the TRS seedling-layer PAR (24-h mean, area-weighted and the brightest FATES patch, MJ/m2/day),
   ! the TRS seedling-depth soil suction (mm) and canopy closure (TRS fields are -1 without -demog_trs).
   character(len=1024), save :: regen_file = ''
   integer, parameter :: nrg = 10
   real(r8), allocatable, save :: regen_acc(:,:,:)     ! (nrg, pft, site)
   integer, allocatable, save :: regen_ngs(:)          ! May-Sep days per site
   integer, save :: regen_year = -1, regen_ndays = 0
   ! Carbon budget diagnostic (<output prefix>_fates_cbudget.csv, with the regen file): per site x PFT x
   ! canopy layer (1 canopy, 2 understory) and year: mean stems (/ha) and biomass, leaf, store, repro C
   ! (kgC/m2); GPP, maintenance respiration, NPP, excess-C efflux and tissue turnover (kgC/m2/yr); N gain
   ! (gN/m2/yr); biomass-weighted shares of days with growth limited by C (cnp_limiter 1) and by N (2);
   ! N store / N store target; May-Sep btran of the cohort's PFT, canopy trim and l2fr (biomass weighted).
   character(len=1024), save :: cb_file = ''
   integer, parameter :: ncb = 21
   real(r8), allocatable, save :: cb_acc(:,:,:,:)      ! (ncb, pft, layer, site)
   ! env RHESSYS_FATES_CB_SMALL (dbh, cm; diagnostic, off when unset): canopy-layer-1 cohorts with dbh below it
   ! are booked as layer 3 in the cbudget file (not in cbmonth, which keeps them in layer 1). Output only.
   real(r8), save :: cb_small_dbh = -1.0_r8
   ! monthly companion (<output prefix>_fates_cbmonth.csv): GPP, maintenance respiration, NPP (kgC/m2/month)
   ! and mean biomass (kgC/m2) per site x PFT x canopy layer x month
   character(len=1024), save :: cbm_file = ''
   real(r8), allocatable, save :: cbm_acc(:,:,:,:,:)   ! (4, month, pft, layer, site)
   integer, allocatable, save :: cbm_nd(:)             ! days per month accumulated
   logical, save :: mort_first = .true.
   integer, save :: mort_year = -1, mort_ndays = 0
   real(r8), allocatable, save :: mort_acc(:,:,:)      ! (nmd, pft, layer)
   ! env RHESSYS_FATES_MORT_PATCHES (comma list of patch IDs): restrict the mort file to those sites (mean over them)
   logical,  allocatable, save :: mort_use(:)
   integer,  save :: mort_nsel = 0
   real(r8), allocatable, save :: term_base(:,:,:)     ! (3 types + imort, pft, layer) cumulative at year start
   ! Annual organ C budget (<output prefix>_fates_alloc.csv), site mean by PFT, organs leaf, fine
   ! root, sapwood, structure, storage, reproduction: production = stock change + turnover + C lost
   ! with dying plants (cohort mortality rates; termination and impact mortality are whole-plant and
   ! reported separately), i.e. the year's growth per organ, comparable with legacy RHESSys allocation.
   character(len=1024), save :: alloc_file = ''
   logical, save :: alloc_first = .true.
   integer, save :: alloc_jul_days = 0
   real(r8), allocatable, save :: alloc_acc(:,:,:)     ! (1 turnover, 2 mortality loss; organ; pft), kgC/m2 summed over sites
   real(r8), allocatable, save :: alloc_stock0(:,:)    ! (organ, pft) stocks at year start, summed over sites
   real(r8), allocatable, save :: alloc_jul(:,:)       ! (1 leaf area, 2 leaf C; pft) July sums over sites and days
   real(r8), allocatable, save :: alloc_nst(:,:)       ! (1 N store, 2 FATES N store target; pft) sums over sites and days (CNP)

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
   character(len=64)   :: envval
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
      alloc_file = logfile(1:i-1)//'_fates_alloc.csv'
      call get_environment_variable('RHESSYS_FATES_REGEN', envval)
      if (len_trim(envval) > 0 .and. trim(envval) /= '0') then
         regen_file = logfile(1:i-1)//'_fates_regen.csv'
         cb_file = logfile(1:i-1)//'_fates_cbudget.csv'
         cbm_file = logfile(1:i-1)//'_fates_cbmonth.csv'
      end if
   else
      mort_file = trim(logfile)//'_mort.csv'
      alloc_file = trim(logfile)//'_alloc.csv'
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
   call set_fates_ctrlparms('seeddisp_cadence',    ival=merge(fates_dispersal_cadence_daily, 0, disp_sel == 1))  ! -demog_disperse
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
   call set_fates_ctrlparms('use_cohort_age_tracking',ival=merge(itrue, ifalse, coage_sel == 1))  ! -demog_coage
   call set_fates_ctrlparms('use_tree_damage',        ival=ifalse)

   ! photosynthesis / physiology model selectors (FATES defaults)
   call set_fates_ctrlparms('photosynth_acclimation',      ival=photosynth_acclim_model_none)
   call set_fates_ctrlparms('use_daylength_factor_switch', ival=merge(itrue, ifalse, dayl_sel == 1))  ! -demog_dayl; bc_in dayl_factor = (dayl/dayl_max)^2
   call set_fates_ctrlparms('stomatal_model',          ival=stomatal_model_sel)  ! Ball-Berry, or Medlyn (-demog_medlyn)
   call set_fates_ctrlparms('stomatal_assim_model',    ival=net_assim_model)
   call set_fates_ctrlparms('maintresp_leaf_model',    ival=lmrmodel_ryan_1991)
   call set_fates_ctrlparms('mort_cstarvation_model',  ival=cstarvation_model_lin)
   call set_fates_ctrlparms('regeneration_model',      ival=regen_model_sel)  ! default, or TRS (-demog_trs)
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
   real(r8), allocatable :: lat_keep(:), lon_keep(:)

   do s = 1, fates%nsites
      out_all(s)%gpp_day = 0; out_all(s)%gpp_night = 0; out_all(s)%resp_m = 0
      out_all(s)%laisun = 0; out_all(s)%laisha = 0; out_all(s)%rssun = 0; out_all(s)%rssha = 0
      out_all(s)%gc_canopy = 0; out_all(s)%swabs_canopy = 0; out_all(s)%btran = 0; out_all(s)%n_sub = 0
   end do

   call set_site_properties(fates%nsites, fates%sites, fates%bc_in)
   if (allocated(inv_idx)) then
      ! per-patch inventories: FATES matches each site to the nearest inventory
      ! by lat/lon; give every site its inventory's placeholder coordinate
      ! (host-written control file, 1-degree grid) for the read only
      allocate(lat_keep(fates%nsites), lon_keep(fates%nsites))
      do s = 1, fates%nsites
         lat_keep(s) = fates%sites(s)%lat; lon_keep(s) = fates%sites(s)%lon
         fates%sites(s)%lat = -80.0_r8 + real(inv_idx(s) / 360, r8)
         fates%sites(s)%lon = real(mod(inv_idx(s), 360), r8)
      end do
   end if
   call init_patches(fates%nsites, fates%sites, fates%bc_in)
   if (allocated(inv_idx)) then
      do s = 1, fates%nsites
         fates%sites(s)%lat = lat_keep(s); fates%sites(s)%lon = lon_keep(s)
      end do
      deallocate(lat_keep, lon_keep, inv_idx)
   end if
   call litter_init_env()                       ! RHESSYS_FATES_LITTER_INIT (soil/litter spin-up carry-over)
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
         ! canopy-area quantities (LAI, SAI, sun/shade LAI, canopy conductance, absorbed
         ! radiation) of a FATES patch are per area covered by crowns; weight them, as
         ! CLM does, by bc_out%canopy_fraction_pa (crown-covered fraction x patch area)
         w = bo%canopy_fraction_pa(ifp)
         out%tlai = out%tlai + w * bo%tlai_pa(ifp)
         out%elai = out%elai + w * bo%elai_pa(ifp)
         out%tsai = out%tsai + w * bo%tsai_pa(ifp)
         out%htop = out%htop + cpatch%area * area_inv * bo%htop_pa(ifp)
         if (with_fluxes) then
            rb = bi%rb_pa(ifp)
            out%laisun = out%laisun + w * bo%laisun_pa(ifp)
            out%laisha = out%laisha + w * bo%laisha_pa(ifp)
            out%rssun  = out%rssun  + cpatch%area * area_inv * bo%rssun_pa(ifp)
            out%rssha  = out%rssha  + cpatch%area * area_inv * bo%rssha_pa(ifp)
            if (bo%laisun_pa(ifp) > 0.0_r8) &
               out%gc_canopy = out%gc_canopy + w * bo%laisun_pa(ifp) / (rb + bo%rssun_pa(ifp))
            if (bo%laisha_pa(ifp) > 0.0_r8) &
               out%gc_canopy = out%gc_canopy + w * bo%laisha_pa(ifp) / (rb + bo%rssha_pa(ifp))
            do ib = 1, VEG_NSWB
               out%swabs_canopy = out%swabs_canopy + w * &
                    (bi%solad_parb(ifp,ib) * bo%fabd_parb(ifp,ib) + &
                     bi%solai_parb(ifp,ib) * bo%fabi_parb(ifp,ib))
            end do
            out%btran = out%btran + cpatch%area * area_inv * bo%btran_pa(ifp)
         end if
         cpatch => cpatch%younger
      end do
   end associate
   out%lai_evg_frac = 0.0_r8
   if (wl > 0.0_r8) out%lai_evg_frac = wl_evg / wl
end subroutine sample_site

!-------------------------------------------------------------------------------
! Diagnostic only (no state change): if env RHESSYS_RAD_DIAG holds this site's
! RHESSys patch ID, print to stderr the daylight shortwave FATES's canopy
! radiation absorbs in the canopy and transmits to the ground (direct->direct,
! direct->diffuse, diffuse->diffuse), crown-covered share weighted as the host
! weights swabs_canopy, uncovered share transmitted in full.
!-------------------------------------------------------------------------------
subroutine rad_ground_diag(s)
   integer, intent(in) :: s
   integer, save :: dpid = -2
   character(len=32) :: ev
   integer :: stat, ifp, ib
   real(r8) :: pf, cf, sd, si, tr, ab, lai, sai
   type(fates_patch_type), pointer :: cpatch
   if (dpid == -2) then
      dpid = -1
      call get_environment_variable('RHESSYS_RAD_DIAG', ev, status=stat)
      if (stat == 0 .and. len_trim(ev) > 0) read(ev, *, iostat=stat) dpid
      if (stat /= 0) dpid = -1
   end if
   if (dpid < 0 .or. fates%sites(s)%h_gid /= dpid) return
   sd = 0; si = 0; tr = 0; ab = 0; lai = 0; sai = 0
   associate(bi => fates%bc_in(s), bo => fates%bc_out(s))
      cpatch => fates%sites(s)%oldest_patch
      do while (associated(cpatch))
         ifp = cpatch%patchno
         pf = cpatch%area * area_inv
         cf = bo%canopy_fraction_pa(ifp)
         lai = lai + cf * bo%elai_pa(ifp)
         sai = sai + cf * bo%esai_pa(ifp)
         do ib = 1, VEG_NSWB
            sd = sd + pf * bi%solad_parb(ifp,ib)
            si = si + pf * bi%solai_parb(ifp,ib)
            tr = tr + cf * (bi%solad_parb(ifp,ib) * (bo%ftdd_parb(ifp,ib) + bo%ftid_parb(ifp,ib)) + &
                            bi%solai_parb(ifp,ib) * bo%ftii_parb(ifp,ib)) + &
                      (pf - cf) * (bi%solad_parb(ifp,ib) + bi%solai_parb(ifp,ib))
            ab = ab + cf * (bi%solad_parb(ifp,ib) * bo%fabd_parb(ifp,ib) + &
                            bi%solai_parb(ifp,ib) * bo%fabi_parb(ifp,ib))
         end do
         cpatch => cpatch%younger
      end do
   end associate
   write(0,'(a,i0,8(1x,f9.4))') 'RADGFATES ', dpid, sd, si, ab, tr, &
        ab / max(sd + si, 1.0e-9_r8), tr / max(sd + si, 1.0e-9_r8), lai, sai
end subroutine rad_ground_diag

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
! -demog_diurnal: the daylight period as n equal flux sub-steps instead of one
! daylight-mean step (photosynthesis, conductance and respiration only; the
! daily demography is unchanged). Step k is centred at hour angle
! h_k = -H + (k - 1/2) 2H/n (H = sunset hour angle):
!   cos(zenith)  cz_k = a + b cos(h_k)                         (bc cz_a, cz_b)
!   shortwave    daylight-mean direct and diffuse x cz_k / mean(cz)  (daylight means kept)
!   temperature  Parton & Logan (1981) day curve tmin + (tmax - tmin) sin(pi m / (D + 2 x 1.86 h)),
!                m = time after sunrise, D = day length, shifted so its mean equals the
!                daylight-mean temperature the single step uses
!   vapour       actual vapour pressure constant through the day; saturation pressure
!                scaled from the host's daylight value by the Magnus ratio
! GPP and maintenance respiration are summed over the sub-steps; out%gc_canopy and
! out%swabs_canopy are time means (diagnostics), and the per-step conductance,
! absorbed shortwave, temperature (C) and VPD (Pa) go to the host transpiration.
!-------------------------------------------------------------------------------
subroutine daylight_diurnal(s, bc, dt_day, out)
   integer,         intent(in)    :: s
   type(veg_bc_t),  intent(in)    :: bc
   real(r8),        intent(in)    :: dt_day
   type(veg_out_t), intent(inout) :: out
   real(r8), parameter :: lag = 1.86_r8 * 3600.0_r8
   integer  :: n, k, ib
   real(r8) :: hsun, ha, cz(VEG_NSUB_MAX), tk(VEG_NSUB_MAX), czm, tm, dt, gpp, resp_m, es_k, ea_k
   real(r8) :: gpp_sum, rm_sum, gc_sum, sw_sum
   n = min(max(nint(bc%n_diurnal), 1), VEG_NSUB_MAX)
   hsun = bc%hhalf
   dt = dt_day / real(n, r8)
   do k = 1, n
      ha = -hsun + (real(k, r8) - 0.5_r8) * 2.0_r8 * hsun / real(n, r8)
      cz(k) = max(bc%cz_a + bc%cz_b * cos(ha), 0.0_r8)
      tk(k) = bc%tmin + (bc%tmax - bc%tmin) * sin(pi_const * (real(k, r8) - 0.5_r8) * dt / (dt_day + 2.0_r8 * lag))
   end do
   czm = sum(cz(1:n)) / real(n, r8)
   tm  = sum(tk(1:n)) / real(n, r8)
   tk(1:n) = tk(1:n) + (bc%t_veg - tm)
   gpp_sum = 0.0_r8; rm_sum = 0.0_r8; gc_sum = 0.0_r8; sw_sum = 0.0_r8
   do k = 1, n
      es_k = bc%esat_tv * esat_magnus(tk(k)) / esat_magnus(bc%t_veg)
      ea_k = min(bc%eair, es_k)
      associate(bi => fates%bc_in(s))
         bi%coszen = min(max(cz(k), 1.0e-3_r8), 1.0_r8)
         do ib = 1, VEG_NSWB
            bi%solad_parb(:,ib) = bc%solad(ib) * cz(k) / max(czm, 1.0e-9_r8)
            bi%solai_parb(:,ib) = bc%solai(ib) * cz(k) / max(czm, 1.0e-9_r8)
         end do
         bi%t_veg_pa(:)   = tk(k)
         bi%tgcm_pa(:)    = tk(k)
         bi%esat_tv_pa(:) = es_k
         bi%eair_pa(:)    = ea_k
      end associate
      call flux_substep(s, dt)
      call sample_site(s, out, .true.)
      call site_carbon(s, gpp, resp_m)
      gpp_sum = gpp_sum + gpp
      rm_sum  = rm_sum + resp_m
      out%sub_gc(k)    = out%gc_canopy
      out%sub_swabs(k) = out%swabs_canopy
      out%sub_t(k)     = tk(k) - tfrz_const
      out%sub_vpd(k)   = max(es_k - ea_k, 0.0_r8)
      gc_sum = gc_sum + out%gc_canopy
      sw_sum = sw_sum + out%swabs_canopy
   end do
   out%n_sub = real(n, r8)
   out%gc_canopy    = gc_sum / real(n, r8)
   out%swabs_canopy = sw_sum / real(n, r8)
   out%gpp_day = gpp_sum
   out%resp_m  = out%resp_m + rm_sum
   ! restore the daylight means in bc_in (the night step and the running means reset what they use)
   associate(bi => fates%bc_in(s))
      bi%coszen = bc%coszen
      do ib = 1, VEG_NSWB
         bi%solad_parb(:,ib) = bc%solad(ib)
         bi%solai_parb(:,ib) = bc%solai(ib)
      end do
      bi%t_veg_pa(:) = bc%t_veg; bi%tgcm_pa(:) = bc%tgcm
      bi%esat_tv_pa(:) = bc%esat_tv; bi%eair_pa(:) = bc%eair
   end associate
contains
   real(r8) function esat_magnus(tkel)
      real(r8), intent(in) :: tkel
      esat_magnus = 611.2_r8 * exp(17.67_r8 * (tkel - 273.15_r8) / (tkel - 273.15_r8 + 243.5_r8))
   end function esat_magnus
end subroutine daylight_diurnal

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
   real(r8) :: amf, exc                        ! -demog_nforage share; cohort excess C
   real(r8) :: dt_day, dt_night, gpp, resp_m
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort

   s = site
   if (s < 1 .or. s > fates%nsites) error stop 'rhessys_veg_do: bad site'
   dt_day   = min(max(bc%dayl, 0.0_r8), 86400.0_r8)
   dt_night = 86400.0_r8 - dt_day
   out%gpp_day = 0; out%gpp_night = 0; out%resp_m = 0
   call sample_site(s, out, .true.)               ! zero/structure before the steps

   ! --- daylight sub-step(s) (bc_in holds the daylight forcing) ---
   out%n_sub = 0.0_r8
   out%sub_gc(:) = 0.0_r8; out%sub_swabs(:) = 0.0_r8; out%sub_t(:) = 0.0_r8; out%sub_vpd(:) = 0.0_r8
   if (dt_day > 0.0_r8 .and. nint(bc%n_diurnal) > 1 .and. bc%hhalf > 0.0_r8) then
      call daylight_diurnal(s, bc, dt_day, out)   ! -demog_diurnal
   else if (dt_day > 0.0_r8) then
      call flux_substep(s, dt_day)
      call sample_site(s, out, .true.)            ! daylight radiation/conductance/btran
      call rad_ground_diag(s)                     ! env RHESSYS_RAD_DIAG=<patch ID>: diagnostic only
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
   if (regen_model_sel == TRS_regeneration) then
      ! TRS seedling-layer light (UpdateFatesRMeansTStep -> SeedlingParPatch) uses the bc_in PAR
      ! and the normalized canopy PAR profiles, both zeroed by the night sub-step: redo the
      ! normalized radiation for the daylight sun and pass the 24-h mean PAR (W/m2)
      associate(bi => fates%bc_in(s))
         bi%coszen = bc%coszen
         do ib = 1, VEG_NSWB
            bi%solad_parb(:,ib) = bc%solad(ib)
            bi%solai_parb(:,ib) = bc%solai(ib)
         end do
         if (bi%coszen > 0.0_r8) call FatesNormalizedCanopyRadiation(fates%sites(s:s), fates%bc_in(s:s), fates%bc_out(s:s))
         do ib = 1, VEG_NSWB
            bi%solad_parb(:,ib) = bc%solad(ib) * dt_day / 86400.0_r8
            bi%solai_parb(:,ib) = bc%solai(ib) * dt_day / 86400.0_r8
         end do
      end associate
   end if
   call UpdateFatesRMeansTStep(fates%sites(s:s), fates%bc_in(s:s), fates%bc_out(s:s))

   ! --- daily demography (Phase 3; CLM dynamics_driv): phenology, growth,
   !     allocation, mortality, recruitment, disturbance, patch/cohort fusion.
   !     Carbon-only PARTEH: nutrients supplied (as CLM with carbon-only). ---
   if (cnp_coupled) then
      ! RHESSys-granted uptake -> cohorts (competitor = cohort), share by capacity
      amf = 0.0_r8
      if (nforage_cost > 0.0_r8) amf = min(1.0_r8, max(bc%n_uptake_am_frac, 0.0_r8))
      call distribute_n_uptake(s, bc%n_uptake_nh4 * (1.0_r8 - amf), bc%n_uptake_no3 * (1.0_r8 - amf), &
           out%n_uptake_recv, bc%n_avail)
      out%n_am_recv = 0.0_r8
      if (amf > 0.0_r8) call distribute_n_excess(s, bc%n_uptake_nh4 * amf, bc%n_uptake_no3 * amf, &
           nforage_mask, out%n_am_recv)
      out%n_mine_recv = 0.0_r8
      if (nmine_cost > 0.0_r8 .and. bc%n_uptake_mine > 0.0_r8) &
         call distribute_n_excess(s, bc%n_uptake_mine, 0.0_r8, nmine_mask, out%n_mine_recv)
      out%n_uptake_recv = out%n_uptake_recv + out%n_am_recv + out%n_mine_recv
      call UnPackNutrientAquisitionBCs(fates%sites(s:s), fates%bc_in(s:s), .false., .true.)
   else
      out%n_uptake_recv = 0.0_r8
      call UnPackNutrientAquisitionBCs(fates%sites(s:s), fates%bc_in(s:s), .true., .true.)
   end if

   ! Phase 6: host-imposed fire / harvest for today (always resets yesterday's)
   call apply_disturbance(s, bc)
   if (cnp_coupled) call mgmt_pre(fates%sites(s), s, nint(bc%dist_kind), bc%thin_frac)
   ! stream-side site: FATES diverts this share of today's fresh leaf litter to the host
   host_leaf_export(:) = 0.0_r8
   host_leaf_export_frac = 0.0_r8
   if (allocated(leaf_exp_frac)) host_leaf_export_frac = leaf_exp_frac(s)
   host_litter_passthrough = litter_direct
   ! today's symbiotic N fixation (accumulated over the flux sub-steps; FATES zeroes it in its daily step)
   out%n_fix = 0.0_r8
   if (cnp_coupled) then
      cpatch => fates%sites(s)%oldest_patch
      do while (associated(cpatch))
         ccohort => cpatch%tallest
         do while (associated(ccohort))
            out%n_fix = out%n_fix + ccohort%sym_nfix_daily * ccohort%n
            ccohort => ccohort%shorter
         end do
         cpatch => cpatch%younger
      end do
      out%n_fix = out%n_fix * area_inv * 1000.0_r8               ! kgN/site/day -> gN/m2/day
   end if
   call ed_ecosystem_dynamics(fates%sites(s), fates%bc_in(s), fates%bc_out(s))
   if (allocated(regen_acc)) call regen_accumulate(s)
   if (allocated(leaf_exp)) leaf_exp(:, s) = host_leaf_export(:) * 1000.0_r8   ! kg -> g /m2/day
   host_leaf_export_frac = 0.0_r8
   host_litter_passthrough = .false.
   if (cnp_coupled) call mgmt_post(fates%sites(s), s, fates%bc_in(s))
   call ed_update_site(fates%sites(s), fates%bc_in(s), fates%bc_out(s), is_restarting=.false.)
   call UpdateFatesAvgSnowDepth(fates%sites(s:s), fates%bc_in(s:s))
   call canopy_summarization(1, fates%sites(s:s), fates%bc_in(s:s))
   call update_hlm_dynamics(1, fates%sites(s:s), [s], fates%bc_out(s:s))
   call sample_site(s, out, .false.)              ! end-of-day canopy structure + stocks
   call nutrient_out(s, out)                      ! capacity, litter to host, total stocks
   ! -demog_nmine / -demog_nforage: today's excess C of the mining / foraging PFTs -> N they
   ! bid for tomorrow (kgC/site/day -> gC/m2/day -> gN/m2/day)
   out%n_mine_demand = 0.0_r8
   out%n_am_demand = 0.0_r8
   if (cnp_coupled .and. (nmine_cost > 0.0_r8 .or. nforage_cost > 0.0_r8)) then
      cpatch => fates%sites(s)%oldest_patch
      do while (associated(cpatch))
         ccohort => cpatch%tallest
         do while (associated(ccohort))
            exc = max(nm_fin(ccohort%resp_excess_hold), 0.0_r8) * ccohort%n
            if (nmine_cost > 0.0_r8 .and. btest(nmine_mask, ccohort%pft - 1)) &
               out%n_mine_demand = out%n_mine_demand + exc
            if (nforage_cost > 0.0_r8 .and. btest(nforage_mask, ccohort%pft - 1)) &
               out%n_am_demand = out%n_am_demand + exc
            ccohort => ccohort%shorter
         end do
         cpatch => cpatch%younger
      end do
      if (nmine_cost > 0.0_r8) out%n_mine_demand = out%n_mine_demand * area_inv * 1000.0_r8 / nmine_cost
      if (nforage_cost > 0.0_r8) out%n_am_demand = out%n_am_demand * area_inv * 1000.0_r8 / nforage_cost
   end if
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
   if (len_trim(regen_file) > 0) call regen_write(int(nsites))
   if (nfun_k(1) > 0.0_r8) then
      if (.not. allocated(nfun_diag)) allocate(nfun_diag(5, nsites))
      nfun_diag(:,:) = 0.0_r8
   end if
   !$omp parallel do schedule(dynamic, 4) default(shared) private(s)
   do s = 1, nsites
      call rhessys_veg_do(s, bcs(s), outs(s))
   end do
   !$omp end parallel do
   call mort_accumulate(int(nsites))
   if (disp_sel == 1) call disperse_seeds(int(nsites))
   if (cnp_coupled) call mgmt_write_day()
   if (allocated(regen_acc)) regen_ndays = regen_ndays + 1
end subroutine rhessys_veg_do_all

!-------------------------------------------------------------------------------
! Recruitment diagnostic (see regen_file): today's terms for site s, after the
! daily demography. Called inside the OpenMP site loop; touches only site s.
!-------------------------------------------------------------------------------
subroutine regen_accumulate(s)
   integer, intent(in) :: s
   type(fates_patch_type), pointer :: cpatch
   type(litter_type),      pointer :: litt
   integer  :: ft
   real(r8) :: w, par, parmax
   logical  :: gs
   gs = (hlm_current_month >= 5 .and. hlm_current_month <= 9)
   parmax = 0.0_r8
   cpatch => fates%sites(s)%oldest_patch
   do while (associated(cpatch))
      litt => cpatch%litter(1)                     ! carbon
      w = cpatch%area * area_inv
      par = -1.0_r8
      if (regen_model_sel == TRS_regeneration) then
         par = cpatch%seedling_layer_par24%GetMean() * 86400.0_r8 * 1.0e-6_r8   ! W/m2 -> MJ/m2/day
         parmax = max(parmax, par)
      end if
      do ft = 1, numpft
         associate(a => regen_acc(:, ft, s))
         a(1) = a(1) + w * (litt%seed_in_local(ft) + litt%seed_in_extern(ft))
         a(2) = a(2) + w * litt%seed_germ_in(ft)
         a(3) = a(3) + w * litt%seed_germ_decay(ft)
         a(5) = a(5) + w * litt%seed(ft)
         a(6) = a(6) + w * litt%seed_germ(ft)
         if (gs) then
            a(7)  = a(7)  + w * par
            a(10) = a(10) + w * min(1.0_r8, cpatch%total_canopy_area / max(cpatch%area, 1.0e-9_r8))
            if (regen_model_sel == TRS_regeneration) then
               a(9) = a(9) + w * cpatch%sdlng_emerg_smp(ft)%p%GetMean()
            else
               a(9) = a(9) - w
            end if
         end if
         end associate
      end do
      cpatch => cpatch%younger
   end do
   call cb_accumulate(s, gs)
   if (s == 1) cbm_nd(hlm_current_month) = cbm_nd(hlm_current_month) + 1
   do ft = 1, numpft
      regen_acc(4, ft, s) = regen_acc(4, ft, s) + fates%sites(s)%recruitment_rate(ft)   ! individuals per site (1 ha)
      if (gs) regen_acc(8, ft, s) = regen_acc(8, ft, s) + merge(parmax, -1.0_r8, regen_model_sel == TRS_regeneration)
   end do
   fates%sites(s)%recruitment_rate(:) = 0.0_r8
   if (gs) regen_ngs(s) = regen_ngs(s) + 1
end subroutine regen_accumulate

!-------------------------------------------------------------------------------
! Carbon budget diagnostic (see cb_file): today's cohort terms for site s.
!-------------------------------------------------------------------------------
subroutine cb_accumulate(s, gs)
   integer, intent(in) :: s
   logical, intent(in) :: gs
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort
   integer  :: ft, il, io, ilc
   real(r8) :: nm2, ctot, tov
   cpatch => fates%sites(s)%oldest_patch
   do while (associated(cpatch))
      ccohort => cpatch%tallest
      do while (associated(ccohort))
         ft = ccohort%pft; il = min(max(ccohort%canopy_layer, 1), 2)
         nm2 = ccohort%n * area_inv
         ilc = il
         if (il == 1 .and. ccohort%dbh < cb_small_dbh) ilc = 3
         ctot = 0.0_r8; tov = 0.0_r8
         do io = 1, n_aorg
            ctot = ctot + ccohort%prt%GetState(aorg(io), carbon12_element)
            tov = tov + cbf(ccohort%prt%GetTurnover(aorg(io), carbon12_element))
         end do
         associate(a => cb_acc(:, ft, ilc, s))
         a(1) = a(1) + nm2 * 1.0e4_r8
         a(2) = a(2) + nm2 * ctot
         a(3) = a(3) + nm2 * ccohort%prt%GetState(leaf_organ, carbon12_element)
         a(4) = a(4) + nm2 * ccohort%prt%GetState(store_organ, carbon12_element)
         a(5) = a(5) + nm2 * ccohort%prt%GetState(repro_organ, carbon12_element)
         associate(b => cbm_acc(:, hlm_current_month, ft, il, s))
         if (.not. ccohort%isnew) then
            b(1) = b(1) + nm2 * cbf(ccohort%gpp_acc_hold) / 365.0_r8
            b(2) = b(2) + nm2 * cbf(ccohort%resp_m_acc_hold) / 365.0_r8
            b(3) = b(3) + nm2 * cbf(ccohort%npp_acc_hold) / 365.0_r8
         end if
         b(4) = b(4) + nm2 * ctot
         end associate
         if (.not. ccohort%isnew) then
            a(6) = a(6) + nm2 * cbf(ccohort%gpp_acc_hold) / 365.0_r8
            a(7) = a(7) + nm2 * cbf(ccohort%resp_m_acc_hold) / 365.0_r8
            a(8) = a(8) + nm2 * cbf(ccohort%npp_acc_hold) / 365.0_r8
            a(21) = a(21) + nm2 * cbf(ccohort%resp_excess_hold)   ! kgC/indiv/day; already removed from npp_acc_hold
         end if
         a(10) = a(10) + nm2 * tov
         if (cnp_coupled) then
            a(9)  = a(9)  + nm2 * cbf(ccohort%daily_c_efflux)
            a(11) = a(11) + nm2 * cbf(ccohort%daily_n_gain) * 1000.0_r8
            if (ccohort%cnp_limiter == 1) a(12) = a(12) + nm2 * ctot
            if (ccohort%cnp_limiter == 2) a(13) = a(13) + nm2 * ctot
            a(14) = a(14) + nm2 * ccohort%prt%GetState(store_organ, nitrogen_element)
            a(15) = a(15) + nm2 * ccohort%prt%GetNutrientTarget(nitrogen_element, store_organ, stoich_growth_min)
         end if
         a(16) = a(16) + nm2 * ctot            ! biomass weight for the shares and the May-Sep means
         if (gs) then
            a(17) = a(17) + nm2 * ctot * cpatch%btran_ft(ft)
            a(18) = a(18) + nm2 * ctot * cbf(ccohort%canopy_trim)
            a(19) = a(19) + nm2 * ctot * cbf(ccohort%l2fr)
            a(20) = a(20) + nm2 * ctot           ! May-Sep biomass weight
         end if
         end associate
         ccohort => ccohort%shorter
      end do
      cpatch => cpatch%younger
   end do
contains
   real(r8) function cbf(x)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(r8), intent(in) :: x
      cbf = 0.0_r8
      if (ieee_is_finite(x)) then
         if (abs(x) < 1.0e30_r8) cbf = x
      end if
   end function cbf
end subroutine cb_accumulate

!-------------------------------------------------------------------------------
! Recruitment diagnostic: allocate on the first call; on the first day of a new
! year write last year's rows (one per site x PFT) and reset. Serial.
!-------------------------------------------------------------------------------
subroutine regen_write(nsites)
   integer, intent(in) :: nsites
   integer :: s, ft, u, il, im, cbs_stat
   real(r8) :: g
   character(len=64) :: cbs_ev
   if (.not. allocated(regen_acc)) then
      allocate(regen_acc(nrg, numpft, nsites), regen_ngs(nsites), cb_acc(ncb, numpft, 3, nsites))
      allocate(cbm_acc(4, 12, numpft, 2, nsites), cbm_nd(12))
      regen_acc = 0.0_r8; regen_ngs = 0; regen_ndays = 0; cb_acc = 0.0_r8; cbm_acc = 0.0_r8; cbm_nd = 0
      call get_environment_variable('RHESSYS_FATES_CB_SMALL', cbs_ev, status=cbs_stat)
      if (cbs_stat == 0 .and. len_trim(cbs_ev) > 0) then
         read(cbs_ev, *, iostat=cbs_stat) cb_small_dbh
         if (cbs_stat /= 0) cb_small_dbh = -1.0_r8
         if (cb_small_dbh > 0.0_r8) write(*,'(a,f6.2,a)') 'FATES host: cbudget layer 3 = canopy layer 1 with dbh < ', &
              cb_small_dbh, ' cm (RHESSYS_FATES_CB_SMALL)'
      end if
      open(newunit=u, file=trim(cbm_file), status='replace', action='write')
      write(u,'(a)') 'year,patchID,pft,layer,month,gpp_kgC_m2,resp_m_kgC_m2,npp_kgC_m2,biomass_kgC_m2'
      close(u)
      open(newunit=u, file=trim(cb_file), status='replace', action='write')
      write(u,'(a)') 'year,patchID,pft,layer,stems_per_ha,biomass_kgC_m2,leafc_kgC_m2,storec_kgC_m2,reproc_kgC_m2,'// &
           'gpp_kgC_m2_yr,resp_m_kgC_m2_yr,npp_kgC_m2_yr,c_efflux_kgC_m2_yr,turnover_kgC_m2_yr,n_gain_gN_m2_yr,'// &
           'frac_C_lim,frac_N_lim,nstore_over_target,gs_btran,gs_canopy_trim,gs_l2fr,resp_excess_kgC_m2_yr'
      close(u)
      regen_year = hlm_current_year
      open(newunit=u, file=trim(regen_file), status='replace', action='write')
      write(u,'(a)') 'year,patchID,pft,seed_rain_kgC_m2_yr,germination_kgC_m2_yr,seedling_mort_kgC_m2_yr,'// &
           'recruits_per_ha_yr,seed_bank_kgC_m2,seedling_pool_kgC_m2,gs_floor_par_MJ_m2_d,'// &
           'gs_floor_par_max_MJ_m2_d,gs_seedling_smp_mm,gs_canopy_closure'
      close(u)
      return
   end if
   if (hlm_current_year == regen_year .or. regen_ndays == 0) return
   open(newunit=u, file=trim(regen_file), status='old', position='append', action='write')
   do s = 1, nsites
      g = real(max(regen_ngs(s), 1), r8)
      do ft = 1, numpft
         associate(a => regen_acc(:, ft, s))
         write(u,'(i0,",",i0,",",i0,4(",",es12.5),2(",",es12.5),4(",",es12.5))') regen_year, fates%sites(s)%h_gid, ft, &
              a(1), a(2), a(3), a(4), a(5) / regen_ndays, a(6) / regen_ndays, a(7) / g, a(8) / g, a(9) / g, a(10) / g
         end associate
      end do
   end do
   close(u)
   open(newunit=u, file=trim(cb_file), status='old', position='append', action='write')
   do s = 1, nsites
      do ft = 1, numpft
         do il = 1, 3
            associate(a => cb_acc(:, ft, il, s))
            if (a(16) <= 0.0_r8) cycle
            write(u,'(i0,",",i0,",",i0,",",i0,18(",",es12.5))') regen_year, fates%sites(s)%h_gid, ft, il, &
                 a(1) / regen_ndays, a(2) / regen_ndays, a(3) / regen_ndays, a(4) / regen_ndays, a(5) / regen_ndays, &
                 a(6), a(7), a(8), a(9), a(10), a(11), a(12) / a(16), a(13) / a(16), a(14) / max(a(15), 1.0e-30_r8), &
                 a(17) / max(a(20), 1.0e-30_r8), a(18) / max(a(20), 1.0e-30_r8), a(19) / max(a(20), 1.0e-30_r8), a(21)
            end associate
         end do
      end do
   end do
   close(u)
   open(newunit=u, file=trim(cbm_file), status='old', position='append', action='write')
   do s = 1, nsites
      do ft = 1, numpft
         do il = 1, 2
            if (sum(cbm_acc(4, :, ft, il, s)) <= 0.0_r8) cycle
            do im = 1, 12
               write(u,'(i0,",",i0,",",i0,",",i0,",",i0,4(",",es12.5))') regen_year, fates%sites(s)%h_gid, ft, il, im, &
                    cbm_acc(1:3, im, ft, il, s), cbm_acc(4, im, ft, il, s) / max(cbm_nd(im), 1)
            end do
         end do
      end do
   end do
   close(u)
   regen_acc = 0.0_r8; regen_ngs = 0; regen_ndays = 0; cb_acc = 0.0_r8; cbm_acc = 0.0_r8; cbm_nd = 0
   regen_year = hlm_current_year
end subroutine regen_write

!-------------------------------------------------------------------------------
! Forest management (rhessys_fates_mgmt.F90): read the driver file (may be
! empty: TEC thinning only), RHESSys patch / zone / hillslope ID per site.
!-------------------------------------------------------------------------------
subroutine rhessys_veg_mgmt_init(c_file, flen, nsites, pid, zid, hid) bind(C, name='rhessys_veg_mgmt_init')
   character(kind=c_char), intent(in) :: c_file(*)
   integer(c_size_t),      intent(in) :: flen
   integer(c_int),         intent(in) :: nsites
   integer(c_int),         intent(in) :: pid(nsites), zid(nsites), hid(nsites)
   character(len=1024) :: fname
   integer :: i
   fname = ''
   do i = 1, min(int(flen), len(fname))
      fname(i:i) = c_file(i)
   end do
   i = index(mort_file, '_fates_mort.csv', back=.true.)
   if (i > 0) then
      call mgmt_init(trim(fname), int(nsites), pid, zid, hid, fates_log_unit, &
           mort_file(1:i-1)//'_fates_mgmt.csv')
   else
      call mgmt_init(trim(fname), int(nsites), pid, zid, hid, fates_log_unit, 'rhessys_fates_mgmt.csv')
   end if
end subroutine rhessys_veg_mgmt_init

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
      if (allocated(mort_use)) then
         if (.not. mort_use(s)) cycle
      end if
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

subroutine alloc_stocks(nsites, stock)
   ! organ C stocks by PFT, kgC per m2 of site, summed over sites
   integer,  intent(in)  :: nsites
   real(r8), intent(out) :: stock(n_aorg, numpft)
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort
   integer :: s, io
   stock = 0.0_r8
   do s = 1, nsites
      cpatch => fates%sites(s)%oldest_patch
      do while (associated(cpatch))
         ccohort => cpatch%tallest
         do while (associated(ccohort))
            do io = 1, n_aorg
               stock(io, ccohort%pft) = stock(io, ccohort%pft) + &
                    ccohort%n * area_inv * ccohort%prt%GetState(aorg(io), carbon12_element)
            end do
            ccohort => ccohort%shorter
         end do
         cpatch => cpatch%younger
      end do
   end do
end subroutine alloc_stocks

subroutine alloc_write(nsites)
   ! one row per PFT: site-mean organ budget over the accumulation period (kgC/m2, not annualised),
   ! NPP and C efflux over the period, whole-plant termination + impact mortality C, July means of
   ! leaf area and leaf C (legacy proj_sla = lai_jul / leafc_jul), and (CNP) the period-mean N store and
   ! FATES's N store target (fates_cnp_nitr_store_ratio x full-leaf N), gN per m2 of site
   integer, intent(in) :: nsites
   real(r8) :: stock1(n_aorg, numpft), cum(4, numpft, 2), r, npp, eff, wpm, ds, p
   integer :: u, ft, io
   character(len=6), parameter :: oname(n_aorg) = ['leaf  ', 'fnrt  ', 'sapw  ', 'struct', 'store ', 'repro ']
   character(len=2048) :: hdr
   if (len_trim(alloc_file) == 0 .or. mort_ndays == 0) return
   call alloc_stocks(nsites, stock1)
   call mort_term_cum(nsites, cum)
   r = 1.0_r8 / real(nsites, r8)
   if (alloc_first) then
      hdr = 'year,ndays,pft'
      do io = 1, n_aorg
         hdr = trim(hdr)//','//trim(oname(io))//'_stock,'//trim(oname(io))//'_dstock,'//trim(oname(io))// &
               '_turnover,'//trim(oname(io))//'_mortloss,'//trim(oname(io))//'_production'
      end do
      hdr = trim(hdr)//',npp,c_efflux,term_imort_c,lai_jul,leafc_jul,nstore_gN,nstore_target_gN'
      open(newunit=u, file=trim(alloc_file), status='replace', action='write')
      write(u,'(a)') trim(hdr)
      alloc_first = .false.
   else
      open(newunit=u, file=trim(alloc_file), status='old', position='append', action='write')
   end if
   do ft = 1, numpft
      npp = (mort_acc(3, ft, 1) + mort_acc(3, ft, 2)) * r
      eff = (mort_acc(17, ft, 1) + mort_acc(17, ft, 2)) * r
      wpm = sum(cum(1:4, ft, 1:2) - term_base(1:4, ft, 1:2)) * r
      write(u,'(i0,",",i0,",",i0)', advance='no') mort_year, mort_ndays, ft
      do io = 1, n_aorg
         ds = (stock1(io, ft) - alloc_stock0(io, ft)) * r
         p  = ds + (alloc_acc(1, io, ft) + alloc_acc(2, io, ft)) * r
         write(u,'(5(",",es12.5))', advance='no') stock1(io, ft) * r, ds, alloc_acc(1, io, ft) * r, &
              alloc_acc(2, io, ft) * r, p
      end do
      write(u,'(7(",",es12.5))') npp, eff, wpm, alloc_jul(1, ft) * r / max(alloc_jul_days, 1), &
           alloc_jul(2, ft) * r / max(alloc_jul_days, 1), 1.0e3_r8 * alloc_nst(1, ft) * r / mort_ndays, &
           1.0e3_r8 * alloc_nst(2, ft) * r / mort_ndays
   end do
   close(u)
end subroutine alloc_write

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
   d  = real(mort_ndays, r8) * real(mort_nsel, r8)    ! site-days: state means
   fd = 365.0_r8 / real(mort_ndays, r8) / real(mort_nsel, r8)   ! flux sums -> per year, site mean
   do ft = 1, numpft
      do il = 1, 2
         associate(a => mort_acc(:, ft, il))
         write(u,'(i0,",",i0,",",i0,",",i0,21(",",es12.5))') mort_year, mort_ndays, ft, il, &
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

subroutine mort_select(nsites)
   ! sites included in the mort file: all, or those whose patch ID is in env RHESSYS_FATES_MORT_PATCHES
   integer, intent(in) :: nsites
   character(len=4096) :: ev
   character(len=32) :: tok
   integer :: s, i, j, n, ios, pid, stat
   allocate(mort_use(nsites))
   call get_environment_variable('RHESSYS_FATES_MORT_PATCHES', ev, status=stat)
   if (stat /= 0 .or. len_trim(ev) == 0) then
      mort_use = .true.
   else
      mort_use = .false.
      i = 1; n = len_trim(ev)
      do while (i <= n)
         j = index(ev(i:n), ',')
         if (j == 0) then
            tok = ev(i:n); i = n + 1
         else
            tok = ev(i:i+j-2); i = i + j
         end if
         read(tok, *, iostat=ios) pid
         if (ios /= 0) cycle
         do s = 1, nsites
            if (fates%sites(s)%h_gid == pid) mort_use(s) = .true.
         end do
      end do
      write(*,'(a,i0,a)') 'FATES host: mort file restricted to ', count(mort_use), ' sites (RHESSYS_FATES_MORT_PATCHES)'
   end if
   mort_nsel = max(count(mort_use), 1)
end subroutine mort_select

subroutine mort_accumulate(nsites)
   integer, intent(in) :: nsites
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort
   real(r8) :: nm2, ctot, leafc, fnrtc, storec, mr
   integer  :: s, ft, il, io
   real(r8) :: cum(4, numpft, 2)

   if (len_trim(mort_file) == 0) return
   if (.not. allocated(mort_use)) call mort_select(nsites)
   if (.not. allocated(mort_acc)) then
      allocate(mort_acc(nmd, numpft, 2), term_base(4, numpft, 2))
      mort_acc = 0.0_r8
      call mort_term_cum(nsites, cum)
      term_base = cum
      mort_year = hlm_current_year
      allocate(alloc_acc(2, n_aorg, numpft), alloc_stock0(n_aorg, numpft), alloc_jul(2, numpft), alloc_nst(2, numpft))
      alloc_acc = 0.0_r8; alloc_jul = 0.0_r8; alloc_jul_days = 0; alloc_nst = 0.0_r8
      call alloc_stocks(nsites, alloc_stock0)
   end if
   if (hlm_current_year /= mort_year) then
      call alloc_write(nsites)          ! before mort_write, which moves term_base on
      alloc_acc = 0.0_r8; alloc_jul = 0.0_r8; alloc_jul_days = 0; alloc_nst = 0.0_r8
      call alloc_stocks(nsites, alloc_stock0)
      call mort_write(nsites)
      mort_acc = 0.0_r8
      mort_ndays = 0
      mort_year = hlm_current_year
   end if
   mort_ndays = mort_ndays + 1
   if (hlm_current_month == 7) alloc_jul_days = alloc_jul_days + 1
   do s = 1, nsites
      if (.not. mort_use(s)) cycle
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
            ! today's recruits have no mortality rates yet (FATES leaves hmort
            ! uninitialised until their first mortality_rates call)
            if (.not. ccohort%isnew) then
            a(4)  = a(4)  + nm2 * ctot * fin(ccohort%bmort)  / 365.0_r8
            a(5)  = a(5)  + nm2 * ctot * fin(ccohort%cmort)  / 365.0_r8
            a(6)  = a(6)  + nm2 * ctot * fin(ccohort%hmort)  / 365.0_r8
            a(7)  = a(7)  + nm2 * ctot * fin(ccohort%frmort) / 365.0_r8
            a(8)  = a(8)  + nm2 * ctot * fin(ccohort%smort)  / 365.0_r8
            a(9)  = a(9)  + nm2 * ctot * fin(ccohort%asmort) / 365.0_r8
            a(10) = a(10) + nm2 * ctot * fin(ccohort%dgmort) / 365.0_r8
            end if
            a(11) = a(11) + nm2 * storec
            a(14) = a(14) + nm2 * ctot
            a(18) = a(18) + nm2 * leafc
            mr = 0.0_r8
            if (.not. ccohort%isnew) mr = (fin(ccohort%bmort) + fin(ccohort%cmort) + fin(ccohort%hmort) + &
                 fin(ccohort%frmort) + fin(ccohort%smort) + fin(ccohort%asmort) + fin(ccohort%dgmort)) / 365.0_r8
            do io = 1, n_aorg
               alloc_acc(1, io, ft) = alloc_acc(1, io, ft) + nm2 * fin(ccohort%prt%GetTurnover(aorg(io), carbon12_element))
               alloc_acc(2, io, ft) = alloc_acc(2, io, ft) + nm2 * ccohort%prt%GetState(aorg(io), carbon12_element) * mr
            end do
            if (cnp_coupled) then
               alloc_nst(1, ft) = alloc_nst(1, ft) + nm2 * ccohort%prt%GetState(store_organ, nitrogen_element)
               alloc_nst(2, ft) = alloc_nst(2, ft) + nm2 * &
                    ccohort%prt%GetNutrientTarget(nitrogen_element, store_organ, stoich_growth_min)
            end if
            if (hlm_current_month == 7) then
               alloc_jul(1, ft) = alloc_jul(1, ft) + ccohort%treelai * ccohort%c_area * area_inv
               alloc_jul(2, ft) = alloc_jul(2, ft) + nm2 * leafc
            end if
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
! Plant N demand mode (see ndemand_mode); called by the host before the run
!-------------------------------------------------------------------------------
subroutine rhessys_veg_set_l2fr_max(x) bind(C, name='rhessys_veg_set_l2fr_max')
   ! upper bound on FATES's dynamic leaf-to-fineroot ratio (host patch variable; <= 0 = none)
   use PRTAllometricCNPMod, only : host_l2fr_max
   real(c_double), intent(in) :: x
   host_l2fr_max = x
end subroutine rhessys_veg_set_l2fr_max

subroutine rhessys_veg_set_nforage(cost, mask) bind(C, name='rhessys_veg_set_nforage')
   ! -demog_nforage: cost (gC/gN) and PFT bit mask of AM foraging; cost <= 0 = off
   real(c_double), intent(in) :: cost
   integer(c_int), intent(in) :: mask
   nforage_cost = cost
   nforage_mask = mask
end subroutine rhessys_veg_set_nforage

subroutine rhessys_veg_set_nfun(k, charge, mask) bind(C, name='rhessys_veg_set_nfun')
   ! -demog_nfun: k = kN_AM, kC_AM, kN_ECM, kC_ECM; charge 1 = cohorts pay N x cost from NPP;
   ! mask = ECM PFTs (bit k-1 = PFT k); k(1) <= 0 = off
   real(c_double), intent(in) :: k(4)
   integer(c_int), intent(in) :: charge, mask
   nfun_k = k
   nfun_charge = charge
   nfun_mask = mask
end subroutine rhessys_veg_set_nfun

subroutine rhessys_veg_nfun_diag(site, v) bind(C, name='rhessys_veg_nfun_diag')
   ! -demog_nfun: today's site terms (g/m2/day): N to AM, N to ECM, N x cost AM, N x cost ECM, C charged
   integer(c_int), intent(in)  :: site
   real(c_double), intent(out) :: v(5)
   v = 0.0_r8
   if (.not. allocated(nfun_diag)) return
   if (site < 1 .or. site > size(nfun_diag, 2)) return
   v = nfun_diag(:, site)
end subroutine rhessys_veg_nfun_diag

subroutine rhessys_veg_set_browse(rate, height, mask) bind(C, name='rhessys_veg_set_browse')
   ! -demog_browse: browse mortality (1/yr) of tree cohorts shorter than height (m) in the PFT bit mask
   ! (FATES host patch variables in EDMortalityFunctionsMod); rate <= 0 = off
   use EDMortalityFunctionsMod, only : host_browse_rate, host_browse_height, host_browse_mask
   real(c_double), intent(in) :: rate, height
   integer(c_int), intent(in) :: mask
   host_browse_rate = rate
   host_browse_height = height
   host_browse_mask = mask
end subroutine rhessys_veg_set_browse

subroutine rhessys_veg_set_suppress(gcrit, dbhmax, npft, rate, minlayer) bind(C, name='rhessys_veg_set_suppress')
   ! -demog_suppress: growth-dependent mortality, rate(pft) (1/yr) x (1 - g/gcrit) for cohorts with dbh < dbhmax
   ! (cm) whose ~1-yr mean dbh growth g < gcrit (cm/yr) (FATES host patch variables in
   ! EDMortalityFunctionsMod) in canopy layer >= minlayer (1 = all, 2 = understory only); gcrit <= 0 = off
   use EDMortalityFunctionsMod, only : host_supp_rate, host_supp_gcrit, host_supp_dbhmax, host_supp_minlayer
   real(c_double), intent(in) :: gcrit, dbhmax
   integer(c_int), intent(in) :: npft, minlayer
   real(c_double), intent(in) :: rate(npft)
   host_supp_rate(:) = 0.0_r8
   host_supp_rate(1:min(npft, size(host_supp_rate))) = rate(1:min(npft, size(host_supp_rate)))
   host_supp_gcrit = gcrit
   host_supp_dbhmax = dbhmax
   host_supp_minlayer = max(1, minlayer)
end subroutine rhessys_veg_set_suppress

subroutine rhessys_veg_set_mort_disturb(site, frac) bind(C, name='rhessys_veg_set_mort_disturb')
   ! -demog_mort_disturb: this site's share of canopy mortality that opens a new (gap) patch
   ! (FATES host patch field ed_site_type%host_mort_disturb_frac; < 0 = fates_mort_disturb_frac)
   integer(c_int), intent(in) :: site
   real(c_double), intent(in) :: frac
   if (site < 1 .or. site > fates%nsites) return
   fates%sites(site)%host_mort_disturb_frac = frac
end subroutine rhessys_veg_set_mort_disturb

subroutine rhessys_veg_set_nmine(cost, mask) bind(C, name='rhessys_veg_set_nmine')
   ! -demog_nmine: cost (gC/gN) and PFT bit mask of C-paid N acquisition; cost <= 0 = off
   real(c_double), intent(in) :: cost
   integer(c_int), intent(in) :: mask
   nmine_cost = cost
   nmine_mask = mask
end subroutine rhessys_veg_set_nmine

subroutine rhessys_veg_set_nvcmax(x) bind(C, name='rhessys_veg_set_nvcmax')
   ! photosynthetic capacity follows leaf N status (host patch variable host_nvcmax_floor in
   ! FatesCohortMod): x = floor on the scaling factor, in (0, 1]; <= 0 = stock FATES
   use FatesCohortMod, only : host_nvcmax_floor
   real(c_double), intent(in) :: x
   if (x > 0.0_r8) then
      host_nvcmax_floor = min(1.0_r8, x)
   else
      host_nvcmax_floor = -1.0_r8
   end if
end subroutine rhessys_veg_set_nvcmax

subroutine rhessys_veg_set_min_n(x) bind(C, name='rhessys_veg_set_min_n')
   ! extra FATES cohort termination density (host patch variable host_min_npm2, plants per m2 of
   ! patch); x in stems/ha, <= 0 = stock FATES
   use EDCohortDynamicsMod, only : host_min_npm2
   real(c_double), intent(in) :: x
   if (x > 0.0_r8) then
      host_min_npm2 = x * 1.0e-4_r8
   else
      host_min_npm2 = -1.0_r8
   end if
end subroutine rhessys_veg_set_min_n

subroutine rhessys_veg_set_stomatal(model) bind(C, name='rhessys_veg_set_stomatal')
   ! stomatal conductance model for set_fates_ctrlparms: 1 = Ball-Berry, 2 = Medlyn; call before rhessys_veg_init
   integer(c_int), intent(in) :: model
   stomatal_model_sel = merge(2, 1, model == 2)
end subroutine rhessys_veg_set_stomatal

subroutine rhessys_veg_set_regen(model) bind(C, name='rhessys_veg_set_regen')
   ! regeneration model for set_fates_ctrlparms: 1 = default, 2 = Tree Recruitment Scheme with
   ! seedling dynamics (Hanbury-Brown et al. 2022); call before rhessys_veg_init
   integer(c_int), intent(in) :: model
   regen_model_sel = merge(TRS_regeneration, default_regeneration, model == 2)
end subroutine rhessys_veg_set_regen

subroutine rhessys_veg_set_coage(on) bind(C, name='rhessys_veg_set_coage')
   ! cohort age tracking for set_fates_ctrlparms (needed for age senescence); call before rhessys_veg_init
   integer(c_int), intent(in) :: on
   coage_sel = merge(1, 0, on == 1)
end subroutine rhessys_veg_set_coage

subroutine rhessys_veg_set_dayl(on) bind(C, name='rhessys_veg_set_dayl')
   ! day-length factor on photosynthetic capacity for set_fates_ctrlparms; call before rhessys_veg_init
   integer(c_int), intent(in) :: on
   dayl_sel = merge(1, 0, on == 1)
end subroutine rhessys_veg_set_dayl

subroutine rhessys_veg_set_disperse(on) bind(C, name='rhessys_veg_set_disperse')
   ! seed dispersal between sites (FATES cadence daily); call before rhessys_veg_init
   integer(c_int), intent(in) :: on
   disp_sel = merge(1, 0, on == 1)
end subroutine rhessys_veg_set_disperse

subroutine rhessys_veg_set_disperse_geom(n, x, y, area) bind(C, name='rhessys_veg_set_disperse_geom')
   ! site centroids (m) and areas (m2) by FATES site index; builds the neighbour lists and kernel shares.
   ! Call after rhessys_veg_init (needs the FATES parameters).
   integer(c_int), intent(in) :: n
   real(c_double), intent(in) :: x(n), y(n), area(n)
   integer  :: i, j, k, ft, nfix
   real(r8) :: d, dmax, a, f, wsum, pi
   real(r8), allocatable :: rad(:)
   pi = 4.0_r8 * atan(1.0_r8)
   if (disp_sel /= 1) return
   dmax = 0.0_r8
   do ft = 1, numpft
      if (.not. (EDPftvarcon_inst%seed_dispersal_fraction(ft) >= 0.0_r8 .and. &
                 EDPftvarcon_inst%seed_dispersal_fraction(ft) <= 1.0_r8)) then
         write(fates_log_unit,*) 'rhessys_fates: -demog_disperse: fates_seed_dispersal_fraction must be set (0-1) ', &
              'for every pft; pft ', ft, ' = ', EDPftvarcon_inst%seed_dispersal_fraction(ft)
         error stop 'rhessys_fates: fates_seed_dispersal_fraction unset'
      end if
      if (EDPftvarcon_inst%seed_dispersal_fraction(ft) > 0.0_r8) then
         a = EDPftvarcon_inst%seed_dispersal_pdf_scale(ft); d = EDPftvarcon_inst%seed_dispersal_max_dist(ft)
         if (.not. (a > 0.0_r8 .and. a < 1.0e6_r8 .and. d > 0.0_r8 .and. d < 1.0e6_r8)) then
            write(fates_log_unit,*) 'rhessys_fates: -demog_disperse: pft ', ft, ' disperses seed but ', &
                 'fates_seed_dispersal_pdf_scale / _max_dist are not set (m): ', a, d
            error stop 'rhessys_fates: seed dispersal parameters missing'
         end if
         dmax = max(dmax, d)
      end if
   end do
   allocate(disp_nnb(n), disp_area(n), disp_tot(3, numpft), rad(n))
   disp_area(:) = area(:); disp_tot = 0.0_r8; disp_nnb = 0
   rad(:) = sqrt(area(:) / pi)
   do i = 1, n
      do j = 1, n
         if (j /= i .and. hypot(x(i) - x(j), y(i) - y(j)) - rad(i) - rad(j) <= dmax) disp_nnb(i) = disp_nnb(i) + 1
      end do
   end do
   disp_maxnb = max(1, maxval(disp_nnb))
   allocate(disp_nb(disp_maxnb, n), disp_w(numpft, disp_maxnb, n))
   disp_nb = 0; disp_w = 0.0_r8; nfix = 0
   do i = 1, n
      k = 0
      do j = 1, n
         d = hypot(x(i) - x(j), y(i) - y(j)) - rad(i) - rad(j)    ! edge gap (m)
         if (j == i .or. d > dmax) cycle
         k = k + 1; disp_nb(k, i) = j
         do ft = 1, numpft
            if (EDPftvarcon_inst%seed_dispersal_fraction(ft) <= 0.0_r8) cycle
            if (d > EDPftvarcon_inst%seed_dispersal_max_dist(ft)) cycle
            a = EDPftvarcon_inst%seed_dispersal_pdf_scale(ft)
            disp_w(ft, k, i) = exp(-max(d, 0.0_r8) / a) * rad(j)
         end do
      end do
      do ft = 1, numpft
         wsum = sum(disp_w(ft, 1:disp_nnb(i), i))
         if (wsum > 0.0_r8) then
            disp_w(ft, 1:disp_nnb(i), i) = disp_w(ft, 1:disp_nnb(i), i) / wsum
         else if (EDPftvarcon_inst%seed_dispersal_fraction(ft) > 0.0_r8) then
            nfix = nfix + 1                                        ! no neighbour: dispersed seed lost
         end if
      end do
   end do
   write(fates_log_unit,'(a,i0,a,f8.1,a,f8.1,a,i0,a)') ' rhessys_fates: seed dispersal: ', n, ' sites, max distance ', dmax, &
        ' m (edge gap), mean ', real(sum(disp_nnb), r8) / real(n, r8), ' neighbours per site; ', nfix, &
        ' site x pft without a neighbour (dispersed seed lost)'
   do ft = 1, numpft
      if (EDPftvarcon_inst%seed_dispersal_fraction(ft) > 0.0_r8) &
         write(fates_log_unit,'(a,i0,a,f6.3,a,f8.1,a,f8.3)') '   pft ', ft, ': dispersed fraction ', &
              EDPftvarcon_inst%seed_dispersal_fraction(ft), ', kernel scale (m) ', EDPftvarcon_inst%seed_dispersal_pdf_scale(ft), &
              ', mean share landing on other vegetated sites ', sum(disp_w(ft, :, :)) / real(n, r8)
   end do
   deallocate(rad)
end subroutine rhessys_veg_set_disperse_geom

!-------------------------------------------------------------------------------
! -demog_disperse: today's dispersed seed (site%seed_out, kgC per FATES site area per day, C only -- host
! patch to SeedUpdate) spread over the neighbours; site%seed_in = what each site receives, used by
! tomorrow's SeedUpdate; seed_out reset. Serial (after the OpenMP site loop). Yearly totals to the log.
!-------------------------------------------------------------------------------
subroutine disperse_seeds(nsites)
   use EDTypesMod, only : area
   integer, intent(in) :: nsites
   integer  :: i, k, j, ft
   real(r8) :: out_kg, got
   real(r8), allocatable :: recv(:)
   if (.not. allocated(disp_w)) return
   if (disp_year < 0) disp_year = hlm_current_year
   if (hlm_current_year /= disp_year) then
      do ft = 1, numpft
         if (disp_tot(1, ft) > 0.0_r8) write(fates_log_unit,'(a,i0,a,i0,3(a,es11.4))') ' rhessys_fates: seed dispersal ', &
              disp_year, ' pft ', ft, ': out ', disp_tot(1, ft), ' kgC, received ', disp_tot(2, ft), ' kgC, lost ', disp_tot(3, ft)
      end do
      disp_tot = 0.0_r8; disp_year = hlm_current_year
   end if
   allocate(recv(nsites))
   do ft = 1, numpft
      recv = 0.0_r8
      do i = 1, nsites
         out_kg = fates%sites(i)%seed_out(ft) / area * disp_area(i)     ! kgC/day leaving site i
         if (out_kg <= 0.0_r8) cycle
         got = 0.0_r8
         do k = 1, disp_nnb(i)
            j = disp_nb(k, i)
            recv(j) = recv(j) + out_kg * disp_w(ft, k, i)
            got = got + out_kg * disp_w(ft, k, i)
         end do
         disp_tot(1, ft) = disp_tot(1, ft) + out_kg
         disp_tot(2, ft) = disp_tot(2, ft) + got
         disp_tot(3, ft) = disp_tot(3, ft) + (out_kg - got)
      end do
      do j = 1, nsites
         fates%sites(j)%seed_in(ft) = recv(j) / disp_area(j) * area          ! kgC per FATES site area per day
         fates%sites(j)%seed_out(ft) = 0.0_r8
      end do
   end do
   deallocate(recv)
end subroutine disperse_seeds

subroutine rhessys_veg_set_ndemand(mode, mult) bind(C, name='rhessys_veg_set_ndemand')
   integer(c_int), intent(in) :: mode
   real(c_double), intent(in) :: mult
   ! bit 1 = growth-based N demand / plant priority (-demog_nneed); bit 2 = store-regulated (-demog_nstore_reg);
   ! bit 3 = need-capped (-demog_ncap); bit 4 = cohort split by need (-demog_nsplit, split only)
   ndemand_mode = iand(int(mode), 7)
   nsplit_need = btest(int(mode), 3)
   nsplit_min  = btest(int(mode), 4)
   ncap_mult = mult
end subroutine rhessys_veg_set_ndemand

!-------------------------------------------------------------------------------
! A cohort's N uptake claim per plant (kgN/plant/s): RD capacity, or (mode bit 1) the
! larger of capacity and potential growth N need; (mode bit 2) times the cohort's N store room;
! (mode bit 3) at most the growth need plus the store refill
!-------------------------------------------------------------------------------
real(r8) function n_claim(c)
   type(fates_cohort_type), pointer, intent(in) :: c
   integer  :: ft
   ft  = c%pft
   n_claim = c%prt%GetState(fnrt_organ, carbon12_element) * &
        (EDPftvarcon_inst%vmax_nh4(ft) + EDPftvarcon_inst%vmax_no3(ft))
   if (iand(ndemand_mode, 1) == 1) n_claim = max(n_claim, n_need(c))
   if (iand(ndemand_mode, 2) == 2) n_claim = n_claim * n_store_room(c)
   if (iand(ndemand_mode, 4) == 4) n_claim = min(n_claim, ncap_mult * (n_need(c) + n_store_deficit(c)))
end function n_claim

! N to refill a cohort's N store to its target over ncap_refill_days (kgN/plant/s)
real(r8) function n_store_deficit(c)
   type(fates_cohort_type), pointer, intent(in) :: c
   real(r8) :: tgt
   n_store_deficit = 0.0_r8
   if (.not. cnp_coupled) return
   tgt = c%prt%GetNutrientTarget(nitrogen_element, store_organ, stoich_growth_min)
   n_store_deficit = max(0.0_r8, tgt - c%prt%GetState(store_organ, nitrogen_element)) / &
        (ncap_refill_days * 86400.0_r8)
end function n_store_deficit

! room left in a cohort's N store, 0-1: max(0, 1 - N store / FATES's N store overflow level)
real(r8) function n_store_room(c)
   type(fates_cohort_type), pointer, intent(in) :: c
   real(r8) :: tgt
   n_store_room = 1.0_r8
   if (.not. cnp_coupled) return
   tgt = c%prt%GetNutrientTarget(nitrogen_element, store_organ, stoich_max)
   if (tgt > 0.0_r8) n_store_room = min(1.0_r8, max(0.0_r8, 1.0_r8 - c%prt%GetState(store_organ, nitrogen_element) / tgt))
end function n_store_room

! potential growth N need per plant (kgN/plant/s): potential C gain for growth
! (annualised GPP - maintenance resp, less growth resp) x tissue N:C at target
! stoichiometry (leaf, fine root, sapwood, structure)
real(r8) function n_need(c)
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   type(fates_cohort_type), pointer, intent(in) :: c
   real(r8) :: gain, cl, cf, cs, ct, nt
   integer  :: ft
   ft  = c%pft
   ! potential C gain for growth (kgC/plant/s): annualised GPP - maintenance resp,
   ! less growth respiration
   gain = 0.0_r8
   if (ieee_is_finite(c%gpp_acc_hold) .and. ieee_is_finite(c%resp_m_acc_hold)) &
        gain = max(c%gpp_acc_hold - c%resp_m_acc_hold, 0.0_r8) * &
               (1.0_r8 - prt_params%grperc(ft)) / (365.0_r8 * 86400.0_r8)
   ! tissue N:C at target stoichiometry (leaf, fine root, sapwood, structure)
   cl = c%prt%GetState(leaf_organ,   carbon12_element)
   cf = c%prt%GetState(fnrt_organ,   carbon12_element)
   cs = c%prt%GetState(sapw_organ,   carbon12_element)
   ct = c%prt%GetState(struct_organ, carbon12_element)
   nt = cl * prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(leaf_organ)) + &
        cf * prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(fnrt_organ)) + &
        cs * prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(sapw_organ)) + &
        ct * prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(struct_organ))
   if (cl + cf + cs + ct > 0.0_r8) then
      n_need = gain * nt / (cl + cf + cs + ct)
   else
      n_need = gain * prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(leaf_organ))
   end if
end function n_need

!-------------------------------------------------------------------------------
! Site potential growth N need (kgN/m2/day), for the host's plant-priority uptake
!-------------------------------------------------------------------------------
subroutine rhessys_veg_site_nneed(site, need) bind(C, name='rhessys_veg_site_nneed')
   integer(c_int), intent(in)  :: site
   real(c_double), intent(out) :: need
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort
   need = 0.0_r8
   if (site < 1 .or. site > fates%nsites) return
   cpatch => fates%sites(site)%oldest_patch
   do while (associated(cpatch))
      ccohort => cpatch%tallest
      do while (associated(ccohort))
         need = need + n_need(ccohort) * ccohort%n
         ccohort => ccohort%shorter
      end do
      cpatch => cpatch%younger
   end do
   need = need * area_inv * 86400.0_r8
end subroutine rhessys_veg_site_nneed

!-------------------------------------------------------------------------------
! Split the site's daily N uptake (gN/m2) over the cohorts (FATES competitors in
! coupled N uptake mode: one per cohort, in oldest-patch / tallest-cohort order)
! in proportion to their RD capacity, fine-root C x (Vmax NH4 + Vmax NO3).
! recv = total actually placed (0 if no cohort has fine roots).
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! -demog_nmine: add the mined N (gN/m2/day, granted by the host yesterday) to the NH4
! uptake of the mining cohorts, shared by their excess C of yesterday (resp_excess_hold
! still holds it: FATES resets it in today's dynamics). If no mining cohort has excess C
! left (cohorts died or fused), share it over all cohorts by fine-root C, so no N is lost.
! recv = N placed (gN/m2/day). Call after distribute_n_uptake (adds to its fluxes).
!-------------------------------------------------------------------------------
real(r8) function nm_fin(x)
   ! x if finite and sane, else 0 (resp_excess_hold is NaN on cohorts created today)
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   real(r8), intent(in) :: x
   nm_fin = 0.0_r8
   if (ieee_is_finite(x)) then
      if (abs(x) < 1.0e30_r8) nm_fin = x
   end if
end function nm_fin

subroutine distribute_n_excess(s, nh4, no3, mask, recv)
   ! NH4 and NO3 (gN/m2/day) added to the uptake of the cohorts whose PFT is in mask, shared by
   ! their excess C of yesterday (resp_excess_hold still holds it: FATES resets it in today's
   ! dynamics). If none of them has excess C left (cohorts died or fused), shared over all
   ! cohorts by fine-root C, so no N is lost. recv = N placed. Call after distribute_n_uptake.
   integer,  intent(in)  :: s, mask
   real(r8), intent(in)  :: nh4, no3
   real(r8), intent(out) :: recv
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort
   real(r8) :: wsum, w
   integer  :: icomp, pass, mode

   recv = 0.0_r8
   associate(bi => fates%bc_in(s))
      do mode = 1, 2
         wsum = 0.0_r8
         do pass = 1, 2
            icomp = 0
            cpatch => fates%sites(s)%oldest_patch
            do while (associated(cpatch))
               ccohort => cpatch%tallest
               do while (associated(ccohort))
                  icomp = icomp + 1
                  if (mode == 1) then
                     w = 0.0_r8
                     if (btest(mask, ccohort%pft - 1)) w = max(nm_fin(ccohort%resp_excess_hold), 0.0_r8) * ccohort%n
                  else
                     w = ccohort%prt%GetState(fnrt_organ, carbon12_element) * ccohort%n
                  end if
                  if (pass == 1) then
                     wsum = wsum + w
                  else if (wsum > 0.0_r8 .and. icomp <= size(bi%plant_nh4_uptake_flux, 1)) then
                     bi%plant_nh4_uptake_flux(icomp,1) = bi%plant_nh4_uptake_flux(icomp,1) + nh4 * w / wsum
                     bi%plant_no3_uptake_flux(icomp,1) = bi%plant_no3_uptake_flux(icomp,1) + no3 * w / wsum
                     recv = recv + (nh4 + no3) * w / wsum
                  end if
                  ccohort => ccohort%shorter
               end do
               cpatch => cpatch%younger
            end do
         end do
         if (wsum > 0.0_r8) exit
      end do
   end associate
end subroutine distribute_n_excess

subroutine distribute_n_uptake(s, nh4, no3, recv, navail)
   integer,  intent(in)  :: s
   real(r8), intent(in)  :: nh4, no3
   real(r8), intent(out) :: recv
   real(r8), intent(in)  :: navail                ! -demog_nfun: root-zone mineral N (gN/m2)
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort
   real(r8) :: wsum, w, cost, ni, wbase_sum
   integer  :: icomp, pass, k
   logical  :: use_min

   recv = 0.0_r8
   wsum = 0.0_r8
   wbase_sum = 0.0_r8
   use_min = nsplit_min
   associate(bi => fates%bc_in(s))
      bi%plant_nh4_uptake_flux(:,:) = 0.0_r8
      bi%plant_no3_uptake_flux(:,:) = 0.0_r8
      do pass = 1, 2
         ! -demog_nsplit 2: if no cohort has a growth N need today (e.g. leaf-off), share by the claim as without
         ! the flag, so the N RHESSys granted is still all placed (no N lost); nfun cost is not used with it
         if (pass == 2 .and. use_min .and. wsum <= 0.0_r8 .and. wbase_sum > 0.0_r8) then
            use_min = .false.; wsum = wbase_sum
         end if
         icomp = 0
         cpatch => fates%sites(s)%oldest_patch
         do while (associated(cpatch))
            ccohort => cpatch%tallest
            do while (associated(ccohort))
               icomp = icomp + 1
               if (ndemand_mode /= 0) then
                  w = n_claim(ccohort) * ccohort%n
               else
                  w = ccohort%prt%GetState(fnrt_organ, carbon12_element) * ccohort%n * &
                      (EDPftvarcon_inst%vmax_nh4(ccohort%pft) + EDPftvarcon_inst%vmax_no3(ccohort%pft))
               end if
               if (nsplit_need) w = max(w, n_need(ccohort) * ccohort%n)
               if (pass == 1) wbase_sum = wbase_sum + w
               if (use_min)     w = min(w, n_need(ccohort) * ccohort%n)
               cost = 0.0_r8
               if (nfun_k(1) > 0.0_r8) then
                  cost = nfun_cost(ccohort, navail)
                  w = w / cost
               end if
               if (pass == 1) then
                  wsum = wsum + w
               else if (wsum > 0.0_r8 .and. icomp <= size(bi%plant_nh4_uptake_flux, 1)) then
                  bi%plant_nh4_uptake_flux(icomp,1) = nh4 * w / wsum
                  bi%plant_no3_uptake_flux(icomp,1) = no3 * w / wsum
                  ni = (nh4 + no3) * w / wsum                    ! gN/m2/day
                  recv = recv + ni
                  if (nfun_k(1) > 0.0_r8 .and. ni > 0.0_r8) then
                     k = 1
                     if (btest(nfun_mask, ccohort%pft - 1)) k = 2
                     nfun_diag(k, s) = nfun_diag(k, s) + ni
                     nfun_diag(k + 2, s) = nfun_diag(k + 2, s) + ni * cost
                     if (nfun_charge == 1 .and. ccohort%n > 0.0_r8) then
                        ! N x cost (gC/m2/day) -> kgC/indiv/day, paid from today's NPP
                        ccohort%resp_m_acc = ccohort%resp_m_acc + ni * cost * 1.0e-3_r8 / (ccohort%n * area_inv)
                        nfun_diag(5, s) = nfun_diag(5, s) + ni * cost
                     end if
                  end if
               end if
               ccohort => ccohort%shorter
            end do
            cpatch => cpatch%younger
         end do
      end do
   end associate
end subroutine distribute_n_uptake

!-------------------------------------------------------------------------------
! -demog_nfun: FUN-style cost of active N uptake (gC/gN) for one cohort, as CLM5
! fun_cost_active: kN / mineral N + kC / root C, with the cohort's fine-root C per
! m2 of its crown as root C and AM or ECM parameters by PFT.
!-------------------------------------------------------------------------------
real(r8) function nfun_cost(c, navail)
   type(fates_cohort_type), intent(in) :: c
   real(r8), intent(in) :: navail
   real(r8) :: rootc, kn, kc
   real(r8), parameter :: big_cost = 1.0e6_r8
   if (btest(nfun_mask, c%pft - 1)) then
      kn = nfun_k(3); kc = nfun_k(4)
   else
      kn = nfun_k(1); kc = nfun_k(2)
   end if
   rootc = 0.0_r8
   if (c%c_area > 0.0_r8 .and. c%n > 0.0_r8) &
      rootc = 1000.0_r8 * c%prt%GetState(fnrt_organ, carbon12_element) / (c%c_area / c%n)
   if (rootc > 1.0e-6_r8 .and. navail > 1.0e-9_r8) then
      nfun_cost = min(kn / navail + kc / rootc, big_cost)
   else
      nfun_cost = big_cost
   end if
end function nfun_cost

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
         if (ndemand_mode /= 0) then
            out%ncap = out%ncap + n_claim(ccohort) * ccohort%n * area_inv * 1000.0_r8 * 86400.0_r8
         else
            out%ncap = out%ncap + ccohort%prt%GetState(fnrt_organ, carbon12_element) * &
                 ccohort%n * area_inv * 1000.0_r8 * 86400.0_r8 * &
                 (EDPftvarcon_inst%vmax_nh4(ccohort%pft) + EDPftvarcon_inst%vmax_no3(ccohort%pft))
         end if
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
   ! FATES's net_root_uptake includes symbiotic fixation (daily_n_gain = NH4 + NO3 + fixation)
   if (cnp_coupled) out%efflux_n = max(out%n_uptake_recv + out%n_fix - &
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
   if (kind == 4) then
      call apply_wmfire(s, bc)
      return
   end if
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
! WMFire fire on site s (dist_kind 4, -demog 3 + -firespread): WMFire decides where
! and how hard it burns; FATES does the fire. In the burned fraction wm_frac of every
! FATES patch, cohort mortality follows WMFire's own rules (compute_family_fire_effects)
! at the mean intensity p = wm_burn of the burned cells:
!   understory (height <= under_h): m_u = (k^p - 1)/(k - 1), k = understory_mort (k = 1: p)
!   overstory  (height >  over_h) : m_o = 1 - 1/(1 + exp(-k1 (U - k2))), U = understory
!       C consumed (kgC/m2: m_u x consumed share x above-ground C of understory cohorts)
!       + litter C consumed per m2 burned (host)
!   in between: linear blend by height. Consumed share of the killed (crown burned):
!   same form with k = consumption. Litter / CWD fuel consumed: WMFire's fractions
!   (leaves 1, twigs + small branches 0.85, large branches + trunks 0.339).
!-------------------------------------------------------------------------------
subroutine apply_wmfire(s, bc)
   integer,        intent(in) :: s
   type(veg_bc_t), intent(in) :: bc
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort
   real(r8) :: p, mu, cu, mo, m, up, under_c, U, hu, ho, ag
   hu = bc%wm_under_h; ho = max(bc%wm_over_h, hu + 1.0e-6_r8)
   p  = min(max(bc%wm_burn, 0.0_r8), 1.0_r8)
   mu = wm_mort(bc%wm_under_mort, p)
   cu = wm_mort(bc%wm_consumption, mu)
   under_c = 0.0_r8
   cpatch => fates%sites(s)%oldest_patch
   do while (associated(cpatch))
      ccohort => cpatch%tallest
      do while (associated(ccohort))
         if (ccohort%height <= hu) then
            ag = prt_params%allom_agb_frac(ccohort%pft)
            under_c = under_c + ccohort%n * area_inv * ( &
                 ccohort%prt%GetState(leaf_organ, carbon12_element) + ag * ( &
                 ccohort%prt%GetState(sapw_organ, carbon12_element) + &
                 ccohort%prt%GetState(struct_organ, carbon12_element)))
         end if
         ccohort => ccohort%shorter
      end do
      cpatch => cpatch%younger
   end do
   U  = under_c * mu * cu + max(bc%wm_litter_consumed, 0.0_r8)
   mo = 1.0_r8 - 1.0_r8 / (1.0_r8 + exp(-(bc%wm_over_k1 * (U - bc%wm_over_k2))))
   cpatch => fates%sites(s)%oldest_patch
   do while (associated(cpatch))
      cpatch%fire = 1
      cpatch%frac_burnt = min(max(bc%wm_frac, 0.0_r8), 0.98_r8)
      cpatch%fuel%frac_burnt(:) = 0.0_r8
      cpatch%fuel%frac_burnt(1:2) = 0.85_r8
      cpatch%fuel%frac_burnt(3:4) = 0.339_r8
      cpatch%fuel%frac_burnt(5:6) = 1.0_r8
      ccohort => cpatch%tallest
      do while (associated(ccohort))
         if (ccohort%height <= hu) then
            m = mu
         else if (ccohort%height > ho) then
            m = mo
         else
            up = (ho - ccohort%height) / (ho - hu)
            m = up * mu + (1.0_r8 - up) * mo
         end if
         ccohort%fire_mort = min(max(m, 0.0_r8), 1.0_r8)
         ccohort%fraction_crown_burned = min(max(wm_mort(bc%wm_consumption, m), 0.0_r8), 1.0_r8)
         ccohort => ccohort%shorter
      end do
      cpatch => cpatch%younger
   end do
contains
   real(r8) function wm_mort(k, x)
      real(r8), intent(in) :: k, x
      if (abs(k - 1.0_r8) < 1.0e-9_r8) then
         wm_mort = k * x
      else if (k > 0.0_r8) then
         wm_mort = (k**x - 1.0_r8) / (k - 1.0_r8)
      else
         wm_mort = x
      end if
   end function wm_mort
end subroutine apply_wmfire

!-------------------------------------------------------------------------------
! Per-patch inventories: inventory index (0-based, order of the host-written
! FATES control file) of every site; used by the next cold start
!-------------------------------------------------------------------------------
subroutine rhessys_veg_inv_assign(nsites, idx) bind(C, name='rhessys_veg_inv_assign')
   integer(c_int), intent(in) :: nsites
   integer(c_int), intent(in) :: idx(nsites)
   if (allocated(inv_idx)) deallocate(inv_idx)
   allocate(inv_idx(nsites))
   inv_idx = idx
end subroutine rhessys_veg_inv_assign

!-------------------------------------------------------------------------------
! Today's FATES litter C to the host by source (gC/m2/day), see lit_src
!-------------------------------------------------------------------------------
subroutine rhessys_veg_litter_src(site, v) bind(C, name='rhessys_veg_litter_src')
   integer(c_int), intent(in)  :: site
   real(c_double), intent(out) :: v(10)
   v = 0.0_r8
   if (.not. allocated(lit_src)) return
   if (site < 1 .or. site > fates%nsites) return
   v = lit_src(:, site)
end subroutine rhessys_veg_litter_src

!-------------------------------------------------------------------------------
! -demog_litter_direct on (flag /= 0): FATES releases all fine litter to the host
! on the day it falls (local FATES patch, host_litter_passthrough)
!-------------------------------------------------------------------------------
subroutine rhessys_veg_set_litter_direct(flag) bind(C, name='rhessys_veg_set_litter_direct')
   integer(c_int), intent(in) :: flag
   litter_direct = (flag /= 0)
end subroutine rhessys_veg_set_litter_direct

!-------------------------------------------------------------------------------
! -demog_litter_direct: today's litter to the host (g/m2/day) by class (labile,
! cellulose, lignin), source (leaf fines, seeds, AG CWD, fine-root fines, BG CWD)
! and element (C, N): v(class + 3 (source - 1) + 15 (element - 1))
!-------------------------------------------------------------------------------
subroutine rhessys_veg_litter_dcmp(site, v) bind(C, name='rhessys_veg_litter_dcmp')
   integer(c_int), intent(in)  :: site
   real(c_double), intent(out) :: v(30)
   v = 0.0_r8
   if (.not. allocated(lit_dc)) return
   if (site < 1 .or. site > fates%nsites) return
   v = reshape(lit_dc(:, :, :, site), [30])
end subroutine rhessys_veg_litter_dcmp

!-------------------------------------------------------------------------------
! Stream-side sites: share of the site's fresh leaf litter that falls into the
! channel (RHESSys stream_litter_frac); FATES takes it out before its litter pools
! (local FATES patch, host_leaf_export). Called by the host before the day's step.
!-------------------------------------------------------------------------------
subroutine rhessys_veg_set_leaf_export(site, frac) bind(C, name='rhessys_veg_set_leaf_export')
   integer(c_int), intent(in) :: site
   real(c_double), intent(in) :: frac
   if (site < 1 .or. site > fates%nsites) return
   if (.not. allocated(leaf_exp_frac)) then
      if (frac <= 0.0_r8) return
      allocate(leaf_exp_frac(fates%nsites), leaf_exp(4, fates%nsites))
      leaf_exp_frac(:) = 0.0_r8
      leaf_exp(:,:) = 0.0_r8
   end if
   leaf_exp_frac(site) = min(max(frac, 0.0_r8), 1.0_r8)
end subroutine rhessys_veg_set_leaf_export

!-------------------------------------------------------------------------------
! Today's leaf litter taken out of FATES for the stream (g/m2/day): C, N, labile C, labile N
!-------------------------------------------------------------------------------
subroutine rhessys_veg_leaf_export(site, v) bind(C, name='rhessys_veg_leaf_export')
   integer(c_int), intent(in)  :: site
   real(c_double), intent(out) :: v(4)
   v = 0.0_r8
   if (.not. allocated(leaf_exp)) return
   if (site < 1 .or. site > fates%nsites) return
   v = leaf_exp(:, site)
end subroutine rhessys_veg_leaf_export

!-------------------------------------------------------------------------------
! FATES CWD fragmentation split into cellulose / lignin (fates_frag_cwd_fcel/flig)
!-------------------------------------------------------------------------------
subroutine rhessys_veg_cwd_frac(v) bind(C, name='rhessys_veg_cwd_frac')
   real(c_double), intent(out) :: v(2)
   v(1) = ED_val_cwd_fcel
   v(2) = ED_val_cwd_flig
end subroutine rhessys_veg_cwd_frac

!-------------------------------------------------------------------------------
! Live leaf C of site s (kgC/m2), WMFire canopy fuel
!-------------------------------------------------------------------------------
subroutine rhessys_veg_leaf_c(site, leafc) bind(C, name='rhessys_veg_leaf_c')
   integer(c_int), intent(in)  :: site
   real(c_double), intent(out) :: leafc
   type(fates_patch_type),  pointer :: cpatch
   type(fates_cohort_type), pointer :: ccohort
   leafc = 0.0_r8
   if (site < 1 .or. site > fates%nsites) return
   cpatch => fates%sites(site)%oldest_patch
   do while (associated(cpatch))
      ccohort => cpatch%tallest
      do while (associated(ccohort))
         leafc = leafc + ccohort%n * ccohort%prt%GetState(leaf_organ, carbon12_element)
         ccohort => ccohort%shorter
      end do
      cpatch => cpatch%younger
   end do
   leafc = leafc * area_inv
end subroutine rhessys_veg_leaf_c

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
   if (.not. allocated(lit_src)) then
      !$omp critical (lit_src_alloc)
      if (.not. allocated(lit_src)) then
         allocate(lit_src(10, fates%nsites)); lit_src = 0.0_r8
      end if
      !$omp end critical (lit_src_alloc)
   end if
   lit_src(:, s) = 0.0_r8
   if (litter_direct .and. .not. allocated(lit_dc)) then
      !$omp critical (lit_dc_alloc)
      if (.not. allocated(lit_dc)) then
         allocate(lit_dc(3, 5, 2, fates%nsites)); lit_dc = 0.0_r8
      end if
      !$omp end critical (lit_dc_alloc)
   end if
   if (allocated(lit_dc)) lit_dc(:, :, :, s) = 0.0_r8
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
         if (allocated(lit_dc)) then          ! -demog_litter_direct: by source and class
            lit_dc(:, 1, el, s) = lit_dc(:, 1, el, s) + litt%leaf_fines_frag(:) * af * 1000.0_r8
            do ipft = 1, numpft
               sd = (litt%seed_decay(ipft) + litt%seed_germ_decay(ipft)) * af * 1000.0_r8
               lit_dc(1, 2, el, s) = lit_dc(1, 2, el, s) + sd * EDPftvarcon_inst%lf_flab(ipft)
               lit_dc(2, 2, el, s) = lit_dc(2, 2, el, s) + sd * EDPftvarcon_inst%lf_fcel(ipft)
               lit_dc(3, 2, el, s) = lit_dc(3, 2, el, s) + sd * EDPftvarcon_inst%lf_flig(ipft)
            end do
            lit_dc(2, 3, el, s) = lit_dc(2, 3, el, s) + sum(litt%ag_cwd_frag(:)) * ED_val_cwd_fcel * af * 1000.0_r8
            lit_dc(3, 3, el, s) = lit_dc(3, 3, el, s) + sum(litt%ag_cwd_frag(:)) * ED_val_cwd_flig * af * 1000.0_r8
            do j = 1, nlev_eff
               lit_dc(:, 4, el, s) = lit_dc(:, 4, el, s) + litt%root_fines_frag(:, j) * af * 1000.0_r8
               lit_dc(2, 5, el, s) = lit_dc(2, 5, el, s) + sum(litt%bg_cwd_frag(:, j)) * ED_val_cwd_fcel * af * 1000.0_r8
               lit_dc(3, 5, el, s) = lit_dc(3, 5, el, s) + sum(litt%bg_cwd_frag(:, j)) * ED_val_cwd_flig * af * 1000.0_r8
            end do
         end if
         if (el == 1) then
            lit_src(1, s) = lit_src(1, s) + sum(litt%leaf_fines_frag(:)) * af * 1000.0_r8
            lit_src(2, s) = lit_src(2, s) + (sum(litt%seed_decay(:)) + sum(litt%seed_germ_decay(:))) * af * 1000.0_r8
            lit_src(3, s) = lit_src(3, s) + sum(litt%ag_cwd_frag(:)) * af * 1000.0_r8
            lit_src(4, s) = lit_src(4, s) + sum(litt%root_fines_frag(:, 1:nlev_eff)) * af * 1000.0_r8
            lit_src(5, s) = lit_src(5, s) + sum(litt%bg_cwd_frag(:, 1:nlev_eff)) * af * 1000.0_r8
            out%leaf_litter_c = out%leaf_litter_c + sum(litt%leaf_fines(:)) * af
            out%ag_cwd_c      = out%ag_cwd_c + sum(litt%ag_cwd(:)) * af
         else
            lit_src(6, s) = lit_src(6, s) + sum(litt%leaf_fines_frag(:)) * af * 1000.0_r8
            lit_src(7, s) = lit_src(7, s) + (sum(litt%seed_decay(:)) + sum(litt%seed_germ_decay(:))) * af * 1000.0_r8
            lit_src(8, s) = lit_src(8, s) + sum(litt%ag_cwd_frag(:)) * af * 1000.0_r8
            lit_src(9, s) = lit_src(9, s) + sum(litt%root_fines_frag(:, 1:nlev_eff)) * af * 1000.0_r8
            lit_src(10, s) = lit_src(10, s) + sum(litt%bg_cwd_frag(:, 1:nlev_eff)) * af * 1000.0_r8
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
   character(len=16) :: ev
   integer :: u, ivar, kind, n, stat
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
   call get_environment_variable('RHESSYS_FATES_LITTER_OUT', ev, status=stat)
   if (stat == 0 .and. len_trim(ev) > 0 .and. trim(ev) /= '0') call litter_write(trim(fname) // '.litter')
   flush(fates_log_unit)
end subroutine rhessys_veg_restart_write

!-------------------------------------------------------------------------------
! FATES litter carry-over for soil / litter equilibration spin-ups (2026-11-05).
! The inventory init ships no litter or CWD (FatesInventoryInitMod), so an
! inventory start has an empty FATES litter layer.
!  RHESSYS_FATES_LITTER_OUT=1: every restart write also writes <restart>.litter,
!    the area-weighted site-mean FATES litter (kg/m2) per element: leaf fines
!    (ndcmpy classes), above-ground CWD (ncwd), then per soil layer fine-root
!    fines (ndcmpy) and below-ground CWD (ncwd); one line per site and element,
!    keyed by RHESSys patch ID.
!  RHESSYS_FATES_LITTER_INIT=<file>: on a cold start, every FATES patch of each
!    listed site gets those litter / CWD pools (replacing the init's zeros); seed
!    pools are not touched. FATES's mass-balance baselines (old_stock, iflux
!    litter) are reset afterwards, so the added mass is a start state, not a flux.
! Unset: no change (bit-identical).
!-------------------------------------------------------------------------------
subroutine litter_write(fname)
   character(len=*), intent(in) :: fname
   type(fates_patch_type), pointer :: cpatch
   type(litter_type),      pointer :: litt
   real(r8), allocatable :: lf(:), ag(:), rf(:,:), bg(:,:)
   real(r8) :: af
   integer  :: u, s, el, ne, nl
   ne = 1
   if (cnp_coupled) ne = 2
   nl = fates%bc_in(1)%nlevsoil
   allocate(lf(ndcmpy), ag(ncwd), rf(ndcmpy, nl), bg(ncwd, nl))
   open(newunit=u, file=fname, status='replace', action='write')
   write(u, '(a,4i6)') 'RHFLIT1', nl, ncwd, ndcmpy, ne
   do s = 1, fates%nsites
      do el = 1, ne
         lf = 0.0_r8; ag = 0.0_r8; rf = 0.0_r8; bg = 0.0_r8
         cpatch => fates%sites(s)%oldest_patch
         do while (associated(cpatch))
            litt => cpatch%litter(el)
            af = cpatch%area * area_inv
            lf = lf + litt%leaf_fines(1:ndcmpy) * af
            ag = ag + litt%ag_cwd(1:ncwd) * af
            rf = rf + litt%root_fines(1:ndcmpy, 1:nl) * af
            bg = bg + litt%bg_cwd(1:ncwd, 1:nl) * af
            cpatch => cpatch%younger
         end do
         write(u, '(i10,i3,1000es16.8)') fates%sites(s)%h_gid, el, lf, ag, rf, bg
      end do
   end do
   close(u)
   deallocate(lf, ag, rf, bg)
   write(fates_log_unit,*) 'rhessys_fates: wrote FATES litter ', trim(fname)
end subroutine litter_write

subroutine litter_init_env()
   type(fates_patch_type), pointer :: cpatch
   type(litter_type),      pointer :: litt
   character(len=1024) :: fname
   character(len=8)    :: tag
   real(r8), allocatable :: lf(:), ag(:), rf(:,:), bg(:,:)
   real(r8) :: biomass_stock, litter_stock, seed_stock
   integer  :: u, s, el, ne, nl, nlf, ncf, ndf, nef, pid, ios, stat, nset
   logical, allocatable :: touched(:)
   call get_environment_variable('RHESSYS_FATES_LITTER_INIT', fname, status=stat)
   if (stat /= 0 .or. len_trim(fname) == 0) return
   ne = 1
   if (cnp_coupled) ne = 2
   nl = fates%bc_in(1)%nlevsoil
   open(newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
   if (ios /= 0) then
      write(*,*) 'FATAL: cannot open RHESSYS_FATES_LITTER_INIT file ', trim(fname)
      error stop 'litter_init_env'
   end if
   read(u, *) tag, nlf, ncf, ndf, nef
   if (tag /= 'RHFLIT1' .or. nlf /= nl .or. ncf /= ncwd .or. ndf /= ndcmpy) then
      write(*,*) 'FATAL: litter file ', trim(fname), ' does not match this run: layers ', nlf, nl, &
           ' cwd ', ncf, ncwd, ' dcmpy ', ndf, ndcmpy
      error stop 'litter_init_env'
   end if
   allocate(lf(ndcmpy), ag(ncwd), rf(ndcmpy, nl), bg(ncwd, nl), touched(fates%nsites))
   touched = .false.; nset = 0
   do
      read(u, *, iostat=ios) pid, el, lf, ag, rf, bg
      if (ios /= 0) exit
      if (el < 1 .or. el > ne) cycle
      do s = 1, fates%nsites
         if (fates%sites(s)%h_gid /= pid) cycle
         cpatch => fates%sites(s)%oldest_patch
         do while (associated(cpatch))
            litt => cpatch%litter(el)
            litt%leaf_fines(1:ndcmpy)      = lf
            litt%ag_cwd(1:ncwd)            = ag
            litt%root_fines(1:ndcmpy, 1:nl) = rf
            litt%bg_cwd(1:ncwd, 1:nl)      = bg
            cpatch => cpatch%younger
         end do
         touched(s) = .true.; nset = nset + 1
      end do
   end do
   close(u)
   do s = 1, fates%nsites
      if (.not. touched(s)) cycle
      do el = 1, ne
         call SiteMassStock(fates%sites(s), el, fates%sites(s)%mass_balance(el)%old_stock, &
              biomass_stock, litter_stock, seed_stock)
         fates%sites(s)%iflux_balance(el)%iflux_liveveg = (biomass_stock + seed_stock) * area_inv
         fates%sites(s)%iflux_balance(el)%iflux_litter  = litter_stock * area_inv
      end do
   end do
   write(fates_log_unit,*) 'rhessys_fates: FATES litter set from ', trim(fname), ' (', nset, &
        ' site-element lines, ', count(touched), ' sites)'
   flush(fates_log_unit)
   deallocate(lf, ag, rf, bg, touched)
end subroutine litter_init_env

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
      out_all(s)%gc_canopy = 0; out_all(s)%swabs_canopy = 0; out_all(s)%btran = 0; out_all(s)%n_sub = 0
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
   real(r8) :: leafn_w, storen_w
   character(len=1024) :: sname
   real(r8) :: tot, bio_c, lit_c, seed_c, bio_n, lit_n, seed_n

   fname = ''
   do i = 1, min(int(flen), len(fname))
      fname(i:i) = c_file(i)
   end do

   ! Per-site mass stocks (<prefix>_fates_stock.csv, written with the cohorts, kg/m2 of the RHESSys patch):
   ! FATES plant (all organs, above + below ground), FATES litter (leaf/fine-root fines and coarse woody debris
   ! not yet fragmented into the RHESSys litter pools) and seed, for C and (CNP) N. Output only.
   i = index(fname, '_fates_cohorts.csv', back=.true.)
   if (i > 0) then
      sname = fname(1:i-1)//'_fates_stock.csv'
      call SiteMassStock(fates%sites(site), 1, tot, bio_c, lit_c, seed_c)
      bio_n = 0.0_r8; lit_n = 0.0_r8; seed_n = 0.0_r8
      if (cnp_coupled) call SiteMassStock(fates%sites(site), 2, tot, bio_n, lit_n, seed_n)
      if (first /= 0) then
         open(newunit=u, file=trim(sname), status='replace', action='write')
         write(u,'(a)') 'date,patchID,plant_c,litter_c,seed_c,plant_n,litter_n,seed_n'
      else
         open(newunit=u, file=trim(sname), status='old', position='append', action='write')
      end if
      write(u,'(i8,",",i0,6(",",es13.6))') ymd, patch_id, bio_c * area_inv, lit_c * area_inv, &
           seed_c * area_inv, bio_n * area_inv, lit_n * area_inv, seed_n * area_inv
      close(u)
   end if
   if (first /= 0) then
      open(newunit=u, file=trim(fname), status='replace', action='write')
      write(u,'(a)') 'date,patchID,fates_patch,patch_age_yr,patch_area_frac,pft,dbh_cm,height_m,'// &
           'stems_per_ha,canopy_layer,leafc_kg_per_plant,crown_area_m2,treelai_m2_per_m2_crown,'// &
           'storec_kg_per_plant,fnrtc_kg_per_plant,status_coh,site_cstatus,site_gdd,site_nchilldays,'// &
           'site_ncolddays,site_leafondate,site_leafoffdate,site_date,sapwc_kg_per_plant,structc_kg_per_plant,'// &
           'leafn_kg_per_plant,storen_kg_per_plant,vcmax25top'
   else
      open(newunit=u, file=trim(fname), status='old', position='append', action='write')
   end if
   ipa = 0
   cpatch => fates%sites(site)%oldest_patch
   do while (associated(cpatch))
      ipa = ipa + 1
      ccohort => cpatch%tallest
      do while (associated(ccohort))
         ! treelai: FATES's leaf area per m2 of this cohort's crown (cohort LAI on the ground =
         ! treelai x crown_area_m2 / ground area)
         ! storec/fnrtc: storage and fine-root C per plant; status_coh 2 leaves on, 1 off, 3 shedding;
         ! site_*: cold-deciduous phenology state (cstatus 2 not cold, 1 cold; dates are model days);
         ! sapwc/structc: sapwood and structural (dead wood) C per plant, above + below ground
         ! leafn/storen: leaf and storage N per plant (0 when not N-coupled); vcmax25top: the cohort's
         ! canopy-top Vcmax25 (umol/m2/s), scaled by leaf N status under -demog_nvcmax
         if (cnp_coupled) then
            leafn_w = ccohort%prt%GetState(leaf_organ, nitrogen_element)
            storen_w = ccohort%prt%GetState(store_organ, nitrogen_element)
         else
            leafn_w = 0.0_r8; storen_w = 0.0_r8
         end if
         write(u,'(i8,",",i0,",",i0,",",f9.3,",",f8.5,",",i0,",",f8.3,",",f8.3,",",es12.5,",",i0,",",es12.5,",",f10.3,",",es12.5,'// &
               '",",es12.5,",",es12.5,",",i0,",",i0,",",f9.2,",",i0,",",i0,",",i0,",",i0,",",i0,",",es12.5,",",es12.5,'// &
               '",",es12.5,",",es12.5,",",f8.3)') &
              ymd, patch_id, ipa, cpatch%age, cpatch%area * area_inv, ccohort%pft, ccohort%dbh, &
              ccohort%height, ccohort%n * 1.0e4_r8 / cpatch%area * (cpatch%area * area_inv), &
              ccohort%canopy_layer, ccohort%prt%GetState(leaf_organ, carbon12_element), ccohort%c_area, &
              ccohort%treelai, ccohort%prt%GetState(store_organ, carbon12_element), &
              ccohort%prt%GetState(fnrt_organ, carbon12_element), ccohort%status_coh, &
              fates%sites(site)%cstatus, fates%sites(site)%grow_deg_days, fates%sites(site)%nchilldays, &
              fates%sites(site)%ncolddays, fates%sites(site)%cleafondate, fates%sites(site)%cleafoffdate, &
              fates%sites(site)%phen_model_date, ccohort%prt%GetState(sapw_organ, carbon12_element), &
              ccohort%prt%GetState(struct_organ, carbon12_element), leafn_w, storen_w, ccohort%vcmax25top
         ccohort => ccohort%shorter
      end do
      cpatch => cpatch%younger
   end do
   close(u)
end subroutine rhessys_veg_dump_cohorts

!-------------------------------------------------------------------------------
subroutine rhessys_veg_clean() bind(C, name=_VEG_CLEAN)
   ! TODO(Phase 5): deallocate site structures (FATES has no fates_clean yet).
   if (allocated(alloc_acc)) call alloc_write(size(fates%sites))
   if (allocated(mort_acc)) call mort_write(size(fates%sites))
   if (cnp_coupled .and. fates_log_open) call mgmt_summary(fates_log_unit)
   if (fates_log_open) then
      write(fates_log_unit,*) 'rhessys_fates: set_bc calls ', n_setbc, ', bc_in check failures ', n_bad_total
      close(fates_log_unit)
   end if
   fates_is_init = .false.
   fates_log_open = .false.
end subroutine rhessys_veg_clean

end module rhessys_fates
