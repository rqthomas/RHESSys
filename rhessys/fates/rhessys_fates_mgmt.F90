!===============================================================================
! rhessys_fates_mgmt.F90 -- forest management on FATES cohorts for RHESSys
! (-demog mode 3; plan: RHESSys_development/FATES_management_module_plan.md).
!
! Re-implements operations of the Vegetation Management (VM) module for FATES
! by Joshua Rady (github.com/JoshuaRady/fates, branch VegetationManagement /
! VM_Next; FATES <= sci.1.53) on FATES sci.1.93.2_api.46.0.0, targeted by
! RHESSys patch / zone / hillslope instead of latitude / longitude.
!
! Removals (harvest, thinning, understory control) go through FATES's own
! logging machinery: the host sets per-cohort logging rates through the hook
! added by fates/patches/rhessys_host_harvest_hook.patch (EDLoggingMortalityMod
! host_harvest_override). Canopy removals are disturbance generating (FATES
! splits off the removed crown area as a new secondary patch, wood products,
! slash, understory damage); understory removals are not (cohort density
! falls in place). Planting builds new cohorts as FATES recruitment does and
! books their mass as a generic input of the FATES mass balance.
!
! Driver file (one event per line, '!' starts a comment, first non-comment
! line is a header and is skipped):
!   <date> <target> <event>(<arg>=<value>, ...)
! date   : YYYY-MM-DD (separators - / .), '*' matches any digit (VM v1.1)
! target : all | patch:<id> | patches:[id,id,...] | patches:@<file> |
!          zone:<id> | zones:[...] | hill:<id> | hills:[...]
! events : harvest_trees (alias harvest_timber), harvest_biomass, thin_proportional,
!          thin_low_perfect (alias thin_low), thin_low_probabilistic, thin_row_low,
!          understory_control, hardwood_control, plant (alias plant_site), replant
!          -- arguments: see parse_event(). Any removal event takes the VM v1.1
!          conditions if_age_min / if_age_max (FATES patch age, yr), if_bai (min
!          basal area of the eligible trees, m2/ha) and if_dbh_max (their mean dbh
!          must be below it, cm); '<event>_if' names are accepted as aliases.
! Several events on one site and day combine (removal fractions compound).
!===============================================================================
module rhessys_fates_mgmt

   use iso_c_binding
   use FatesConstantsMod,      only : r8 => fates_r8, itrue, ifalse, nearzero
   use FatesConstantsMod,      only : ihard_season_decid, ihard_stress_decid, isemi_stress_decid
   use FatesConstantsMod,      only : leaves_on, leaves_off, ievergreen
   use EDTypesMod,             only : ed_site_type, area_inv, init_recruit_trim
   use EDTypesMod,             only : phen_cstat_nevercold, phen_cstat_iscold
   use FatesPatchMod,          only : fates_patch_type
   use FatesCohortMod,         only : fates_cohort_type
   use FatesInterfaceTypesMod, only : bc_in_type, numpft, nleafage, hlm_parteh_mode
   use FatesInterfaceTypesMod, only : hlm_current_year, hlm_current_month, hlm_current_day
   use EDLoggingMortalityMod,  only : host_harvest_override, host_harvest_active, &
                                      host_coll_under_frac, host_logging_trigger, host_agb_harvest
   use EDCohortDynamicsMod,    only : create_cohort, InitPRTObject
   use EDPftvarcon,            only : EDPftvarcon_inst
   use PRTParametersMod,       only : prt_params
   use PRTGenericMod,          only : prt_vartypes, SetState, StorageNutrientTarget
   use PRTGenericMod,          only : num_elements, element_list
   use PRTGenericMod,          only : carbon12_element, nitrogen_element, phosphorus_element
   use PRTGenericMod,          only : carbon_only, carbon_nitrogen_phosphorus
   use PRTGenericMod,          only : leaf_organ, fnrt_organ, sapw_organ, store_organ, &
                                      repro_organ, struct_organ
   use FatesAllometryMod,      only : h2d_allom, h_allom, bleaf, bfineroot, bsap_allom, &
                                      bagw_allom, bbgw_allom, bdead_allom, bstore_allom

   implicit none
   private

   public :: mgmt_init, mgmt_pre, mgmt_post, mgmt_write_day, mgmt_is_on, mgmt_summary

   integer, parameter :: op_harvest = 1, op_thin_prop = 2, op_thin_low = 3, &
                         op_understory = 4, op_plant = 5, op_thin_prob = 6, op_thin_row = 7, &
                         op_hardwood = 8, op_replant = 9
   character(len=22), parameter :: op_name(9) = [character(len=22) :: 'harvest_trees', &
        'thin_proportional', 'thin_low_perfect', 'understory_control', 'plant', &
        'thin_low_probabilistic', 'thin_row_low', 'hardwood_control', 'replant']
   integer, parameter :: prof_bole = 1, prof_in_place = 2, prof_agb = 3
   ! VM thin_low_probabilistic model (loblolly pine thinning trials, Rady): removal
   ! probability 1/(1 + exp(-k (dbh - mean dbh - m))), k = steepness, m solved for the goal
   real(r8), parameter :: vm_prob_steepness = -0.4428_r8
   ! VM replant: patches recently disturbed by people and bare of trees
   real(r8), parameter :: vm_replant_age_max = 5.0_r8
   integer, parameter :: where_all = 1, where_youngest = 2
   integer, parameter :: maxpft_mg = 16
   real(r8), parameter :: big = 1.0e30_r8
   real(r8), parameter :: pi_mg = 3.14159265358979323846_r8

   type :: mg_event_t
      character(len=16) :: pat(3) = ''              ! year, month, day patterns
      integer  :: op = 0
      logical  :: all_sites = .false.
      logical, allocatable :: site_mask(:)
      logical  :: pft_mask(maxpft_mg) = .true.
      real(r8) :: fraction = -1.0_r8, final_ba = -1.0_r8, final_sd = -1.0_r8
      real(r8) :: dbh_min = 0.0_r8, dbh_max = big, ht_min = 0.0_r8, ht_max = big
      integer  :: profile = prof_bole
      real(r8) :: under_damage = -1.0_r8              ! < 0: FATES logging_coll_under_frac
      real(r8) :: efficiency = 1.0_r8
      real(r8) :: density(maxpft_mg) = 0.0_r8         ! plant: stems/ha per PFT
      real(r8) :: p_dbh = -1.0_r8, p_height = -1.0_r8
      integer  :: where = where_all
      real(r8) :: row_fraction = -1.0_r8                ! thin_row_low
      real(r8) :: steepness = vm_prob_steepness         ! thin_low_probabilistic
      real(r8) :: if_age_min = -1.0_r8, if_age_max = -1.0_r8   ! conditions (patch age, yr)
      real(r8) :: if_bai = -1.0_r8, if_dbh_max = -1.0_r8
      real(r8) :: bare_dbh = 2.0_r8                     ! replant: "bare" = no tree >= this dbh (cm)
      logical  :: pfts_given = .false.
      character(len=200) :: text = ''
   end type mg_event_t

   type :: cptr_t
      type(fates_cohort_type), pointer :: p => null()
   end type cptr_t

   type(mg_event_t), allocatable, save :: ev(:)
   integer, save :: nev = 0
   logical, save :: mg_on = .false.
   integer, save :: nsites_mg = 0
   integer, allocatable, save :: site_pid(:)
   integer, save :: logu = -1
   character(len=1024), save :: diag_file = ''
   logical, save :: diag_first = .true.

   ! per-site, per-day diagnostics (each site written only by the thread running it)
   logical,  allocatable, save :: d_on(:)
   character(len=96), allocatable, save :: d_ops(:)
   real(r8), allocatable, save :: d_ba0(:), d_ba1(:), d_sd0(:), d_sd1(:), d_c0(:), d_c1(:)
   real(r8), allocatable, save :: d_prod_c(:), d_prod_n(:), d_plant_n(:), d_plant_c(:), d_plant_stems(:)
   ! run totals (site means, written at the end)
   integer,  save :: n_site_events = 0
   real(r8), save :: tot_prod_c = 0.0_r8, tot_removed_c = 0.0_r8, tot_plant_c = 0.0_r8

   ! per-site selection for the hook: cohort address -> logging rates (threadprivate:
   ! one site is processed start-to-end by one thread)
   integer, parameter :: maxsel = 4000
   type(c_ptr), save :: sel_ptr(maxsel)
   real(r8),    save :: sel_dir(maxsel), sel_col(maxsel)
   integer,     save :: nsel = 0
   logical,     save :: sel_agb = .false.          ! a whole-AGB harvest selected today
   !$omp threadprivate(sel_ptr, sel_dir, sel_col, nsel, sel_agb)

contains

!-------------------------------------------------------------------------------
logical function mgmt_is_on()
   mgmt_is_on = mg_on
end function mgmt_is_on

!-------------------------------------------------------------------------------
! Read the driver file. pid/zid/hid: RHESSys patch, zone and hillslope IDs of
! each FATES site (site order). Always installs the hook (inactive until used).
!-------------------------------------------------------------------------------
subroutine mgmt_init(file, nsites, pid, zid, hid, log_unit, diagfile)
   character(len=*), intent(in) :: file, diagfile
   integer,          intent(in) :: nsites, pid(:), zid(:), hid(:), log_unit
   integer :: u, ios, n, k, nbad
   character(len=4096) :: line
   character(len=64)   :: tok_date
   character(len=1024) :: tok_target
   character(len=4000) :: tok_event
   logical :: header_seen
   type(mg_event_t), allocatable :: tmp(:)

   host_harvest_override => mgmt_hook
   logu = log_unit
   nsites_mg = nsites
   allocate(site_pid(nsites)); site_pid = pid(1:nsites)
   allocate(d_on(nsites), d_ops(nsites), d_ba0(nsites), d_ba1(nsites), d_sd0(nsites), &
            d_sd1(nsites), d_c0(nsites), d_c1(nsites), d_prod_c(nsites), d_prod_n(nsites), &
            d_plant_n(nsites), d_plant_c(nsites), d_plant_stems(nsites))
   d_on = .false.
   diag_file = diagfile
   if (len_trim(file) == 0) return

   open(newunit=u, file=trim(file), status='old', action='read', iostat=ios)
   if (ios /= 0) call mg_fail('cannot open management driver file '//trim(file))
   allocate(ev(16)); nev = 0; header_seen = .false.
   do
      read(u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      k = index(line, '!')
      if (k > 0) line = line(1:k-1)
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (.not. header_seen) then          ! first non-comment line = header
         header_seen = .true.
         cycle
      end if
      call split3(line, tok_date, tok_target, tok_event)
      if (nev == size(ev)) then
         allocate(tmp(2*size(ev))); tmp(1:nev) = ev(1:nev)
         call move_alloc(tmp, ev)
      end if
      nev = nev + 1
      ev(nev)%text = trim(tok_date)//' '//trim(tok_target)//' '//trim(tok_event)
      call parse_date(tok_date, ev(nev))
      allocate(ev(nev)%site_mask(nsites)); ev(nev)%site_mask = .false.
      call parse_target(tok_target, nsites, pid, zid, hid, ev(nev), nbad)
      if (nbad > 0) write(logu,*) 'rhessys_fates_mgmt: WARNING ', nbad, &
           ' target IDs match no vegetated RHESSys patch in: ', trim(ev(nev)%text)
      call parse_event(tok_event, ev(nev))
   end do
   close(u)
   mg_on = (nev > 0)
   write(logu,*) 'rhessys_fates_mgmt: ', nev, ' management events from ', trim(file)
   do n = 1, nev
      write(logu,*) '  event ', n, ': ', trim(ev(n)%text), '  -> sites ', &
           merge(nsites, count(ev(n)%site_mask), ev(n)%all_sites)
   end do
   flush(logu)
end subroutine mgmt_init

!-------------------------------------------------------------------------------
! Before ed_ecosystem_dynamics for site s: find today's events, build the
! per-cohort removal selection, arm the FATES hook. tec_kind/tec_thin: the
! Phase 6 TEC resolution (3 = partial removal -> proportional thinning by leaf
! habit, [1] evergreen, [2] deciduous).
!-------------------------------------------------------------------------------
subroutine mgmt_pre(site, s, tec_kind, tec_thin)
   type(ed_site_type), intent(inout), target :: site
   integer,  intent(in) :: s, tec_kind
   real(r8), intent(in) :: tec_thin(2)
   type(mg_event_t) :: tev
   integer :: n, ft
   logical :: any_rem, any_ev
   real(r8) :: cu

   nsel = 0
   sel_agb = .false.
   host_harvest_active = .false.
   host_coll_under_frac = -1.0_r8
   host_agb_harvest = .false.
   d_on(s) = .false.
   any_rem = .false.; any_ev = .false.; cu = -1.0_r8
   d_ops(s) = ''

   ! Phase 6 TEC partial removal -> thin_proportional per leaf habit
   if (tec_kind == 3) then
      do ft = 1, min(numpft, maxpft_mg)
         tev%op = op_thin_prop
         tev%pft_mask = .false.
         tev%pft_mask(ft) = .true.
         if (prt_params%phen_leaf_habit(ft) == ievergreen) then
            tev%fraction = tec_thin(1)
         else
            tev%fraction = tec_thin(2)
         end if
         if (tev%fraction > 0.0_r8) then
            if (.not. any_ev) call stand_state(site, s, .true.)
            any_ev = .true.; any_rem = .true.
            call select_removal(site, tev)
         end if
      end do
      if (any_ev) d_ops(s) = 'tec_thin'
   end if

   if (mg_on) then
      do n = 1, nev
         if (.not. event_today(ev(n), s)) cycle
         if (.not. any_ev) call stand_state(site, s, .true.)
         any_ev = .true.
         if (len_trim(d_ops(s)) + len_trim(op_name(ev(n)%op)) < 94) then
            if (len_trim(d_ops(s)) > 0) then
               d_ops(s) = trim(d_ops(s))//' '//trim(op_name(ev(n)%op))
            else
               d_ops(s) = trim(op_name(ev(n)%op))
            end if
         end if
         if (ev(n)%op == op_plant .or. ev(n)%op == op_replant) cycle   ! after dynamics
         if (host_logging_trigger) then                     ! TEC clear-cut today wins
            write(logu,*) 'rhessys_fates_mgmt: site ', s, ' TEC harvest today, skipping ', &
                 trim(ev(n)%text)
            cycle
         end if
         any_rem = .true.
         if (ev(n)%under_damage >= 0.0_r8) cu = max(cu, ev(n)%under_damage)
         call select_removal(site, ev(n))
      end do
   end if

   if (any_ev) d_on(s) = .true.
   if (any_rem .and. .not. host_logging_trigger) then
      host_harvest_active  = .true.
      host_logging_trigger = .true.          ! FATES logging_time for this site today
      host_coll_under_frac = cu
      host_agb_harvest     = sel_agb         ! whole-AGB export for directly logged trees
   end if
end subroutine mgmt_pre

!-------------------------------------------------------------------------------
! After ed_ecosystem_dynamics (before ed_update_site): planting; disarm hook;
! diagnostics.
!-------------------------------------------------------------------------------
subroutine mgmt_post(site, s, bc_in)
   type(ed_site_type), intent(inout), target :: site
   integer,          intent(in) :: s
   type(bc_in_type), intent(in) :: bc_in
   integer :: n

   host_harvest_active = .false.
   host_coll_under_frac = -1.0_r8
   host_agb_harvest = .false.
   sel_agb = .false.
   nsel = 0
   if (.not. d_on(s)) return
   d_prod_c(s) = sum(site%mass_balance(1)%wood_product_harvest(:)) * area_inv
   d_prod_n(s) = 0.0_r8
   if (num_elements >= 2) d_prod_n(s) = sum(site%mass_balance(2)%wood_product_harvest(:)) * area_inv
   call stand_state(site, s, .false.)                   ! after removals, before planting
   d_plant_n(s) = 0.0_r8; d_plant_c(s) = 0.0_r8; d_plant_stems(s) = 0.0_r8
   if (mg_on) then
      do n = 1, nev
         if (ev(n)%op /= op_plant .and. ev(n)%op /= op_replant) cycle
         if (.not. event_today(ev(n), s)) cycle
         call plant(site, s, bc_in, ev(n))
      end do
   end if
end subroutine mgmt_post

!-------------------------------------------------------------------------------
! Serial, after all sites: append today's managed sites to <prefix>_fates_mgmt.csv
!-------------------------------------------------------------------------------
subroutine mgmt_write_day()
   integer :: s, u
   if (.not. allocated(d_on)) return
   if (.not. any(d_on)) return
   if (len_trim(diag_file) == 0) return
   if (diag_first) then
      open(newunit=u, file=trim(diag_file), status='replace', action='write')
      write(u,'(a)') 'date,patchID,events,stems_per_ha_before,stems_per_ha_after,'// &
           'basal_area_m2_ha_before,basal_area_m2_ha_after,live_c_kg_m2_before,live_c_kg_m2_after,'// &
           'wood_product_c_kg_m2,wood_product_n_kg_m2,planted_stems_per_ha,planted_c_kg_m2,planted_n_kg_m2'
      diag_first = .false.
   else
      open(newunit=u, file=trim(diag_file), status='old', position='append', action='write')
   end if
   do s = 1, nsites_mg
      if (.not. d_on(s)) cycle
      write(u,'(i4.4,i2.2,i2.2,",",i0,",",a,11(",",es13.6))') hlm_current_year, hlm_current_month, &
           hlm_current_day, site_pid(s), trim(d_ops(s)), d_sd0(s), d_sd1(s), d_ba0(s), d_ba1(s), &
           d_c0(s), d_c1(s), d_prod_c(s), d_prod_n(s), d_plant_stems(s), d_plant_c(s), d_plant_n(s)
      n_site_events = n_site_events + 1
      tot_prod_c    = tot_prod_c + d_prod_c(s)
      tot_removed_c = tot_removed_c + (d_c0(s) - d_c1(s))
      tot_plant_c   = tot_plant_c + d_plant_c(s)
      d_on(s) = .false.
   end do
   close(u)
end subroutine mgmt_write_day

subroutine mgmt_summary(u)
   integer, intent(in) :: u
   if (.not. allocated(d_on)) return
   write(u,*) 'rhessys_fates_mgmt: managed site-days ', n_site_events, &
        ' sum over site-days (kgC/m2): live C removed ', tot_removed_c, ' wood products ', &
        tot_prod_c, ' planted ', tot_plant_c
end subroutine mgmt_summary

!===============================================================================
! The hook FATES calls (EDLoggingMortalityMod host_harvest_override) for every
! cohort of a site with host_harvest_active: today's rates for this cohort.
!===============================================================================
subroutine mgmt_hook(site, cohort, lmort_direct, lmort_collateral, lmort_infra, l_degrad)
   type(ed_site_type),      intent(inout), target :: site
   type(fates_cohort_type), intent(inout), target :: cohort
   real(r8),                intent(inout) :: lmort_direct, lmort_collateral
   real(r8),                intent(inout) :: lmort_infra, l_degrad
   integer :: k
   type(c_ptr) :: p
   lmort_direct = 0.0_r8; lmort_collateral = 0.0_r8
   lmort_infra  = 0.0_r8; l_degrad = 0.0_r8
   p = c_loc(cohort)
   do k = 1, nsel
      if (c_associated(sel_ptr(k), p)) then
         lmort_direct     = sel_dir(k)
         lmort_collateral = sel_col(k)
         return
      end if
   end do
end subroutine mgmt_hook

!===============================================================================
! Selection
!===============================================================================
subroutine add_sel(c, f, profile)
   type(fates_cohort_type), pointer, intent(in) :: c
   real(r8), intent(in) :: f
   integer,  intent(in) :: profile
   integer :: k
   real(r8) :: ff, keep
   type(c_ptr) :: p
   ff = min(max(f, 0.0_r8), 1.0_r8)
   if (ff <= 0.0_r8) return
   p = c_loc(c)
   do k = 1, nsel
      if (c_associated(sel_ptr(k), p)) exit
   end do
   if (k > nsel) then
      if (nsel == maxsel) call mg_fail('selection list full (maxsel)')
      nsel = nsel + 1; k = nsel
      sel_ptr(k) = p; sel_dir(k) = 0.0_r8; sel_col(k) = 0.0_r8
   end if
   keep = 1.0_r8 - sel_dir(k) - sel_col(k)               ! survivors so far
   if (profile == prof_in_place) then
      sel_col(k) = sel_col(k) + keep * ff
   else
      sel_dir(k) = sel_dir(k) + keep * ff
      if (profile == prof_agb) sel_agb = .true.
   end if
end subroutine add_sel

logical function tree_ok(c, e, use_dbh_max, use_ht)
   type(fates_cohort_type), pointer, intent(in) :: c
   type(mg_event_t), intent(in) :: e
   logical, intent(in) :: use_dbh_max, use_ht
   tree_ok = .false.
   if (prt_params%woody(c%pft) /= itrue) return
   if (c%pft > maxpft_mg) return
   if (.not. e%pft_mask(c%pft)) return
   if (c%dbh < e%dbh_min) return
   if (use_dbh_max .and. c%dbh >= e%dbh_max) return
   if (use_ht .and. (c%height < e%ht_min .or. c%height >= e%ht_max)) return
   tree_ok = .true.
end function tree_ok

subroutine select_removal(site, e)
   type(ed_site_type), intent(inout), target :: site
   type(mg_event_t), intent(in) :: e
   type(fates_patch_type),  pointer :: cp
   type(fates_cohort_type), pointer :: c
   type(cptr_t), allocatable :: cl(:)
   real(r8), allocatable :: f(:)
   real(r8) :: ba, sd, goal, frow, mean_dbh
   integer  :: nc, i, gkind

   cp => site%oldest_patch
   do while (associated(cp))
      if (.not. patch_ok(cp, e)) then          ! VM v1.1 conditions (per FATES patch)
         cp => cp%younger
         cycle
      end if
      select case (e%op)

      case (op_harvest)
         c => cp%tallest
         do while (associated(c))
            if (tree_ok(c, e, .true., .true.)) then
               if (e%fraction >= 0.0_r8) then
                  call add_sel(c, e%fraction, e%profile)
               else
                  call add_sel(c, 1.0_r8, e%profile)
               end if
            end if
            c => c%shorter
         end do

      case (op_understory)
         c => cp%tallest
         do while (associated(c))
            if (c%canopy_layer > 1 .and. c%pft <= maxpft_mg) then
               if (e%pft_mask(c%pft) .and. c%dbh < e%dbh_max) &
                    call add_sel(c, e%efficiency, prof_in_place)
            end if
            c => c%shorter
         end do

      case (op_hardwood)
         ! VM hardwood_control: kill 'hardwood' PFTs (all sizes, all layers) in place.
         ! Default PFTs: woody, not evergreen (VM used broadleaf PFT indices of its
         ! 12-PFT file); pfts= overrides.
         c => cp%tallest
         do while (associated(c))
            if (prt_params%woody(c%pft) == itrue .and. c%pft <= maxpft_mg) then
               if (hard_pft(c%pft, e)) call add_sel(c, e%efficiency, prof_in_place)
            end if
            c => c%shorter
         end do

      case (op_thin_prop, op_thin_low, op_thin_prob, op_thin_row)
         call eligible(cp, e, cl, nc, ba, sd, mean_dbh)
         if (nc > 0) then
            allocate(f(nc)); f = 0.0_r8
            ! goal: 1 = basal area (m2/ha), 2 = stems (/ha)
            if (e%final_ba >= 0.0_r8) then
               gkind = 1; goal = e%final_ba
            else if (e%final_sd >= 0.0_r8) then
               gkind = 2; goal = e%final_sd
            else
               gkind = 2; goal = sd * (1.0_r8 - max(min(e%fraction, 1.0_r8), 0.0_r8))
            end if
            select case (e%op)
            case (op_thin_prop)
               if (gkind == 1 .and. ba > 0.0_r8) then
                  f = max(1.0_r8 - goal / ba, 0.0_r8)
               else if (sd > 0.0_r8) then
                  f = max(1.0_r8 - goal / sd, 0.0_r8)
               end if
            case (op_thin_low)
               call thin_low_fracs(cl, nc, cp%area, gkind, goal, f)
            case (op_thin_prob)
               call thin_prob_fracs(cl, nc, cp%area, gkind, goal, mean_dbh, e%steepness, f)
            case (op_thin_row)
               ! VM thin_row_low: only if above the goal; remove every 1/X row (same
               ! fraction of every eligible tree), then thin from below to the goal
               if ((gkind == 1 .and. ba > goal) .or. (gkind == 2 .and. sd > goal)) then
                  frow = max(min(e%row_fraction, 1.0_r8), 0.0_r8)
                  call thin_low_fracs(cl, nc, cp%area, gkind, goal, f, keep=1.0_r8 - frow)
                  f = frow + (1.0_r8 - frow) * f
               end if
            end select
            do i = 1, nc
               call add_sel(cl(i)%p, f(i), e%profile)
            end do
            deallocate(f)
         end if
         if (allocated(cl)) deallocate(cl)
      end select
      cp => cp%younger
   end do
end subroutine select_removal

! eligible trees of a FATES patch (pfts, dbh >= dbh_min), sorted by dbh (ascending);
! basal area (m2/ha), stems (/ha) and stem-weighted mean dbh (cm) of them
subroutine eligible(cp, e, cl, nc, ba, sd, mean_dbh)
   type(fates_patch_type), pointer, intent(in) :: cp
   type(mg_event_t), intent(in) :: e
   type(cptr_t), allocatable, intent(out) :: cl(:)
   integer,  intent(out) :: nc
   real(r8), intent(out) :: ba, sd, mean_dbh
   type(fates_cohort_type), pointer :: c, tp
   integer :: i, j
   nc = 0; ba = 0.0_r8; sd = 0.0_r8; mean_dbh = 0.0_r8
   c => cp%tallest
   do while (associated(c))
      if (tree_ok(c, e, .false., .false.)) nc = nc + 1
      c => c%shorter
   end do
   allocate(cl(max(nc, 1)))
   if (nc == 0 .or. cp%area <= 0.0_r8) then
      nc = 0
      return
   end if
   i = 0
   c => cp%tallest
   do while (associated(c))
      if (tree_ok(c, e, .false., .false.)) then
         i = i + 1; cl(i)%p => c
         ba = ba + c%n * pi_mg * (c%dbh / 200.0_r8)**2
         sd = sd + c%n
         mean_dbh = mean_dbh + c%n * c%dbh
      end if
      c => c%shorter
   end do
   if (sd > 0.0_r8) mean_dbh = mean_dbh / sd
   ba = ba / cp%area * 1.0e4_r8; sd = sd / cp%area * 1.0e4_r8
   do i = 2, nc                                    ! insertion sort by dbh
      j = i
      do while (j > 1)
         if (cl(j-1)%p%dbh <= cl(j)%p%dbh) exit
         tp => cl(j-1)%p; cl(j-1)%p => cl(j)%p; cl(j)%p => tp
         j = j - 1
      end do
   end do
end subroutine eligible

! 'perfect' low thinning: remove the smallest trees first until the stand (after
! an earlier removal leaving the fraction keep of every tree) reaches the goal
subroutine thin_low_fracs(cl, nc, area, gkind, goal, f, keep)
   type(cptr_t), intent(in) :: cl(:)
   integer,  intent(in) :: nc, gkind
   real(r8), intent(in) :: area, goal
   real(r8), intent(out) :: f(:)
   real(r8), intent(in), optional :: keep
   real(r8) :: kp, tot, rem, w, contrib
   integer :: i
   kp = 1.0_r8
   if (present(keep)) kp = keep
   f = 0.0_r8
   tot = 0.0_r8
   do i = 1, nc
      tot = tot + kp * wgt(cl(i)%p, area, gkind)
   end do
   rem = max(tot - goal, 0.0_r8)
   do i = 1, nc
      if (rem <= 0.0_r8) exit
      contrib = kp * wgt(cl(i)%p, area, gkind)
      if (contrib <= 0.0_r8) cycle
      w = min(1.0_r8, rem / contrib)
      f(i) = w
      rem = rem - w * contrib
   end do
end subroutine thin_low_fracs

! VM thin_low_probabilistic: removal probability p = 1/(1 + exp(-k (dbh - mean - m)))
! (k < 0: smaller trees more likely); the midpoint m is solved by bisection so the
! remaining basal area / stems equal the goal (VM: linear+binary search, 0.1 trees/ha)
subroutine thin_prob_fracs(cl, nc, area, gkind, goal, mean_dbh, k, f)
   type(cptr_t), intent(in) :: cl(:)
   integer,  intent(in) :: nc, gkind
   real(r8), intent(in) :: area, goal, mean_dbh, k
   real(r8), intent(out) :: f(:)
   real(r8) :: lo, hi, m, st
   integer :: it
   f = 0.0_r8
   if (state(1.0e3_r8) > goal) then              ! even m -> +inf barely reaches it
      m = 1.0e3_r8
   else if (state(-1.0e3_r8) <= goal) then       ! nothing to remove
      return
   else
      lo = -1.0e3_r8; hi = 1.0e3_r8
      do it = 1, 200                              ! remaining state decreases with m
         m = 0.5_r8 * (lo + hi)
         st = state(m)
         if (st > goal) then
            lo = m
         else
            hi = m
         end if
         if (hi - lo < 1.0e-10_r8) exit
      end do
      m = 0.5_r8 * (lo + hi)
   end if
   do it = 1, nc
      f(it) = prob(cl(it)%p%dbh, m)
   end do
contains
   real(r8) function prob(d, mm)
      real(r8), intent(in) :: d, mm
      prob = 1.0_r8 / (1.0_r8 + exp(min(max(-k * (d - mean_dbh - mm), -500.0_r8), 500.0_r8)))
   end function prob
   real(r8) function state(mm)
      real(r8), intent(in) :: mm
      integer :: ii
      state = 0.0_r8
      do ii = 1, nc
         state = state + (1.0_r8 - prob(cl(ii)%p%dbh, mm)) * wgt(cl(ii)%p, area, gkind)
      end do
   end function state
end subroutine thin_prob_fracs

! a cohort's basal area (m2/ha, gkind 1) or stems (/ha, gkind 2) per ha of its FATES patch
real(r8) function wgt(c, area, gkind)
   type(fates_cohort_type), pointer, intent(in) :: c
   real(r8), intent(in) :: area
   integer,  intent(in) :: gkind
   wgt = c%n / area * 1.0e4_r8
   if (gkind == 1) wgt = wgt * pi_mg * (c%dbh / 200.0_r8)**2
end function wgt

logical function hard_pft(ft, e)
   integer, intent(in) :: ft
   type(mg_event_t), intent(in) :: e
   if (e%pfts_given) then
      hard_pft = e%pft_mask(ft)
   else
      hard_pft = (prt_params%phen_leaf_habit(ft) /= ievergreen)
   end if
end function hard_pft

! VM v1.1 conditions on a FATES patch for an event (all optional, < 0 = unset)
logical function patch_ok(cp, e)
   type(fates_patch_type), pointer, intent(in) :: cp
   type(mg_event_t), intent(in) :: e
   type(cptr_t), allocatable :: cl(:)
   integer  :: nc
   real(r8) :: ba, sd, md
   patch_ok = .true.
   if (e%if_age_min >= 0.0_r8 .and. cp%age < e%if_age_min) patch_ok = .false.
   if (e%if_age_max >= 0.0_r8 .and. cp%age > e%if_age_max) patch_ok = .false.
   if (.not. patch_ok) return
   if (e%if_bai >= 0.0_r8 .or. e%if_dbh_max >= 0.0_r8) then
      call eligible(cp, e, cl, nc, ba, sd, md)
      if (e%if_bai >= 0.0_r8 .and. ba < e%if_bai) patch_ok = .false.
      if (e%if_dbh_max >= 0.0_r8 .and. nc > 0 .and. md >= e%if_dbh_max) patch_ok = .false.
      if (allocated(cl)) deallocate(cl)
   end if
end function patch_ok

! VM replant: patch recently disturbed by people (age since anthropogenic disturbance
! in [if_age_min (default 0), if_age_max (default 5 yr)]) and bare of trees (no woody
! cohort with dbh >= bare_dbh; VM: no tree cohort at all, which FATES's own recruits
! would always defeat)
logical function replant_ok(cp, e)
   type(fates_patch_type), pointer, intent(in) :: cp
   type(mg_event_t), intent(in) :: e
   type(fates_cohort_type), pointer :: c
   real(r8) :: amin, amax
   amin = max(e%if_age_min, 0.0_r8)
   amax = vm_replant_age_max
   if (e%if_age_max >= 0.0_r8) amax = e%if_age_max
   replant_ok = (cp%age_since_anthro_disturbance >= amin .and. &
                 cp%age_since_anthro_disturbance <= amax)
   if (.not. replant_ok) return
   c => cp%tallest
   do while (associated(c))
      if (prt_params%woody(c%pft) == itrue .and. c%dbh >= e%bare_dbh) then
         replant_ok = .false.
         return
      end if
      c => c%shorter
   end do
end function replant_ok

!-------------------------------------------------------------------------------
! Stand state of site s (trees, all FATES patches, per ha / per m2 of site)
!-------------------------------------------------------------------------------
subroutine stand_state(site, s, before)
   type(ed_site_type), intent(inout), target :: site
   integer, intent(in) :: s
   logical, intent(in) :: before
   type(fates_patch_type),  pointer :: cp
   type(fates_cohort_type), pointer :: c
   real(r8) :: ba, sd, cl
   ba = 0.0_r8; sd = 0.0_r8; cl = 0.0_r8
   cp => site%oldest_patch
   do while (associated(cp))
      c => cp%tallest
      do while (associated(c))
         cl = cl + c%n * (c%prt%GetState(leaf_organ, carbon12_element) + &
              c%prt%GetState(fnrt_organ, carbon12_element) + c%prt%GetState(sapw_organ, carbon12_element) + &
              c%prt%GetState(store_organ, carbon12_element) + c%prt%GetState(struct_organ, carbon12_element) + &
              c%prt%GetState(repro_organ, carbon12_element))
         if (prt_params%woody(c%pft) == itrue) then
            ba = ba + c%n * pi_mg * (c%dbh / 200.0_r8)**2
            sd = sd + c%n
         end if
         c => c%shorter
      end do
      cp => cp%younger
   end do
   ! site area = 1/area_inv m2
   if (before) then
      d_ba0(s) = ba * area_inv * 1.0e4_r8; d_sd0(s) = sd * area_inv * 1.0e4_r8; d_c0(s) = cl * area_inv
      d_ba1(s) = d_ba0(s); d_sd1(s) = d_sd0(s); d_c1(s) = d_c0(s)
      d_prod_c(s) = 0.0_r8; d_prod_n(s) = 0.0_r8
      d_plant_n(s) = 0.0_r8; d_plant_c(s) = 0.0_r8; d_plant_stems(s) = 0.0_r8
   else
      d_ba1(s) = ba * area_inv * 1.0e4_r8; d_sd1(s) = sd * area_inv * 1.0e4_r8; d_c1(s) = cl * area_inv
   end if
end subroutine stand_state

!===============================================================================
! Planting (VM plant / plant_site): new cohorts built as FATES recruitment
! builds them (EDPhysiologyMod recruitment), at the requested size and
! density, mass booked as flux_generic_in of the FATES mass balance.
!===============================================================================
subroutine plant(site, s, bc_in, e)
   type(ed_site_type), intent(inout), target :: site
   integer,          intent(in) :: s
   type(bc_in_type), intent(in) :: bc_in
   type(mg_event_t), intent(in) :: e
   type(fates_patch_type), pointer :: cp, young
   class(prt_vartypes), pointer :: prt
   integer  :: ft, el, element_id, leaf_status, iage
   real(r8) :: height, dbh, nn, l2fr, efleaf, effnrt, efstem
   real(r8) :: c_leaf, c_fnrt, c_sapw, a_sapw, c_agw, c_bgw, c_struct, c_store
   real(r8) :: m_leaf, m_fnrt, m_sapw, m_struct, m_store, m_tot

   young => null()
   cp => site%oldest_patch
   do while (associated(cp))
      if (.not. associated(young)) then
         young => cp
      else if (cp%age < young%age) then
         young => cp
      end if
      cp => cp%younger
   end do

   do ft = 1, min(numpft, maxpft_mg)
      if (e%density(ft) <= 0.0_r8) cycle
      if (e%p_dbh > 0.0_r8) then
         dbh = e%p_dbh
         call h_allom(dbh, ft, height)
      else
         height = EDPftvarcon_inst%hgt_min(ft)
         if (e%p_height > 0.0_r8) height = e%p_height
         call h2d_allom(height, ft, dbh)
      end if
      efleaf = 1.0_r8; effnrt = 1.0_r8; efstem = 1.0_r8; leaf_status = leaves_on
      select case (prt_params%phen_leaf_habit(ft))
      case (ihard_season_decid)
         if (site%cstatus == phen_cstat_nevercold .or. site%cstatus == phen_cstat_iscold) then
            efleaf = 0.0_r8
            effnrt = 1.0_r8 - prt_params%phen_fnrt_drop_fraction(ft)
            efstem = 1.0_r8 - prt_params%phen_stem_drop_fraction(ft)
            leaf_status = leaves_off
         end if
      case (ihard_stress_decid, isemi_stress_decid)
         efleaf = site%elong_factor(ft)
         effnrt = 1.0_r8 - (1.0_r8 - efleaf) * prt_params%phen_fnrt_drop_fraction(ft)
         efstem = 1.0_r8 - (1.0_r8 - efleaf) * prt_params%phen_stem_drop_fraction(ft)
         if (efleaf <= 0.0_r8) leaf_status = leaves_off
      end select

      cp => site%oldest_patch
      do while (associated(cp))
         if (e%op == op_replant) then
            if (.not. replant_ok(cp, e)) then
               cp => cp%younger
               cycle
            end if
         else
            if (e%where == where_youngest .and. .not. associated(cp, young)) then
               cp => cp%younger
               cycle
            end if
            if (.not. patch_ok(cp, e)) then
               cp => cp%younger
               cycle
            end if
         end if
         nn = e%density(ft) * 1.0e-4_r8 * cp%area                 ! plants in this FATES patch
         if (nn > 0.0_r8) then
            l2fr = site%rec_l2fr(ft, cp%NCL_p)
            call bleaf(dbh, ft, 1, init_recruit_trim, efleaf, c_leaf)
            call bfineroot(dbh, ft, init_recruit_trim, l2fr, effnrt, c_fnrt)
            call bsap_allom(dbh, ft, 1, init_recruit_trim, efstem, a_sapw, c_sapw)
            call bagw_allom(dbh, ft, 1, efstem, c_agw)
            call bbgw_allom(dbh, ft, efstem, c_bgw)
            call bdead_allom(c_agw, c_bgw, c_sapw, ft, c_struct)
            call bstore_allom(dbh, ft, 1, init_recruit_trim, c_store)
            prt => null()
            call InitPRTObject(prt)
            do el = 1, num_elements
               element_id = element_list(el)
               select case (element_id)
               case (carbon12_element)
                  m_struct = c_struct; m_leaf = c_leaf; m_fnrt = c_fnrt
                  m_sapw = c_sapw; m_store = c_store
               case (nitrogen_element)
                  m_struct = c_struct * prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(struct_organ))
                  m_leaf   = c_leaf   * prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(leaf_organ))
                  m_fnrt   = c_fnrt   * prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(fnrt_organ))
                  m_sapw   = c_sapw   * prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(sapw_organ))
                  m_store  = StorageNutrientTarget(ft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)
               case (phosphorus_element)
                  m_struct = c_struct * prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(struct_organ))
                  m_leaf   = c_leaf   * prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(leaf_organ))
                  m_fnrt   = c_fnrt   * prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(fnrt_organ))
                  m_sapw   = c_sapw   * prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(sapw_organ))
                  m_store  = StorageNutrientTarget(ft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)
               end select
               select case (hlm_parteh_mode)
               case (carbon_only, carbon_nitrogen_phosphorus)
                  call SetState(prt, leaf_organ, element_id, m_leaf, 1)
                  do iage = 2, nleafage
                     call SetState(prt, leaf_organ, element_id, 0.0_r8, iage)
                  end do
                  call SetState(prt, fnrt_organ,   element_id, m_fnrt)
                  call SetState(prt, sapw_organ,   element_id, m_sapw)
                  call SetState(prt, store_organ,  element_id, m_store)
                  call SetState(prt, struct_organ, element_id, m_struct)
                  call SetState(prt, repro_organ,  element_id, 0.0_r8)
               case default
                  call mg_fail('planting: unsupported PARTEH mode')
               end select
               m_tot = nn * (m_struct + m_leaf + m_fnrt + m_sapw + m_store)
               site%mass_balance(el)%flux_generic_in = site%mass_balance(el)%flux_generic_in + m_tot
               if (element_id == carbon12_element) d_plant_c(s) = d_plant_c(s) + m_tot * area_inv
               if (element_id == nitrogen_element) d_plant_n(s) = d_plant_n(s) + m_tot * area_inv
            end do
            call prt%CheckInitialConditions()
            call create_cohort(site, cp, ft, nn, height, 0.0_r8, dbh, prt, efleaf, effnrt, &
                 efstem, leaf_status, 1, init_recruit_trim, 0.0_r8, cp%NCL_p, 1, site%spread, bc_in)
            d_plant_stems(s) = d_plant_stems(s) + nn * area_inv * 1.0e4_r8
         end if
         cp => cp%younger
      end do
   end do
end subroutine plant

!===============================================================================
! Driver-file parsing
!===============================================================================
subroutine split3(line, a, b, rest)
   character(len=*), intent(in)  :: line
   character(len=*), intent(out) :: a, b, rest
   character(len=len(line)) :: t
   integer :: k
   t = adjustl(line)
   k = index(t, ' ')
   a = t(1:k-1); t = adjustl(t(k:))
   k = index(t, ' ')
   b = t(1:k-1); rest = adjustl(t(k:))
   if (len_trim(a) == 0 .or. len_trim(b) == 0 .or. len_trim(rest) == 0) &
        call mg_fail('need <date> <target> <event> in: '//trim(line))
end subroutine split3

subroutine parse_date(tok, e)
   character(len=*), intent(in) :: tok
   type(mg_event_t), intent(inout) :: e
   integer :: i, f, st
   f = 1; st = 1
   do i = 1, len_trim(tok) + 1
      if (i > len_trim(tok)) then
         e%pat(f) = tok(st:i-1)
      else if (tok(i:i) == '-' .or. tok(i:i) == '/' .or. tok(i:i) == '.') then
         if (f == 3) call mg_fail('bad date '//trim(tok))
         e%pat(f) = tok(st:i-1); f = f + 1; st = i + 1
      end if
   end do
   if (f /= 3) call mg_fail('bad date (need Y-M-D) '//trim(tok))
   do f = 1, 3
      if (len_trim(e%pat(f)) == 0 .or. verify(trim(e%pat(f)), '0123456789*') /= 0) &
           call mg_fail('bad date '//trim(tok))
   end do
end subroutine parse_date

logical function field_match(pat, val)
   character(len=*), intent(in) :: pat
   integer, intent(in) :: val
   character(len=16) :: v, w
   integer :: n, i, iv, ios
   n = len_trim(pat)
   if (index(pat(1:n), '*') == 0) then
      read(pat(1:n), *, iostat=ios) iv
      field_match = (ios == 0 .and. iv == val)
      return
   end if
   write(w, '(i16.16)') val
   v = w(17-n:16)                                  ! zero-padded to the pattern width
   field_match = .true.
   do i = 1, n
      if (pat(i:i) /= '*' .and. pat(i:i) /= v(i:i)) field_match = .false.
   end do
   if (val >= 10**min(n, 9)) field_match = .false.  ! value wider than the pattern
end function field_match

logical function event_today(e, s)
   type(mg_event_t), intent(in) :: e
   integer, intent(in) :: s
   event_today = .false.
   if (.not. (e%all_sites .or. e%site_mask(s))) return
   if (.not. field_match(e%pat(1), hlm_current_year)) return
   if (.not. field_match(e%pat(2), hlm_current_month)) return
   if (.not. field_match(e%pat(3), hlm_current_day)) return
   event_today = .true.
end function event_today

subroutine parse_target(tok, nsites, pid, zid, hid, e, nbad)
   character(len=*), intent(in) :: tok
   integer, intent(in) :: nsites, pid(:), zid(:), hid(:)
   type(mg_event_t), intent(inout) :: e
   integer, intent(out) :: nbad
   character(len=16) :: kind
   character(len=len(tok)) :: val
   integer, allocatable :: ids(:)
   integer :: k, i, nid
   nbad = 0
   if (trim(tok) == 'all') then
      e%all_sites = .true.
      return
   end if
   k = index(tok, ':')
   if (k == 0) call mg_fail('bad target '//trim(tok))
   kind = tok(1:k-1); val = tok(k+1:)
   if (val(1:1) == '@') then
      call read_id_file(trim(val(2:)), ids, nid)
   else
      call parse_int_list(val, ids, nid)
   end if
   do i = 1, nid
      select case (trim(kind))
      case ('patch', 'patches')
         where (pid(1:nsites) == ids(i)) e%site_mask = .true.
         if (.not. any(pid(1:nsites) == ids(i))) nbad = nbad + 1
      case ('zone', 'zones')
         where (zid(1:nsites) == ids(i)) e%site_mask = .true.
         if (.not. any(zid(1:nsites) == ids(i))) nbad = nbad + 1
      case ('hill', 'hills')
         where (hid(1:nsites) == ids(i)) e%site_mask = .true.
         if (.not. any(hid(1:nsites) == ids(i))) nbad = nbad + 1
      case default
         call mg_fail('bad target kind '//trim(kind))
      end select
   end do
end subroutine parse_target

subroutine read_id_file(fn, ids, nid)
   character(len=*), intent(in) :: fn
   integer, allocatable, intent(out) :: ids(:)
   integer, intent(out) :: nid
   integer :: u, ios, v
   integer, allocatable :: tmp(:)
   open(newunit=u, file=fn, status='old', action='read', iostat=ios)
   if (ios /= 0) call mg_fail('cannot open target ID file '//trim(fn))
   allocate(ids(64)); nid = 0
   do
      read(u, *, iostat=ios) v
      if (ios /= 0) exit
      if (nid == size(ids)) then
         allocate(tmp(2*size(ids))); tmp(1:nid) = ids(1:nid); call move_alloc(tmp, ids)
      end if
      nid = nid + 1; ids(nid) = v
   end do
   close(u)
end subroutine read_id_file

subroutine parse_int_list(val, ids, nid)
   character(len=*), intent(in) :: val
   integer, allocatable, intent(out) :: ids(:)
   integer, intent(out) :: nid
   real(r8), allocatable :: r(:)
   integer :: i
   call parse_real_list(val, r)
   nid = size(r)
   allocate(ids(nid))
   do i = 1, nid
      ids(i) = nint(r(i))
   end do
end subroutine parse_int_list

! "3.5" or "[1, 2, 3]" -> real array
subroutine parse_real_list(val, r)
   character(len=*), intent(in) :: val
   real(r8), allocatable, intent(out) :: r(:)
   character(len=len(val)) :: t
   integer :: i, n, st, ios, k
   t = adjustl(val)
   k = len_trim(t)
   if (t(1:1) == '[') then
      if (t(k:k) /= ']') call mg_fail('unclosed list '//trim(val))
      t = t(2:k-1)
   end if
   n = 1
   do i = 1, len_trim(t)
      if (t(i:i) == ',') n = n + 1
   end do
   allocate(r(n))
   st = 1; k = 0
   do i = 1, len_trim(t) + 1
      if (i > len_trim(t) .or. t(min(i, len(t)):min(i, len(t))) == ',') then
         k = k + 1
         read(t(st:i-1), *, iostat=ios) r(k)
         if (ios /= 0) call mg_fail('bad number in '//trim(val))
         st = i + 1
      end if
   end do
end subroutine parse_real_list

! name(key=value, key=[a,b], key=word)
subroutine parse_event(spec, e)
   character(len=*), intent(in) :: spec
   type(mg_event_t), intent(inout) :: e
   character(len=64) :: name, key
   character(len=len(spec)) :: args, item, val
   integer :: k1, k2, i, st, depth, ft, nl
   real(r8), allocatable :: r(:)
   logical :: got_pfts, got_density

   k1 = index(spec, '('); k2 = index(spec, ')', back=.true.)
   if (k1 < 2 .or. k2 < k1) call mg_fail('bad event '//trim(spec))
   name = adjustl(spec(1:k1-1))
   nl = len_trim(name)
   if (nl > 3) then                                   ! VM v1.1 '<event>_if' aliases
      if (name(nl-2:nl) == '_if') name = name(1:nl-3)
   end if
   select case (trim(name))
   case ('harvest_trees', 'harvest_timber')
      e%op = op_harvest
   case ('harvest_biomass')
      e%op = op_harvest; e%profile = prof_agb
   case ('thin_proportional')
      e%op = op_thin_prop
   case ('thin_low_perfect', 'thin_low')
      e%op = op_thin_low
   case ('thin_low_probabilistic')
      e%op = op_thin_prob
   case ('thin_row_low')
      e%op = op_thin_row
   case ('understory_control')
      e%op = op_understory
   case ('hardwood_control')
      e%op = op_hardwood
   case ('plant', 'plant_site')
      e%op = op_plant
   case ('replant')
      e%op = op_replant
   case default
      call mg_fail('unknown management event '//trim(name))
   end select
   args = spec(k1+1:k2-1)
   got_pfts = .false.; got_density = .false.
   st = 1; depth = 0
   do i = 1, len_trim(args) + 1
      if (i <= len_trim(args)) then
         if (args(i:i) == '[') depth = depth + 1
         if (args(i:i) == ']') depth = depth - 1
         if (.not. (args(i:i) == ',' .and. depth == 0)) cycle
      end if
      item = adjustl(args(st:i-1)); st = i + 1
      if (len_trim(item) == 0) cycle
      k1 = index(item, '=')
      if (k1 < 2) call mg_fail('argument needs name=value: '//trim(item)//' in '//trim(spec))
      key = adjustl(item(1:k1-1)); val = adjustl(item(k1+1:))
      select case (trim(key))
      case ('pfts', 'pft')
         call parse_real_list(val, r)
         e%pft_mask = .false.
         do ft = 1, size(r)
            if (nint(r(ft)) < 1 .or. nint(r(ft)) > maxpft_mg) call mg_fail('bad pft in '//trim(spec))
            e%pft_mask(nint(r(ft))) = .true.
         end do
         got_pfts = .true.
         e%pfts_given = .true.
      case ('thin_fraction', 'harvest_fraction', 'fraction')
         call one(val, e%fraction)
      case ('final_basal_area')
         call one(val, e%final_ba)
      case ('final_stem_density')
         call one(val, e%final_sd)
      case ('dbh_min')
         call one(val, e%dbh_min)
      case ('dbh_max')
         call one(val, e%dbh_max)
      case ('ht_min')
         call one(val, e%ht_min)
      case ('ht_max')
         call one(val, e%ht_max)
      case ('efficiency')
         call one(val, e%efficiency)
      case ('understory_damage')
         call one(val, e%under_damage)
      case ('row_fraction')
         call one(val, e%row_fraction)
      case ('steepness')
         call one(val, e%steepness)
      case ('if_age_min', 'delay')
         call one(val, e%if_age_min)
      case ('if_age_max')
         call one(val, e%if_age_max)
      case ('if_bai')
         call one(val, e%if_bai)
      case ('if_dbh_max')
         call one(val, e%if_dbh_max)
      case ('bare_dbh')
         call one(val, e%bare_dbh)
      case ('method')
         continue                                         ! VM: mow / herbicide are equivalent
      case ('flux_profile', 'profile')
         select case (trim(val))
         case ('bole_harvest', 'bole', 'harvest')
            e%profile = prof_bole
         case ('in_place', 'kill')
            e%profile = prof_in_place
         case ('agb_harvest', 'biomass', 'agb')
            e%profile = prof_agb
         case default
            call mg_fail('unknown flux_profile '//trim(val))
         end select
      case ('density')
         call parse_real_list(val, r)
         got_density = .true.
         if (size(r) == 1) then
            e%density = r(1)                               ! applied to selected pfts below
         else
            e%density = 0.0_r8
            e%density(1:min(size(r), maxpft_mg)) = r(1:min(size(r), maxpft_mg))
            e%density = -e%density                         ! by position in pfts list: fix below
         end if
      case ('dbh')
         call one(val, e%p_dbh)
      case ('height')
         call one(val, e%p_height)
      case ('where')
         select case (trim(val))
         case ('all', 'everywhere')
            e%where = where_all
         case ('youngest', 'youngest_patch')
            e%where = where_youngest
         case default
            call mg_fail('unknown where '//trim(val))
         end select
      case default
         call mg_fail('unknown argument '//trim(key)//' in '//trim(spec))
      end select
   end do

   if (e%op == op_plant .or. e%op == op_replant) then
      if (.not. got_density) call mg_fail('plant needs density= in '//trim(spec))
      if (.not. got_pfts) call mg_fail('plant needs pfts= in '//trim(spec))
      call fix_plant_density(spec, e)
   end if
   if ((e%op == op_thin_prop .or. e%op == op_thin_low .or. e%op == op_thin_prob) .and. &
        e%fraction < 0.0_r8 .and. e%final_ba < 0.0_r8 .and. e%final_sd < 0.0_r8) &
        call mg_fail('thinning needs thin_fraction, final_basal_area or final_stem_density: '//trim(spec))
   ! herbicide-type control kills only its targets: no understory damage under the
   ! killed crowns unless asked for (FATES's logging gap would kill everything there)
   if (e%op == op_hardwood .and. e%under_damage < 0.0_r8) e%under_damage = 0.0_r8
   if (e%op == op_thin_row) then
      if (e%row_fraction <= 0.0_r8 .or. e%row_fraction > 1.0_r8) &
           call mg_fail('thin_row_low needs row_fraction in (0,1]: '//trim(spec))
      if (e%final_ba < 0.0_r8 .and. e%final_sd < 0.0_r8) &
           call mg_fail('thin_row_low needs final_basal_area or final_stem_density: '//trim(spec))
   end if
contains
   subroutine one(v, x)
      character(len=*), intent(in) :: v
      real(r8), intent(out) :: x
      real(r8), allocatable :: rr(:)
      call parse_real_list(v, rr)
      x = rr(1)
   end subroutine one
end subroutine parse_event

! plant density: one value -> every selected PFT; a list -> by position in pfts=[...]
subroutine fix_plant_density(spec, e)
   character(len=*), intent(in) :: spec
   type(mg_event_t), intent(inout) :: e
   real(r8) :: d(maxpft_mg)
   integer :: ft, k
   if (all(e%density >= 0.0_r8)) then
      do ft = 1, maxpft_mg
         if (.not. e%pft_mask(ft)) e%density(ft) = 0.0_r8
      end do
      return
   end if
   d = -e%density; e%density = 0.0_r8; k = 0
   do ft = 1, maxpft_mg
      if (e%pft_mask(ft)) then
         k = k + 1
         if (k <= maxpft_mg) e%density(ft) = d(k)
      end if
   end do
   if (count(d > 0.0_r8) > count(e%pft_mask)) call mg_fail('more densities than pfts in '//trim(spec))
end subroutine fix_plant_density

subroutine mg_fail(msg)
   character(len=*), intent(in) :: msg
   if (logu /= -1) then
      write(logu,*) 'rhessys_fates_mgmt: FATAL: ', trim(msg)
      flush(logu)
   end if
   write(0,*) 'rhessys_fates_mgmt: FATAL: ', trim(msg)
   error stop 'rhessys_fates_mgmt'
end subroutine mg_fail

end module rhessys_fates_mgmt
