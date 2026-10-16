# RHESSys extensions: user's guide

This guide covers what has been added to RHESSys in this fork (`github.com/rqthomas/RHESSys`) since upstream commit `90abf784` (2026-01-20). It is written for people who already run RHESSys and want to use the new options.

Status as of 2026-10-14. Development notes, tests and calibration results are in `RHESSys_development/README.md` (outside this repo).

---

## 1. Ground rules

- **Opt-in.** Every new process is off by default, through a command-line flag or a def-file parameter whose default reproduces the old behaviour. With everything off, output and state files are byte-identical to the build before each feature was added. Each feature was checked that way.
- **Bug fixes are not opt-in.** Several upstream bugs were fixed (section 12). These change results compared with upstream RHESSys even when you turn nothing on.
- **Restarts are exact.** A run restarted from a state file reproduces a continuous run bit for bit (section 3). The one exception is the FATES backend, which reproduces it to round-off.

---

## 2. Building

The binary is built in the `ccr-rhessys` docker image (`RHESSys_development/ccr_rhessys_dwh/.devcontainer`):

```sh
cd rhessys
make clobber && make                 # legacy build -> rhessys7.5
```

- Use `make clobber`, not `make clean`. `make clean` also deletes generated parser code that needs flex/bison.
- Build serially. `make -j` has a race.

| Build option | Command | What it adds |
|---|---|---|
| FATES vegetation backend | `make fates=T PGM=rhessys7.5_fates` | `-demog*` options (section 10). Needs `libfates` built from `RHESSys_development/fates` (`cmake -S . -B build`); override `FATES_DIR` to move it. |
| WMFire | `make wmfire=T` (or `make fates=T wmfire=T`) | Fire spread. Build `lib/libwmfire.so` first (see `rhessys/fates/WMFIRE.md`); run with `LD_LIBRARY_PATH=<repo>/lib`. |
| N budget debug | `make CMD_OPTS=-DN_BUDGET` | Daily N mass-balance checks, printed to the log. |
| Stream heat budget debug | `make CMD_OPTS=-DSTREAM_HEAT_BUDGET` | Per-reach heat-closure check for `-strtemp`. |

`rhessys/Testing/bit_identity_test.sh` runs the byte-identity regression test against a blessed reference snapshot.

---

## 3. Exact restarts and editing state files

`output_current_state` (a tec event) writes a state file holding every value carried from one day to the next, at full precision (`%.17g`). This includes:
- phenology timing;
- running means;
- decomposition scalars;
- saved sunlit/shade splits;
- N pools;
- the previous day's soil temperature.

Optional blocks are written after the basins:

| Block | Written when | Contents (one line per reach) |
|---|---|---|
| `stream_routing_state` | `-str` | `reach_ID initial_flow previous_lateral_input water_depth reservoir_storage Qout previous_Qin` |
| `stream_benthic_state` | `-strbgc` | `reach_ID benthic_POC benthic_PON wc_NO3 wc_NH4 wc_DON wc_DOC wc_POC wc_PON` |
| `stream_temperature_state` | `-strtemp` | `reach_ID T_water T_gw_mem T_sub_mem` |

The blocks can come in any order and any of them can be absent.

**Restart convention.** A state written by the event `Y M D 1 output_current_state` is restarted with `-st Y M D 1`.

**Older state files** still load. Missing values fall back to their cold-start values, with a warning for each missing block.

**Editing a state file** (e.g. for data assimilation or a manual harvest). Values are tag-keyed text, so you can edit them directly. For any value, `-9999` means "recompute this from the other pools on read". Use it like this:

- **Vegetation C pools edited:** set `epv.day_livestem_turnover`, `epv.day_livecroot_turnover`, `epv.day_deadleaf_turnover`, `phen.leaflitfallc`, `phen.frootlitfallc` and `phen.leaflitfallc_wstress` to `-9999`.
- **Litter or soil C pool edited:** set its N pool (`litter_ns.litr2n`, `soil_ns.soil1n`, …) to `-9999`. It then resets to the fixed C:N.
- **LAI:** LAI is not a state. Edit `cs.leafc` (and `ns.leafn`).

---

## 4. Command-line options (quick reference)

| Option | Needs | What it does | Section |
|---|---|---|---|
| `-strbgc <file>` | `-str`, `-g` | In-stream storage and processing of POC/PON/DOM/N | 5.3 |
| `-strtemp <file>` | `-str` | Stream water temperature per reach | 5.4 |
| `-spinmode <ndep_mult> <years>` | | Accelerated spin-up: N deposition × `ndep_mult` and no erosion detachment for the first `years` | 9.2 |
| `-demog <params.json> [1\|2\|3]` | `-g`, FATES build | FATES vegetation backend | 10 |
| `-demog_inv <ctrl file>` | `-demog` | FATES initial stand from an inventory | 10 |
| `-demog_restart` | `-demog` | Restart FATES from `<world>.fates` + `<world>.fates_host` | 10 |
| `-demog_mgmt <driver file>` | `-demog … 3` | Forest management on FATES cohorts | 10 |
| `-demog_nneed` | `-demog` | FATES plant N demand from potential growth, not only root capacity | 10 |
| `-demog_wood_resp <f>` | `-demog … 3` | Fraction (0–1) of FATES CWD fragmentation C respired by wood decay before it reaches RHESSys litter | 10 |

`-str` is now recognised as an option by `valid_option`. Before, a value placed after it could be misread as an argument.

---

## 5. Streams

### 5.1 New stream-routing output

`<prefix>_streamrouting.daily` gains columns, in this order:

```
day month year reachID Qout lateralinput Qin waterdepth reservoir.store
[-g]       NO3_out NH4_out DON_out DOC_out
           sediment_out
[-g]       lateral_NO3 lateral_NH4 lateral_DON lateral_DOC
           lateral_sediment
           POC_out PON_out lateral_POC lateral_PON lateral_POC_labile lateral_PON_labile
           lateral_POC_sed lateral_PON_sed
[-strbgc]  23 processing columns (5.3)
[-strtemp] 17 temperature columns (5.4)
```

Units:
- Flows: m³/day. The `*_out` columns include upstream inputs, and the `lateral_*` columns cover only the reach's own lateral-input patches.
- Loads: kg/day.
- `lateral_POC/PON` includes the eroded-soil part, `lateral_POC_sed/PON_sed` (section 6).

**Nutrients and DOM in reaches.** Patch and hillslope-groundwater NO3, NH4, DON and DOC loads are routed through the reach network with the water (with `-g`).

**Groundwater in routing.** Hillslope groundwater discharge (`gw.Qout`) is now routed through the kinematic wave. Before, it appeared in basin streamflow but not in stream routing.

### 5.2 Leaf litter into streams

| Parameter (file) | Default | Meaning |
|---|---|---|
| `stream_litter_overhang` (landuse def) | 0 (off) | Each reach's stream-side (lateral-input) patches send (overhang × channel area / their total area) of the day's leaf litterfall straight into the reach as POC/PON. The litter pools are debited, and the labile (litr1) part is tracked. |

Without `-strbgc` this POC passes through the network on the day it falls.

### 5.3 In-stream storage and processing: `-strbgc <file>`

Each reach is a well-mixed water column with volume V = alfa·Q^0.6·L, plus a streambed store. Daily processes:
- suspended POC/PON settles to the bed;
- bed POC/PON fragments to DOC/DON and is entrained at high flow;
- DOM is decomposed with microbial stoichiometry (CUE, microbial C:N);
- the microbes immobilize N or mineralize it to NH4;
- NO3 is denitrified;
- optionally, attached microbes decompose part of the bed material in place.

Rates scale with fT = θ^(T−20), using the `-strtemp` water temperature if it is on, otherwise max(2.99 + 0.72·Tair, 0). Each process is solved exactly over the day. The `*_out` columns become the processed loads.

| Parameter | Default | Units | Meaning |
|---|---|---|---|
| `stream_theta` | 1.047 | – | temperature coefficient |
| `stream_v_dep` | 10 | m/day | deposition velocity of suspended POC/PON |
| `stream_k_frag` | 0.01 | 1/day at 20 °C | bed POC/PON fragmentation to DOC/DON |
| `stream_e_max` | 0.5 | 1/day | entrainment rate of the bed store at bankfull |
| `stream_export_p` | 2 | – | exponent on Q/Q_bankfull for entrainment |
| `stream_vf_DOM` | 0.1 | m/day | DOM decomposition uptake velocity |
| `stream_CUE` | 0.3 | – | microbial carbon use efficiency |
| `stream_CN_mic` | 10 | kgC/kgN | microbial biomass C:N |
| `stream_vf_N` | 2 | m/day | NH4+NO3 supply velocity for immobilization |
| `stream_vf_denit` | 0.1 | m/day | denitrification uptake velocity |
| `stream_f_mic` | 0 | – | share of bed breakdown done in place by attached microbes |
| `stream_bgc_lake_reaches` | none | file | reach IDs inside a lake/reservoir: loads pass through unprocessed |

Extra output columns:
```
benthic_POC benthic_PON POC_deposit POC_entrain Frag_C Frag_N DOM_dec_C N_mineral N_immob
N_limit stream_CO2 stream_denitrif V_water wc_DOC wc_DON wc_NO3 wc_NH4 wc_POC wc_PON
Bed_mic_C Bed_N_mineral Bed_N_immob Bed_N_limit
```

Bankfull discharge comes from the stream-table geometry (Manning). An example parameter file is `ccr_rhessys_dwh/defs/stream_bgc.def`.

### 5.4 Stream temperature: `-strtemp <file>`

Gives a daily water temperature for every reach. It is meant for lake-model inflows (GLM-AED, ELCOM) and for `-strbgc`. Design notes: `RHESSys_development/stream_temperature_design.md`.

**Sources.** Each source water has its own temperature, mixed by flow:

| Source (RHESSys flux) | Temperature |
|---|---|
| hillslope groundwater (`base_flow`) | running mean of air temperature over `stream_tau_gw` days + `stream_dT_gw` (or a fixed `stream_T_gw_const`) |
| patch subsurface flow | running mean over `stream_tau_sub` days |
| return flow | the day's rain (one third of the way from dew point to air temperature) and snowmelt (0 °C) first; any extra is soil water |
| upstream reaches | their daily mean outflow temperature |

The running means are driven by max(Tair, 0), and by 0 °C on snow-covered patches.

**Heat budget.** Each reach is a well-mixed water column exchanging heat at its surface:
- shortwave, reduced by shade;
- longwave from the sky and from the canopy;
- emitted longwave;
- evaporation (Brady–Graves–Geyer wind function);
- sensible heat (Bowen ratio);
- optional bed exchange toward the groundwater temperature.

The energy balance is linearized and the day is solved exactly, so it is stable for any reach volume. Water does not go below 0 °C; there is no ice model.

**Shading from LAI.** Shading uses the reach's stream-side patches:

```
S_v = min(1, c_overhang · H_canopy / surface width) · (1 − exp(−k · PAI))
S   = 1 − (1 − S_topo)(1 − S_v)
```

- PAI is the area-weighted Σ cover × `epv.proj_pai` (leaf plus stem area) of the strata.
- k is the strata `ext_coef`.
- H_canopy is the strata height.

S cuts shortwave, swaps sky longwave for canopy longwave, and shelters the wind. Daily LAI is used, so leaf-out, harvest, fire or an edited `leafc` change stream temperature the next day.

| Parameter | Default | Meaning |
|---|---|---|
| `stream_albedo_w`, `stream_emiss_w`, `stream_emiss_veg` | 0.08, 0.97, 0.97 | water albedo, water and canopy emissivity |
| `stream_k_shade` | −1 | extinction coefficient for shading; −1 = strata `ext_coef` |
| `stream_SAI` | 0 | plant area (m²/m²) added to `proj_pai`, e.g. to raise leaf-off shade |
| `stream_c_overhang` | 1 | shaded width share = min(1, c·H/w) |
| `stream_c_wind_shelter` | 0.5 | wind over water = U·(1 − c·S_v) |
| `stream_wind_fn_mult` | 1 | multiplier on the evaporation wind function |
| `stream_k_bed` | 0 | bed heat exchange, kJ/m²/day/K (86.4 = 1 W/m²/K) |
| `stream_tau_gw`, `stream_dT_gw`, `stream_T_gw_const` | 365 d, 1 °C, −999 (off) | groundwater temperature |
| `stream_tau_sub` | 60 d | soil-water temperature |
| `stream_T_init`, `stream_T_gw_init` | 10, 12.5 °C | cold-start values (no state block) |
| `stream_V_min` | 1 m³ | smaller reaches use their steady-state temperature |
| `stream_T_snow_insulate` | 1 | 0 °C driver under snow |
| `stream_temp_exchange` | 1 | 0 = mixing only (testing) |
| `stream_temp_lake_reaches` | none | file of reach IDs inside a reservoir (mixing only) |
| `stream_temp_reach_file` | none | per-reach `reach_ID S_topo shade_override c_overhang`; `shade_override` < 0 uses the LAI model, `c_overhang` < 0 uses the default |

Extra output columns:

| Columns | Units | Meaning |
|---|---|---|
| `T_mean T_end` | °C | daily mean outflow temperature; end-of-day temperature (the state) |
| `T_lat T_up` | °C | flow-weighted lateral and upstream inflow temperatures; −9999 = no flow |
| `T_gw T_sub` | °C | source temperatures used that day |
| `shade pai_rip` | 0–1; m²/m² | shade; riparian plant area index |
| `H_sw H_lw H_evap H_sens H_bed H_ice` | kJ/m²/day | surface fluxes, + = into the water |
| `Q_lat_gw Q_lat_sub Q_lat_ret` | m³/day | lateral inflow by pathway |

**Which column to use as a lake inflow temperature:**
- for a reach entering a lake: `T_mean`;
- for hillslope water draining straight into the lake: `T_lat`.

An example parameter file is `ccr_rhessys_dwh/defs/stream_temp.def`. **Not yet calibrated:** leaf-off shade from RHESSys stem area is probably low.

---

## 6. Erosion, sediment and particulate organic matter

Sediment is detached by rain splash on each patch. Overland flow carries it, limited by a transport capacity; what is not carried redeposits. In the streams it passes through conservatively (or settles, with `-strbgc`).

```
rainsplash (kg/m²/day) = 0.029 · rain_mm · soil_erodibility_K · C · [splash_saturated_factor if saturated to the surface]
C  = cover_and_management_C                       if in (0, 1]
   = exp(−cover_canopy_k · LAI) · exp(−cover_litter_a · litter C)   otherwise
Tc = sediment_transport_capacity_c · Qout^1.5 · sin(slope)^sediment_transport_slope_exp · exp(−sediment_transport_litter_b · litter C)
```

All parameters are in the soil def:

| Parameter | Default | Meaning |
|---|---|---|
| `soil_erodibility_K` | 0 (off) | erodibility in the splash equation |
| `cover_and_management_C` | 0 (unset) | fixed USLE-style C factor; unset = from LAI and litter |
| `cover_canopy_k` | 0.5 | steepness of the LAI part of C (0.5 = old hard-coded value) |
| `cover_litter_a` | 0 | m²/kg C, litter-cover part of C (above-ground litr1–4) |
| `splash_saturated_factor` | 1 | splash multiplier on days saturated to the surface (water films absorb drop impact) |
| `sediment_transport_capacity_c` | 0 (off) | overland-flow transport capacity coefficient |
| `sediment_transport_slope_exp` | 0 | slope exponent on capacity |
| `sediment_transport_litter_b` | 0 | litter reduction of capacity |
| `sediment_OC_mixing_mass` | 0 | kg soil/m² (e.g. 130 = 0.1 m × 1300 kg/m³): sediment takes the fraction rainsplash/M of each soil pool, so its organic content follows the soil's |
| `sediment_OC_frac` | 0 | fixed kg C per kg sediment (alternative to the mixing mass) |

- Eroded soil organic C/N travels with the sediment (patch states `surface_sedC/N`) and reaches the streams as `lateral_POC_sed/PON_sed`.
- Sediment not exported by the next detachment day redeposits, and its C/N returns to soil1–4.
- Under the FATES backend, FATES's unfragmented dead leaves count as litter cover.

---

## 7. Soil water, groundwater and nitrogen

All of these are soil def parameters:

| Parameter | Default | Meaning |
|---|---|---|
| `soil_percolation_rate`, `soil_percolation_exp` | 0 (off), 1 | the saturated zone leaks q = rate × (saturated thickness/soil depth)^exp × Ksat_vertical (m/day; at most half the saturated water) to the hillslope groundwater store, carrying saturated-zone nitrate and soil DOM. Gives deep soils an exit at the soil bottom. |
| `gw_recharge_through_soil` | 0 | 1 = surface solutes in recharge water enter the soil first, and recharge carries soil-water NO3 to groundwater |
| `gw_NO3_recharge_frac` | 1 | multiplier on the soil-water NO3 concentration carried to groundwater (with the option above) |
| `gw_DOM_recharge_frac` | 0 | fraction of the soil-water DOC/DON concentration carried with recharge (0 = groundwater gets only surface DOM) |
| `gw_DOC_recharge_max` | −1 (off) | cap on recharge DOC, mg C/L (DON in proportion) |
| `soil_riparian` | 0 | 1 marks the soil as riparian for `-gwtoriparian` (soil ID 42 is still riparian as before) |
| `nitrif_proportion` | 1 | scales the maximum nitrification rate |
| `denitrif_proportion` | 1 | scales denitrification |

**Behaviour changes, not optional:**
- NH4 deposition now goes into soil NH4 (`sminn`), not the surface pool. Before, it could reach groundwater without passing through the soil.
- With `-gwtoriparian`, solutes passed to riparian patches are no longer multiplied by the hourly water depth. That bug made most groundwater N and DOC disappear.

---

## 8. Vegetation physiology and mortality

All of these are stratum (veg) def parameters:

| Parameter | Default | Meaning |
|---|---|---|
| `epc.reserve_max_frac` | 0 (off) | N-limited surplus photosynthate is kept in `cpool` up to this fraction of plant biomass C. The reserve pays maintenance respiration, which prevents carbon-starvation die-offs from short deficits. |
| `epc.starv_mort_scale` | 0 (off) | graded starvation mortality, applied daily: killed share = min(1, (cpool_mort_fract·B − cpool) / (scale·live biomass C)). Replaces the all-or-nothing kill when cpool < 0. |
| `epc.hydr_mort_days0` | 0 (off) | hydraulic-stress mortality starts after this many stress days |
| `epc.hydr_mort_days_span` | 300 | stress days from onset to full rate |
| `epc.hydr_mort_psi` | psi_close | predawn soil water potential threshold (MPa) |
| `epc.hydr_mort_recovery` | 1 | stress days removed per non-stress day |
| `epc.hydr_mort_rate` | 1 | full mortality rate (1/yr) |

**How the stress index works.** It goes up by 1 on each leaf-on, unfrozen day with predawn psi ≤ the threshold, and down by `recovery` on other such days. It is saved in the state file when it is not zero. Daily mortality is rate × min(1, (index − days0)/span) / 365.

Dead plants reach litter and CWD gradually, and any carbon debt is taken from the CWD that the dying plants add.

---

## 9. Climate, radiation, spin-up and disturbance

### 9.1 Climate and radiation

| Item | Where | Notes |
|---|---|---|
| Daily CO2 sequence | base station: list `CO2` as a non-critical daily sequence, file `<prefix>.CO2` (ppm) | Feeds Farquhar photosynthesis and stomatal conductance. |
| `atm_CO2` | zone def | Constant CO2 when there is no sequence. **Upstream bug fixed:** the key was looked up as `"atm_CO2 "` (trailing space), so the def value was never read and every run used 322 ppm. |
| `radiation_slope_aspect` | zone def, default 0 | 1 = scale supplied (flat-surface) Kdown_direct/diffuse by each zone's computed sloped/flat clear-sky ratio (sun geometry, slope, aspect, horizons). Upstream applied flat radiation unchanged to every slope. |
| `-climrepeat` | – | **Bug fixed:** each repeat of the climate record started one day late, so long spin-ups drifted out of phase with the calendar (about 6 weeks after 1000 years). |

### 9.2 Spin-up mode: `-spinmode <ndep_mult> <years>`

For the first `<years>` of the run, NO3 and NH4 deposition are multiplied by `<ndep_mult>` and erosion detachment is off (redeposition still runs). This replaces making separate high-deposition climate files.

### 9.3 New tec event: `redefine_world_thin_fire_remain` (thintyp 4)

This event works like `redefine_world_thin_remain` (thintyp 1), except that killed trees' above-ground non-structural C and N go to CWD with the dead stems instead of to labile litter. Those pools are cpool, npool, leaf and stem stores and transfers, gresp and retransn.

Use it for the "killed, left on site" part of a fire. With thintyp 1, a large `cpool` produced an unrealistic jump in labile litter.

```
2026 8 1 1 redefine_world_thin_harvest       # consumed by fire
2026 8 1 2 redefine_world_thin_fire_remain   # killed, left as CWD
```

---

## 10. FATES vegetation backend (`-demog`)

The real FATES model (cohorts, light competition, demography) can replace RHESSys's plant module. RHESSys keeps hydrology, soil biogeochemistry, routing, radiation and aerodynamics. It needs a FATES build (section 2) and `-g`.

### 10.1 Modes

| Mode | Command | What FATES does |
|---|---|---|
| 1 (default) | `-demog <params.json>` | Read-only: FATES runs alongside; RHESSys results are byte-identical to legacy |
| 2 | `-demog <params.json> 2` | FATES transpiration and canopy structure (LAI, SAI, height, sunlit share) drive RHESSys water and energy; legacy plant C/N still runs |
| 3 | `-demog <params.json> 3` | Full coupling: FATES also takes N from RHESSys soil mineral N and returns litter/CWD to RHESSys litter; legacy plant C/N is off |

### 10.1a When to use mode 2

Mode 2 couples FATES to RHESSys **water and energy only**. Each day:
- FATES sets the canopy that RHESSys's radiation, interception, aerodynamics and evaporation see (LAI split by leaf habit, SAI, height, sunlit share);
- RHESSys scales its transpiration demand to FATES transpiration, then applies its own soil-water supply limits and accounting. The basin water balance closes exactly.

The legacy RHESSys plants keep running underneath, with growth, allocation, phenology, litterfall, mortality, turnover and plant N uptake. They supply all of the soil carbon and nitrogen cycle. FATES runs carbon-only with nutrients supplied, so it is never N-limited, and its litter does not go to RHESSys.

Mode 2 is the right choice for:
- **Attributing water effects to vegetation structure.** You can ask how streamflow, ET and soil moisture change when the canopy is a FATES size- and age-structured stand instead of the legacy strata, while soil BGC and stream chemistry keep coming from the same legacy processes. Comparing mode 2 with a legacy run isolates the hydrologic effect of FATES. Comparing mode 2 with mode 3 isolates the effect of N coupling and FATES litter.
- **Calibrating FATES canopy physiology before turning on N.** You can check FATES GPP, transpiration, LAI seasonality and canopy height against flux, ET or LAI data, and tune the JSON, without N limitation and FATES litter chemistry also changing. This is the step between mode 1 (read-only) and mode 3.
- **Hydrology studies where nutrients are not the question.** Examples: water yield over succession, under drought, or with demographic shifts simulated by FATES.
- **Keeping a calibrated soil BGC and stream chemistry setup.** Mode 3 replaces the litter inputs (different chemistry, an extra fragmentation stage, FATES CWD rates), so soil C/N and stream DOC/NO3 need recalibrating after a FATES spin-up. Mode 2 leaves them on the legacy calibration.

Use mode 3 instead when any of these apply:
- **Carbon, nitrogen or stream chemistry should follow the FATES stand.** In mode 2, litter, soil N, plant N uptake and stream NO3/DOC come from the legacy shadow plants. Their carbon pools are not the FATES stand, so C and N are not consistent with the canopy that hydrology sees.
- **N limitation matters.** FATES GPP in mode 2 is unconstrained by N.
- **You need disturbance or management.**
  - TEC thinning, harvest and fire events (`redefine_world_thin_*`) act on the legacy strata in mode 2. The next day FATES structure overwrites the canopy, so the event has essentially no hydrologic effect.
  - `-demog_mgmt` needs mode 3, and so does the WMFire–FATES coupling.
- **You report plant C.** RHESSys's standard plant-C outputs show the legacy shadow, not FATES. FATES plant C is in `<prefix>_fates_*.csv`, in every mode.

Cost: modes 2 and 3 cost about the same, since nearly all the time is FATES itself. Measured for CCR (1,425 sites, `fates_params_ccr2pft_fast*.json`) on one thread, before sites ran in parallel:
- modes 1–2: about 170 s per simulated year;
- mode 3: about 240 s per simulated year.

With OpenMP (12 threads), mode 3 takes about 45 s per year, against about 5 s for legacy.

### 10.2 Inputs

- **Parameter file (JSON).** For CCR, use `RHESSys_development/fates_inputs/fates_params_ccr2pft_fast.json` (modes 1–2) or `fates_params_ccr2pft_fast_cnp.json` (mode 3). Both have 2 PFTs (needleleaf evergreen, broadleaf cold-deciduous) and coarser cohort fusion, about 3.6× faster.
- **Initial stand:** `-demog_inv <ctrl>` takes ED2-style inventory files.

  The control file is one header line, then `type latitude longitude pss-path css-path`:
  ```
  type latitude longitude pss-name css-name
  1 37.36 0.0 /abs/path/stand.pss /abs/path/stand.css
  ```
  The patch file (`.pss`):
  ```
  time patch trk age area
  1850 p1 1 100.0 1.0
  ```
  The cohort file (`.css`):
  ```
  time patch dbh height pft nplant
  1850 p1 30.0 0.0 2 0.020
  ```
  - Paths are absolute and at most 256 characters.
  - `dbh` is in cm. `height` 0 means it comes from allometry.
  - `nplant` is stems/m².

  Example: `fates_inputs/ccr_static/`.
- **Sites.** Only vegetated patches become FATES sites; water and reservoir patches are skipped.

### 10.3 Restart

At `output_current_state`, RHESSys also writes `<state>.fates` and `<state>.fates_host`. Read them back with `-w <state> -demog_restart`.

FATES restarts match to round-off (about 1e-10 relative), not bit for bit.

### 10.4 Disturbance and management

- **TEC events in mode 3.**
  - Removal ≥ 0.99 alone → clear-cut harvest through FATES logging: stems → wood product (exported), slash → litter/CWD.
  - Removal + kill on the same day → fire with uniform severity.
  - Removal < 0.99 alone → proportional thinning by leaf habit.
- **`-demog_mgmt <driver>`.** Dated thinning, harvest, understory/hardwood control, planting and replanting on chosen patches, zones or hillslopes. Events can carry conditions (stand age, basal area, mean dbh). Output: `<prefix>_fates_mgmt.csv`. Full syntax: `rhessys/fates/MANAGEMENT.md`.
- **WMFire with FATES (`-demog … 3 -firespread <cell size>`).** FATES fuels go into WMFire, and WMFire's height-based mortality is applied per cohort. See `rhessys/fates/WMFIRE.md`. The fire regime is not calibrated for CCR.

### 10.5 Outputs and performance

| File | Contents |
|---|---|
| `<prefix>_fates_basin.csv` | daily area-weighted FATES vs legacy: GPP, transpiration, LAI, btran, canopy conductance, absorbed shortwave (columns vary by phase/mode) |
| `<prefix>_fates_flux.csv` | daily fluxes for one diagnostic patch (choose it with env `RHESSYS_FATES_DIAG_PATCH`) |
| `<prefix>_fates_bc.csv` | soil boundary conditions for the diagnostic patch |
| `<prefix>_fates.log` | FATES log, warnings, run totals |
| `<prefix>_fates_mgmt.csv` | management events (with `-demog_mgmt`) |

- **Plant C in mode 3:** RHESSys's standard plant-C outputs show frozen legacy values. FATES plant C/N is in the `_fates_*.csv` files.
- **Threads:** FATES sites run in parallel (OpenMP), so set `OMP_NUM_THREADS`. RHESSys's printed "time cost" is CPU time summed over threads, so time runs by wall clock.

### 10.6 Known limitations

- FATES GPP is about 1.7–1.9× legacy with default parameters (not calibrated).
- Soil BGC needs recalibration after a FATES spin-up, because FATES litter chemistry and CWD rates differ from legacy.
- P is prescribed (RHESSys has no P cycle).

Design and test history: `RHESSys_development/ED_into_RHESSys_design.md`.

---

## 11. Other output changes

- **Basin outputs.** `<prefix>_basin.yearly`/`.monthly` gain `streamflow_NH4 streamflow_DON streamflow_DOC` (stream patches + groundwater, g/m²) and `N_loss` (denitrification + stream NO3 + NH4 + DON). `_grow_basin.daily` gains `N_loss`, `streamflow_POC` and `streamflow_PON`. The existing `DON`/`DOC` columns are net patch DOM outflow, not stream export.
- **Yearly label.** A yearly row is printed at the start of the next year and labelled with the printing year. The row labelled 2026 therefore holds calendar 2025.
- **Mislabelled column.** The patch output column labelled `SED` is `snowpack.energy_deficit`.

---

## 12. Bug fixes that change results compared with upstream (no flag needed)

| Fix | Effect |
|---|---|
| N mass balance in `update_decomp_root.c` | immobilisation shortfall created ~3 g N/m²/yr |
| `mr_deficit` repayment (`compute_potential_N_uptake*.c`) | double-counted ~50 g C/m²/yr into cpool |
| `update_mortality`: npool zeroed when cpool < 0 | deleted ~0.26 g N/m²/yr |
| CWD fragmentation moved N at a fixed C:N below the C:N of incoming dead wood | CWD N went negative, creating ~0.06–0.09 g N/m²/yr of litter N |
| `-climrepeat` off-by-one | climate drifted 1 day per repeat cycle |
| Kinematic-wave alfa used slope^+0.5 instead of slope^−0.5 | channel areas and volumes ~slope^0.6 too small; routing attenuation too weak |
| `atm_CO2` key never matched | def value was ignored (322 ppm always) |
| NH4 deposition to the surface pool | now to soil NH4 |
| `-gwtoriparian` solute transfer | most groundwater N/DOC vanished |
| Stream routing ignored hillslope groundwater | gw now routed through reaches |
| Out-of-bounds write for reaches with no neighbouring hillslope | fixed in `construct_stream_routing_topology.c` |

**Consequence:** spin-ups made with upstream RHESSys, or with earlier builds of this fork, are not in equilibrium with the current code. Spin up again before comparing.
