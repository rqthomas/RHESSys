# WMFire with FATES (`-demog 3` + `-firespread`)

WMFire (`RHESSys/FIRE/`, a C++ library) spreads fire across a raster grid once a month. RHESSys then applies the fire's effects.

In mode 3, the vegetation it acts on is FATES:
- **Where and how hard it burns:** WMFire decides.
- **What it burns:** FATES fuels.
- **Who dies:** FATES cohorts, by height, following WMFire's own severity rules.

Legacy-only runs (no `-demog`) are unchanged.

## Build (`ccr-rhessys` docker image; Boost headers are installed there)

```
cd RHESSys/FIRE
g++ -O2 -fPIC -c RanNums.cpp -o /tmp/RanNums.o
g++ -O2 -fPIC -c WMFire.cpp -o /tmp/WMFire.o
g++ -shared -fPIC -o ../lib/libwmfire.so /tmp/RanNums.o /tmp/WMFire.o
cd ../rhessys
make clobber PGM=__none__
make fates=T wmfire=T PGM=rhessys7.5_fates_fire      # or: make wmfire=T (legacy)
```

Running needs `LD_LIBRARY_PATH=<repo>/RHESSys/lib`.

## Inputs (CCR set-up, `ed_bitid/`)

| What | Where | How made |
|---|---|---|
| Fire grids `ccr.dem/.patch/.zone/.hillslope` (369 × 373, 26.465 m, 65,403 cells) | `ed_bitid/wmfire/grid/` | `RHESSysPreprocessing::write_fire_grids('ccr', template_ccr_v2.txt, spatial_data)`; it writes into the map folder, so move the files out afterwards |
| Fire parameters | `benchmark/defs/fire_ccr.def` | `n_rows 369`, `n_cols 373`, `calc_fire_effects 1`, `fire_write 1`, `ran_seed 12345`; everything else is WMFire default |
| World header | `benchmark/hdr_fire/ccr/ccr_fire.hdr` | `ccr.hdr` plus `1 num_fire_default_files` / def path / `1 num_fire_grid_prefix` / grid prefix (absolute path), after the stratum defaults |
| World file | `ed_bitid/wmfire/start_fire.state` | seed state + `1 fire_parm_ID` in every patch (mandatory with `-firespread`) |

Command-line additions: `-firespread 26.46504`, the cell size in metres.

WMFire writes `FireSizes0.txt` (cells burned, year, month, wind, direction, ignitions) into the run directory.

## Coupling (mode 3)

1. **Fuels.**
   - Live fuel: FATES leaf C.
   - Litter fuel: RHESSys litter pools plus the FATES standing dead leaves.
   - CWD fuel: FATES above-ground CWD.
   - Legacy strata are not read for FATES sites.
2. **One fire per patch per month.** WMFire's per-cell burns are summed for each patch family:
   - burned fraction `f` = Σ share of the burned cells;
   - intensity `b` = mean burn of those cells.

   The upstream code instead applies effects once per burned cell. With CCR's ~43-cell patches, that strips the litter of the whole patch for each burned cell.
3. **Severity in FATES.** In the burned fraction `f` (FATES `frac_burnt`, capped at 0.98), each cohort gets WMFire's mortality at intensity `b`:
   - understory (h ≤ 4 m): `(k^b − 1)/(k − 1)`, with `k = understory_mort`;
   - overstory (h > 6 m): `1 − 1/(1 + exp(−k1 (U − k2)))`, where `U` = understory C consumed + litter C consumed;
   - between 4 and 6 m: a linear blend.

   Crown consumption uses the same form with `k = consumption`. FATES fuel consumption uses WMFire's fractions: leaves 1, twigs and small branches 0.85, large branches and trunks 0.339. FATES does the patch split, killed trees → litter/CWD, and burned C/N → atmosphere. The burned C/N is in the host N budget.
4. **RHESSys bookkeeping.**
   - RHESSys litter and surface-soil losses (WMFire's fractions) are scaled by `f`.
   - Legacy strata are left untouched (inert shadow).
   - TEC disturbance events on the same day take precedence.

## Test at CCR (1850–1856, same seed)

| Run | Fires | Effect |
|---|---|---|
| Legacy WMFire | 12 fires; 3 large (86 %, 28 %, 48 % of the watershed) | CWD −80 %, litter −50…60 %, soil3 −20 %; plant C up (less N immobilisation) |
| Mode 3, before coupling | fires happen | FATES untouched; only RHESSys litter burns; legacy strata modified |
| **Mode 3, coupled** | 7,887 site fires, 21,286 ha burned (watershed 4,580 ha) | understory stems −75…85 %, trees ≥ 20 cm unchanged; FATES litter + CWD 2.1 vs 3.3 kgC/m²; LAI 4.4 vs 6.0; 1,646 gC/m² + 19.9 gN/m² burned; N budget 1e-17; FATES mass check never fails |

Regression gates:
- legacy (no fire) 32/32 byte-identical;
- legacy WMFire reproduces its pre-coupling run (33/33);
- the coupled binary without `-firespread` matches Phase 6 (34/34).

## Caveats

- **The fire regime is not calibrated for CCR.** With WMFire's default parameters (Sierra Nevada), CCR burns several times its area in 6 years: a fire rotation of about 1–4 years, against decades to centuries for Piedmont oak–pine. Calibrating `mean_ign`, `ignition_tmin`, the moisture / load / deficit parameters (`moisture_k*`, `load_k*`, `ign_def_mod`) against a target fire-return interval comes next.
- **Simplifications:**
  - `U` uses FATES understory above-ground C per m² of site; legacy uses C per m² of understory cover.
  - The intermediate story is not separated in `U`.
- **Upstream WMFire issues noticed** (left as is in the legacy path):
  - effects applied per burned cell (above);
  - `if (multiscale_flag = 1)` assignments (`execute_firespread_event.c`);
  - an unreachable intermediate-height branch in `compute_family_fire_effects.c`.
