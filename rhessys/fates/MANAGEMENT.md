# Forest management on FATES cohorts (`-demog_mgmt`)

Applies thinning, harvest, understory control and planting to the FATES cohorts of chosen RHESSys patches.

- **Requirements:** coupled mode 3 (`-demog <params.json> 3`).
- **Code:**
  - `fates/rhessys_fates_mgmt.F90` (management module);
  - one hook patch to FATES: `fates/patches/rhessys_host_patches.patch`.
- **Origin:** re-implements operations of J. Rady's Vegetation Management (VM) module for FATES (`github.com/JoshuaRady/fates`), targeted by RHESSys IDs instead of lat/lon.

```
rhessys7.5_fates ... -demog fates_params.json 3 -demog_mgmt my_events.txt
```

## Driver file

- One event per line.
- `!` starts a comment.
- The **first non-comment line is a header** and is skipped.
- Fields are separated by spaces. The target must not contain spaces.

```
Date        Target            Event
1955-03-15  patch:896         thin_low_perfect(final_basal_area=18, dbh_min=5)
****-11-01  hill:64           understory_control(efficiency=0.8)
1980-08-01  zone:896          harvest_trees(harvest_fraction=1.0)
1981-04-01  zone:896          plant(pfts=[1], density=1200, height=1.0, where=youngest)
1990-07-01  patches:@ids.txt  thin_proportional(thin_fraction=0.3, flux_profile=in_place)
```

**Date**
- Format `YYYY-MM-DD`; the separators `-`, `/` and `.` are all accepted.
- `*` matches any digit, so `****-01-01` fires every 1 January.

**Target**
- `all`;
- `patch:<id>` or `patches:[id,id,...]`;
- `patches:@<file>`: a file of IDs separated by whitespace; the path is relative to the run directory;
- `zone:<id>` or `zones:[...]`;
- `hill:<id>` or `hills:[...]`.

IDs that match no vegetated RHESSys patch give a warning in `<prefix>_fates.log`.

**Events**

Arguments are `name=value`, in any order. Lists use `[a,b]`.

| Event (VM alias) | Arguments (default) | Effect |
|---|---|---|
| `harvest_trees` (`harvest_timber`) | `pfts` (all), `dbh_min` (0), `dbh_max` (∞), `ht_min`, `ht_max` (m), `harvest_fraction` (1), `flux_profile`, `understory_damage` | removes that fraction of every tree in the PFT / size window |
| `thin_proportional` | `pfts`, `dbh_min` (0), one of `thin_fraction`, `final_basal_area` (m²/ha), `final_stem_density` (/ha); `flux_profile`, `understory_damage` | same fraction from every eligible tree; goals are per FATES patch |
| `thin_low_perfect` (`thin_low`) | as `thin_proportional` | removes the smallest eligible trees first until the goal is met (`thin_fraction` = fraction of stems) |
| `understory_control` | `efficiency` (1), `pfts`, `dbh_max`; `method` accepted, ignored | kills that fraction of understory (canopy layer > 1) plants in place; no new FATES patch |
| `plant` (`plant_site`) | `pfts` (required), `density` (stems/ha; one value or one per PFT), `dbh` (cm) or `height` (m; default PFT recruit height), `where` = `all` / `youngest` | adds cohorts in every FATES patch, or only the youngest |
| `harvest_biomass` | as `harvest_trees` | whole above-ground biomass harvest: trunks **and branches** of the removed trees are exported (VM `agb_harvest`); same as `flux_profile=agb_harvest` on any removal |
| `thin_low_probabilistic` | as `thin_proportional`; `steepness` (−0.4428) | VM's realistic low thinning: removal probability `1/(1+exp(−k(dbh − mean dbh − m)))`, with `m` solved so the goal is met exactly; mostly small trees, but some large ones too (VM fit to loblolly-pine thinning trials) |
| `thin_row_low` | `row_fraction` (required), `final_basal_area` or `final_stem_density`, `pfts`, `dbh_min` | only if the FATES patch is above the goal: removes the row fraction of every tree, then thins from below to the goal |
| `hardwood_control` | `efficiency` (1), `pfts` (default: woody PFTs that are not evergreen) | kills that fraction of hardwoods of all sizes in place; understory damage defaults to 0 (herbicide) |
| `replant` | as `plant`, plus `if_age_min` (alias `delay`, 0), `if_age_max` (5 yr), `bare_dbh` (2 cm) | plants only in FATES patches recently disturbed by people (time since anthropogenic disturbance in [if_age_min, if_age_max]) that are bare of trees, i.e. have no woody cohort with dbh ≥ `bare_dbh`. VM required no tree cohort at all, which FATES's own recruits would always defeat |

**Conditions (VM v1.1).** Every event, including `plant`, accepts conditions that are checked on each FATES patch. The event acts only on patches that pass. The `<event>_if` names (e.g. `thin_low_probabilistic_if`) are accepted as aliases.
- `if_age_min`, `if_age_max`: FATES patch age (yr);
- `if_bai`: minimum basal area of the eligible trees (m²/ha);
- `if_dbh_max`: their stem-weighted mean dbh must be below this (cm).

Shared arguments:
- **`flux_profile`:**
  - `bole_harvest` (default): trunks are exported as wood product (`fates_landuse_logging_export_frac`), and the rest is left as slash.
  - `in_place`: everything stays on site as litter and CWD.
- `agb_harvest`: trunks and branches are exported (VM).
- **`understory_damage`:** the fraction of understory plants killed under removed canopy crowns. The default is the FATES parameter `fates_landuse_logging_coll_under_frac`.

**`dbh_min` and stem goals.** The default `dbh_min` is 0, so stem-based goals count seedlings. Use `dbh_min=5` (or similar) for thinning by merchantable stems; VM's documentation warns about the same thing.

## How it works

**Removals.** Harvest, thinning and understory control use FATES's own logging path, with per-cohort rates set by the host:
- **Canopy trees:** the removed crown area becomes a new secondary FATES patch. Wood products, slash and understory damage follow, all booked by FATES.
- **Understory trees:** they die in place, with no new patch.

**Planting.** New cohorts are built as FATES builds recruits. Their C/N/P is booked as a generic input of the FATES mass balance, and the host N budget picks it up.

**Combining events.**
- Several events on one site and day combine; removal fractions compound.
- A RHESSys TEC clear-cut on the same day takes precedence over removals in the driver file.

**RHESSys TEC thinning events** (`redefine_world_thin_*`):
- removal < 0.99 with no kill → `thin_proportional` by leaf habit;
- removal ≥ 0.99 → clear-cut harvest (Phase 6);
- removal + kill → fire (Phase 6).

**Restarts.** Events depend only on the date, so a restarted run fires them exactly as an uninterrupted run would.

## Output

`<prefix>_fates_mgmt.csv` has one row per managed site and day:
- the events;
- stems/ha and basal area (m²/ha) before and after;
- live C before and after (kg/m²);
- wood-product C and N exported (kg/m²);
- stems planted, with their C and N (kg/m²).

The run total is in `<prefix>_fates.log`.

## VM coverage

All VM v1.0/v1.1 operations are ported, with these differences:
- **Removed:** `patch_fraction`, partial-patch operations. FATES's own disturbance machinery decides the disturbed area.
- **Defaults are by leaf habit:** hardwood PFTs default to "woody, not evergreen", because VM's index list belonged to its 12-PFT file.
- **Replant's bareness test** is relaxed to `bare_dbh` (see `replant` above).
