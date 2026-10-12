#!/usr/bin/env bash
#==============================================================================
# bit_identity_test.sh -- off-mode (legacy) bit-identity regression
#
# Part of the ED-into-RHESSys / FATES-backend work (design doc 6.5, Phase 0
# exit criterion #4). PROVES that adding the pluggable vegetation-backend
# scaffolding (demog_flag / prime_veg / rhessys_veg.* ) does NOT change legacy
# behavior: a benchmark run with NO -demog must be BYTE-IDENTICAL to a
# reference snapshot captured from the pre-change code.
#
# Workflow:
#   1. Build+run the UNMODIFIED code, then:  ./bit_identity_test.sh --bless
#      (records the reference snapshot).
#   2. Apply the scaffolding, rebuild, then:  ./bit_identity_test.sh
#      (reruns + diffs; exits non-zero on ANY byte difference).
#
# Config via environment (or edit the CONFIG block below):
#   RHESSYS    path to the rhessys binary          (default ../rhessys7.5)
#   RUN_DIR    dir to run the benchmark from        (default ./benchmark)
#   RUN_ARGS   rhessys args (legacy; NO -demog)     (REQUIRED)
#   OUT_SUBDIR output subdir the run writes into     (default out)
#   REF_DIR    reference snapshot dir               (default ./reference)
#   SEED_DIR   optional dir whose contents are copied into OUT_SUBDIR before
#              each run (e.g. the worldfile, so output_current_state dumps
#              land in the compared dir)
#==============================================================================
set -euo pipefail

#------------------------------- CONFIG ---------------------------------------
RHESSYS="${RHESSYS:-../rhessys7.5}"
RUN_DIR="${RUN_DIR:-./benchmark}"
OUT_SUBDIR="${OUT_SUBDIR:-out}"
REF_DIR="${REF_DIR:-./reference}"
SEED_DIR="${SEED_DIR:-}"
# Example (edit to your benchmark; MUST be legacy mode -- no -demog):
#   RUN_ARGS="-st 2000 1 1 1 -ed 2003 1 1 1 -b worldfiles/world \
#             -t tecfiles/tec.full -r flowtables/sub.flow flowtables/surf.flow \
#             -pre out/bench -g"
RUN_ARGS="${RUN_ARGS:-}"
#------------------------------------------------------------------------------

die() { echo "ERROR: $*" >&2; exit 2; }

[ -n "$RUN_ARGS" ] || die "RUN_ARGS not set (benchmark rhessys args; legacy, no -demog)."
[ -x "$RHESSYS" ]  || die "rhessys binary not found/executable: $RHESSYS"
[ -d "$RUN_DIR" ]  || die "RUN_DIR does not exist: $RUN_DIR"
case " $RUN_ARGS " in
  *" -demog "*) die "RUN_ARGS must NOT contain -demog -- this is the LEGACY baseline." ;;
esac

# Absolute path to the binary so we can cd into RUN_DIR to run.
RHESSYS_ABS="$(cd "$(dirname "$RHESSYS")" && pwd)/$(basename "$RHESSYS")"
MODE="${1:-compare}"

run_benchmark() {
  echo ">> running benchmark: $RHESSYS_ABS $RUN_ARGS   (cwd=$RUN_DIR)"
  rm -rf "${RUN_DIR:?}/${OUT_SUBDIR}"
  mkdir -p "${RUN_DIR}/${OUT_SUBDIR}"
  [ -z "$SEED_DIR" ] || cp -R "${SEED_DIR}/." "${RUN_DIR}/${OUT_SUBDIR}/"
  ( cd "$RUN_DIR" && "$RHESSYS_ABS" $RUN_ARGS ) || die "rhessys run failed"
}

case "$MODE" in
  --bless)
    run_benchmark
    rm -rf "$REF_DIR"; mkdir -p "$REF_DIR"
    cp -R "${RUN_DIR}/${OUT_SUBDIR}/." "$REF_DIR/"
    echo "BLESSED: reference snapshot saved to $REF_DIR"
    ;;

  compare)
    [ -d "$REF_DIR" ] || die "no reference at $REF_DIR -- run './bit_identity_test.sh --bless' first."
    run_benchmark
    echo ">> diffing ${RUN_DIR}/${OUT_SUBDIR} against reference $REF_DIR"
    if diff -r --brief "$REF_DIR" "${RUN_DIR}/${OUT_SUBDIR}"; then
      echo "PASS: off-mode output is BYTE-IDENTICAL to the reference."
      exit 0
    else
      echo "FAIL: off-mode output DIFFERS from the reference (legacy path was perturbed)." >&2
      echo "      Investigate before proceeding -- Phase 0 requires bit-identity." >&2
      exit 1
    fi
    ;;

  *)
    echo "usage: $0 [--bless | compare]" >&2
    exit 2
    ;;
esac
