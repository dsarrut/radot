#!/usr/bin/env zsh
# Run the OpenGATE reference simulations of p142 (VRT radionuclides) with opengate_jobs, locally.
#
#   validation/opengate/p142/run_references.zsh [number_of_jobs]     (default 12)
#
# Needs opengate and the data in p142/data/ (python validation/opengate/download_data.py).
# For each radionuclide: builds the electron and photon campaigns (create_campaigns.py, bremsstrahlung
# included), runs each one in parallel jobs (split, run, merge), then merges both parts into
# p142/data/new/output_<radionuclide>_vrt/, with the opengate version.
# Ac225 (analog mode) is not run here.
#
# Activities (more activity = more statistics; results are normalised per MBq), measured on a laptop:
#   Y90 30 MBq ~3.2 CPU-h, Lu177 10 MBq ~1.2 CPU-h, I131 5 MBq ~2.3 CPU-h: ~35 min with 12 jobs.
set -euo pipefail

JOBS="${1:-12}"
P142_DIR="${0:A:h}"
cd "$P142_DIR/data"

# radionuclide and activity (Bq); the e-factor is the one of the p142 configurations (1)
run() {
  local name="$1" activity="$2"
  local r="${name:l}"
  echo "=== $name: $activity Bq, $JOBS jobs"
  python "$P142_DIR/create_campaigns.py" --radionuclide "$name" --activity "$activity"
  for part in e gamma_tle; do
    opengate_jobs run "campaigns/$r/$part" -j "$JOBS"
  done
  opengate_dose_rate --merge-only "campaigns/$r/e/output" "campaigns/$r/gamma_tle/output" \
      --e-factor 1 -o "new/output_${r}_vrt"
  cp "campaigns/$r/opengate_version.json" "new/output_${r}_vrt/"
}

run Y90 30e6
run Lu177 10e6
run I131 5e6

echo "Done: $P142_DIR/data/new/"
