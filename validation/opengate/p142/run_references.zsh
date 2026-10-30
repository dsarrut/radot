#!/usr/bin/env zsh
# Run the OpenGATE reference simulations of p142 (VRT radionuclides) with opengate_jobs, locally.
#
#   validation/opengate/p142/run_references.zsh [number_of_jobs]     (default 12)
#
# Needs opengate and the data in p142/data/ (python validation/opengate/download_data.py).
# For each radionuclide, create_campaigns.py builds the electron and photon campaigns (bremsstrahlung
# included), runs each one in parallel jobs (split, run, merge), then assembles
# p142/data/new/output_<radionuclide>_vrt/ (merged images, e/ and gamma_tle/, reference.json with the
# activity, opengate version).
# Ac225 (analog mode) is not run here.
#
# Activities (more activity = more statistics; results are normalised per MBq), measured on a laptop:
#   Y90 30 MBq ~3.2 CPU-h, Lu177 10 MBq ~1.2 CPU-h, I131 5 MBq ~2.3 CPU-h: ~35 min with 12 jobs.
set -euo pipefail

JOBS="${1:-12}"
P142_DIR="${0:A:h}"
cd "$P142_DIR/data"

for name activity in Y90 30e6 Lu177 10e6 I131 5e6; do
  echo "=== $name: $activity Bq, $JOBS jobs"
  python "$P142_DIR/create_campaigns.py" --radionuclide "$name" --activity "$activity" \
      --run -j "$JOBS" -o "new/output_${name:l}_vrt"
done

echo "Done: $P142_DIR/data/new/"
