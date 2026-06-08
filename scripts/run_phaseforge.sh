#!/usr/bin/env bash
set -euo pipefail

source /home/kk/.venvs/phaseforge-grace-orb/bin/activate
export PATH="/home/kk/.venvs/phaseforge-grace-orb/bin:$PATH"

export PHASEFORGE_ORB_DEVICE="${PHASEFORGE_ORB_DEVICE:-cuda}"
export PHASEFORGE_ORB_COMPILE="${PHASEFORGE_ORB_COMPILE:-false}"

cd "$(dirname "$0")"

SPECIES="${PHASEFORGE_SPECIES:-}"
if [[ -z "$SPECIES" ]]; then
  if [[ -s species.in ]]; then
    SPECIES="$(tr -d '[:space:]' < species.in)"
  else
    echo "Set PHASEFORGE_SPECIES, e.g. PHASEFORGE_SPECIES=Er,Fe,Ti" >&2
    exit 2
  fi
fi
echo "$SPECIES" > species.in

LEVEL="${PHASEFORGE_LEVEL:-1}"
MLIP="${PHASEFORGE_MLIP:-ORB}"
MODEL="${PHASEFORGE_MODEL:-orb-v3-conservative-inf-omat}"
LATTICES="${PHASEFORGE_LATTICES:-BCC_A2,FCC_A1,HCP_A3}"

IFS=',' read -ra lattice_list <<< "$LATTICES"

for phase in "${lattice_list[@]}"; do
  echo "Calculating ${phase} phase"

  # sqs2tdb intentionally needs the copy command twice: the first call creates
  # the phase/species files, the second populates SQS directories.
  sqs2tdb -cp -sp="$SPECIES" -l="$phase" -lv="$LEVEL"
  sqs2tdb -cp -sp="$SPECIES" -l="$phase" -lv="$LEVEL"

  while IFS= read -r str_file; do
    structure_dir="$(dirname "$str_file")"
    echo "  MLIP relaxing ${structure_dir}"
    (
      cd "$structure_dir"
      if [[ -s energy ]]; then
        echo "    energy exists; skipping"
        exit 0
      fi
      runstruct_vasp -nr >/dev/null 2>&1 || true
      if [[ ! -s POSCAR ]]; then
        cellcvrt -c -sig=9 < str.out > cellcvrt.tmp
        python /home/kk/bin/atat_to_poscar.py cellcvrt.tmp POSCAR
        rm -f cellcvrt.tmp
      fi
      MLIPrelax -mlip="$MLIP" -model="$MODEL"
      extract_MLIP
    )
  done < <(find "$phase" -mindepth 2 -maxdepth 2 -name str.out | sort)

  (
    cd "$phase"
    if [[ "$phase" == "BCC_A2" || "$phase" == "FCC_A1" || "$phase" == "HCP_A3" || "$phase" == "LIQUID" ]]; then
      printf '1,0\n2,0\n' > terms.in
    else
      printf '1,0:1,0\n2,0:1,0\n' > terms.in
    fi
    sqs2tdb -fit
  )
done

sqs2tdb -tdb -oc
