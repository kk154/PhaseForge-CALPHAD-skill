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
TARGET_TEMPERATURE="${PHASEFORGE_TEMPERATURE_K:-}"
LIQUID_TEMPERATURES="${PHASEFORGE_LIQUID_TEMPERATURES:-}"

contains_liquid=false
IFS=',' read -ra lattice_list <<< "$LATTICES"
for lattice in "${lattice_list[@]}"; do
  if [[ "${lattice^^}" == "LIQUID" ]]; then
    contains_liquid=true
  fi
done

if [[ "$contains_liquid" == "true" && -z "$TARGET_TEMPERATURE" ]]; then
  echo "PHASEFORGE_TEMPERATURE_K is required when PHASEFORGE_LATTICES includes LIQUID." >&2
  exit 2
fi

if [[ "$contains_liquid" == "true" && -z "$LIQUID_TEMPERATURES" ]]; then
  LIQUID_TEMPERATURES="$(python - "$TARGET_TEMPERATURE" <<'PY'
import sys
t = float(sys.argv[1])
temps = [max(1.0, t - 100.0), t, t + 100.0]
print(",".join(f"{value:g}" for value in temps))
PY
)"
fi

temperature_dir_name() {
  python - "$1" <<'PY'
import sys
value = float(sys.argv[1])
print("T_" + f"{value:g}".replace(".", "p") + "K")
PY
}

for phase in "${lattice_list[@]}"; do
  echo "Calculating ${phase} phase"

  # sqs2tdb intentionally needs the copy command twice: the first call creates
  # the phase/species files, the second populates SQS directories.
  sqs2tdb -cp -sp="$SPECIES" -l="$phase" -lv="$LEVEL"
  sqs2tdb -cp -sp="$SPECIES" -l="$phase" -lv="$LEVEL"

  while IFS= read -r str_file; do
    structure_dir="$(dirname "$str_file")"
    echo "  Processing ${structure_dir}"
    (
      cd "$structure_dir"
      if [[ -s energy && "${phase^^}" != "LIQUID" ]]; then
        echo "    energy exists; skipping"
        exit 0
      fi
      runstruct_vasp -nr >/dev/null 2>&1 || true
      if [[ ! -s POSCAR ]]; then
        cellcvrt -c -sig=9 < str.out > cellcvrt.tmp
        python /home/kk/bin/atat_to_poscar.py cellcvrt.tmp POSCAR
        rm -f cellcvrt.tmp
      fi
      if [[ "${phase^^}" == "LIQUID" ]]; then
        IFS=',' read -ra liquid_temp_list <<< "$LIQUID_TEMPERATURES"
        for temp in "${liquid_temp_list[@]}"; do
          temp_dir="$(temperature_dir_name "$temp")"
          mkdir -p "$temp_dir"
          cp POSCAR "$temp_dir/POSCAR"
          (
            cd "$temp_dir"
            if [[ -s liquid_md_summary.json ]]; then
              echo "    liquid MD exists for ${temp} K; skipping"
              exit 0
            fi
            if [[ "${PHASEFORGE_LIQUID_NO_RUN:-false}" =~ ^(1|true|yes|on)$ ]]; then
              python ../../../run_liquid_lammps.py --poscar POSCAR --temperature "$temp" --mlip "$MLIP" --model "$MODEL" --no-run
            else
              python ../../../run_liquid_lammps.py --poscar POSCAR --temperature "$temp" --mlip "$MLIP" --model "$MODEL"
            fi
          )
        done
        target_dir="$(temperature_dir_name "$TARGET_TEMPERATURE")"
        if [[ ! -s "$target_dir/energy" ]]; then
          echo "Missing target-temperature liquid energy: ${structure_dir}/${target_dir}/energy" >&2
          exit 3
        fi
        cp "$target_dir/energy" energy
      else
        MLIPrelax -mlip="$MLIP" -model="$MODEL"
        extract_MLIP
      fi
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

if [[ "$contains_liquid" == "true" ]]; then
  python collect_liquid_md_energies.py --root LIQUID --output liquid_md_energies.csv
fi

sqs2tdb -tdb -oc
