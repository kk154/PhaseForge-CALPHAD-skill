---
name: mat-phaseforge-calphad
description: Compute calculation-driven ternary isothermal phase diagrams for any user-specified chemical system and temperature using PhaseForge, ORB/MLIP structure energetics, literature/binary-phase-diagram compound discovery, Materials Project structure downloads, local CIFs, pycalphad TDB construction, phonopy vibrational free-energy and formation-vibrational-free-energy corrections for the relaxed structures, high-resolution grid sampling, and Thermo-Calc-style colored straight-boundary plotting. Use when asked to identify binary/ternary compounds in a system, download their structures from MP, calculate or refine a CALPHAD-like computed phase diagram workflow, or plot a computed isothermal section without fitting to experiment.
---

# PhaseForge CALPHAD Phase Diagrams

Use this skill to calculate an isothermal ternary phase diagram for a user-given system and temperature. Treat Er-Fe-Ti only as the provenance of the workflow, not as a default system. The user must provide or imply the component list, temperature, available structures/CIFs, and target plotting style.

Do not fit, tune, or optimize parameters against an experimental phase diagram unless the user explicitly asks for a separate inverse-modeling task. Experimental figures may be used for visual comparison only when the user allows it.

The default final diagram must use phonon-corrected Gibbs energies. After ORB relaxation, calculate phonopy vibrational free energies for the relaxed compound structures and the elemental reference structures, convert them to formation vibrational free energies, apply those corrections to the TDB, and sample/plot the corrected TDB. A 0 K ORB-only TDB is an intermediate diagnostic, not the final result, unless the user explicitly asks to skip phonons or a phonon calculation is infeasible.

## Inputs To Establish

- Components: exactly three real components for the ternary section, for example `ER,FE,TI` or `NI,AL,CR`.
- Temperature: Kelvin preferred; convert Celsius by `T_K = T_C + 273.15`.
- PhaseForge level and lattices: usually `PHASEFORGE_LEVEL=2` and `BCC_A2,FCC_A1,HCP_A3,FCC_L12,...` as relevant.
- Structures: local CIFs first; query Materials Project only for missing phases.
- Elemental phonon reference structures: one relaxed structure for each component, with phase labels such as `A_REF`, `B_REF`, `C_REF`.
- Compound and solution-phase models: derive from structures and chemistry of the target system. Do not reuse another system's sublattice model.
- Plot resolution: use `step=0.005` for detailed line-style Thermo-Calc plots when runtime permits.

Never store Materials Project API keys in the skill, scripts, docs, logs, or committed files. Use an environment variable such as `MP_API_KEY`. When the user gives a new chemical system, browse current literature and phase-diagram sources before choosing compounds; cite sources in the response.

## Bundled Scripts

Scripts live in `scripts/`; copy them into the run directory and patch per-system scripts there when needed.

- `run_phaseforge.sh`: generic PhaseForge SQS solid-solution energetics. Set `PHASEFORGE_SPECIES=A,B,C`.
- `calc_cif_compound_energies.py`: ORB relaxation for local CIF compounds. This script may need minor phase-name ordering edits if the current filename-to-phase naming is not appropriate.
- `download_mp_structures.py`: download CIFs from Materials Project for candidate compounds found from literature/phase-diagram review.
- `merge_compound_energies.py`: merge compound energy tables.
- `build_compound_tdb.py`: add stoichiometric compound phases using a per-system reference CSV.
- `prepare_phonon_targets.py`: create `phonon_targets.csv` and `phonon_references.csv` from relaxed compounds and elemental reference structures.
- `calc_vibrational_free_energy.py`: phonopy/ORB vibrational free energies for a per-system target CSV.
- `make_vibrational_corrections.py`: convert vibrational free energies to formation vibrational corrections using a per-system element-reference CSV.
- `apply_thermal_corrections.py`: apply finite-temperature Gibbs shifts from a correction CSV to a TDB.
- `sample_ternary_isotherm.py`: generic pycalphad ternary grid sampler.
- `plot_thermocalc_style_ternary.py`: generic Thermo-Calc-style colored ternary plotter with optional straight fitted boundaries.

## Workflow

1. Create a run directory and copy scripts.

```powershell
$Skill = 'C:\Users\15461\.codex\skills\mat-phaseforge-calphad'
$Run = 'E:\codex项目\CALPHAD\phaseforge_runs\<SYSTEM>_<TEMP>'
New-Item -ItemType Directory -Force -Path $Run | Out-Null
Copy-Item -LiteralPath "$Skill\scripts\*" -Destination $Run -Force
```

2. Discover binary and ternary compounds before calculating compound energies.

For a new system `A-B-C`, first build `candidate_compounds.csv` from literature, assessed binary phase diagrams, ternary phase-diagram papers, handbooks, and database cross-checks. Browse the web for current information and cite the sources used. Include all known stable and important metastable binary and ternary intermetallic/compound phases that may affect the isothermal section; include solution phases or homogeneity ranges as notes for later sublattice modeling.

Recommended source order:

- Assessed binary phase diagrams for `A-B`, `A-C`, and `B-C`.
- Ternary phase-diagram or isothermal-section papers for `A-B-C`.
- ICSD/COD/Materials Project formula checks when literature structures are unclear.
- Materials Project hull phases for each binary and the ternary as a coverage check, not as the sole source of phase knowledge.

Candidate table format:

```csv
phase,formula,system,kind,material_id,source,notes
P_AB2,AB2,A-B,binary,,<citation or URL>,prototype or phase range if known
P_ABC,A1B1C1,A-B-C,ternary,mp-0000,<citation or URL>,use material_id when known
```

Then download structures from MP where possible:

```powershell
wsl bash -lc "cd /mnt/e/codex项目/CALPHAD/phaseforge_runs/<SYSTEM>_<TEMP> && source /home/kk/.venvs/phaseforge-grace-orb/bin/activate && export MP_API_KEY=\$MP_API_KEY && python download_mp_structures.py --candidates candidate_compounds.csv --output-dir mp_structures --manifest mp_structure_manifest.csv --e-above-hull-max 0.10"
```

If MP has no structure for an experimentally reported compound, keep it in `candidate_compounds.csv`, note `status=not_found` in the manifest, and search COD/ICSD/literature CIFs or construct a prototype-substituted structure. Do not silently drop literature phases because MP lacks them.

3. Run PhaseForge/ORB solid-solution energetics.

```powershell
wsl bash -lc "cd /mnt/e/codex项目/CALPHAD/phaseforge_runs/<SYSTEM>_<TEMP> && source /home/kk/.venvs/phaseforge-grace-orb/bin/activate && PHASEFORGE_SPECIES=A,B,C PHASEFORGE_LEVEL=2 PHASEFORGE_ORB_DEVICE=cuda PHASEFORGE_LATTICES=BCC_A2,FCC_A1,HCP_A3 bash run_phaseforge.sh"
```

4. Relax compound CIFs with ORB. Use CIFs supplied by the user first, then MP-downloaded CIFs, then manually sourced/prototype structures. Document every structure origin in the manifest.

```powershell
wsl bash -lc "cd /mnt/e/codex项目/CALPHAD/phaseforge_runs/<SYSTEM>_<TEMP> && source /home/kk/.venvs/phaseforge-grace-orb/bin/activate && python calc_cif_compound_energies.py --components A,B,C --cif-dir mp_structures --output-dir compound_orb --max-orderings 10"
```

5. Build the compound TDB.

Before running `build_compound_tdb.py`, inspect the base TDB and the PhaseForge endpoint energies for the current system. Create `element_references.csv`:

```csv
element,stable_ref_func,orb_reference_eV_atom
A,FUNC0000,-0.000000
B,FUNC0001,-0.000000
C,FUNC0002,-0.000000
```

`stable_ref_func` must map each element to its stable-element reference function in the base TDB, and `orb_reference_eV_atom` must come from the current PhaseForge/ORB run. Never reuse values from another system.

Then run:

```powershell
wsl bash -lc "cd /mnt/e/codex项目/CALPHAD/phaseforge_runs/<SYSTEM>_<TEMP> && source /home/kk/.venvs/phaseforge-grace-orb/bin/activate && python build_compound_tdb.py --base-tdb <BASE>.tdb --components A,B,C --references element_references.csv --energies compound_orb/compound_energies.csv --temperature <T_K> --output <SYSTEM>_orb_compounds.tdb"
```

6. Add solution phases only when justified by structure or chemistry.

For a target system with mixed-sublattice compounds or known homogeneity ranges, create a system-specific endmember script in the run directory. Derive the sublattice model from CIF mixed occupancies, crystallographic sites, or chemically meaningful substitutions. Calibrate interaction terms only from computed endmember/compound energies unless the user explicitly requests experimental fitting.

7. Calculate phonon free energies and correct the TDB before sampling.

Use ideal configurational entropy from CIF partial occupancies where present. Then calculate vibrational free energies for the relaxed structures that are actually used in the TDB. Include all stoichiometric compound phases that may appear in the final diagram and one elemental reference structure per component. Do not treat phonons as optional for the final result unless the user explicitly asks to skip them or the calculation fails after reasonable attempts.

First create `element_phonon_references.csv` from current-system elemental reference structures:

```csv
element,phase,path
A,A_REF,path/to/A/CONTCAR
B,B_REF,path/to/B/CONTCAR
C,C_REF,path/to/C/CONTCAR
```

Generate the phonon target and reference tables from the ORB-relaxed compound structures:

```powershell
wsl bash -lc "cd /mnt/e/codex项目/CALPHAD/phaseforge_runs/<SYSTEM>_<TEMP> && source /home/kk/.venvs/phaseforge-grace-orb/bin/activate && python prepare_phonon_targets.py --compound-energies compound_orb/compound_energies.csv --reference-structures element_phonon_references.csv --targets-output phonon_targets.csv --references-output phonon_references.csv"
```

Run at the target temperature:

```powershell
wsl bash -lc "cd /mnt/e/codex项目/CALPHAD/phaseforge_runs/<SYSTEM>_<TEMP> && source /home/kk/.venvs/phaseforge-grace-orb/bin/activate && python calc_vibrational_free_energy.py --targets phonon_targets.csv --temperature <T_K> --mesh 12,12,12 --supercell 2,2,2 --output-dir phonon_free_energy"
```

Use `--supercell 2,2,2` as the rigorous default for production phase diagrams. For very large structures, `--supercell 1,1,1` can be used only as a fast approximation, and the final response must report that approximation clearly.

Convert to formation vibrational corrections, where `correction_eV_atom = F_vib(phase) - sum(x_i F_vib(element_i))`, then apply it:

```powershell
wsl bash -lc "cd /mnt/e/codex项目/CALPHAD/phaseforge_runs/<SYSTEM>_<TEMP> && source /home/kk/.venvs/phaseforge-grace-orb/bin/activate && python make_vibrational_corrections.py --fvib phonon_free_energy/vibrational_free_energies.csv --references phonon_references.csv --output thermal_corrections_from_phonons.csv"
wsl bash -lc "cd /mnt/e/codex项目/CALPHAD/phaseforge_runs/<SYSTEM>_<TEMP> && source /home/kk/.venvs/phaseforge-grace-orb/bin/activate && python apply_thermal_corrections.py --base-tdb <SYSTEM>_orb_compounds.tdb --corrections thermal_corrections_from_phonons.csv --output <SYSTEM>_orb_phonon_<T_K>K.tdb"
```

The corrected `<SYSTEM>_orb_phonon_<T_K>K.tdb` is the default final TDB for the phase diagram. Compare it with the uncorrected TDB only as a diagnostic.

8. Sample the ternary section at high resolution.

```powershell
wsl bash -lc "cd /mnt/e/codex项目/CALPHAD/phaseforge_runs/<SYSTEM>_<TEMP> && source /home/kk/.venvs/phaseforge-grace-orb/bin/activate && python sample_ternary_isotherm.py --tdb <SYSTEM>_orb_phonon_<T_K>K.tdb --components A,B,C --temperature <T_K> --step 0.005 --csv <SYSTEM>_<T_K>K_phonon_step0005_grid.csv"
```

9. Plot in Thermo-Calc style with colored regions and straight fitted phase boundaries.

```powershell
wsl bash -lc "cd /mnt/e/codex项目/CALPHAD/phaseforge_runs/<SYSTEM>_<TEMP> && source /home/kk/.venvs/phaseforge-grace-orb/bin/activate && python plot_thermocalc_style_ternary.py --grid <SYSTEM>_<T_K>K_phonon_step0005_grid.csv --components A,B,C --temperature '<T_K> K' --output <SYSTEM>_<T_K>K_phonon_thermocalc_style.png --label-column assemblage --max-internal-labels 26 --straight-boundaries --min-straight-boundary-points 25"
```

Use `--markers markers.csv` for important computed compounds. The marker CSV columns are `label,x_a,x_b,x_c`. Use `--phase-aliases aliases.csv` for cleaner display labels with columns `phase,label`.

## Validation Checklist

- Confirm the component order is consistent across PhaseForge, TDB construction, grid sampling, and plotting.
- Confirm `candidate_compounds.csv` covers all three binary subsystems and the ternary literature/MP cross-check.
- Confirm `mp_structure_manifest.csv` records MP IDs, formulas, hull energies, CIF paths, and not-found compounds.
- Confirm all hard-coded element references, endpoint ORB energies, phase names, and sublattice models are for the current system.
- Confirm `phonon_targets.csv` includes elemental references and every relaxed compound phase used in the corrected TDB.
- Confirm `thermal_corrections_from_phonons.csv` contains formation vibrational free-energy corrections for all compound phases that require correction.
- Check that the final grid CSV contains `i`, `j`, three `x_<element>` columns, `assemblage`, `dominant`, and `plot_label`.
- Check phonon result summaries for imaginary modes and report any adopted approximation.
- Confirm the final plot has colored assemblage regions, a complete legend, readable internal labels, and straight black boundary segments when requested.
- Report final phonon-corrected `.tdb`, grid `.csv`, `.png`, and `.pdf` paths.
