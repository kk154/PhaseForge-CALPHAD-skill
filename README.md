# PhaseForge CALPHAD Skill

中文 | [English](#english)

## 简介

`mat-phaseforge-calphad` 是一个 Codex skill，用于基于 PhaseForge、ORB/MLIP 结构能量、Materials Project 结构检索、pycalphad TDB 构建、声子自由能校正、SQS/LAMMPS 液相建模和三元等温截面采样，生成计算驱动的 CALPHAD-like 三元相图工作流。

它适用于：

- 为任意三元体系建立候选二元/三元化合物清单；
- 下载或整理 Materials Project、本地 CIF、文献原型结构；
- 计算固溶体和化合物的 MLIP/ORB 能量；
- 对 relaxed 结构计算声子自由能和形成振动自由能校正；
- 对目标温度下可能出现的液相进行 SQS/LAMMPS 多温度焓计算；
- 构建包含化合物相、声子校正固相和 practical CALPHAD 液相的 TDB；
- 使用 pycalphad 采样高分辨率等温截面；
- 绘制接近 Thermo-Calc 风格的彩色三元相图。

> 该 skill 默认是“计算驱动”工作流，不会拟合实验相图。只有在用户明确要求反演或参数拟合时，才应使用实验数据进行拟合。

## 仓库结构

```text
.
├── SKILL.md
├── agents/
│   └── openai.yaml
├── scripts/
│   ├── run_phaseforge.sh
│   ├── run_liquid_lammps.py
│   ├── collect_liquid_md_energies.py
│   ├── fit_liquid_tdb.py
│   ├── download_mp_structures.py
│   ├── calc_cif_compound_energies.py
│   ├── merge_compound_energies.py
│   ├── build_compound_tdb.py
│   ├── prepare_phonon_targets.py
│   ├── calc_vibrational_free_energy.py
│   ├── make_vibrational_corrections.py
│   ├── apply_thermal_corrections.py
│   ├── sample_ternary_isotherm.py
│   └── plot_thermocalc_style_ternary.py
└── LICENSE
```

## 安装到 Codex

将本仓库克隆到 Codex skills 目录：

```powershell
git clone https://github.com/kk154/PhaseForge-CALPHAD-skill.git C:\Users\15461\.codex\skills\mat-phaseforge-calphad
```

如果你的 Codex skills 目录不同，请把目标路径替换成你的本地 skills 根目录。

## 典型用法

在 Codex 中请求类似任务：

```text
使用 mat-phaseforge-calphad skill，帮我计算 Ni-Al-Cr 在 1000 K 的三元等温相图。
```

skill 会引导 Codex 建立运行目录、复制脚本、检索候选化合物、运行 PhaseForge/ORB、对 relaxed 固相结构计算声子自由能、对 LIQUID SQS 结构运行 LAMMPS 多温度液相焓计算、构建最终 TDB、采样相图并绘图。具体流程和命令模板见 `SKILL.md`。

## 依赖

该 skill 的脚本通常需要以下环境：

- PhaseForge；
- ORB/MLIP 推理环境；
- Python 3；
- pycalphad；
- pymatgen 和 mp-api；
- phonopy；
- 支持 MLIP/GNNP 的 LAMMPS 液相计算环境；
- matplotlib、numpy、pandas、scipy；
- 可选 GPU/CUDA 环境。

实际依赖取决于你执行的步骤。例如，只绘图和采样不需要运行 ORB；下载 Materials Project 结构需要 `MP_API_KEY`。

## 许可证

本仓库使用 MIT License。详见 `LICENSE`。

---

## English

`mat-phaseforge-calphad` is a Codex skill for calculation-driven CALPHAD-like ternary isothermal phase-diagram workflows. It combines PhaseForge, ORB/MLIP structure energetics, Materials Project structure retrieval, pycalphad TDB generation, phonon free-energy corrections, SQS/LAMMPS liquid modeling, high-resolution ternary grid sampling, and Thermo-Calc-style plotting.

It is useful for:

- building binary and ternary candidate-compound lists for a user-specified ternary system;
- collecting structures from Materials Project, local CIF files, and literature prototypes;
- computing solid-solution and compound energetics with MLIP/ORB workflows;
- calculating phonon free energies and formation vibrational free-energy corrections for relaxed structures;
- computing multi-temperature liquid enthalpies for LIQUID SQS structures with LAMMPS;
- constructing TDB files with stoichiometric compounds, phonon-corrected solids, and practical CALPHAD LIQUID terms;
- sampling ternary isothermal sections with pycalphad;
- plotting colored ternary phase diagrams with readable labels and fitted straight boundaries.

> This skill is intended for a computation-first workflow. It does not fit or tune parameters against experimental phase diagrams unless the user explicitly asks for a separate inverse-modeling task.

## Repository Layout

```text
.
├── SKILL.md
├── agents/
│   └── openai.yaml
├── scripts/
│   ├── run_phaseforge.sh
│   ├── run_liquid_lammps.py
│   ├── collect_liquid_md_energies.py
│   ├── fit_liquid_tdb.py
│   ├── download_mp_structures.py
│   ├── calc_cif_compound_energies.py
│   ├── merge_compound_energies.py
│   ├── build_compound_tdb.py
│   ├── prepare_phonon_targets.py
│   ├── calc_vibrational_free_energy.py
│   ├── make_vibrational_corrections.py
│   ├── apply_thermal_corrections.py
│   ├── sample_ternary_isotherm.py
│   └── plot_thermocalc_style_ternary.py
└── LICENSE
```

## Install For Codex

Clone this repository into your Codex skills directory:

```powershell
git clone https://github.com/kk154/PhaseForge-CALPHAD-skill.git C:\Users\15461\.codex\skills\mat-phaseforge-calphad
```

If your Codex skills root is different, replace the destination path with your local skills directory.

## Typical Usage

Ask Codex for a task such as:

```text
Use the mat-phaseforge-calphad skill to calculate a Ni-Al-Cr ternary isothermal phase diagram at 1000 K.
```

The skill guides Codex through run-directory setup, script copying, compound discovery, PhaseForge/ORB calculations, phonon free-energy calculations for relaxed solids, multi-temperature LAMMPS liquid enthalpy calculations for LIQUID SQS structures, final TDB construction, ternary sampling, and plotting. See `SKILL.md` for the complete workflow and command templates.

## Dependencies

The bundled scripts commonly rely on:

- PhaseForge;
- an ORB/MLIP inference environment;
- Python 3;
- pycalphad;
- pymatgen and mp-api;
- phonopy;
- LAMMPS with MLIP/GNNP support for liquid calculations;
- matplotlib, numpy, pandas, scipy;
- optional GPU/CUDA support.

The exact dependencies depend on the selected workflow steps. For example, plotting and pycalphad sampling do not require ORB execution, while Materials Project downloads require `MP_API_KEY`.

## License

This repository is released under the MIT License. See `LICENSE`.
