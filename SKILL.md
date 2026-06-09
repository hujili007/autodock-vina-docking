---
name: autodock-vina-docking
description: Use when running AutoDock Vina protein-ligand docking from a target/PDB ID and PubChem CID, especially when results need Excel summaries, best-pose files, and PyMOL figure a/b rendering.
---

# AutoDock Vina Docking

Use this skill when the user wants to run a learning-oriented but complete AutoDock Vina docking workflow for one protein target and one small molecule, or for a CSV batch of protein-ligand pairs. Single-run inputs are:

- Target name, for example `S100A8`
- PDB ID or a local receptor PDB, for example `5HLO`
- Compound name, for example `Epicatechin`
- PubChem CID, for example `72276`

## Prerequisites

This skill expects a conda environment with AutoDock Vina, Meeko, Open Babel, RDKit, ProDy, PyMOL open source, and Pillow. The recommended environment name is `docking`.

If the environment does not exist, create it from the project root:

```bash
bash autodock-vina-docking/scripts/setup_conda_env.sh
```

To use a different environment name:

```bash
bash autodock-vina-docking/scripts/setup_conda_env.sh my-docking-env
```

Manual equivalent:

```bash
conda create -n docking -c conda-forge python=3.11 meeko openbabel vina rdkit prody pymol-open-source pillow
```

If the `docking` environment already exists and only PyMOL rendering is missing:

```bash
conda install -n docking -c conda-forge pymol-open-source pillow
```

Run commands with the environment's Python, for example:

```bash
conda run -n docking python autodock-vina-docking/scripts/run_pipeline.py --help
```

## Main Command

## Batch Command

Use a CSV with exactly these standard headers:

```csv
compound_name,pubchem_cid,target,pdbid
Epicatechin,72276,S100A8,5HLO
```

Run from the project root:

```bash
conda run -n docking python autodock-vina-docking/scripts/run_batch.py \
  --input-csv ligands.csv
```

Batch mode creates:

```text
runs/batch_<YYYYMMDD_HHMMSS>/
```

Outputs include:

- `batch_manifest.json`
- `batch_results_summary.csv`
- Shared receptor and box files under `receptors/<PDBID>/`
- Per-row docking outputs under `dockings/<row_index>_<PDBID>_<CID>/`

Rows with errors are recorded in `batch_results_summary.csv`; later rows continue running.

## Single Command

Run from the project root:

```bash
conda run -n docking python autodock-vina-docking/scripts/run_pipeline.py \
  --target S100A8 \
  --pdb-id 5HLO \
  --compound-name Epicatechin \
  --pubchem-cid 72276
```

Each run creates an isolated directory:

```text
runs/<PDBID>_<CID>_<YYYYMMDD_HHMMSS>/
```

The project root summary workbook is appended or created at:

```text
docking_results_summary.xlsx
```

## Outputs

The per-run `manifest.json` is the source of truth for output paths. A successful full run should include:

- Original receptor PDB
- Protein-only receptor PDB
- Receptor PDBQT
- Ligand SDF
- Ligand PDBQT
- `vina_box_config.txt`
- Docked multi-pose PDBQT
- Best-pose PDBQT
- Best-pose PDB
- `figure_a_overall_pose.png`
- `figure_b_binding_site.png`
- `pymol_docking_session.pse`

The figure PNGs are also copied to the project root for quick inspection.

## Workflow Rules

- Prefer PubChem CID as the ligand structure key. Use the compound name mainly for logs, summaries, and filenames.
- Prefer PubChem 3D SDF. If unavailable, fall back to PubChem 2D SDF plus local RDKit 3D generation.
- Use ligand-based docking box only when the receptor has a credible non-additive co-crystal ligand.
- If no credible ligand is found, continue with blind docking and write a warning to the config, logs, summary, and manifest.
- If PyMOL rendering fails, preserve docking outputs and record the rendering error in `manifest.json`.

## References

Read these only as needed:

- `references/workflow.md` for the step-by-step workflow.
- `references/docking_box.md` for box selection behavior and caveats.
- `references/pubchem_ligand.md` for ligand preparation details.
- `references/pymol_figures.md` for figure a and figure b rendering recipes.
