# AutoDock Vina Docking Skill

English | [中文](./README.zh.md)

AutoDock Vina Docking Skill is a reproducible, script-based workflow for protein-ligand docking. It is designed for learning, small-scale docking experiments, and agent-assisted use in Codex or Claude.

The workflow accepts a protein target/PDB ID and a small molecule name/PubChem CID, then prepares receptor and ligand files, chooses a docking box, runs AutoDock Vina, extracts the best pose, renders PyMOL figures, and writes result summaries.

## Features

- Single protein-ligand docking from command-line arguments.
- CSV batch docking with one or many rows.
- Reuses receptor and docking box preparation for repeated `target+pdbid` pairs in batch mode.
- Downloads ligand structures from PubChem by CID.
- Prefers PubChem 3D SDF and falls back to PubChem 2D SDF plus RDKit 3D generation.
- Generates receptor and ligand PDBQT files.
- Runs AutoDock Vina through the Python API.
- Extracts the best pose as both PDBQT and PDB.
- Generates PyMOL figure a and figure b style images.
- Writes `docking_results_summary.xlsx`, per-run `manifest.json`, and batch summaries.

## Directory Structure

```text
autodock-vina-docking/
├── SKILL.md
├── README.md
├── README.zh.md
├── references/
│   ├── docking_box.md
│   ├── pubchem_ligand.md
│   ├── pymol_figures.md
│   └── workflow.md
└── scripts/
    ├── setup_conda_env.sh
    ├── run_pipeline.py
    ├── run_batch.py
    ├── PrepareReceptor.py
    ├── IdentifyBindingSite.py
    ├── PrepareLigand.py
    ├── RunDocking.py
    ├── AnalyzeDocking.py
    ├── render_pymol_figures.py
    ├── docking_pipeline.py
    └── test_docking_pipeline.py
```

## Prerequisites

Install Miniconda, Anaconda, or Mambaforge first. Then create the recommended `docking` environment from the project root:

```bash
bash scripts/setup_conda_env.sh
```

Manual equivalent:

```bash
conda create -n docking -c conda-forge python=3.11 meeko openbabel vina rdkit prody pymol-open-source pillow
```

If you already have a `docking` environment and only need PyMOL figure rendering:

```bash
conda install -n docking -c conda-forge pymol-open-source pillow
```

Verify the environment:

```bash
conda run -n docking python scripts/run_pipeline.py --help
conda run -n docking python scripts/run_batch.py --help
```

## Single Docking Example

Run from the `autodock-vina-docking/` directory:

```bash
conda run -n docking python scripts/run_pipeline.py \
  --target S100A8 \
  --pdb-id 5HLO \
  --compound-name Epicatechin \
  --pubchem-cid 72276
```

This creates:

```text
runs/<PDBID>_<CID>_<YYYYMMDD_HHMMSS>/
```

The project-level workbook is appended or created:

```text
docking_results_summary.xlsx
```

## Batch CSV Docking Example

Create `ligands.csv` with exactly these headers:

```csv
compound_name,pubchem_cid,target,pdbid
Epicatechin,72276,S100A8,5HLO
DL-Tryptophan,1148,S100A8,5HLO
```

Run:

```bash
conda run -n docking python scripts/run_batch.py --input-csv ligands.csv
```

Batch mode creates:

```text
runs/batch_<YYYYMMDD_HHMMSS>/
```

For repeated `target+pdbid` pairs, receptor preparation and docking box detection are performed once and reused. Each CSV row still gets its own docking directory.

## Output Files

Single-run and per-row batch outputs include:

- Original receptor PDB.
- Protein-only receptor PDB.
- Receptor PDBQT.
- Ligand SDF.
- Ligand PDBQT.
- `vina_box_config.txt`.
- Docked multi-pose PDBQT.
- Best-pose PDBQT.
- Best-pose PDB.
- `figure_a_overall_pose.png`.
- `figure_b_binding_site.png`.
- `pymol_docking_session.pse`.
- `manifest.json`.

Batch-level outputs include:

- `batch_manifest.json`: input CSV, row counts, receptor cache files, per-row manifests, and errors.
- `batch_results_summary.csv`: success or failed status for every row.
- `docking_results_summary.xlsx`: successful docking rows only.

## Using as a Codex or Claude Skill

This repository is structured as an agent skill. To use it with Codex, place the `autodock-vina-docking/` directory under your Codex skills directory, for example:

```text
~/.codex/skills/autodock-vina-docking/
```

Then start a new Codex session and ask naturally, for example:

```text
Use the autodock-vina-docking skill to dock S100A8 / 5HLO with DL-Tryptophan / 1148.
```

For batch docking, provide a CSV file and ask the agent to run all rows.

## Notes on Interpretation

Docking scores and poses are computational predictions. They are useful for hypothesis generation, workflow learning, and prioritization, but they are not experimental evidence of binding or biological activity.

When no credible co-crystal ligand is detected, the workflow uses blind docking and writes a warning to the config, logs, summaries, and manifests. Blind docking can produce false-positive poses and should be interpreted cautiously.

## Acknowledgements

This workflow builds on open-source scientific software, including AutoDock Vina, RDKit, Meeko, Open Babel, ProDy, PubChem PUG REST, and PyMOL open source. Please cite the relevant tools and databases when using this workflow in research outputs.

## License

This project is licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/). Non-commercial sharing and adaptation are welcome, with attribution.
