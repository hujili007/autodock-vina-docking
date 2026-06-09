# AutoDock Vina Workflow

This skill is organized as single-run and batch entrypoints plus reusable step scripts. Use `scripts/run_pipeline.py` for one ligand and `scripts/run_batch.py` for CSV batches.

## One-Command Pipeline

```bash
python autodock-vina-docking/scripts/run_pipeline.py \
  --target S100A8 \
  --pdb-id 5HLO \
  --compound-name Epicatechin \
  --pubchem-cid 72276
```

The pipeline creates `runs/<PDBID>_<CID>_<timestamp>/` and writes `manifest.json` inside that directory.

## CSV Batch Pipeline

Batch CSV files must use these headers:

```csv
compound_name,pubchem_cid,target,pdbid
Epicatechin,72276,S100A8,5HLO
```

Run:

```bash
conda run -n docking python autodock-vina-docking/scripts/run_batch.py --input-csv ligands.csv
```

Batch mode creates `runs/batch_<timestamp>/`. Receptor preparation and binding box detection are reused for identical `target+pdbid` pairs under `receptors/<PDBID>/`; each ligand row gets its own directory under `dockings/<row_index>_<PDBID>_<CID>/`.

Batch outputs:

- `batch_manifest.json` with input CSV, row counts, receptor cache files, per-row manifests, and errors.
- `batch_results_summary.csv` with success or failed status for every row.
- `docking_results_summary.xlsx`, appended only for successful docking rows.

## Stages

1. Receptor preparation
   - Downloads the original PDB when a PDB ID is provided.
   - Writes `<PDBID>_original.pdb`.
   - Writes protein-only `<PDBID>_protein.pdb`.
   - Converts the protein receptor to `<PDBID>_receptor.pdbqt`.

2. Binding box selection
   - Reads the original PDB, not the cleaned protein-only receptor.
   - Uses a credible co-crystal ligand when present.
   - Falls back to blind docking when no credible ligand is present.
   - Writes `vina_box_config.txt`.

3. Ligand preparation
   - Uses PubChem CID as the structure source.
   - Downloads a PubChem 3D SDF when available.
   - Falls back to PubChem 2D SDF plus RDKit 3D generation.
   - Converts to ligand PDBQT.

4. Docking
   - Runs AutoDock Vina with receptor PDBQT, ligand PDBQT, and `vina_box_config.txt`.
   - Writes multi-pose docked PDBQT.
   - Writes a JSON summary and appends `docking_results_summary.xlsx`.

5. Post-processing and figures
   - Extracts `MODEL 1` as the best-pose PDBQT.
   - Converts the best-pose PDBQT to PDB for visualization.
   - Renders figure a, figure b, and a PyMOL session when PyMOL is available.

## Failure Behavior

The pipeline should fail clearly for missing receptor, invalid PubChem CID, invalid box config, or docking execution errors. PyMOL rendering is non-blocking: if figure generation fails, docking outputs are kept and `manifest.json` records the error.

In batch mode, row-level failures are recorded and later rows continue. If receptor preparation fails for a shared `target+pdbid`, all rows using that receptor are marked failed while other receptors continue.
