import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

from docking_pipeline import (
    BatchDockingRecord,
    compact_timestamp,
    convert_pdbqt_to_pdb,
    extract_best_pose_pdbqt,
    log_step,
    parse_batch_csv,
    receptor_cache_key,
    safe_stem,
    write_batch_results_csv,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run AutoDock Vina docking for every row in a CSV file."
    )
    parser.add_argument("--input-csv", required=True, help="CSV with compound_name,pubchem_cid,target,pdbid.")
    parser.add_argument("--runs-dir", default="runs", help="Directory for batch outputs.")
    parser.add_argument("--results-xlsx", default="docking_results_summary.xlsx", help="Workbook summary path.")
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--n-poses", type=int, default=9)
    parser.add_argument("--skip-pymol", action="store_true", help="Skip figure rendering.")
    return parser.parse_args()


def run_step(command, cwd):
    log_step("Running: " + " ".join(str(part) for part in command))
    subprocess.run(command, cwd=cwd, check=True)


def prepare_receptor(record: BatchDockingRecord, receptor_dir: Path, scripts_dir: Path, python: str):
    output_prefix = record.pdbid.upper()
    receptor_dir.mkdir(parents=True, exist_ok=True)
    run_step(
        [
            python,
            str(scripts_dir / "PrepareReceptor.py"),
            "--pdb-id",
            output_prefix,
            "--target",
            record.target,
            "--output-prefix",
            output_prefix,
        ],
        receptor_dir,
    )
    run_step(
        [
            python,
            str(scripts_dir / "IdentifyBindingSite.py"),
            "--pdb",
            f"{output_prefix}_original.pdb",
            "--output",
            "vina_box_config.txt",
        ],
        receptor_dir,
    )
    return {
        "original_pdb": str(receptor_dir / f"{output_prefix}_original.pdb"),
        "protein_pdb": str(receptor_dir / f"{output_prefix}_protein.pdb"),
        "receptor_pdbqt": str(receptor_dir / f"{output_prefix}_receptor.pdbqt"),
        "box_config": str(receptor_dir / "vina_box_config.txt"),
    }


def symlink_or_copy(source: Path, dest: Path):
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    try:
        dest.symlink_to(source)
    except OSError:
        shutil.copyfile(source, dest)


def run_docking_record(
    record: BatchDockingRecord,
    run_dir: Path,
    receptor_files: dict,
    scripts_dir: Path,
    python: str,
    project_root: Path,
    results_xlsx: Path,
    exhaustiveness: int,
    n_poses: int,
    skip_pymol: bool,
) -> dict:
    run_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = record.pdbid.upper()
    ligand_prefix = safe_stem(record.compound_name, record.pubchem_cid)

    symlink_or_copy(Path(receptor_files["original_pdb"]), run_dir / f"{output_prefix}_original.pdb")
    symlink_or_copy(Path(receptor_files["receptor_pdbqt"]), run_dir / f"{output_prefix}_receptor.pdbqt")
    symlink_or_copy(Path(receptor_files["box_config"]), run_dir / "vina_box_config.txt")

    run_step(
        [
            python,
            str(scripts_dir / "PrepareLigand.py"),
            "--compound-name",
            record.compound_name,
            "--pubchem-cid",
            record.pubchem_cid,
            "--output-prefix",
            ligand_prefix,
        ],
        run_dir,
    )

    docked_pose = f"{ligand_prefix}_{output_prefix}_docked.pdbqt"
    best_pose = f"{ligand_prefix}_{output_prefix}_best_pose.pdbqt"
    best_pose_pdb = f"{ligand_prefix}_{output_prefix}_best_pose.pdb"
    summary_json = f"{ligand_prefix}_{output_prefix}_summary.json"
    run_step(
        [
            python,
            str(scripts_dir / "RunDocking.py"),
            "--receptor",
            f"{output_prefix}_receptor.pdbqt",
            "--ligand",
            f"{ligand_prefix}.pdbqt",
            "--box-config",
            "vina_box_config.txt",
            "--out",
            docked_pose,
            "--summary",
            summary_json,
            "--results-xlsx",
            str(results_xlsx),
            "--compound-name",
            record.compound_name,
            "--target",
            record.target,
            "--pdb-id",
            output_prefix,
            "--pubchem-cid",
            record.pubchem_cid,
            "--exhaustiveness",
            str(exhaustiveness),
            "--n-poses",
            str(n_poses),
        ],
        run_dir,
    )

    extract_best_pose_pdbqt(run_dir / docked_pose, run_dir / best_pose)
    convert_pdbqt_to_pdb(run_dir / best_pose, run_dir / best_pose_pdb)

    figure_outputs = {}
    pymol_error = None
    if not skip_pymol:
        try:
            run_step(
                [
                    python,
                    str(scripts_dir / "render_pymol_figures.py"),
                    "--protein",
                    f"{output_prefix}_original.pdb",
                    "--ligand",
                    best_pose_pdb,
                    "--figure-a",
                    "figure_a_overall_pose.png",
                    "--figure-b",
                    "figure_b_binding_site.png",
                    "--session",
                    "pymol_docking_session.pse",
                ],
                run_dir,
            )
            figure_outputs = {
                "figure_a": str(run_dir / "figure_a_overall_pose.png"),
                "figure_b": str(run_dir / "figure_b_binding_site.png"),
                "pymol_session": str(run_dir / "pymol_docking_session.pse"),
            }
        except Exception as exc:
            pymol_error = str(exc)
            print(f"[Warning] PyMOL rendering failed, docking results were kept: {pymol_error}")

    summary_path = run_dir / summary_json
    summary = json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.exists() else {}
    pymol_metadata_path = run_dir / "pymol_figures_metadata.json"
    pymol_metadata = (
        json.loads(pymol_metadata_path.read_text(encoding="utf-8"))
        if pymol_metadata_path.exists()
        else None
    )
    warnings = []
    if summary.get("box_warning"):
        warnings.append(summary["box_warning"])
    if pymol_error:
        warnings.append(f"PyMOL rendering failed: {pymol_error}")

    manifest = {
        "row_index": record.row_index,
        "target": record.target,
        "pdb_id": output_prefix,
        "compound_name": record.compound_name,
        "pubchem_cid": record.pubchem_cid,
        "run_dir": str(run_dir),
        "results_xlsx": str(project_root / results_xlsx),
        "best_affinity_kcal_per_mol": summary.get("best_affinity_kcal_per_mol"),
        "warnings": warnings,
        "files": {
            "original_pdb": str(run_dir / f"{output_prefix}_original.pdb"),
            "receptor_pdbqt": str(run_dir / f"{output_prefix}_receptor.pdbqt"),
            "ligand_sdf": str(run_dir / f"{ligand_prefix}.sdf"),
            "ligand_pdbqt": str(run_dir / f"{ligand_prefix}.pdbqt"),
            "box_config": str(run_dir / "vina_box_config.txt"),
            "docked_pose_pdbqt": str(run_dir / docked_pose),
            "best_pose_pdbqt": str(run_dir / best_pose),
            "best_pose_pdb": str(run_dir / best_pose_pdb),
            **figure_outputs,
        },
        "summary": summary,
        "pymol_metadata": pymol_metadata,
        "pymol_error": pymol_error,
    }
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def main():
    args = parse_args()
    project_root = Path.cwd()
    scripts_dir = Path(__file__).resolve().parent
    python = sys.executable
    records = parse_batch_csv(Path(args.input_csv))
    timestamp = compact_timestamp()
    batch_dir = project_root / args.runs_dir / f"batch_{timestamp}"
    receptors_dir = batch_dir / "receptors"
    dockings_dir = batch_dir / "dockings"
    batch_dir.mkdir(parents=True, exist_ok=True)

    receptor_results = {}
    receptor_errors = {}
    batch_rows = []
    row_results = []

    for key in dict.fromkeys(receptor_cache_key(record) for record in records):
        representative = next(record for record in records if receptor_cache_key(record) == key)
        receptor_dir = receptors_dir / representative.pdbid.upper()
        try:
            receptor_results[key] = prepare_receptor(representative, receptor_dir, scripts_dir, python)
        except Exception as exc:
            receptor_errors[key] = str(exc)
            print(f"[Warning] Receptor preparation failed for {key}: {exc}")

    for record in records:
        key = receptor_cache_key(record)
        run_dir = dockings_dir / f"{record.row_index:03d}_{record.pdbid.upper()}_{record.pubchem_cid}"
        if key in receptor_errors:
            error = receptor_errors[key]
            batch_rows.append(
                {
                    "status": "failed",
                    "error": error,
                    "compound_name": record.compound_name,
                    "pubchem_cid": record.pubchem_cid,
                    "target": record.target,
                    "pdbid": record.pdbid.upper(),
                    "run_dir": str(run_dir),
                    "best_affinity_kcal_per_mol": "",
                }
            )
            row_results.append({"row_index": record.row_index, "status": "failed", "error": error})
            continue

        try:
            manifest = run_docking_record(
                record,
                run_dir,
                receptor_results[key],
                scripts_dir,
                python,
                project_root,
                project_root / args.results_xlsx,
                args.exhaustiveness,
                args.n_poses,
                args.skip_pymol,
            )
            batch_rows.append(
                {
                    "status": "success",
                    "error": "",
                    "compound_name": record.compound_name,
                    "pubchem_cid": record.pubchem_cid,
                    "target": record.target,
                    "pdbid": record.pdbid.upper(),
                    "run_dir": str(run_dir),
                    "best_affinity_kcal_per_mol": manifest.get("best_affinity_kcal_per_mol", ""),
                }
            )
            row_results.append(
                {
                    "row_index": record.row_index,
                    "status": "success",
                    "manifest": str(run_dir / "manifest.json"),
                }
            )
        except Exception as exc:
            error = str(exc)
            print(f"[Warning] Docking failed for row {record.row_index}: {error}")
            batch_rows.append(
                {
                    "status": "failed",
                    "error": error,
                    "compound_name": record.compound_name,
                    "pubchem_cid": record.pubchem_cid,
                    "target": record.target,
                    "pdbid": record.pdbid.upper(),
                    "run_dir": str(run_dir),
                    "best_affinity_kcal_per_mol": "",
                }
            )
            row_results.append({"row_index": record.row_index, "status": "failed", "error": error})

    write_batch_results_csv(batch_dir / "batch_results_summary.csv", batch_rows)
    success_count = sum(1 for row in batch_rows if row["status"] == "success")
    failed_count = sum(1 for row in batch_rows if row["status"] == "failed")
    batch_manifest = {
        "input_csv": str(Path(args.input_csv).resolve()),
        "batch_dir": str(batch_dir),
        "total_rows": len(records),
        "success_count": success_count,
        "failed_count": failed_count,
        "results_xlsx": str(project_root / args.results_xlsx),
        "receptors": {
            f"{target}|{pdbid}": value
            for (target, pdbid), value in receptor_results.items()
        },
        "receptor_errors": {
            f"{target}|{pdbid}": value
            for (target, pdbid), value in receptor_errors.items()
        },
        "rows": row_results,
    }
    (batch_dir / "batch_manifest.json").write_text(
        json.dumps(batch_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("Batch pipeline complete.")
    print(f"Batch manifest: {batch_dir / 'batch_manifest.json'}")
    print(f"Batch summary: {batch_dir / 'batch_results_summary.csv'}")
    print(f"Results workbook: {project_root / args.results_xlsx}")
    print(f"Success: {success_count}; Failed: {failed_count}")


if __name__ == "__main__":
    main()
