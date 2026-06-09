import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

from docking_pipeline import compact_timestamp, log_step, run_directory_name, safe_stem


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the complete AutoDock Vina docking workflow."
    )
    parser.add_argument("--target", required=True, help="Target name, e.g. S100A8.")
    parser.add_argument("--pdb-id", required=True, help="PDB ID, e.g. 5HLO.")
    parser.add_argument("--compound-name", required=True, help="Compound name, e.g. Epicatechin.")
    parser.add_argument("--pubchem-cid", required=True, help="PubChem CID, e.g. 72276.")
    parser.add_argument("--runs-dir", default="runs", help="Directory for per-run outputs.")
    parser.add_argument("--results-xlsx", default="docking_results_summary.xlsx", help="Workbook summary path.")
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--n-poses", type=int, default=9)
    parser.add_argument("--skip-pymol", action="store_true", help="Skip figure rendering.")
    return parser.parse_args()


def run_step(command, cwd):
    log_step("Running: " + " ".join(str(part) for part in command))
    subprocess.run(command, cwd=cwd, check=True)


def main():
    args = parse_args()
    project_root = Path.cwd()
    scripts_dir = Path(__file__).resolve().parent
    timestamp = compact_timestamp()
    run_dir = project_root / args.runs_dir / run_directory_name(args.pdb_id, args.pubchem_cid, timestamp)
    run_dir.mkdir(parents=True, exist_ok=True)

    python = sys.executable
    output_prefix = args.pdb_id.upper()
    ligand_prefix = safe_stem(args.compound_name, args.pubchem_cid)

    log_step(f"Run directory: {run_dir}")
    run_step(
        [
            python,
            str(scripts_dir / "PrepareReceptor.py"),
            "--pdb-id",
            args.pdb_id,
            "--target",
            args.target,
            "--output-prefix",
            output_prefix,
        ],
        run_dir,
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
        run_dir,
    )
    run_step(
        [
            python,
            str(scripts_dir / "PrepareLigand.py"),
            "--compound-name",
            args.compound_name,
            "--pubchem-cid",
            args.pubchem_cid,
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
            str(project_root / args.results_xlsx),
            "--compound-name",
            args.compound_name,
            "--target",
            args.target,
            "--pdb-id",
            args.pdb_id.upper(),
            "--pubchem-cid",
            args.pubchem_cid,
            "--exhaustiveness",
            str(args.exhaustiveness),
            "--n-poses",
            str(args.n_poses),
        ],
        run_dir,
    )

    from docking_pipeline import convert_pdbqt_to_pdb, extract_best_pose_pdbqt

    extract_best_pose_pdbqt(run_dir / docked_pose, run_dir / best_pose)
    convert_pdbqt_to_pdb(run_dir / best_pose, run_dir / best_pose_pdb)

    figure_outputs = {}
    pymol_error = None
    if not args.skip_pymol:
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
        "target": args.target,
        "pdb_id": args.pdb_id.upper(),
        "compound_name": args.compound_name,
        "pubchem_cid": args.pubchem_cid,
        "run_dir": str(run_dir),
        "results_xlsx": str(project_root / args.results_xlsx),
        "best_affinity_kcal_per_mol": summary.get("best_affinity_kcal_per_mol"),
        "warnings": warnings,
        "files": {
            "original_pdb": str(run_dir / f"{output_prefix}_original.pdb"),
            "protein_pdb": str(run_dir / f"{output_prefix}_protein.pdb"),
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

    # Convenience copies for users browsing the project root.
    for key in ("figure_a", "figure_b"):
        if key in figure_outputs:
            shutil.copyfile(figure_outputs[key], project_root / Path(figure_outputs[key]).name)

    print("Pipeline complete.")
    print(f"Run manifest: {run_dir / 'manifest.json'}")
    print(f"Results workbook: {project_root / args.results_xlsx}")


if __name__ == "__main__":
    main()
