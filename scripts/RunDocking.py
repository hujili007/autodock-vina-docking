import argparse
from pathlib import Path

from vina import Vina

from docking_pipeline import append_docking_result_xlsx, formatted_timestamp, log_step, read_vina_config, write_summary_json


def parse_args():
    parser = argparse.ArgumentParser(
        description="Step 4: run AutoDock Vina docking with the Python API."
    )
    parser.add_argument("--receptor", required=True, help="Prepared receptor PDBQT.")
    parser.add_argument("--ligand", required=True, help="Prepared ligand PDBQT.")
    parser.add_argument("--box-config", default="vina_box_config.txt", help="Vina box config from Step 2.")
    parser.add_argument("--out", default="docked_pose.pdbqt", help="Output docked pose PDBQT.")
    parser.add_argument("--summary", default="docking_summary.json", help="Output docking summary JSON.")
    parser.add_argument("--results-xlsx", default="docking_results_summary.xlsx", help="Excel summary workbook to create or append.")
    parser.add_argument("--compound-name", default="", help="Compound name for the Excel summary.")
    parser.add_argument("--target", default="", help="Target name for the Excel summary.")
    parser.add_argument("--pdb-id", default="", help="PDB ID for the Excel summary.")
    parser.add_argument("--pubchem-cid", default="", help="PubChem CID for the Excel summary.")
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--n-poses", type=int, default=9)
    return parser.parse_args()


def main():
    args = parse_args()
    receptor = Path(args.receptor)
    ligand = Path(args.ligand)
    box_config = Path(args.box_config)
    for path, label in [(receptor, "receptor"), (ligand, "ligand"), (box_config, "box config")]:
        if not path.exists():
            raise FileNotFoundError(f"Missing {label} file: {path}")

    log_step("Step 4 runs AutoDock Vina docking.")
    log_step(f"Receptor: {receptor}")
    log_step(f"Ligand: {ligand}")
    log_step(f"Box config: {box_config}")

    config = read_vina_config(box_config)
    center = [float(config["center_x"]), float(config["center_y"]), float(config["center_z"])]
    box_size = [float(config["size_x"]), float(config["size_y"]), float(config["size_z"])]
    log_step(f"Docking center: {center}")
    log_step(f"Docking box size: {box_size}")

    vina = Vina(sf_name="vina")
    vina.set_receptor(str(receptor))
    vina.set_ligand_from_file(str(ligand))
    vina.compute_vina_maps(center=center, box_size=box_size)

    vina.dock(exhaustiveness=args.exhaustiveness, n_poses=args.n_poses)
    vina.write_poses(args.out, n_poses=args.n_poses, overwrite=True)
    energies = vina.energies(n_poses=args.n_poses)
    best_affinity = float(energies[0][0]) if len(energies) else None

    summary = {
        "receptor": str(receptor),
        "ligand": str(ligand),
        "box_config": str(box_config),
        "box_mode": config.get("mode", "unknown"),
        "box_source": config.get("source", "unknown"),
        "box_warning": config.get("warning"),
        "center": center,
        "box_size": box_size,
        "exhaustiveness": args.exhaustiveness,
        "n_poses": args.n_poses,
        "pose_file": args.out,
        "best_affinity_kcal_per_mol": best_affinity,
        "compound_name": args.compound_name,
        "target": args.target,
        "pdb_id": args.pdb_id,
        "pubchem_cid": args.pubchem_cid,
    }
    write_summary_json(Path(args.summary), summary)

    result_row = {
        "Compound": args.compound_name,
        "Target": args.target,
        "PDBID": args.pdb_id,
        "CID": args.pubchem_cid,
        "Affinity(kcal/mol)": round(best_affinity, 3) if best_affinity is not None else "",
        "Timestamp": formatted_timestamp(),
        "Center_X": center[0],
        "Center_Y": center[1],
        "Center_Z": center[2],
        "Size_X": box_size[0],
        "Size_Y": box_size[1],
        "Size_Z": box_size[2],
    }
    append_docking_result_xlsx(Path(args.results_xlsx), result_row)

    print(f"Best affinity: {best_affinity} kcal/mol")
    print(f"Docked poses: {args.out}")
    print(f"Docking summary: {args.summary}")
    print(f"Excel results summary: {args.results_xlsx}")


if __name__ == "__main__":
    main()
