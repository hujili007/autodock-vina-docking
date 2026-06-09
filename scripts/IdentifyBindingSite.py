import argparse
from pathlib import Path

import prody

from docking_pipeline import (
    BoxResult,
    compute_box_from_coords,
    ligand_candidates_from_selection,
    log_step,
    write_vina_config,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Step 2: identify a docking box for AutoDock Vina."
    )
    parser.add_argument("--pdb", required=True, help="Original PDB file with ligands/additives preserved.")
    parser.add_argument("--output", default="vina_box_config.txt", help="Output Vina box config.")
    parser.add_argument("--padding", type=float, default=5.0, help="Padding around a co-crystallized ligand.")
    parser.add_argument("--blind-padding", type=float, default=5.0, help="Padding around all protein atoms.")
    parser.add_argument("--large-box-threshold", type=float, default=40.0, help="Warn if any box axis exceeds this size.")
    return parser.parse_args()


def main():
    args = parse_args()
    pdb_path = Path(args.pdb)
    if not pdb_path.exists():
        raise FileNotFoundError(f"Original PDB file does not exist: {pdb_path}")

    log_step("Step 2 determines the docking search box.")
    log_step(f"Reading original PDB: {pdb_path}")
    structure = prody.parsePDB(str(pdb_path))

    hetero = structure.select("hetero and not water")
    candidates = ligand_candidates_from_selection(hetero)
    if candidates:
        ligand = candidates[0]
        center, size = compute_box_from_coords(ligand.coords, args.padding)
        source = f"co-crystallized ligand candidate {ligand.label}"
        result = BoxResult(mode="ligand_based", center=center, size=size, source=source)
        log_step(f"Selected ligand candidate: {ligand.label} ({ligand.heavy_atoms} heavy atoms).")
    else:
        protein = structure.select("protein")
        if protein is None:
            raise ValueError(f"No protein atoms found in {pdb_path}")
        center, size = compute_box_from_coords(protein.getCoords(), args.blind_padding)
        warning = (
            "No credible drug-like co-crystallized ligand was found. "
            "Using blind docking over all protein chains; interpret results cautiously."
        )
        result = BoxResult(
            mode="blind_docking",
            center=center,
            size=size,
            source="all protein chains",
            warning=warning,
        )
        print(f"[Warning] {warning}")

    if (result.size > args.large_box_threshold).any():
        large_warning = (
            f"At least one box axis is larger than {args.large_box_threshold:.1f} A. "
            "Vina may be slower and may produce more false-positive poses."
        )
        print(f"[Warning] {large_warning}")
        result.warning = f"{result.warning}; {large_warning}" if result.warning else large_warning

    write_vina_config(result, Path(args.output))
    print(f"Binding site mode: {result.mode}")
    print(f"Binding site source: {result.source}")
    print(f"Center: {result.center}")
    print(f"Size: {result.size}")
    print(f"Vina box config: {args.output}")


if __name__ == "__main__":
    main()
