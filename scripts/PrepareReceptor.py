import argparse
import shutil
from pathlib import Path

import prody

from docking_pipeline import env_script, log_step, run_command


def parse_args():
    parser = argparse.ArgumentParser(
        description="Step 1: prepare a protein receptor for AutoDock Vina."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pdb-id", help="PDB ID to download, for example 5HLO.")
    group.add_argument("--pdb-file", help="Local original PDB file.")
    parser.add_argument("--target", default="target", help="Target name for teaching logs.")
    parser.add_argument("--output-prefix", help="Output prefix. Defaults to PDB ID or input stem.")
    return parser.parse_args()


def main():
    args = parse_args()
    output_prefix = args.output_prefix

    log_step("Step 1 prepares a protein-only receptor and converts it to PDBQT.")
    log_step(f"Target label: {args.target}")

    if args.pdb_id:
        pdb_id = args.pdb_id.upper()
        output_prefix = output_prefix or pdb_id
        log_step(f"Downloading original PDB structure for {pdb_id}.")
        fetched_path = prody.fetchPDB(pdb_id, compressed=False)
        if not fetched_path:
            raise RuntimeError(f"Could not download PDB ID {pdb_id}.")
        original_pdb = Path(str(fetched_path))
    else:
        original_pdb = Path(args.pdb_file)
        if not original_pdb.exists():
            raise FileNotFoundError(f"Input PDB file does not exist: {original_pdb}")
        output_prefix = output_prefix or original_pdb.stem

    canonical_original = Path(f"{output_prefix}_original.pdb")
    if original_pdb.resolve() != canonical_original.resolve():
        shutil.copyfile(original_pdb, canonical_original)
        log_step(f"Saved a stable copy of the original PDB: {canonical_original}")
    else:
        log_step(f"Using original PDB: {canonical_original}")

    log_step("Extracting protein atoms only for receptor preparation.")
    structure = prody.parsePDB(str(canonical_original))
    protein = structure.select("protein")
    if protein is None:
        raise ValueError(f"No protein atoms found in {canonical_original}")

    protein_pdb = Path(f"{output_prefix}_protein.pdb")
    prody.writePDB(str(protein_pdb), protein)
    print(f"Protein atoms: {protein.numAtoms()}")
    print(f"Clean protein PDB: {protein_pdb}")

    receptor_prefix = f"{output_prefix}_receptor"
    receptor_pdbqt = Path(f"{receptor_prefix}.pdbqt")
    meeko_receptor_script = env_script("mk_prepare_receptor.py")
    if not meeko_receptor_script.exists():
        raise FileNotFoundError(f"Meeko receptor script not found: {meeko_receptor_script}")

    run_command(
        [
            str(meeko_receptor_script),
            "-i",
            str(protein_pdb),
            "-o",
            receptor_prefix,
            "-p",
        ],
        "Converting protein PDB to receptor PDBQT with Meeko.",
    )
    print(f"Receptor PDBQT: {receptor_pdbqt}")


if __name__ == "__main__":
    main()
