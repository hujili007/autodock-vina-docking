import argparse
from pathlib import Path

from docking_pipeline import (
    download_file,
    env_script,
    generate_3d_sdf_from_2d,
    log_step,
    pubchem_2d_sdf_url,
    pubchem_3d_sdf_url,
    run_command,
    safe_stem,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Step 3: prepare a PubChem ligand for AutoDock Vina."
    )
    parser.add_argument("--compound-name", required=True, help="Compound name for logs and filenames.")
    parser.add_argument("--pubchem-cid", required=True, help="PubChem compound CID.")
    parser.add_argument("--output-prefix", help="Output prefix. Defaults to sanitized name plus CID.")
    return parser.parse_args()


def main():
    args = parse_args()
    prefix = args.output_prefix or safe_stem(args.compound_name, args.pubchem_cid)
    ligand_sdf = Path(f"{prefix}.sdf")
    ligand_2d_sdf = Path(f"{prefix}_2d.sdf")
    ligand_pdbqt = Path(f"{prefix}.pdbqt")

    log_step("Step 3 prepares the ligand structure and converts it to PDBQT.")
    log_step(f"Compound: {args.compound_name}")
    log_step(f"PubChem CID: {args.pubchem_cid}")

    try:
        url = pubchem_3d_sdf_url(args.pubchem_cid)
        log_step(f"Downloading PubChem 3D SDF: {url}")
        download_file(url, ligand_sdf)
        log_step(f"Saved PubChem 3D SDF: {ligand_sdf}")
    except RuntimeError as exc:
        print(f"[Warning] PubChem 3D SDF was not available: {exc}")
        url = pubchem_2d_sdf_url(args.pubchem_cid)
        log_step(f"Downloading PubChem 2D SDF for local 3D generation: {url}")
        download_file(url, ligand_2d_sdf)
        generate_3d_sdf_from_2d(ligand_2d_sdf, ligand_sdf)
        log_step(f"Generated local 3D SDF with RDKit: {ligand_sdf}")

    meeko_ligand_script = env_script("mk_prepare_ligand.py")
    if not meeko_ligand_script.exists():
        raise FileNotFoundError(f"Meeko ligand script not found: {meeko_ligand_script}")

    run_command(
        [
            str(meeko_ligand_script),
            "-i",
            str(ligand_sdf),
            "-o",
            str(ligand_pdbqt),
        ],
        "Converting ligand SDF to ligand PDBQT with Meeko.",
    )
    print(f"Ligand SDF: {ligand_sdf}")
    print(f"Ligand PDBQT: {ligand_pdbqt}")


if __name__ == "__main__":
    main()
