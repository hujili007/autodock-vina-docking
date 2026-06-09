import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Render docking figure A and B with PyMOL.")
    parser.add_argument("--protein", required=True, help="Original protein PDB file.")
    parser.add_argument("--ligand", required=True, help="Best pose ligand PDBQT/PDB file.")
    parser.add_argument("--figure-a", default="figure_a_overall_pose.png")
    parser.add_argument("--figure-b", default="figure_b_binding_site.png")
    parser.add_argument("--session", default="pymol_docking_session.pse")
    parser.add_argument("--metadata", default="pymol_figures_metadata.json")
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=1200)
    return parser.parse_args()


def find_pymol_executable() -> str:
    candidates = [
        os.environ.get("PYMOL_EXE"),
        str(Path(sys.executable).with_name("pymol")),
        shutil.which("pymol"),
        "/Applications/PyMOL.app/Contents/bin/pymol",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return candidate
    raise RuntimeError(
        "PyMOL executable was not found. Set PYMOL_EXE or install PyMOL.app."
    )


def overlay_binding_site_labels(figure_b: Path, metadata_path: Path) -> None:
    try:
        from PIL import Image, ImageDraw, ImageFont
    except ImportError:
        print("[Warning] Pillow is not available; skipping 2D residue label overlay.")
        return

    figure_b = Path(figure_b)
    metadata_path = Path(metadata_path)
    if not figure_b.exists() or not metadata_path.exists():
        return

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    labels = metadata.get("labeled_residues", [])
    if not labels:
        return

    image = Image.open(figure_b).convert("RGBA")
    draw = ImageDraw.Draw(image)
    width, height = image.size
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 34)
    except OSError:
        font = ImageFont.load_default()

    slots = [
        ((0.09, 0.22), (0.38, 0.39)),
        ((0.72, 0.23), (0.59, 0.39)),
        ((0.09, 0.50), (0.39, 0.56)),
        ((0.73, 0.53), (0.60, 0.57)),
        ((0.43, 0.84), (0.49, 0.67)),
    ]
    for label, (text_xy, anchor_xy) in zip(labels[: len(slots)], slots):
        text_x = int(width * text_xy[0])
        text_y = int(height * text_xy[1])
        anchor_x = int(width * anchor_xy[0])
        anchor_y = int(height * anchor_xy[1])
        bbox = draw.textbbox((text_x, text_y), label, font=font)
        pad = 4
        draw.rectangle(
            (bbox[0] - pad, bbox[1] - pad, bbox[2] + pad, bbox[3] + pad),
            fill=(255, 255, 255, 210),
        )
        line_start_x = bbox[2] + 4 if text_x < anchor_x else bbox[0] - 4
        line_start_y = (bbox[1] + bbox[3]) // 2
        draw.line((line_start_x, line_start_y, anchor_x, anchor_y), fill=(80, 80, 80, 210), width=2)
        draw.text((text_x, text_y), label, fill=(0, 0, 0, 255), font=font)

    image.convert("RGB").save(figure_b)


def render_script(args, protein: Path, ligand: Path) -> str:
    return textwrap.dedent(
        f"""
        import json
        import math
        from pathlib import Path
        from pymol import cmd

        protein = {json.dumps(str(protein.resolve()))}
        ligand = {json.dumps(str(ligand.resolve()))}
        figure_a = {json.dumps(str(Path(args.figure_a).resolve()))}
        figure_b = {json.dumps(str(Path(args.figure_b).resolve()))}
        session = {json.dumps(str(Path(args.session).resolve()))}
        metadata_path = {json.dumps(str(Path(args.metadata).resolve()))}
        width = {args.width}
        height = {args.height}

        cmd.load(protein, "protein_original")
        cmd.load(ligand, "best_pose")
        if cmd.count_atoms("protein_original") == 0:
            raise RuntimeError("Protein file loaded with zero atoms.")
        if cmd.count_atoms("best_pose") == 0:
            raise RuntimeError("Ligand file loaded with zero atoms.")

        cmd.hide("everything")
        cmd.bg_color("white")
        cmd.set_color("icam_ligand", [0.45, 0.78, 0.55])
        cmd.set_color("icam_ligand_dark", [0.24, 0.58, 0.36])
        cmd.set_color("icam_cartoon_pink", [1.00, 0.70, 0.82])
        cmd.set_color("icam_residue_pink", [1.00, 0.56, 0.72])
        cmd.set_color("icam_surface", [0.93, 0.93, 0.92])
        cmd.create("protein_cartoon", "protein_original and polymer.protein")
        cmd.create("protein_surface", "protein_original and polymer.protein")
        cmd.disable("protein_original")

        cmd.show("cartoon", "protein_cartoon")
        cmd.color("icam_cartoon_pink", "protein_cartoon")
        cmd.set("cartoon_transparency", 0.0, "protein_cartoon")
        cmd.show("surface", "protein_surface")
        cmd.color("icam_surface", "protein_surface")
        cmd.set("transparency", 0.72, "protein_surface")
        cmd.show("sticks", "best_pose")
        cmd.color("icam_ligand", "best_pose")
        cmd.set("stick_radius", 0.32, "best_pose")
        cmd.hide("nonbonded")
        cmd.hide("spheres")
        cmd.hide("labels")
        cmd.hide("dashes")
        cmd.set("antialias", 2)
        cmd.set("ambient", 0.45)
        cmd.set("specular", 0.2)
        cmd.set("depth_cue", 0)
        cmd.set("ray_opaque_background", 1)
        cmd.orient("protein_cartoon")
        cmd.zoom("protein_cartoon or best_pose", 4)
        cmd.png(figure_a, width=width, height=height, dpi=150)

        cmd.hide("surface", "protein_surface")
        cmd.show("cartoon", "protein_cartoon")
        cmd.color("gray85", "protein_cartoon")
        cmd.set("cartoon_transparency", 0.0, "protein_cartoon")
        cmd.select("pocket_residues", "byres (protein_cartoon within 4.0 of best_pose)")
        pocket_atom_count = cmd.count_atoms("pocket_residues")
        cmd.show("sticks", "pocket_residues")
        cmd.color("icam_residue_pink", "pocket_residues")
        cmd.show("sticks", "best_pose")
        cmd.color("icam_ligand", "best_pose")
        cmd.set("stick_radius", 0.34, "best_pose")
        cmd.delete("polar_contacts")
        cmd.distance("polar_contacts", "best_pose", "pocket_residues", cutoff=3.5, mode=2, label=0)
        polar_contacts_present = "polar_contacts" in cmd.get_names("objects")
        cmd.color("yellow", "polar_contacts")
        cmd.set("dash_width", 3.0, "polar_contacts")
        cmd.set("dash_gap", 0.25, "polar_contacts")
        cmd.hide("labels")
        ligand_atoms = cmd.get_model("best_pose").atom
        ligand_center = [
            sum(atom.coord[i] for atom in ligand_atoms) / len(ligand_atoms)
            for i in range(3)
        ]
        label_candidates = []
        seen_residues = set()
        for atom in cmd.get_model("pocket_residues and name CA").atom:
            key = (atom.chain, atom.resi, atom.resn)
            if key in seen_residues:
                continue
            seen_residues.add(key)
            distance_to_ligand = math.sqrt(sum((atom.coord[i] - ligand_center[i]) ** 2 for i in range(3)))
            label_candidates.append((distance_to_ligand, atom))
        label_candidates.sort(key=lambda item: item[0])
        labeled_residue_count = min(5, len(label_candidates))
        labeled_residues = []
        for idx, (_, atom) in enumerate(label_candidates[:labeled_residue_count]):
            labeled_residues.append("%s-%s" % (atom.resn, atom.resi))
        cmd.zoom("best_pose or pocket_residues", 7)
        cmd.png(figure_b, width=width, height=height, dpi=150)
        cmd.save(session)

        metadata = {{
            "figure_a": str(Path(figure_a).resolve()),
            "figure_b": str(Path(figure_b).resolve()),
            "session": str(Path(session).resolve()),
            "pocket_atom_count": pocket_atom_count,
            "labeled_residue_count": labeled_residue_count,
            "labeled_residues": labeled_residues,
            "polar_contacts_object_present": polar_contacts_present,
            "polar_contacts_count": None if polar_contacts_present else 0,
        }}
        Path(metadata_path).write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\\n", encoding="utf-8")
        print(json.dumps(metadata, indent=2, sort_keys=True))
        cmd.quit()
        """
    ).strip()


def main():
    args = parse_args()
    protein = Path(args.protein)
    ligand = Path(args.ligand)
    if not protein.exists():
        raise FileNotFoundError(f"Protein file does not exist: {protein}")
    if not ligand.exists():
        raise FileNotFoundError(f"Ligand file does not exist: {ligand}")

    pymol_exe = find_pymol_executable()
    with tempfile.TemporaryDirectory() as tmp:
        script_path = Path(tmp) / "render_docking_figures.py"
        script_path.write_text(render_script(args, protein, ligand), encoding="utf-8")
        subprocess.run([pymol_exe, "-cq", str(script_path)], check=True)
    overlay_binding_site_labels(Path(args.figure_b), Path(args.metadata))


if __name__ == "__main__":
    main()
