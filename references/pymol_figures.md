# PyMOL Figure Recipes

This skill renders two publication-style draft figures automatically. They are intended as editable starting points, not final journal art.

## Figure A: Overall Docking Pose

Output:

```text
figure_a_overall_pose.png
```

Visual style:

- Protein cartoon in light pink.
- Protein molecular surface in semi-transparent gray.
- Best ligand pose in green sticks.
- White background.
- Orthoscopic camera.

This corresponds to an overall receptor view where the ligand position is visible inside or near the predicted pocket.

## Figure B: Binding-Site Close-Up

Output:

```text
figure_b_binding_site.png
```

Visual style:

- Protein cartoon in light gray.
- Ligand best pose in green sticks.
- Residues within 4 Angstrom of ligand in pink sticks.
- Residue labels on nearby alpha carbons.
- Polar contacts shown as yellow dashed distance objects when PyMOL detects them.

If no polar contacts are detected, the figure is still written and `polar_contacts_count = 0` is recorded in the PyMOL render summary.

## Session File

Output:

```text
pymol_docking_session.pse
```

Open this file in the PyMOL GUI to manually adjust labels, camera angle, residue selections, or contact cutoffs.
