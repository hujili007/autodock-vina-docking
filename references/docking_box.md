# Docking Box Selection

AutoDock Vina requires a search box. This skill chooses the box from the original receptor PDB so that co-crystal ligands are still visible during site detection.

## Ligand-Based Box

Use ligand-based docking when the PDB contains a credible non-protein, non-solvent, non-additive ligand. The script filters common crystallographic additives such as water, salts, glycerol, acetate, sulfate, and metal ions.

The ligand-based box is centered on the selected ligand and padded around its coordinate bounds. This is preferred when the co-crystal ligand marks a biologically plausible binding pocket.

## Blind Docking Box

Use blind docking when no credible co-crystal ligand is found. The box is centered on all protein coordinates and padded to cover the receptor.

Blind docking is allowed for learning and exploration, but it is less reliable than docking into a known or experimentally supported pocket. The script prints a warning and writes the warning into `vina_box_config.txt` and downstream summaries.

## Config Format

`vina_box_config.txt` stores both Vina fields and comments:

```text
# mode = blind_docking
# source = all protein atoms
# warning = No credible co-crystal ligand was found; using blind docking.
center_x = 9.692
center_y = 18.069
center_z = -73.708
size_x = 62.900
size_y = 62.200
size_z = 100.000
```

The required numeric fields are `center_x`, `center_y`, `center_z`, `size_x`, `size_y`, and `size_z`.
