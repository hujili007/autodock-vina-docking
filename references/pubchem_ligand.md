# PubChem Ligand Preparation

The primary ligand input is a PubChem CID. The compound name is used for readable logs, filenames, and the Excel summary, but the CID is treated as the structure identifier.

## Preferred Path

1. Download PubChem 3D SDF:

```text
https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/<CID>/SDF?record_type=3d
```

2. Read the SDF with RDKit.
3. Add hydrogens when needed.
4. Convert the ligand to PDBQT with Meeko.

## Fallback Path

If PubChem 3D SDF is unavailable:

1. Download PubChem 2D SDF:

```text
https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/<CID>/SDF?record_type=2d
```

2. Use RDKit ETKDG to generate a 3D conformer.
3. Optimize geometry with MMFF when possible, otherwise UFF.
4. Convert to PDBQT.

## Why Not SMILES First

SMILES is useful, but it does not contain 3D coordinates. For this learning workflow, PubChem CID plus PubChem 3D SDF is more direct and more reproducible. SMILES can be added later as an optional user input or as a last-resort fallback.

## Output Names

Ligand files use a safe stem built from compound name and CID, for example:

```text
epicatechin_72276.sdf
epicatechin_72276.pdbqt
```
