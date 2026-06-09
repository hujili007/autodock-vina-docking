from __future__ import annotations

import json
import re
import subprocess
import sys
import urllib.error
import urllib.request
import zipfile
import csv
from dataclasses import dataclass
from datetime import datetime
from html import escape
from pathlib import Path
from typing import Iterable
from xml.etree import ElementTree

import numpy as np


COMMON_NON_DRUG_RESNAMES = {
    "ACE",
    "ACT",
    "ACN",
    "ARS",
    "BME",
    "CAC",
    "CA",
    "CD",
    "CL",
    "CO",
    "DMS",
    "DMSO",
    "EDO",
    "GOL",
    "HOH",
    "IOD",
    "K",
    "MG",
    "MN",
    "NA",
    "NH4",
    "NI",
    "NO3",
    "PEG",
    "PO4",
    "SO4",
    "ZN",
}

DOCKING_RESULTS_HEADERS = [
    "Compound",
    "Target",
    "PDBID",
    "CID",
    "Affinity(kcal/mol)",
    "Timestamp",
    "Center_X",
    "Center_Y",
    "Center_Z",
    "Size_X",
    "Size_Y",
    "Size_Z",
]

BATCH_CSV_REQUIRED_HEADERS = ["compound_name", "pubchem_cid", "target", "pdbid"]
BATCH_RESULTS_HEADERS = [
    "status",
    "error",
    "compound_name",
    "pubchem_cid",
    "target",
    "pdbid",
    "run_dir",
    "best_affinity_kcal_per_mol",
]


@dataclass
class BoxResult:
    mode: str
    center: np.ndarray
    size: np.ndarray
    source: str
    warning: str | None = None


@dataclass
class LigandCandidate:
    resname: str
    chain: str
    resnum: int
    coords: np.ndarray
    heavy_atoms: int

    @property
    def label(self) -> str:
        chain = self.chain or "_"
        return f"{self.resname}:{chain}:{self.resnum}"


@dataclass(frozen=True)
class BatchDockingRecord:
    row_index: int
    compound_name: str
    pubchem_cid: str
    target: str
    pdbid: str


def log_step(message: str) -> None:
    print(f"[Docking step] {message}")


def safe_stem(name: str, cid: str | int | None = None) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", name.strip().lower()).strip("_")
    cleaned = cleaned or "compound"
    if cid is None:
        return cleaned
    return f"{cleaned}_{cid}"


def parse_batch_csv(path: Path) -> list[BatchDockingRecord]:
    path = Path(path)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        missing = [header for header in BATCH_CSV_REQUIRED_HEADERS if header not in fieldnames]
        if missing:
            raise ValueError(
                "Batch CSV is missing required header(s): "
                + ", ".join(missing)
                + ". Required headers: "
                + ", ".join(BATCH_CSV_REQUIRED_HEADERS)
            )

        records: list[BatchDockingRecord] = []
        for row_number, row in enumerate(reader, start=2):
            values = {key: (row.get(key) or "").strip() for key in BATCH_CSV_REQUIRED_HEADERS}
            if not any(values.values()):
                continue
            missing_values = [key for key, value in values.items() if not value]
            if missing_values:
                raise ValueError(
                    f"Batch CSV row {row_number} is missing value(s): {', '.join(missing_values)}"
                )
            records.append(
                BatchDockingRecord(
                    row_index=len(records) + 1,
                    compound_name=values["compound_name"],
                    pubchem_cid=values["pubchem_cid"],
                    target=values["target"],
                    pdbid=values["pdbid"].upper(),
                )
            )
    if not records:
        raise ValueError(f"Batch CSV contains no docking rows: {path}")
    return records


def receptor_cache_key(record: BatchDockingRecord) -> tuple[str, str]:
    return (record.target, record.pdbid.upper())


def write_batch_results_csv(path: Path, rows: list[dict]) -> None:
    path = Path(path)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=BATCH_RESULTS_HEADERS)
        writer.writeheader()
        for row in rows:
            writer.writerow({header: row.get(header, "") for header in BATCH_RESULTS_HEADERS})


def compute_box_from_coords(coords: np.ndarray, padding: float) -> tuple[np.ndarray, np.ndarray]:
    if coords.size == 0:
        raise ValueError("No coordinates were provided for box calculation.")
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    center = (mins + maxs) / 2.0
    size = (maxs - mins) + (2.0 * padding)
    return center, size


def candidate_is_drug_like(resname: str, heavy_atoms: int, min_heavy_atoms: int = 6) -> bool:
    return resname.upper() not in COMMON_NON_DRUG_RESNAMES and heavy_atoms >= min_heavy_atoms


def ligand_candidates_from_selection(selection, min_heavy_atoms: int = 6) -> list[LigandCandidate]:
    if selection is None:
        return []

    grouped: dict[tuple[str, str, int], list[int]] = {}
    resnames = selection.getResnames()
    chains = selection.getChids()
    resnums = selection.getResnums()
    coords = selection.getCoords()
    elements = selection.getElements()

    for idx, (resname, chain, resnum, element) in enumerate(zip(resnames, chains, resnums, elements)):
        if str(element).upper() == "H":
            continue
        key = (str(resname).upper(), str(chain), int(resnum))
        grouped.setdefault(key, []).append(idx)

    candidates: list[LigandCandidate] = []
    for (resname, chain, resnum), indices in grouped.items():
        if candidate_is_drug_like(resname, len(indices), min_heavy_atoms=min_heavy_atoms):
            candidates.append(
                LigandCandidate(
                    resname=resname,
                    chain=chain,
                    resnum=resnum,
                    coords=coords[indices],
                    heavy_atoms=len(indices),
                )
            )
    return sorted(candidates, key=lambda item: item.heavy_atoms, reverse=True)


def write_vina_config(result: BoxResult, output_path: Path) -> None:
    output_path = Path(output_path)
    lines = [
        f"# mode = {result.mode}",
        f"# source = {result.source}",
    ]
    if result.warning:
        lines.append(f"# warning = {result.warning}")
    lines.extend(
        [
            f"center_x = {result.center[0]:.3f}",
            f"center_y = {result.center[1]:.3f}",
            f"center_z = {result.center[2]:.3f}",
            f"size_x = {result.size[0]:.3f}",
            f"size_y = {result.size[1]:.3f}",
            f"size_z = {result.size[2]:.3f}",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_vina_config(config_path: Path) -> dict[str, float | str]:
    values: dict[str, float | str] = {}
    for line in Path(config_path).read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            comment = stripped[1:].strip()
            if "=" in comment:
                key, value = comment.split("=", 1)
                values[key.strip()] = value.strip()
            continue
        if "=" not in stripped:
            continue
        key, value = stripped.split("=", 1)
        values[key.strip()] = float(value.strip())

    required = {"center_x", "center_y", "center_z", "size_x", "size_y", "size_z"}
    missing = sorted(required.difference(values))
    if missing:
        raise ValueError(f"Vina config is missing required fields: {', '.join(missing)}")
    return values


def pubchem_3d_sdf_url(cid: str | int) -> str:
    return f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"


def pubchem_2d_sdf_url(cid: str | int) -> str:
    return f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=2d"


def download_file(url: str, output_path: Path, timeout: int = 60) -> None:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            data = response.read()
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"Download failed with HTTP {exc.code}: {url}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Download failed: {url} ({exc.reason})") from exc
    Path(output_path).write_bytes(data)


def generate_3d_sdf_from_2d(input_sdf: Path, output_sdf: Path) -> None:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    supplier = Chem.SDMolSupplier(str(input_sdf), removeHs=False)
    mol = next((item for item in supplier if item is not None), None)
    if mol is None:
        raise RuntimeError(f"RDKit could not read a molecule from {input_sdf}")

    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result != 0:
        raise RuntimeError("RDKit failed to generate a 3D conformer from the 2D structure.")

    try:
        status = AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        status = AllChem.UFFOptimizeMolecule(mol)
    if status not in (0, 1):
        print("[Warning] 3D geometry optimization did not fully converge.")

    writer = Chem.SDWriter(str(output_sdf))
    writer.write(mol)
    writer.close()


def run_command(command: Iterable[str], description: str) -> None:
    command = list(command)
    log_step(description)
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(f"Command not found: {command[0]}") from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"Command failed with exit code {exc.returncode}: {' '.join(command)}") from exc


def env_script(name: str) -> Path:
    return Path(sys.executable).with_name(name)


def write_summary_json(path: Path, payload: dict) -> None:
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def formatted_timestamp(now: datetime | None = None) -> str:
    now = now or datetime.now()
    return f"{now.month}/{now.day}/{now.year} {now.hour}:{now.minute:02d}"


def compact_timestamp(now: datetime | None = None) -> str:
    now = now or datetime.now()
    return now.strftime("%Y%m%d_%H%M%S")


def run_directory_name(pdb_id: str, pubchem_cid: str | int, timestamp: str | None = None) -> str:
    stamp = timestamp or compact_timestamp()
    pdb_part = re.sub(r"[^A-Za-z0-9]+", "_", str(pdb_id).upper()).strip("_") or "PDB"
    cid_part = re.sub(r"[^A-Za-z0-9]+", "_", str(pubchem_cid)).strip("_") or "CID"
    return f"{pdb_part}_{cid_part}_{stamp}"


def extract_best_pose_pdbqt(source_path: Path, output_path: Path) -> Path:
    source_path = Path(source_path)
    output_path = Path(output_path)
    lines = source_path.read_text(encoding="utf-8").splitlines()
    keep: list[str] = []
    in_first_model = False
    saw_model = False

    for line in lines:
        if line.startswith("MODEL"):
            saw_model = True
            in_first_model = line.split()[-1] == "1"
            if in_first_model:
                keep.append(line)
            continue
        if line.startswith("ENDMDL"):
            if in_first_model:
                keep.append(line)
                break
            in_first_model = False
            continue
        if in_first_model or not saw_model:
            keep.append(line)

    atom_lines = [line for line in keep if line.startswith(("ATOM", "HETATM"))]
    if not atom_lines:
        raise RuntimeError(f"No atoms were found in the best pose extracted from {source_path}")

    output_path.write_text("\n".join(keep) + "\n", encoding="utf-8")
    return output_path


def convert_pdbqt_to_pdb(source_path: Path, output_path: Path) -> Path:
    source_path = Path(source_path)
    output_path = Path(output_path)
    pdb_lines: list[str] = []

    for line in source_path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            pdb_lines.append(line[:66].rstrip())
        elif line.startswith("TER"):
            pdb_lines.append("TER")

    if not pdb_lines:
        raise RuntimeError(f"No atom records were found while converting {source_path} to PDB.")

    pdb_lines.append("END")
    output_path.write_text("\n".join(pdb_lines) + "\n", encoding="utf-8")
    return output_path


def _xlsx_col_name(index: int) -> str:
    name = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        name = chr(65 + remainder) + name
    return name


def _is_number(value) -> bool:
    return isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, bool)


def _cell_xml(row_idx: int, col_idx: int, value, style: int | None = None) -> str:
    cell_ref = f"{_xlsx_col_name(col_idx)}{row_idx}"
    style_attr = f' s="{style}"' if style is not None else ""
    if value is None:
        return f'<c r="{cell_ref}"{style_attr}/>'
    if _is_number(value):
        return f'<c r="{cell_ref}"{style_attr}><v>{float(value):.6g}</v></c>'
    text = escape(str(value))
    return f'<c r="{cell_ref}" t="inlineStr"{style_attr}><is><t>{text}</t></is></c>'


def _sheet_xml(rows: list[list]) -> str:
    max_row = len(rows)
    max_col = len(DOCKING_RESULTS_HEADERS)
    dimension = f"A1:{_xlsx_col_name(max_col)}{max(max_row, 1)}"
    column_widths = [24, 14, 12, 12, 18, 20, 12, 12, 12, 12, 12, 12]
    cols_xml = "".join(
        f'<col min="{idx}" max="{idx}" width="{width}" customWidth="1"/>'
        for idx, width in enumerate(column_widths, start=1)
    )
    row_xml = []
    for row_idx, row in enumerate(rows, start=1):
        style = 1 if row_idx == 1 else 2
        cells = "".join(
            _cell_xml(row_idx, col_idx, row[col_idx - 1] if col_idx <= len(row) else "", style)
            for col_idx in range(1, max_col + 1)
        )
        row_xml.append(f'<row r="{row_idx}">{cells}</row>')
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f'<dimension ref="{dimension}"/>'
        '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" '
        'activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>'
        '<sheetFormatPr defaultRowHeight="15"/>'
        f'<cols>{cols_xml}</cols>'
        f'<sheetData>{"".join(row_xml)}</sheetData>'
        '<autoFilter ref="A1:L1"/>'
        '</worksheet>'
    )


def _styles_xml() -> str:
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="2">'
        '<font><sz val="11"/><name val="Calibri"/></font>'
        '<font><b/><sz val="11"/><name val="Calibri"/></font>'
        '</fonts>'
        '<fills count="3">'
        '<fill><patternFill patternType="none"/></fill>'
        '<fill><patternFill patternType="gray125"/></fill>'
        '<fill><patternFill patternType="solid"><fgColor rgb="FFFFFF00"/><bgColor indexed="64"/></patternFill></fill>'
        '</fills>'
        '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
        '<cellXfs count="3">'
        '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
        '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
        '<xf numFmtId="0" fontId="0" fillId="2" borderId="0" xfId="0" applyFill="1"/>'
        '</cellXfs>'
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        '</styleSheet>'
    )


def _read_existing_xlsx_rows(path: Path) -> list[list]:
    if not Path(path).exists():
        return [DOCKING_RESULTS_HEADERS]

    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(path) as workbook:
        root = ElementTree.fromstring(workbook.read("xl/worksheets/sheet1.xml"))

    rows: list[list] = []
    for row_node in root.findall(".//m:sheetData/m:row", ns):
        values = [""] * len(DOCKING_RESULTS_HEADERS)
        for cell in row_node.findall("m:c", ns):
            ref = cell.attrib.get("r", "")
            col_letters = re.sub(r"\d+", "", ref)
            col_idx = 0
            for char in col_letters:
                col_idx = col_idx * 26 + (ord(char) - 64)
            if not 1 <= col_idx <= len(values):
                continue
            if cell.attrib.get("t") == "inlineStr":
                text_node = cell.find("m:is/m:t", ns)
                value = text_node.text if text_node is not None else ""
            else:
                value_node = cell.find("m:v", ns)
                value = value_node.text if value_node is not None else ""
            values[col_idx - 1] = value
        rows.append(values)

    if not rows:
        return [DOCKING_RESULTS_HEADERS]
    return rows


def _write_simple_xlsx(path: Path, rows: list[list]) -> None:
    path = Path(path)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as workbook:
        workbook.writestr(
            "[Content_Types].xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
            '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>'
            '<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>'
            '</Types>',
        )
        workbook.writestr(
            "_rels/.rels",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
            '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>'
            '<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>'
            '</Relationships>',
        )
        workbook.writestr(
            "xl/workbook.xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
            '<sheets><sheet name="Docking Results" sheetId="1" r:id="rId1"/></sheets>'
            '</workbook>',
        )
        workbook.writestr(
            "xl/_rels/workbook.xml.rels",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>'
            '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
            '</Relationships>',
        )
        workbook.writestr("xl/worksheets/sheet1.xml", _sheet_xml(rows))
        workbook.writestr("xl/styles.xml", _styles_xml())
        workbook.writestr(
            "docProps/core.xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
            'xmlns:dc="http://purl.org/dc/elements/1.1/">'
            '<dc:title>Docking Results Summary</dc:title>'
            '</cp:coreProperties>',
        )
        workbook.writestr(
            "docProps/app.xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" '
            'xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">'
            '<Application>Docking Pipeline</Application>'
            '</Properties>',
        )


def append_docking_result_xlsx(path: Path, row: dict) -> None:
    rows = _read_existing_xlsx_rows(path)
    normalized_row = [row.get(header, "") for header in DOCKING_RESULTS_HEADERS]
    rows.append(normalized_row)
    _write_simple_xlsx(Path(path), rows)
