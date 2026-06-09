import tempfile
import unittest
import zipfile
from pathlib import Path

import numpy as np

from docking_pipeline import (
    BoxResult,
    BatchDockingRecord,
    candidate_is_drug_like,
    compute_box_from_coords,
    convert_pdbqt_to_pdb,
    extract_best_pose_pdbqt,
    parse_batch_csv,
    pubchem_3d_sdf_url,
    read_vina_config,
    receptor_cache_key,
    run_directory_name,
    safe_stem,
    append_docking_result_xlsx,
    write_batch_results_csv,
    write_vina_config,
)


class DockingPipelineTests(unittest.TestCase):
    def test_compute_box_from_coords_uses_bounds_and_padding(self):
        coords = np.array([
            [0.0, 1.0, 2.0],
            [10.0, 5.0, 8.0],
        ])

        center, size = compute_box_from_coords(coords, padding=2.0)

        np.testing.assert_allclose(center, [5.0, 3.0, 5.0])
        np.testing.assert_allclose(size, [14.0, 8.0, 10.0])

    def test_vina_config_round_trip_preserves_mode_and_values(self):
        result = BoxResult(
            mode="blind_docking",
            center=np.array([1.1, 2.2, 3.3]),
            size=np.array([20.0, 21.0, 22.0]),
            source="all protein chains",
            warning="large search box",
        )

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "vina_box_config.txt"
            write_vina_config(result, path)
            parsed = read_vina_config(path)

        self.assertEqual(parsed["center_x"], 1.1)
        self.assertEqual(parsed["center_y"], 2.2)
        self.assertEqual(parsed["center_z"], 3.3)
        self.assertEqual(parsed["size_x"], 20.0)
        self.assertEqual(parsed["size_y"], 21.0)
        self.assertEqual(parsed["size_z"], 22.0)
        self.assertEqual(parsed["mode"], "blind_docking")
        self.assertEqual(parsed["warning"], "large search box")

    def test_candidate_filter_rejects_common_crystal_additives(self):
        self.assertFalse(candidate_is_drug_like("ACT", heavy_atoms=4))
        self.assertFalse(candidate_is_drug_like("CAC", heavy_atoms=5))
        self.assertFalse(candidate_is_drug_like("ZN", heavy_atoms=1))
        self.assertTrue(candidate_is_drug_like("LIG", heavy_atoms=18))

    def test_pubchem_url_uses_cid_3d_sdf_endpoint(self):
        self.assertEqual(
            pubchem_3d_sdf_url("72276"),
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/72276/SDF?record_type=3d",
        )

    def test_safe_stem_keeps_name_and_cid_readable(self):
        self.assertEqual(safe_stem("(-)-Epicatechin", "72276"), "epicatechin_72276")

    def test_append_docking_result_xlsx_writes_expected_summary_columns(self):
        row = {
            "Compound": "Epicatechin",
            "Target": "S100A8",
            "PDBID": "5HLO",
            "CID": "72276",
            "Affinity(kcal/mol)": -6.727,
            "Timestamp": "6/3/2026 9:11",
            "Center_X": 9.692,
            "Center_Y": 18.069,
            "Center_Z": -73.708,
            "Size_X": 62.9,
            "Size_Y": 62.2,
            "Size_Z": 100.0,
        }

        with tempfile.TemporaryDirectory() as tmp:
            xlsx_path = Path(tmp) / "docking_results_summary.xlsx"
            append_docking_result_xlsx(xlsx_path, row)

            with zipfile.ZipFile(xlsx_path) as workbook:
                sheet_xml = workbook.read("xl/worksheets/sheet1.xml").decode("utf-8")

        self.assertIn("Compound", sheet_xml)
        self.assertIn("Epicatechin", sheet_xml)
        self.assertIn("S100A8", sheet_xml)
        self.assertIn("5HLO", sheet_xml)
        self.assertIn("72276", sheet_xml)
        self.assertIn("-6.727", sheet_xml)

    def test_extract_best_pose_pdbqt_keeps_first_model_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "poses.pdbqt"
            output = Path(tmp) / "best.pdbqt"
            source.write_text(
                "MODEL 1\n"
                "ATOM      1  C   UNL     1       1.000   2.000   3.000\n"
                "ENDMDL\n"
                "MODEL 2\n"
                "ATOM      2  C   UNL     1       9.000   9.000   9.000\n"
                "ENDMDL\n",
                encoding="utf-8",
            )

            extracted = extract_best_pose_pdbqt(source, output)

            self.assertEqual(extracted, output)
            best_text = output.read_text(encoding="utf-8")
            self.assertIn("MODEL 1", best_text)
            self.assertIn("1.000", best_text)
            self.assertNotIn("MODEL 2", best_text)
            self.assertNotIn("9.000", best_text)

    def test_run_directory_name_uses_target_cid_and_timestamp(self):
        self.assertEqual(
            run_directory_name("5HLO", "72276", "20260608_093000"),
            "5HLO_72276_20260608_093000",
        )

    def test_convert_pdbqt_to_pdb_keeps_atom_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "ligand.pdbqt"
            output = Path(tmp) / "ligand.pdb"
            source.write_text(
                "MODEL 1\n"
                "ATOM      1  C   UNL     1       1.000   2.000   3.000  1.00  0.00     0.000 C\n"
                "ENDMDL\n",
                encoding="utf-8",
            )

            convert_pdbqt_to_pdb(source, output)

            pdb_text = output.read_text(encoding="utf-8")
            self.assertIn("ATOM      1  C   UNL", pdb_text)
            self.assertIn("END", pdb_text)
            self.assertNotIn("MODEL", pdb_text)

    def test_parse_batch_csv_accepts_standard_headers(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ligands.csv"
            path.write_text(
                "compound_name,pubchem_cid,target,pdbid\n"
                "Epicatechin,72276,S100A8,5HLO\n",
                encoding="utf-8",
            )

            records = parse_batch_csv(path)

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].compound_name, "Epicatechin")
        self.assertEqual(records[0].pubchem_cid, "72276")
        self.assertEqual(records[0].target, "S100A8")
        self.assertEqual(records[0].pdbid, "5HLO")

    def test_parse_batch_csv_requires_all_standard_headers(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ligands.csv"
            path.write_text(
                "compound_name,pubchem_cid,target\n"
                "Epicatechin,72276,S100A8\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "pdbid"):
                parse_batch_csv(path)

    def test_receptor_cache_key_deduplicates_same_target_pdbid(self):
        records = [
            BatchDockingRecord(1, "Epicatechin", "72276", "S100A8", "5HLO"),
            BatchDockingRecord(2, "Quercetin", "5280343", "S100A8", "5hlo"),
        ]

        keys = {receptor_cache_key(record) for record in records}

        self.assertEqual(keys, {("S100A8", "5HLO")})

    def test_write_batch_results_csv_records_success_and_failure(self):
        rows = [
            {
                "status": "success",
                "error": "",
                "compound_name": "Epicatechin",
                "pubchem_cid": "72276",
                "target": "S100A8",
                "pdbid": "5HLO",
                "run_dir": "/tmp/run",
                "best_affinity_kcal_per_mol": -6.1,
            },
            {
                "status": "failed",
                "error": "invalid cid",
                "compound_name": "Bad",
                "pubchem_cid": "0",
                "target": "S100A8",
                "pdbid": "5HLO",
                "run_dir": "",
                "best_affinity_kcal_per_mol": "",
            },
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "batch_results_summary.csv"
            write_batch_results_csv(path, rows)
            text = path.read_text(encoding="utf-8")

        self.assertIn("status,error,compound_name,pubchem_cid,target,pdbid,run_dir,best_affinity_kcal_per_mol", text)
        self.assertIn("success", text)
        self.assertIn("failed", text)


if __name__ == "__main__":
    unittest.main()
