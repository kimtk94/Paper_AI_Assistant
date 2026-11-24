import builtins
import copy
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import multiomics_ingest as ingest


class MultiomicsIngestTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.data_dir = self.root / "store"
        self.example = json.loads(
            (ROOT / "examples" / "multiomics_records.example.json").read_text(encoding="utf-8")
        )

    def write_input(self, payload, name="input.json"):
        path = self.root / name
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def test_same_second_backfill_preserves_both_snapshots_and_manifests(self):
        first = copy.deepcopy(self.example)
        second = copy.deepcopy(self.example)
        first["papers"][0]["paper_id"] = "first"
        second["papers"][0]["paper_id"] = "second"
        with patch.object(ingest, "_utc_now", return_value="20261124T120000Z"):
            results = [
                ingest.ingest_bundle(self.write_input(bundle, f"{i}.json"), self.data_dir)
                for i, bundle in enumerate((first, second))
            ]

        self.assertNotEqual(results[0]["run_id"], results[1]["run_id"])
        self.assertEqual(results[0]["ingested_at"], results[1]["ingested_at"])
        snapshots = list((self.data_dir / "snapshots").glob("bundle_*.json"))
        manifests = list((self.data_dir / "manifests").glob("manifest_*.json"))
        self.assertEqual(len(snapshots), 2)
        self.assertEqual(len(manifests), 2)
        self.assertEqual(
            {json.loads(p.read_text())["papers"][0]["paper_id"] for p in snapshots},
            {"first", "second"},
        )
        rows = (self.data_dir / "paper_records.jsonl").read_text().splitlines()
        self.assertEqual(len(rows), len(first["papers"]) + len(second["papers"]))
        for result in results:
            manifest = self.data_dir / "manifests" / f"manifest_{result['run_id']}.json"
            self.assertEqual(json.loads(manifest.read_text()), result)
            self.assertTrue(Path(result["snapshot"]).exists())

    def test_duplicate_run_id_cannot_overwrite_or_append(self):
        path = self.write_input(self.example)
        with (
            patch.object(ingest, "_utc_now", return_value="20261124T120000Z"),
            patch.object(ingest, "uuid4", return_value=SimpleNamespace(hex="collision")),
        ):
            ingest.ingest_bundle(path, self.data_dir)
            before = {str(p): p.read_bytes() for p in self.data_dir.rglob("*") if p.is_file()}
            with self.assertRaises(FileExistsError):
                ingest.ingest_bundle(path, self.data_dir)
        after = {str(p): p.read_bytes() for p in self.data_dir.rglob("*") if p.is_file()}
        self.assertEqual(before, after)

    def test_invalid_envelopes_fail_before_writing(self):
        for payload in (
            [], None, {}, {"count": 1, "results": [{"paper_id": "lost"}]},
            {"papers": {}}, {"ideas": ["not an object"]},
            {"papers": [], "unexpected": []},
        ):
            with self.subTest(payload=payload):
                with self.assertRaises(ValueError):
                    ingest.ingest_bundle(self.write_input(payload), self.data_dir)
                self.assertFalse(self.data_dir.exists())

    @unittest.skipUnless(importlib.util.find_spec("pydantic"), "requires optional Pydantic")
    def test_record_validation_errors_are_not_silently_accepted(self):
        invalid = copy.deepcopy(self.example)
        invalid["papers"][0]["year"] = 2500
        for strict in (False, True):
            with self.subTest(strict=strict):
                with self.assertRaises(ValueError) as error:
                    ingest.ingest_bundle(
                        self.write_input(invalid), self.data_dir, strict_validation=strict
                    )
                self.assertIn("year", str(error.exception))
                self.assertFalse(self.data_dir.exists())

    def missing_dependency(self, name, globals=None, locals=None, fromlist=(), level=0):
        if name == "multiomics_models":
            raise ModuleNotFoundError("No module named 'pydantic'", name="pydantic")
        return self.original_import(name, globals, locals, fromlist, level)

    def test_missing_optional_dependency_warns_in_default_mode(self):
        self.original_import = builtins.__import__
        with patch("builtins.__import__", side_effect=self.missing_dependency):
            with self.assertWarnsRegex(RuntimeWarning, "only the bundle envelope"):
                result = ingest.ingest_bundle(self.write_input(self.example), self.data_dir)
        self.assertEqual(result["papers_appended"], len(self.example["papers"]))

    def test_strict_mode_requires_dependency_before_writing(self):
        self.original_import = builtins.__import__
        with patch("builtins.__import__", side_effect=self.missing_dependency):
            with self.assertRaisesRegex(RuntimeError, "requires Pydantic"):
                ingest.ingest_bundle(
                    self.write_input(self.example), self.data_dir, strict_validation=True
                )
        self.assertFalse(self.data_dir.exists())

    def test_invalid_envelope_is_rejected_even_without_pydantic(self):
        self.original_import = builtins.__import__
        with patch("builtins.__import__", side_effect=self.missing_dependency):
            with self.assertRaises(ValueError):
                ingest.ingest_bundle(self.write_input({"results": []}), self.data_dir)
        self.assertFalse(self.data_dir.exists())

    def test_unrelated_missing_module_does_not_enable_fallback(self):
        original_import = builtins.__import__

        def broken_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "multiomics_models":
                raise ModuleNotFoundError("Broken dependency", name="pydantic_core")
            return original_import(name, globals, locals, fromlist, level)

        with patch("builtins.__import__", side_effect=broken_import):
            with self.assertRaises(ModuleNotFoundError):
                ingest.ingest_bundle(self.write_input(self.example), self.data_dir)
        self.assertFalse(self.data_dir.exists())


if __name__ == "__main__":
    unittest.main()
