import importlib.util
import json
import unittest
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]
SPEC = importlib.util.spec_from_file_location("packet_valid", BASE / "validate.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


class PacketTests(unittest.TestCase):
    def test_valid_existing_packet(self):
        self.assertEqual(M.validate(ROOT), [])

    def test_missing_section_fails(self):
        problems = M.inspect_text("# Example\n" + "content " * 150, "proposal")
        self.assertTrue(any("missing section" in issue for issue in problems))

    def test_email_address_disallowed(self):
        problems = M.inspect_text("# Mail\n" + "text " * 80 + "person@example.com", "email")
        self.assertTrue(any("private" in issue for issue in problems))

    def test_crm_term_disallowed(self):
        problems = M.inspect_text("# Report\n" + "text " * 150 + " SFDC ", "proposal")
        self.assertTrue(any("private" in issue for issue in problems))

    def test_unique_six_entries(self):
        entries = json.loads((BASE / "manifest.json").read_text(encoding="utf-8"))["proposals"]
        self.assertEqual(len(entries), 6)
        self.assertEqual(len({e["proposal"] for e in entries}), 6)
        self.assertEqual(len({e["email"] for e in entries}), 6)
