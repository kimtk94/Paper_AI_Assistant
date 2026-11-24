"""Offline regression tests for official supplementary archive retry and identity gates."""
import io
import unittest
import urllib.error
import zipfile
from unittest.mock import patch

from projects.hs_spatial_gse294009 import original_supplement_methods_audit as audit


def official_zip():
    data = io.BytesIO()
    with zipfile.ZipFile(data, "w") as z:
        z.writestr(audit.SHEET, b"%PDF-1.4\nminimal fixture")
    return data.getvalue()


class Response:
    def __init__(self, data):
        self.data = data
        self.status = 200
        self.headers = {"content-type": "application/zip"}

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def read(self, limit):
        return self.data[:limit]


class ArchiveRetryTests(unittest.TestCase):
    @patch.object(audit.time, "sleep")
    @patch.object(audit.urllib.request, "urlopen")
    def test_transient_503_then_exact_source_succeeds(self, open_url, sleep):
        open_url.side_effect = [
            urllib.error.HTTPError(audit.ARCHIVE, 503, "busy", {}, None),
            Response(official_zip()),
        ]
        data, meta, entries = audit.download()
        self.assertTrue(data.startswith(b"%PDF-"))
        self.assertEqual([x["http"] for x in meta["attempts"]], [503, 200])
        self.assertEqual(len(entries), 1)
        sleep.assert_called_once()

    @patch.object(audit.time, "sleep")
    @patch.object(audit.urllib.request, "urlopen")
    def test_permanent_403_fails_immediately(self, open_url, sleep):
        open_url.side_effect = urllib.error.HTTPError(audit.ARCHIVE, 403, "blocked", {}, None)
        with self.assertRaisesRegex(RuntimeError, '403'):
            audit.download()
        self.assertEqual(open_url.call_count, 1)
        sleep.assert_not_called()

    @patch.object(audit.time, "sleep")
    @patch.object(audit.urllib.request, "urlopen")
    def test_exhausted_transient_source_fails(self, open_url, sleep):
        open_url.side_effect = urllib.error.HTTPError(audit.ARCHIVE, 503, "busy", {}, None)
        with self.assertRaisesRegex(RuntimeError, '503'):
            audit.download()
        self.assertEqual(open_url.call_count, audit.MAX_ATTEMPTS)
        self.assertEqual(sleep.call_count, audit.MAX_ATTEMPTS-1)

    @patch.object(audit.urllib.request, "urlopen")
    def test_wrong_archive_member_is_not_accepted(self, open_url):
        out = io.BytesIO()
        with zipfile.ZipFile(out, "w") as z:
            z.writestr("unverified-supplement.pdf", b"%PDF-1.4")
        open_url.return_value = Response(out.getvalue())
        with self.assertRaisesRegex(ValueError, 'expected one official supplement PDF'):
            audit.download()


if __name__ == "__main__":
    unittest.main()
