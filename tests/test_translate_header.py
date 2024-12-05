# This file is part of astro_metadata_translator.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the LICENSE file at the top-level directory of this distribution
# for details of code ownership.
#
# Use of this source code is governed by a 3-clause BSD-style
# license that can be found in the LICENSE file.

import io
import logging
import os.path
import unittest

from astro_metadata_translator.bin.translate import translate_or_dump_headers

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATA = os.path.join(TESTDIR, "data")


class TestTranslateHeader(unittest.TestCase):
    """Test that the astrometadata translate and dump logic works."""

    def _readlines(self, stream):
        """Return the lines written to the stream.

        Parameters
        ----------
        stream : `io.StringIO`
            The stream to read.

        Returns
        -------
        lines : `list` of `str`
            The lines contained in the stream.
        """
        stream.seek(0)
        return [ll.rstrip() for ll in stream.readlines()]

    def assert_ok_fail(
        self, okay: list[str], failed: list[str], stdout: list[str], expected: tuple[int, int]
    ):
        """Check that we have the expected numbers of successes and
        failures.

        Parameters
        ----------
        okay : `list` [`str`]
            List of successful translations.
        failed : `list` [`str`]
            List of failed translations.
        stdout : `list` [`str`]
            Additional information output.
        expected : `tuple` [`int`, `int`]
            Expected length of ``okay`` and ``failed``.
        """
        # Form message to issue if the test fails.
        newline = "\n"  # f-string can not accept \ in string.
        msg = f"""Converted successfully:
{newline.join(okay)}
Failed conversions:
{newline.join(failed)}
Standard output:
{newline.join(stdout)}
"""
        self.assertEqual((len(okay), len(failed)), expected, msg=msg)

    def test_translate_header(self):
        """Translate some header files."""
        with io.StringIO() as out:
            with self.assertLogs(level=logging.INFO) as cm:
                okay, failed = translate_or_dump_headers(
                    [TESTDATA],
                    r"^fitsheader.*yaml$",
                    0,
                    False,
                    outstream=out,
                    output_mode="none",
                )
            output = self._readlines(out)
            self.assertEqual(output, [])
            lines = [r.getMessage() for r in cm.records]
            self.assertEqual(len(lines), 11)
            self.assertTrue(lines[0].startswith("Analyzing"), f"Line: '{lines[0]}'")

        self.assert_ok_fail(okay, failed, output, (11, 0))

    def test_translate_header_table(self):
        """Translate some header files with table output."""
        with io.StringIO() as out:
            with self.assertLogs(level=logging.WARNING) as cm:
                logging.getLogger().warning("False warning")
                okay, failed = translate_or_dump_headers(
                    [TESTDATA],
                    r"^fitsheader.*yaml$",
                    0,
                    False,
                    outstream=out,
                )
                output = self._readlines(out)
                self.assertIn("ObsId", output[0])
                self.assertTrue(output[2].startswith("-------"))
                self.assertEqual(len(output), 14)  # 3 header lines for QTable.
                self.assertEqual(len(cm.output), 1)  # Should only have the warning this test made.

        self.assert_ok_fail(okay, failed, output, (11, 0))

    def test_translate_bad_header_table(self):
        """Translate a header that has a bad translation in the table."""
        with io.StringIO() as out:
            with self.assertLogs(level=logging.WARNING):
                okay, failed = translate_or_dump_headers(
                    [TESTDATA], "^bad-megaprime.yaml$", 0, False, outstream=out, output_mode="table"
                )
                output = self._readlines(out)
                self.assertIn(" -- ", output[3])  # String masked value.
                self.assertIn(" ———", output[3])  # Quantity masked value.

    def test_translate_header_fails(self):
        """Translate some header files that fail."""
        with io.StringIO() as out:
            with self.assertLogs(level=logging.INFO) as cm:
                okay, failed = translate_or_dump_headers(
                    [TESTDATA], r"^.*yaml$", 0, False, outstream=out, output_mode="none"
                )

            out_lines = self._readlines(out)
            self.assertEqual(len(out_lines), len(failed))
            self.assertTrue(out_lines[0].startswith("Failure processing"), f"Line: '{out_lines[0]}'")
            self.assertIn("not a mapping", out_lines[0], f"Line: '{out_lines[0]}'")

            err_lines = [r.getMessage() for r in cm.records]
            # Filter out warnings.
            analyzed = [e for e in err_lines if e.startswith("Analyzing")]
            self.assertEqual(len(analyzed), 15)  # The number of files analyzed
            self.assertTrue(err_lines[0].startswith("Analyzing"), f"Line: '{err_lines[0]}'")

        self.assert_ok_fail(okay, failed, out_lines, (12, 3))

    def test_translate_header_traceback(self):
        """Translate some header files that fail and trigger traceback."""
        with io.StringIO() as out:
            with self.assertLogs(level=logging.INFO) as cm:
                okay, failed = translate_or_dump_headers(
                    [TESTDATA], r"^.*yaml$", 0, True, outstream=out, output_mode="none"
                )

            out_lines = self._readlines(out)
            self.assertGreaterEqual(len(out_lines), 22, "\n".join(out_lines))
            self.assertTrue(out_lines[0].startswith("Traceback"), f"Line '{out_lines[0]}'")

            lines = [r.getMessage() for r in cm.records]
            self.assertGreaterEqual(len(lines), 13, "\n".join(lines))
            self.assertTrue(lines[0].startswith("Analyzing"), f"Line: '{lines[0]}'")

        self.assert_ok_fail(okay, failed, out_lines, (12, 3))

    def test_translate_header_dump(self):
        """Check that a header is dumped."""
        with io.StringIO() as out:
            with self.assertLogs(level=logging.INFO) as cm:
                okay, failed = translate_or_dump_headers(
                    [os.path.join(TESTDATA, "fitsheader-decam.yaml")],
                    r"^fitsheader.*yaml$",
                    0,
                    False,
                    outstream=out,
                    output_mode="yaml",
                )

            out_lines = self._readlines(out)
            # Look for a DECam header in the output
            header = "\n".join(out_lines)
            self.assertIn("DTINSTRU", header)

            lines = [r.getMessage() for r in cm.records]
            self.assertEqual(len(lines), 1)
            self.assertTrue(lines[0], "Analyzing tests/data/fitsheader-decam.yaml...")

        self.assert_ok_fail(okay, failed, out_lines, (1, 0))

    def test_translate_header_loud(self):
        """Check that ObservationInfo content is displayed."""
        with io.StringIO() as out:
            with self.assertLogs(level=logging.INFO) as cm:
                okay, failed = translate_or_dump_headers(
                    [os.path.join(TESTDATA, "fitsheader-decam.yaml")],
                    r"^fitsheader.*yaml$",
                    0,
                    False,
                    outstream=out,
                    output_mode="verbose",
                )

            out_lines = self._readlines(out)
            # Look for the translated DECam header in the output
            self.assertEqual(out_lines[2], "datetime_begin: 2013-09-01T06:02:55.754")

            lines = [r.getMessage() for r in cm.records]
            self.assertEqual(len(lines), 1)
            self.assertTrue(lines[0], "Analyzing tests/data/fitsheader-decam.yaml...")

        self.assert_ok_fail(okay, failed, out_lines, (1, 0))


if __name__ == "__main__":
    unittest.main()
