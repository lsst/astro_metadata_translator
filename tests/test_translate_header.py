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

    def test_translate_header(self):
        """Translate some header files."""
        with io.StringIO() as out:
            with io.StringIO() as err:
                okay, failed = translate_or_dump_headers(
                    [TESTDATA],
                    r"^fitsheader.*yaml$",
                    0,
                    False,
                    outstream=out,
                    errstream=err,
                    output_mode="none",
                )
                self.assertEqual(self._readlines(out), [])
                lines = self._readlines(err)
                self.assertEqual(len(lines), 10)
                self.assertTrue(lines[0].startswith("Analyzing"), f"Line: '{lines[0]}'")

        self.assertEqual(len(okay), 10)
        self.assertEqual(len(failed), 0)

    def test_translate_header_table(self):
        """Translate some header files with table output."""
        with io.StringIO() as out:
            with io.StringIO() as err:
                okay, failed = translate_or_dump_headers(
                    [TESTDATA], r"^fitsheader.*yaml$", 0, False, outstream=out, errstream=err
                )
                output = self._readlines(out)
                self.assertTrue(output[0].startswith("ObsId"))
                self.assertTrue(output[1].startswith("-------"))
                self.assertEqual(len(output), 12)
                errlines = self._readlines(err)
                self.assertEqual(len(errlines), 0)

        self.assertEqual(len(okay), 10)
        self.assertEqual(len(failed), 0)

    def test_translate_header_fails(self):
        """Translate some header files that fail."""
        with io.StringIO() as out:
            with io.StringIO() as err:
                okay, failed = translate_or_dump_headers(
                    [TESTDATA], r"^.*yaml$", 0, False, outstream=out, errstream=err, output_mode="none"
                )

                out_lines = self._readlines(out)
                self.assertEqual(len(out_lines), len(failed))
                self.assertTrue(out_lines[0].startswith("Failure processing"), f"Line: '{out_lines[0]}'")
                self.assertIn("not a mapping", out_lines[0], f"Line: '{out_lines[0]}'")

                err_lines = self._readlines(err)
                self.assertEqual(len(err_lines), 13)  # The number of files analyzed
                self.assertTrue(err_lines[0].startswith("Analyzing"), f"Line: '{err_lines[0]}'")

        # Form message to issue if the test fails.
        newline = "\n"  # f-string can not accept \ in string.
        msg = f"""Converted successfully:
{newline.join(okay)}
Failed conversions:
{newline.join(failed)}
Standard output:
{newline.join(out_lines)}
"""
        self.assertEqual((len(okay), len(failed)), (10, 3), msg=msg)

    def test_translate_header_traceback(self):
        """Translate some header files that fail and trigger traceback."""
        with io.StringIO() as out:
            with io.StringIO() as err:
                okay, failed = translate_or_dump_headers(
                    [TESTDATA], r"^.*yaml$", 0, True, outstream=out, errstream=err, output_mode="none"
                )

                lines = self._readlines(out)
                self.assertGreaterEqual(len(lines), 22, "\n".join(lines))
                self.assertTrue(lines[0].startswith("Traceback"), f"Line '{lines[0]}'")

                lines = self._readlines(err)
                self.assertGreaterEqual(len(lines), 13, "\n".join(lines))
                self.assertTrue(lines[0].startswith("Analyzing"), f"Line: '{lines[0]}'")

        self.assertEqual(len(okay), 10)
        self.assertEqual(len(failed), 3)

    def test_translate_header_dump(self):
        """Check that a header is dumped."""
        with io.StringIO() as out:
            with io.StringIO() as err:
                okay, failed = translate_or_dump_headers(
                    [os.path.join(TESTDATA, "fitsheader-decam.yaml")],
                    r"^fitsheader.*yaml$",
                    0,
                    False,
                    outstream=out,
                    errstream=err,
                    output_mode="yaml",
                )

                lines = self._readlines(out)
                # Look for a DECam header in the output
                header = "\n".join(lines)
                self.assertIn("DTINSTRU", header)

                lines = self._readlines(err)
                self.assertEqual(len(lines), 1)
                self.assertTrue(lines[0], "Analyzing tests/data/fitsheader-decam.yaml...")

        self.assertEqual(len(okay), 1)
        self.assertEqual(len(failed), 0)

    def test_translate_header_loud(self):
        """Check that ObservationInfo content is displayed."""
        with io.StringIO() as out:
            with io.StringIO() as err:
                okay, failed = translate_or_dump_headers(
                    [os.path.join(TESTDATA, "fitsheader-decam.yaml")],
                    r"^fitsheader.*yaml$",
                    0,
                    False,
                    outstream=out,
                    errstream=err,
                    output_mode="verbose",
                )

                lines = self._readlines(out)
                # Look for the translated DECam header in the output
                self.assertEqual(lines[2], "datetime_begin: 2013-09-01T06:02:55.754")

                lines = self._readlines(err)
                self.assertEqual(len(lines), 1)
                self.assertTrue(lines[0], "Analyzing tests/data/fitsheader-decam.yaml...")

        self.assertEqual(len(okay), 1)
        self.assertEqual(len(failed), 0)


if __name__ == "__main__":
    unittest.main()
