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

from astro_metadata_translator.bin.translateheader import process_files

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATA = os.path.join(TESTDIR, "data")


class TestTranslateHeader(unittest.TestCase):
    """Test that the translate_header CLI logic works."""

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
                okay, failed = process_files(
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
                okay, failed = process_files(
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
                okay, failed = process_files(
                    [TESTDATA], r"^.*yaml$", 0, False, outstream=out, errstream=err, output_mode="none"
                )

                lines = self._readlines(out)
                self.assertEqual(len(lines), len(failed))
                self.assertTrue(lines[0].startswith("Failure processing"), f"Line: '{lines[0]}'")
                self.assertIn("not a mapping", lines[0], f"Line: '{lines[0]}'")

                lines = self._readlines(err)
                self.assertEqual(len(lines), 13)
                self.assertTrue(lines[0].startswith("Analyzing"), f"Line: '{lines[0]}'")

        self.assertEqual(len(okay), 10)
        self.assertEqual(len(failed), 3)

    def test_translate_header_traceback(self):
        """Translate some header files that fail and trigger traceback"""
        with io.StringIO() as out:
            with io.StringIO() as err:
                okay, failed = process_files(
                    [TESTDATA], r"^.*yaml$", 0, True, outstream=out, errstream=err, output_mode="none"
                )

                lines = self._readlines(out)
                self.assertEqual(len(lines), 22)
                self.assertTrue(lines[0].startswith("Traceback"), f"Line '{lines[0]}'")

                lines = self._readlines(err)
                self.assertEqual(len(lines), 13)
                self.assertTrue(lines[0].startswith("Analyzing"), f"Line: '{lines[0]}'")

        self.assertEqual(len(okay), 10)
        self.assertEqual(len(failed), 3)

    def test_translate_header_dump(self):
        """Check that a header is dumped"""
        with io.StringIO() as out:
            with io.StringIO() as err:
                okay, failed = process_files(
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
        """Check that ObservationInfo content is displayed"""
        with io.StringIO() as out:
            with io.StringIO() as err:
                okay, failed = process_files(
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
