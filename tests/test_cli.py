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

import logging
import os.path
import tempfile
import traceback
import unittest

import click.testing
import yaml

from astro_metadata_translator import ObservationGroup
from astro_metadata_translator.cli.astrometadata import main as astrometadata
from astro_metadata_translator.indexing import read_index, read_sidecar

TESTDIR = os.path.abspath(os.path.dirname(__file__))


def click_result_msg(result: click.testing.Result) -> str:
    """Get a standard assert message from a click result.

    Parameters
    ----------
    result : click.testing.Result
        The result object returned from `click.testing.CliRunner.invoke`.

    Returns
    -------
    msg : `str`
        The message string.
    """
    msg = f"""\noutput: {result.output}\nexception: {result.exception}"""
    if result.exception:
        msg += f"""\ntraceback: {"".join(traceback.format_tb(result.exception.__traceback__))}"""
    return msg


class TestCLI(unittest.TestCase):
    """Test astrometadata command line."""

    def setUp(self) -> None:
        self.runner = click.testing.CliRunner()

    def test_package_load(self) -> None:
        """Test that -p option fails properly."""
        # For unknown reasons the click runner does not capture the log
        # message from this so use the unittest capture instead.
        with self.assertLogs(level=logging.WARNING) as cm:
            result = self.runner.invoke(
                astrometadata,
                [
                    "--packages",
                    "lsst.not.found",
                    "dump",
                    os.path.join(TESTDIR, "data", "fitsheader-hsc-HSCA04090107.json"),
                ],
            )
        self.assertEqual(result.exit_code, 0, click_result_msg(result))
        self.assertIn("Failed to import translator module", cm.output[0], click_result_msg(result))

    def test_translate(self) -> None:
        """Test astrometadata translate."""
        result = self.runner.invoke(
            astrometadata, ["translate", os.path.join(TESTDIR, "data", "fitsheader-decam.yaml")]
        )
        self.assertEqual(result.exit_code, 0, click_result_msg(result))
        self.assertIn("instrument: DECam", result.output, click_result_msg(result))

        result = self.runner.invoke(
            astrometadata,
            ["--traceback", "translate", os.path.join(TESTDIR, "data", "bad-sidecar.json")],
        )
        self.assertEqual(result.exit_code, 1, click_result_msg(result))
        self.assertIn("ValueError", result.output, click_result_msg(result))

    def test_dump(self) -> None:
        """Test that YAML is dumped."""
        # Use warning log level to hide the INFO message.
        result = self.runner.invoke(
            astrometadata,
            [
                "--log-level",
                "WARNING",
                "dump",
                os.path.join(TESTDIR, "data", "fitsheader-hsc-HSCA04090107.json"),
            ],
        )
        self.assertEqual(result.exit_code, 0, click_result_msg(result))
        parsed = yaml.safe_load(result.output)
        self.assertEqual(parsed["T_EFMN32"], 50)

    def test_write_sidecar(self):
        """Write a simple sidecar file."""
        # Sidecar files are written next to file.
        with self.runner.isolated_filesystem():
            # Test that we set exit status to 1 if we find no files.
            result = self.runner.invoke(astrometadata, ["write-sidecar", "-c", "metadata", os.path.curdir])
            self.assertEqual(result.exit_code, 1, click_result_msg(result))
            self.assertIn("Found no files", result.output, click_result_msg(result))

            # Test that we can write a sidecar file.
            infile = "fitsheader-decam.yaml"
            os.symlink(os.path.join(TESTDIR, "data", infile), infile)
            result = self.runner.invoke(
                astrometadata, ["write-sidecar", "-c", "metadata", "--regex", r".*\.yaml", os.path.curdir]
            )
            self.assertEqual(result.exit_code, 0, click_result_msg(result))

            jfile = os.path.splitext(infile)[0] + ".json"
            content = read_sidecar(jfile)
            self.assertEqual(content["PROGRAM"], "supernova")

            # Test that we set bad status if a file can not be translated.
            os.symlink(os.path.join(TESTDIR, "data", "small.fits"), "bad.fits")
            result = self.runner.invoke(astrometadata, ["write-sidecar", os.path.curdir])
            self.assertEqual(result.exit_code, 1, click_result_msg(result))
            self.assertIn("No files processed", result.output, click_result_msg(result))

    def test_write_index(self):
        """Write a simple index file."""
        with tempfile.NamedTemporaryFile(suffix=".json") as temp:
            temp.close()
            result = self.runner.invoke(
                astrometadata,
                [
                    "write-index",
                    "-o",
                    temp.name,
                    "--regex",
                    r"fitsheader.*\.yaml",
                    os.path.join(TESTDIR, "data"),
                ],
            )
            self.assertEqual(result.exit_code, 0, click_result_msg(result))
            self.assertNotIn("failed header extraction", result.output, msg="\n\n" + result.output)

            ind = read_index(temp.name)
            self.assertIsInstance(ind, ObservationGroup)
            self.assertGreaterEqual(len(ind), 10)


if __name__ == "__main__":
    unittest.main()
