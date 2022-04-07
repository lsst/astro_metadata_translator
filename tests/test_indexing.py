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
import json
import os
import unittest

from astro_metadata_translator import ObservationGroup, ObservationInfo
from astro_metadata_translator.file_helpers import read_file_info
from astro_metadata_translator.indexing import index_files, process_index_data, process_sidecar_data

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATA = os.path.join(TESTDIR, "data")


class IndexingTestCase(unittest.TestCase):
    def test_indexing(self):
        """Test that we can index two headers"""
        files = ["fitsheader-hsc-HSCA04090107.yaml", "fitsheader-hsc.yaml"]
        files = [os.path.join(TESTDATA, f) for f in files]

        # Index the translated metadata
        index, okay, failed = index_files(files, None, 1, None, "translated")
        self.assertEqual(set(files), set(okay))
        self.assertEqual(failed, [])

        self.assertIn("__COMMON__", index)
        self.assertIn("instrument", index["__COMMON__"])

        # Convert to an ObservationGroup. Filenames are now stored in the
        # corresponding ObservationInfo.
        obs_group = process_index_data(index)
        self.assertIsInstance(obs_group, ObservationGroup)
        self.assertEqual(len(obs_group), 2)
        self.assertEqual(obs_group[0].instrument, "HSC")
        self.assertEqual(set([obs_group[0].filename, obs_group[1].filename]), set(files))

        metadata = process_index_data(index, force_metadata=True)
        self.assertEqual(len(metadata), 2)
        self.assertEqual(metadata[files[0]]["instrument"], "HSC")

        # Index the native FITS headers
        index, okay, failed = index_files(files, None, 1, None, "metadata")
        self.assertEqual(set(files), set(okay))
        self.assertEqual(failed, [])

        # Check that common entries have been factored out
        self.assertIn("__COMMON__", index)
        self.assertIn("TELESCOP", index["__COMMON__"])
        self.assertIn("INSTRUME", index[files[0]])
        self.assertNotIn("INSTRUME", index[files[1]])
        self.assertNotIn("TELESCOP", index[files[0]])

        # Convert back to a dict indexed by filename and check that
        # common has been put back properly.
        metadata = process_index_data(index)
        self.assertEqual(len(metadata), 2)
        self.assertEqual(metadata[files[0]]["INSTRUME"], "Hyper Suprime-Cam")
        self.assertEqual(metadata[files[0]]["TELESCOP"], index["__COMMON__"]["TELESCOP"])
        self.assertEqual(metadata[files[1]]["TELESCOP"], index["__COMMON__"]["TELESCOP"])

    def test_file_reading(self):
        """Test the low-level file reader"""

        # First with a real header (but YAML)
        file = os.path.join(TESTDATA, "fitsheader-hsc-HSCA04090107.yaml")
        info = read_file_info(file, 1, None, "metadata", content_type="simple")
        self.assertEqual(info["PROP-ID"], "o15426")

        info = read_file_info(file, 1, None, "translated", content_type="native")
        self.assertIsInstance(info, ObservationInfo)
        self.assertEqual(info.instrument, "HSC")

        info = read_file_info(file, 1, None, "translated", content_type="simple")
        self.assertIsInstance(info, dict)
        self.assertEqual(info["instrument"], "HSC")

        json_str = read_file_info(file, 1, None, "translated", content_type="json")
        self.assertIsInstance(json_str, str)
        info = json.loads(json_str)
        self.assertEqual(info["instrument"], "HSC")

        processed = process_sidecar_data(info)
        self.assertIsInstance(processed, ObservationInfo)
        self.assertEqual(processed.instrument, "HSC")

        processed = process_sidecar_data(info, force_metadata=True)
        self.assertIsInstance(processed, dict)
        self.assertEqual(processed["instrument"], "HSC")

        json_str = read_file_info(file, 1, None, "metadata", content_type="json")
        self.assertIsInstance(json_str, str)
        info = json.loads(json_str)
        self.assertEqual(info["PROP-ID"], "o15426")

        processed = process_sidecar_data(info)
        self.assertEqual(processed["PROP-ID"], info["PROP-ID"])

        # Read a small fits file
        fits_file = os.path.join(TESTDATA, "small.fits")
        info = read_file_info(fits_file, 0, None, "metadata", content_type="native")
        self.assertEqual(info["FILTER"], "r")

        # The fits file won't translate
        with self.assertRaises(ValueError):
            read_file_info(fits_file, 0, None, "obsInfo")

        with self.assertRaises(ValueError):
            read_file_info(file, 1, None, "unknown")

        with self.assertRaises(FileNotFoundError):
            read_file_info("notthere.not", 1)

        with io.StringIO() as err:
            info = read_file_info("notthere.not", 1, print_trace=False, errstream=err)
            err.seek(0)
            self.assertEqual(err.readlines()[0], "Unable to open file notthere.not\n")

        # Now read a file that can not be translated and should trigger
        # different errors
        bad_file = os.path.join(TESTDATA, "corrections", "SCUBA_test-20000101_00002.yaml")

        with self.assertRaises(ValueError):
            read_file_info(bad_file, 1)

        with io.StringIO() as out:
            with io.StringIO() as err:
                info = read_file_info(bad_file, 1, print_trace=False, errstream=err, outstream=out)
                out.seek(0)
                lines = out.readlines()
                self.assertEqual(len(lines), 1)
                self.assertIn("ValueError", lines[0])

        with io.StringIO() as out:
            with io.StringIO() as err:
                info = read_file_info(bad_file, 1, print_trace=True, errstream=err, outstream=out)
                out.seek(0)
                lines = out.readlines()
                self.assertGreater(len(lines), 4)
                self.assertIn("ValueError", lines[-1])


if __name__ == "__main__":
    unittest.main()
