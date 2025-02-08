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
import logging
import os
import tempfile
import unittest

from astro_metadata_translator import ObservationGroup, ObservationInfo
from astro_metadata_translator.file_helpers import read_file_info
from astro_metadata_translator.indexing import (
    index_files,
    process_index_data,
    process_sidecar_data,
    read_index,
    read_sidecar,
)

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATA = os.path.join(TESTDIR, "data")


class IndexingTestCase(unittest.TestCase):
    """Test indexing and sidecar functionality."""

    def test_indexing(self):
        """Test that we can index two headers."""
        files = ["fitsheader-hsc-HSCA04090107.yaml", "fitsheader-hsc.yaml"]
        files = [os.path.join(TESTDATA, f) for f in files]

        # Index the translated metadata
        index, okay, failed = index_files(files, None, 1, None, "translated")
        self.assertEqual(set(files), set(okay))
        self.assertEqual(failed, [])

        self.assertIn("__COMMON__", index)
        self.assertIn("instrument", index["__COMMON__"])

        # Write the index and check we can read it back.
        with tempfile.NamedTemporaryFile(suffix=".json", mode="w+") as temp:
            print(json.dumps(index), file=temp)
            temp.flush()
            externally_processed = read_index(temp.name)

        # Convert to an ObservationGroup. Filenames are now stored in the
        # corresponding ObservationInfo.
        obs_group = process_index_data(index)
        self.assertIsInstance(obs_group, ObservationGroup)
        self.assertEqual(len(obs_group), 2)
        self.assertEqual(obs_group[0].instrument, "HSC")
        self.assertEqual({obs_group[0].filename, obs_group[1].filename}, set(files))
        self.assertEqual(externally_processed, obs_group)

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
        """Test the low-level file reader."""
        # First with a real header (but YAML)
        file = os.path.join(TESTDATA, "fitsheader-hsc-HSCA04090107.yaml")
        info = read_file_info(file, 1, None, "metadata", content_type="simple")
        self.assertEqual(info["PROP-ID"], "o15426")

        # With metadata sidecar.
        json_file = os.path.splitext(file)[0] + ".json"
        json_info = read_sidecar(json_file)
        # Need to remove the COMMENT fields to avoid confusion between
        # PropertyList and the fallback code with multiple entries.
        json_info.pop("COMMENT", None)
        dict_info = dict(info)  # it may be a PropertyList
        dict_info.pop("COMMENT", None)
        self.assertEqual(json_info, dict_info)

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

        with self.assertLogs(level=logging.WARNING) as cm:
            info = read_file_info("notthere.not", 1, print_trace=False)

        self.assertIn("Unable to open file notthere.not", "\n".join(cm.output))

        # Now read a file that can not be translated and should trigger
        # different errors
        bad_file = os.path.join(TESTDATA, "corrections", "SCUBA_test-20000101_00002.yaml")

        with self.assertLogs(level="DEBUG") as cm:
            with self.assertRaises(ValueError):
                read_file_info(bad_file, 1)
        self.assertIn("Unable to determine translator class", "\n".join(cm.output))

        with io.StringIO() as out:
            info = read_file_info(bad_file, 1, print_trace=False, outstream=out)
            out.seek(0)
            lines = out.readlines()
            self.assertEqual(len(lines), 1)
            self.assertIn("ValueError", lines[0])

        with io.StringIO() as out:
            info = read_file_info(bad_file, 1, print_trace=True, outstream=out)
            out.seek(0)
            lines = out.readlines()
            self.assertGreater(len(lines), 4)
            self.assertIn("ValueError", lines[-1])

        # A sidecar file that is not a dict.
        not_dict = os.path.join(TESTDATA, "bad-sidecar.json")
        with self.assertRaises(ValueError):
            read_sidecar(not_dict)

        with self.assertRaises(ValueError):
            read_index(not_dict)

        # index file that is not JSON.
        with self.assertRaises(ValueError):
            read_index(bad_file)

    def test_obs_info_sidecar(self):
        """Test reading of older files with missing content."""
        # First with a real header (but YAML)
        file = os.path.join(TESTDATA, "fitsheader-hsc.yaml")
        info = read_file_info(file, 1, None, "translated", content_type="native")
        self.assertIsInstance(info, ObservationInfo)
        self.assertEqual(info.instrument, "HSC")

        # With translated metadata sidecar that lacks the group_counter_start.
        json_file = os.path.splitext(file)[0] + ".json"
        json_info = read_sidecar(json_file)
        self.assertIsInstance(json_info, ObservationInfo)
        self.assertEqual(json_info, info)


if __name__ == "__main__":
    unittest.main()
