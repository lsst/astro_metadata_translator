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

import os.path
import unittest

from astro_metadata_translator import ObservationGroup
from astro_metadata_translator.serialize import group_to_fits, info_to_fits
from astro_metadata_translator.tests import read_test_file

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class ObservationGroupTestCase(unittest.TestCase):
    """Test grouping of observations."""

    datadir = os.path.join(TESTDIR, "data")

    def setUp(self):
        self.decam_files = (
            "fitsheader-decam.yaml",
            "fitsheader-decam-0160496.yaml",
            "fitsheader-decam-calexp-0412037_10.yaml",
        )
        self.hsc_files = ("fitsheader-hsc-HSCA04090107.yaml", "fitsheader-hsc.yaml")

    def _files_to_headers(self, files):
        return [read_test_file(os.path.join(self.datadir, f)) for f in files]

    def test_groups(self):
        headers = self._files_to_headers(self.decam_files)

        obs_group = ObservationGroup(headers)
        self.assertEqual(len(obs_group), 3)
        self.assertEqual(
            str(obs_group),
            "[(DECam, 2013-09-01T06:02:55.754), (DECam, 2012-12-11T22:06:32.859),"
            " (DECam, 2015-02-20T00:47:21.127)]",
        )

        sorted_group = ObservationGroup(sorted(obs_group))
        self.assertIsInstance(sorted_group, ObservationGroup)
        self.assertEqual(len(sorted_group), 3)
        self.assertEqual(sorted_group[0], obs_group[1])

        self.assertNotEqual(obs_group, sorted_group)
        obs_group.sort()
        self.assertEqual(obs_group, sorted_group)
        obs_group.reverse()
        self.assertEqual(obs_group[0], sorted_group[-1])

        newest = obs_group.newest()
        oldest = obs_group.oldest()
        self.assertEqual(newest, sorted_group[-1])
        self.assertEqual(oldest, sorted_group[0])

        self.assertLess(oldest, newest)
        self.assertGreater(newest, oldest)

        self.assertNotEqual(oldest, obs_group)
        self.assertNotEqual(obs_group, oldest)

        # Add some headers and check that sorting still works
        obs_group.extend(self._files_to_headers(self.hsc_files))
        self.assertEqual(len(obs_group), 5)
        self.assertEqual(obs_group.newest(), obs_group[3])

        instruments = obs_group.property_values("instrument")
        self.assertEqual(instruments, {"HSC", "DECam"})

        # Check that simplified form round trips
        self.assertEqual(ObservationGroup.from_simple(obs_group.to_simple()), obs_group)

    def test_fits_group(self):
        headers = self._files_to_headers(self.decam_files)

        obs_group = ObservationGroup(headers)
        cards, comments = group_to_fits(obs_group)

        expected = {
            "INSTRUME": "DECam",
            "TIMESYS": "TAI",
            "DATE-OBS": "2012-12-11T22:07:07.859",
            "MJD-OBS": 56272.92161874134,
            "DATE-BEG": "2012-12-11T22:07:07.859",
            "MJD-BEG": 56272.92161874134,
            "DATE-END": "2015-02-20T00:50:11.000",
            "MJD-END": 57073.034849537034,
            "DATE-AVG": "2014-01-15T23:28:39.430",
            "MJD-AVG": 56672.97823413919,
        }
        self.assertEqual(cards, expected)

    def test_fits_info(self):
        header = self._files_to_headers(self.decam_files)[0]
        obs_group = ObservationGroup([header])
        cards, comments = info_to_fits(obs_group[0])

        expected = {
            "INSTRUME": "DECam",
            "TIMESYS": "TAI",
            "MJD-AVG": 56536.25417681625,
            "MJD-END": 56536.25591435185,
            "MJD-OBS": 56536.25243928065,
            "MJD-BEG": 56536.25243928065,
            "DATE-OBS": "2013-09-01T06:03:30.754",
            "DATE-BEG": "2013-09-01T06:03:30.754",
            "DATE-AVG": "2013-09-01T06:06:00.877",
            "DATE-END": "2013-09-01T06:08:31.000",
        }
        self.assertEqual(cards, expected)


if __name__ == "__main__":
    unittest.main()
