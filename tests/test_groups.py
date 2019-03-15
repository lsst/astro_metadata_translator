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

import unittest
import os.path

from astro_metadata_translator.tests import read_test_file
from astro_metadata_translator import ObservationGroup
from astro_metadata_translator.serialize import group_to_fits

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class ObservationGroupTestCase(unittest.TestCase):
    datadir = os.path.join(TESTDIR, "data")

    def setUp(self):
        self.decam_files = ("fitsheader-decam.yaml",
                            "fitsheader-decam-0160496.yaml",
                            "fitsheader-decam-calexp-0412037_10.yaml")
        self.hsc_files = ("fitsheader-hsc-HSCA04090107.yaml",
                          "fitsheader-hsc.yaml")

    def _files_to_headers(self, files):
        return [read_test_file(os.path.join(self.datadir, f)) for f in files]

    def test_groups(self):
        headers = self._files_to_headers(self.decam_files)

        obs_group = ObservationGroup(headers)
        self.assertEqual(len(obs_group), 3)

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

        # Add some headers and check that sorting still works
        obs_group.extend(self._files_to_headers(self.hsc_files))
        self.assertEqual(len(obs_group), 5)
        self.assertEqual(obs_group.newest(), obs_group[3])

        instruments = obs_group.property_values("instrument")
        self.assertEqual(instruments, {"HSC", "DECam"})

    def test_fits_group(self):
        headers = self._files_to_headers(self.decam_files)

        obs_group = ObservationGroup(headers)
        cards, comments = group_to_fits(obs_group)

        expected = {'INSTRUME': 'DECam',
                    'TIMESYS': 'TAI',
                    'DATE-OBS': '2012-12-11T22:07:07.859',
                    'MJD-OBS': 56272.92161874134,
                    'DATE-END': '2015-02-20T00:50:11.000',
                    'MJD-END': 57073.034849537034,
                    'DATE-AVG': '2014-01-15T23:28:39.430',
                    'MJD-AVG': 56672.97823413919}
        self.assertEqual(cards, expected)


if __name__ == "__main__":
    unittest.main()
