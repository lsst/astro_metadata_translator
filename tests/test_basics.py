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

import astropy.time
import astropy.units as u

import astro_metadata_translator
from astro_metadata_translator import ObservationInfo, makeObservationInfo


class BasicTestCase(unittest.TestCase):
    def test_basic(self):
        version = astro_metadata_translator.__version__
        self.assertIsNotNone(version)

    def test_obsinfo(self):
        """Test construction of ObservationInfo without header."""
        obsinfo = makeObservationInfo(boresight_airmass=1.5, tracking_radec=None)
        self.assertIsInstance(obsinfo, ObservationInfo)
        self.assertIsNone(obsinfo.tracking_radec)
        self.assertAlmostEqual(obsinfo.boresight_airmass, 1.5)
        self.assertIsNone(obsinfo.observation_id)
        self.assertEqual(obsinfo.cards_used, set())
        self.assertEqual(obsinfo.stripped_header(), {})

        with self.assertRaises(TypeError):
            ObservationInfo.makeObservationInfo(boresight_airmass=1.5, observation_id=5)

        with self.assertRaises(KeyError):
            obsinfo = ObservationInfo.makeObservationInfo(unrecognized=1.5, keys="unknown")

    def test_simple(self):
        """Test that we can simplify an ObservationInfo."""

        reference = dict(
            boresight_airmass=1.5,
            focus_z=1.0 * u.mm,
            temperature=15 * u.deg_C,
            observation_type="bias",
            exposure_time=5 * u.ks,
            detector_num=32,
            datetime_begin=astropy.time.Time("2021-02-15T12:00:00", format="isot", scale="utc"),
        )

        obsinfo = makeObservationInfo(**reference)
        simple = obsinfo.to_simple()
        newinfo = ObservationInfo.from_simple(simple)
        self.assertEqual(obsinfo, newinfo)

        via_json = ObservationInfo.from_json(newinfo.to_json())
        self.assertEqual(via_json, obsinfo)


if __name__ == "__main__":
    unittest.main()
