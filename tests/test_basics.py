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

import astro_metadata_translator
from astro_metadata_translator import ObservationInfo, makeObservationInfo


class BasicTestCase(unittest.TestCase):

    def test_basic(self):
        version = astro_metadata_translator.__version__
        self.assertIsNotNone(version)

    def test_obsinfo(self):
        """Test construction of ObservationInfo without header."""
        obsinfo = makeObservationInfo(boresight_airmass=1.5)
        self.assertIsInstance(obsinfo, ObservationInfo)
        self.assertAlmostEqual(obsinfo.boresight_airmass, 1.5)
        self.assertIsNone(obsinfo.observation_id)
        self.assertEqual(obsinfo.cards_used, ())
        self.assertEqual(obsinfo.stripped_header(), {})

        with self.assertRaises(TypeError):
            ObservationInfo.makeObservationInfo(boresight_airmass=1.5, observation_id=5)

        with self.assertRaises(KeyError):
            obsinfo = ObservationInfo.makeObservationInfo(unrecognized=1.5, keys="unknown")


if __name__ == "__main__":
    unittest.main()
