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
import astropy.units as u

from astro_metadata_translator.tests import MetadataAssertHelper

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class HscTestCase(unittest.TestCase, MetadataAssertHelper):
    datadir = os.path.join(TESTDIR, "data")

    def test_hsc_translator(self):
        test_data = (("fitsheader-hsc.yaml",
                      dict(telescope="Subaru",
                           instrument="HSC",
                           boresight_rotation_coord="sky",
                           dark_time=30.0*u.s,
                           detector_exposure_id=180804850,
                           detector_name="120",
                           detector_num=50,
                           detector_serial="120",
                           exposure_id=904024,
                           exposure_time=30.0*u.s,
                           object="STRIPE82L",
                           observation_id="HSCA90402400",
                           observation_type="science",
                           physical_filter="HSC-I",
                           pressure=621.7*u.hPa,
                           relative_humidity=33.1,
                           science_program="o13015",
                           temperature=272.35*u.K,
                           visit_id=904024,
                           )),
                     ("fitsheader-hsc-HSCA04090107.yaml",
                      dict(telescope="Subaru",
                           instrument="HSC",
                           boresight_rotation_coord="sky",
                           dark_time=150.0*u.s,
                           detector_exposure_id=8180037,
                           detector_name="061",
                           detector_num=37,
                           detector_serial="061",
                           exposure_id=40900,
                           exposure_time=150.0*u.s,
                           object="SSP-Wide",
                           observation_id="HSCA04090000",
                           observation_type="science",
                           physical_filter="HSC-R",
                           pressure=625.4*u.hPa,
                           relative_humidity=8.6,
                           science_program="o15426",
                           temperature=278.35*u.K,
                           visit_id=40900,
                           )),
                     )
        for file, expected in test_data:
            with self.subTest(f"Testing {file}"):
                self.assertObservationInfoFromYaml(file, dir=self.datadir, **expected)

    def test_suprimecam_translator(self):
        # In this case the airmass is average during observation
        # but it looks like ALTITUDE is from a different time so loosen amdelta
        test_data = (("fitsheader-suprimecam-CORR40535770.yaml",
                      dict(telescope="Subaru",
                           instrument="SuprimeCam",
                           boresight_rotation_coord="unknown",
                           dark_time=200.0*u.s,
                           detector_exposure_id=535770,
                           detector_name="w67c1",
                           detector_num=0,
                           detector_serial="w67c1",
                           exposure_id=53577,
                           exposure_time=200.0*u.s,
                           object="Ecliptic Deep Field",
                           observation_id="SUPE00535770",
                           observation_type="science",
                           physical_filter="W-S-R+",
                           pressure=621.5*u.hPa,
                           relative_humidity=4.9,
                           science_program="o07222",
                           temperature=273.15*u.K,
                           visit_id=53577,
                           wcs_params=dict(amdelta=0.015))),
                     )
        for file, expected in test_data:
            with self.subTest(f"Testing {file}"):
                self.assertObservationInfoFromYaml(file, dir=self.datadir, **expected)


if __name__ == "__main__":
    unittest.main()
