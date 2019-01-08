# This file is part of astro_metadata_translator.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os.path
import unittest
import astropy.units as u

from astro_metadata_translator.tests import read_test_file, MetadataAssertHelper
from astro_metadata_translator import ObservationInfo

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class MegaPrimeTestCase(unittest.TestCase, MetadataAssertHelper):
    datadir = os.path.join(TESTDIR, "data")

    def test_megaprime_translator(self):
        test_data = (("fitsheader-megaprime.yaml",
                      dict(telescope="CFHT 3.6m",
                           instrument="MegaPrime",
                           boresight_rotation_coord="unknown",
                           dark_time=615.0*u.s,
                           detector_exposure_id=37398350,
                           detector_name="Aubin",
                           detector_num=2,
                           detector_serial="8352-15-3",
                           exposure_id=1038843,
                           exposure_time=615.037*u.s,
                           object="w2.+2+2",
                           observation_id="1038843",
                           observation_type="science",
                           physical_filter="i.MP9702",
                           pressure=617.65*u.hPa,
                           relative_humidity=39.77,
                           science_program="08BL05",
                           temperature=0.9*u.deg_C,
                           visit_id=1038843,
                           )),
                     ("fitsheader-megaprime-calexp-849375-14.yaml",
                      dict(telescope="CFHT 3.6m",
                           instrument="MegaPrime",
                           boresight_rotation_coord="unknown",
                           dark_time=300.0*u.s,
                           detector_exposure_id=30577599,
                           detector_name="Andre",
                           detector_num=99,
                           detector_serial="8434-13-5",
                           exposure_id=849375,
                           exposure_time=300.202*u.s,
                           object="D3",
                           observation_id="849375",
                           observation_type="science",
                           physical_filter="r",
                           pressure=615.79*u.hPa,
                           relative_humidity=15.76,
                           science_program="06AL01",
                           temperature=1.75*u.deg_C,
                           visit_id=849375,
                           )),
                     )
        for file, expected in test_data:
            with self.subTest(f"Testing {file}"):
                self.assertObservationInfoFromYaml(file, self.datadir, **expected)

    def test_megaprime_stripping(self):
        header = read_test_file("fitsheader-megaprime.yaml", dir=self.datadir)
        v1 = ObservationInfo(header)

        # Check that headers have been removed
        new_hdr = v1.stripped_header()
        self.assertNotIn("INSTRUME", new_hdr)
        self.assertNotIn("TELESCOP", new_hdr)
        self.assertIn("CCD", new_hdr)


if __name__ == "__main__":
    unittest.main()
