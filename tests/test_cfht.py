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

import unittest
import astropy.units as u

from helper import read_test_file, MetadataAssertHelper
from astro_metadata_translator import ObservationInfo


class MegaPrimeTestCase(unittest.TestCase, MetadataAssertHelper):

    def test_megaprime_translator(self):
        test_data = (("fitsheader-megaprime.yaml",
                      dict(telescope="CFHT 3.6m",
                           instrument="MegaPrime",
                           boresight_rotation_coord="unknown",
                           dark_time=615.0*u.s,
                           detector_exposure_id=37398350,
                           detector_name="8352-15-3",
                           detector_num=2,
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
                     )
        for file, expected in test_data:
            self.assertObservationInfoFromYaml(file, **expected)

    def test_megaprime_stripping(self):
        header = read_test_file("fitsheader-megaprime.yaml")
        v1 = ObservationInfo(header)

        # Check that headers have been removed
        new_hdr = v1.stripped_header()
        self.assertNotIn("INSTRUME", new_hdr)
        self.assertNotIn("TELESCOP", new_hdr)
        self.assertIn("CCD", new_hdr)


if __name__ == "__main__":
    unittest.main()
