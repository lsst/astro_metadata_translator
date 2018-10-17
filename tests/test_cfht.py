# This file is part of obs_metadata.
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
import astropy.units.cds as cds

from helper import readTestFile, MetadataAssertHelper
from lsst.obs.metadata import ObservationInfo


class MegaPrimeTestCase(unittest.TestCase, MetadataAssertHelper):

    def testMegaPrimeTranslator(self):
        test_data = (("fitsheader-megaprime.yaml",
                      dict(telescope="CFHT 3.6m",
                           instrument="MegaPrime",
                           boresight_rotation_coord="unknown",
                           dark_time=615.0,
                           detector_exposure_id=37398350,
                           detector_name="8352-15-3",
                           detector_num=2,
                           exposure=1038843,
                           exposure_time=615.037,
                           object="w2.+2+2",
                           obsid="1038843",
                           obstype="science",
                           physical_filter="i.MP9702",
                           pressure=617.65*cds.mmHg,
                           relative_humidity=39.77,
                           science_program="08BL05",
                           temperature=0.9*u.deg_C,
                           visit=1038843,
                           )),
                     )
        for file, expected in test_data:
            self.assertObservationInfoFromYaml(file, **expected)

    def testMegaPrimeStripping(self):
        header = readTestFile("fitsheader-megaprime.yaml")
        v1 = ObservationInfo(header)

        # Check that headers have been removed
        newHdr = v1.strippedHeader()
        self.assertNotIn("INSTRUME", newHdr)
        self.assertNotIn("TELESCOP", newHdr)
        self.assertIn("CCD", newHdr)


if __name__ == "__main__":
    unittest.main()
