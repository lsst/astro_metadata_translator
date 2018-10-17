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

from helper import MetadataAssertHelper


class DecamTestCase(unittest.TestCase, MetadataAssertHelper):

    def test_decam_translator(self):
        test_data = (("fitsheader-decam.yaml",
                      dict(telescope="CTIO 4.0-m telescope",
                           instrument="DECam",
                           boresight_rotation_coord="unknown",
                           dark_time=201.15662,
                           detector_exposure_id=22938825,
                           detector_name="S1",
                           detector_num=25,
                           exposure=229388,
                           exposure_time=200.0,
                           object="DES supernova hex SN-S1 tiling 22",
                           obsid="ct4m20130901t060255",
                           obstype="science",
                           physical_filter="z DECam SDSS c0004 9260.0 1520.0",
                           pressure=779.0*u.hPa,
                           relative_humidity=23.0,
                           science_program="2012B-0001",
                           temperature=11.9*u.deg_C,
                           visit=229388,
                           wcs_params=dict(max_sep=1.5))),
                     )
        for file, expected in test_data:
            self.assertObservationInfoFromYaml(file, **expected)


if __name__ == "__main__":
    unittest.main()
