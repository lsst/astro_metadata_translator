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
from astropy.time import Time

from lsst.obs.metadata import FitsTranslator, VisitInfo


class TestTranslator(FitsTranslator):

    # Needs a name to be registered
    name = "TestTranslator"

    # Indicate the instrument this class understands
    supportedInstrument = "SCUBA_test"

    # Some new mappings, including an override
    _trivialMap = {"foobar": "BAZ",
                   "telescope": "TELCODE"}

    _constMap = {"format": "HDF5"}


class BasicTestCase(unittest.TestCase):

    def setUp(self):
        # Known simple header
        self.header = {"TELESCOP": "JCMT",
                       "TELCODE": "LSST",
                       "INSTRUME": "SCUBA_test",
                       "DATE-OBS": "2000-01-01T01:00:01.500",
                       "DATE-END": "2000-01-01T02:00:01.500",
                       "BAZ": "bar"}

    def testBasicManualTranslation(self):

        header = self.header

        # Treat the header as standard FITS
        self.assertFalse(FitsTranslator.canTranslate(header))
        self.assertEqual(FitsTranslator.to_telescope(header), "JCMT")
        self.assertEqual(FitsTranslator.to_instrument(header), "SCUBA_test")
        self.assertEqual(FitsTranslator.to_datetime_begin(header),
                         Time(header["DATE-OBS"], format="isot"))

        # Use the special test translator instead
        self.assertTrue(TestTranslator.canTranslate(header))
        self.assertEqual(TestTranslator.to_telescope(header), "LSST")
        self.assertEqual(TestTranslator.to_instrument(header), "SCUBA_test")
        self.assertEqual(TestTranslator.to_format(header), "HDF5")
        self.assertEqual(TestTranslator.to_foobar(header), "bar")

    def testBasicTranslator(self):
        header = self.header

        # Specify a translation class
        v1 = VisitInfo(header, translator=TestTranslator)
        self.assertEqual(v1.instrument, "SCUBA_test")
        self.assertEqual(v1.telescope, "LSST")

        # Now automated class
        v1 = VisitInfo(header)
        self.assertEqual(v1.instrument, "SCUBA_test")
        self.assertEqual(v1.telescope, "LSST")
        print(v1.__dict__)


if __name__ == "__main__":
    unittest.main()
