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

from lsst.obs.metadata import FitsTranslator, VisitInfo


class TestTranslator(FitsTranslator):

    # Needs a name to be registered
    name = "TestTranslator"

    # Indicate the instrument this class understands
    supportedInstrument = "SCUBA_test"

    # Some new mappings, including an override
    _unitMap = {"foobar": "BAZ",
                "telescope": "TELCODE"}

    _constMap = {"format": "HDF5"}


class BasicTestCase(unittest.TestCase):

    def setUp(self):
        # Known simple header
        self.header = {"TELESCOP": "JCMT",
                       "TELCODE": "LSST",
                       "INSTRUME": "SCUBA_test",
                       "BAZ": "bar"}

    def testBasicManualTranslation(self):

        header = self.header

        # Treat the header as standard FITS
        self.assertFalse(FitsTranslator.canTranslate(header))
        self.assertEqual(FitsTranslator.toTelescope(header), "JCMT")
        self.assertEqual(FitsTranslator.toInstrument(header), "SCUBA_test")

        # Use the special test translator instead
        self.assertTrue(TestTranslator.canTranslate(header))
        self.assertEqual(TestTranslator.toTelescope(header), "LSST")
        self.assertEqual(TestTranslator.toInstrument(header), "SCUBA_test")
        self.assertEqual(TestTranslator.toFormat(header), "HDF5")
        self.assertEqual(TestTranslator.toFoobar(header), "bar")

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


if __name__ == "__main__":
    unittest.main()
