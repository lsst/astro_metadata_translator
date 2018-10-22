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
from astropy.time import Time

from astro_metadata_translator import FitsTranslator, StubTranslator, ObservationInfo


class InstrumentTestTranslator(FitsTranslator, StubTranslator):
    """Simple FITS-like translator to test the infrastructure"""

    # Needs a name to be registered
    name = "TestTranslator"

    # Indicate the instrument this class understands
    supported_instrument = "SCUBA_test"

    # Some new mappings, including an override
    _trivial_map = {"foobar": "BAZ",
                    "telescope": "TELCODE",
                    "observation_id": "OBSID"}

    _const_map = {"format": "HDF5"}


class TranslatorTestCase(unittest.TestCase):

    def setUp(self):
        # Known simple header
        self.header = {"TELESCOP": "JCMT",
                       "TELCODE": "LSST",
                       "INSTRUME": "SCUBA_test",
                       "DATE-OBS": "2000-01-01T01:00:01.500",
                       "DATE-END": "2000-01-01T02:00:01.500",
                       "OBSGEO-X": "-5464588.84421314",
                       "OBSGEO-Y": "-2493000.19137644",
                       "OBSGEO-Z": "2150653.35350771",
                       "OBSID": "20000101_00002",
                       "BAZ": "bar"}

    def test_manual_translation(self):

        header = self.header
        translator = FitsTranslator(header)

        # Treat the header as standard FITS
        self.assertFalse(FitsTranslator.can_translate(header))
        self.assertEqual(translator.to_telescope(), "JCMT")
        self.assertEqual(translator.to_instrument(), "SCUBA_test")
        self.assertEqual(translator.to_datetime_begin(),
                         Time(header["DATE-OBS"], format="isot"))

        # Use the special test translator instead
        translator = InstrumentTestTranslator(header)
        self.assertTrue(InstrumentTestTranslator.can_translate(header))
        self.assertEqual(translator.to_telescope(), "LSST")
        self.assertEqual(translator.to_instrument(), "SCUBA_test")
        self.assertEqual(translator.to_format(), "HDF5")
        self.assertEqual(translator.to_foobar(), "bar")

    def test_translator(self):
        header = self.header

        # Specify a translation class
        with self.assertWarns(UserWarning):
            # Since the translator is incomplete it should issue warnings
            v1 = ObservationInfo(header, translator_class=InstrumentTestTranslator)
        self.assertEqual(v1.instrument, "SCUBA_test")
        self.assertEqual(v1.telescope, "LSST")

        # Now automated class
        with self.assertWarns(UserWarning):
            # Since the translator is incomplete it should issue warnings
            v1 = ObservationInfo(header)
        self.assertEqual(v1.instrument, "SCUBA_test")
        self.assertEqual(v1.telescope, "LSST")

        location = v1.location.to_geodetic()
        self.assertAlmostEqual(location.height.to("m").to_value(), 4123.0, places=1)

        # Check that headers have been removed
        new_hdr = v1.stripped_header()
        self.assertNotIn("INSTRUME", new_hdr)
        self.assertNotIn("OBSGEO-X", new_hdr)
        self.assertIn("TELESCOP", new_hdr)

        # Check the list of cards that were used
        used = v1.cards_used
        self.assertIn("INSTRUME", used)
        self.assertIn("OBSGEO-Y", used)
        self.assertNotIn("TELESCOP", used)


if __name__ == "__main__":
    unittest.main()
