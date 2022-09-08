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

from astropy.time import Time

from astro_metadata_translator import FitsTranslator, ObservationInfo, StubTranslator

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class InstrumentTestTranslator(FitsTranslator, StubTranslator):
    """Simple FITS-like translator to test the infrastructure"""

    # Needs a name to be registered
    name = "TestTranslator"

    # Indicate the instrument this class understands
    supported_instrument = "SCUBA_test"

    # Some new mappings, including an override
    _trivial_map = {
        "telescope": "TELCODE",
        "exposure_id": "EXPID",
        "relative_humidity": "HUMIDITY",
        "detector_name": "DETNAME",
        "observation_id": "OBSID",
    }

    # Add translator method to test joining
    def to_physical_filter(self):
        return self._join_keyword_values(["DETNAME", "HUMIDITY"], delim="_")


class MissingMethodsTranslator(FitsTranslator):
    """Translator class that does not implement all the methods."""

    pass


class TranslatorTestCase(unittest.TestCase):
    def setUp(self):
        # Known simple header
        self.header = {
            "TELESCOP": "JCMT",
            "TELCODE": "LSST",
            "INSTRUME": "SCUBA_test",
            "DATE-OBS": "2000-01-01T01:00:01.500",
            "DATE-END": "2000-01-01T02:00:01.500",
            "OBSGEO-X": "-5464588.84421314",
            "OBSGEO-Y": "-2493000.19137644",
            "OBSGEO-Z": "2150653.35350771",
            "OBSID": "20000101_00002",
            "EXPID": "22",  # Should cast to a number
            "DETNAME": 76,  # Should cast to a string
            "HUMIDITY": "55",  # Should cast to a float
            "BAZ": "bar",
        }

    def test_manual_translation(self):

        header = self.header
        translator = FitsTranslator(header)

        # Treat the header as standard FITS
        self.assertFalse(FitsTranslator.can_translate(header))
        self.assertEqual(translator.to_telescope(), "JCMT")
        self.assertEqual(translator.to_instrument(), "SCUBA_test")
        self.assertEqual(translator.to_datetime_begin(), Time(header["DATE-OBS"], format="isot"))

        # This class will issue warnings
        with self.assertLogs("astro_metadata_translator") as cm:

            class InstrumentTestTranslatorExtras(InstrumentTestTranslator):
                """Version of InstrumentTestTranslator with unexpected
                fields."""

                name = "InstrumentTestTranslatorExtras"
                _trivial_map = {"foobar": "BAZ"}
                _const_map = {"format": "HDF5"}

        self.assertIn("Unexpected trivial", cm.output[0])
        self.assertIn("Unexpected constant", cm.output[1])

        # Use the special test translator instead
        translator = InstrumentTestTranslatorExtras(header)
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
        self.assertEqual(v1.exposure_id, 22)
        self.assertIsInstance(v1.exposure_id, int)
        self.assertEqual(v1.detector_name, "76")
        self.assertEqual(v1.relative_humidity, 55.0)
        self.assertIsInstance(v1.relative_humidity, float)
        self.assertEqual(v1.physical_filter, "76_55")

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

        # Stringification
        summary = str(v1)
        self.assertIn("datetime_begin", summary)

        # Create with a subset of properties
        v2 = ObservationInfo(
            header,
            translator_class=InstrumentTestTranslator,
            subset={"telescope", "datetime_begin", "exposure_group"},
        )

        self.assertEqual(v2.telescope, v1.telescope)
        self.assertEqual(v2.datetime_begin, v2.datetime_begin)
        self.assertIsNone(v2.datetime_end)
        self.assertIsNone(v2.location)
        self.assertIsNone(v2.observation_id)

    def test_corrections(self):
        """Apply corrections before translation."""
        header = self.header

        # Specify a translation class
        with self.assertWarns(UserWarning):
            # Since the translator is incomplete it should issue warnings
            v1 = ObservationInfo(
                header,
                translator_class=InstrumentTestTranslator,
                search_path=[os.path.join(TESTDIR, "data", "corrections")],
            )

        # These values should match the expected translation
        self.assertEqual(v1.instrument, "SCUBA_test")
        self.assertEqual(v1.detector_name, "76")
        self.assertEqual(v1.relative_humidity, 55.0)
        self.assertIsInstance(v1.relative_humidity, float)
        self.assertEqual(v1.physical_filter, "76_55")

        # These two should be the "corrected" values
        self.assertEqual(v1.telescope, "AuxTel")
        self.assertEqual(v1.exposure_id, 42)

    def test_failures(self):
        header = {}

        with self.assertRaises(TypeError):
            ObservationInfo(header, translator_class=ObservationInfo)

        with self.assertRaises(ValueError):
            ObservationInfo(
                header, translator_class=InstrumentTestTranslator, subset={"definitely_not_known"}
            )

        with self.assertRaises(ValueError):
            ObservationInfo(
                header, translator_class=InstrumentTestTranslator, required={"definitely_not_known"}
            )

        with self.assertLogs("astro_metadata_translator"):
            with self.assertWarns(UserWarning):
                ObservationInfo(header, translator_class=InstrumentTestTranslator, pedantic=False)

        with self.assertRaises(KeyError):
            with self.assertWarns(UserWarning):
                ObservationInfo(header, translator_class=InstrumentTestTranslator, pedantic=True)

        # Pass in header where the key does exist but the value is bad.
        bad_header = {
            "OBSGEO-X": "not-float",
            "OBSGEO-Y": "not-float",
            "OBSGEO-Z": "not-float",
            "TELESCOP": "JCMT",
            "TELCODE": "LSST",
            "INSTRUME": "SCUBA_test",
        }

        with self.assertLogs("astro_metadata_translator"):
            with self.assertWarns(UserWarning):
                ObservationInfo(bad_header, translator_class=InstrumentTestTranslator, pedantic=False)

        with self.assertLogs("astro_metadata_translator"):
            with self.assertWarns(UserWarning):
                ObservationInfo(
                    header, translator_class=InstrumentTestTranslator, pedantic=False, filename="testfile1"
                )

        with self.assertRaises(KeyError):
            with self.assertWarns(UserWarning):
                ObservationInfo(
                    header, translator_class=InstrumentTestTranslator, pedantic=True, filename="testfile2"
                )

        with self.assertRaises(NotImplementedError):
            with self.assertLogs("astro_metadata_translator", level="WARN"):
                ObservationInfo(header, translator_class=MissingMethodsTranslator)

        with self.assertRaises(KeyError):
            with self.assertWarns(UserWarning):
                with self.assertLogs("astro_metadata_translator", level="WARN"):
                    ObservationInfo(
                        header,
                        translator_class=InstrumentTestTranslator,
                        pedantic=False,
                        required={"boresight_airmass"},
                    )


if __name__ == "__main__":
    unittest.main()
