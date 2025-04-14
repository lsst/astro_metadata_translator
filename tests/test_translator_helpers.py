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

from astropy.coordinates import EarthLocation
from astropy.time import Time

from astro_metadata_translator import StubTranslator
from astro_metadata_translator.translators.helpers import altaz_from_degree_headers


class HelperTranslator(StubTranslator):
    """Base class for testing shadowing."""

    def to_instrument(self):
        return "BaseInstrument"

    def to_observation_type(self):
        return "flat"

    def to_observation_id(self):
        return "OBSID"

    def to_location(self):
        return EarthLocation.from_geodetic(0.0, 0.0)


class ScienceTranslator(HelperTranslator):
    """Translator corresponding to science observation."""

    def to_observation_type(self):
        return "science"


class HelperTestCase(unittest.TestCase):
    """Test translator helpers."""

    def assert_azel(self, azel, az, alt):
        self.assertAlmostEqual(float(azel.az.to_value("deg")), az)
        self.assertAlmostEqual(float(azel.alt.to_value("deg")), alt)

    def test_altaz(self):
        """Test altaz extraction."""
        translator = HelperTranslator({})
        now = Time.now()

        value = altaz_from_degree_headers(translator, [("ELSTART", "AZSTART")], now)
        self.assertIsNone(value)

        translator._header["AZSTART"] = -400.0
        translator._header["ELSTART"] = 45.0

        value = altaz_from_degree_headers(translator, [("ELSTART", "AZSTART")], now)
        self.assertIsNone(value)

        for alt, az in ((80.0, 300.0), (-1.0, 200.0)):
            translator._header["AZSTART"] = az
            translator._header["ELSTART"] = alt
            value = altaz_from_degree_headers(translator, [("ELSTART", "AZSTART")], now, min_alt=-5.0)
            self.assert_azel(value, az, alt)

        translator._header["AZSTART"] = 22.0
        translator._header["ELSTART"] = 34.0
        value = altaz_from_degree_headers(translator, [("AZ1", "EL1"), ("ELSTART", "AZSTART")], now)
        self.assert_azel(value, 22.0, 34.0)

        value = altaz_from_degree_headers(translator, [("ELSTART", "AZSTART")], now, is_zd={"ELSTART"})
        self.assert_azel(value, 22.0, 56.0)

        translator._header["ELSTART"] = -1.0
        with self.assertLogs():
            value = altaz_from_degree_headers(translator, [("ELSTART", "AZSTART")], now, min_alt=-0.5)
        self.assert_azel(value, 22.0, -0.5)

        translator._header["ELSTART"] = 90.2
        with self.assertLogs():
            value = altaz_from_degree_headers(translator, [("ELSTART", "AZSTART")], now, max_alt=95.0)
        self.assert_azel(value, 22.0, 90.0)

    def test_science_altaz(self):
        """Test science altaz."""
        translator = ScienceTranslator({})
        now = Time.now()

        with self.assertRaises(KeyError):
            altaz_from_degree_headers(translator, [("ELSTART", "AZSTART")], now)


if __name__ == "__main__":
    unittest.main()
