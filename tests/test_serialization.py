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

"""Test that archival ObservationInfo JSON files can be read back.

The JSON files in ``tests/data`` are stored long term, so we need to
guarantee that future versions of the code can still deserialize them.
"""

import os
import unittest
from collections.abc import Mapping
from typing import Any

import astropy.coordinates
import astropy.time
import astropy.units as u
import pydantic

from astro_metadata_translator import (
    ObservationInfo,
    PropertyDefinition,
    StubTranslator,
)

TESTDIR = os.path.abspath(os.path.dirname(__file__))
DATADIR = os.path.join(TESTDIR, "data")


# Auto-registers via StubTranslator.__init_subclass__ when this module is
# imported, so the "_translator": "fixture" reference in the archival JSON
# resolves and its extension definitions are available.
class _FixtureTranslator(StubTranslator):
    name = "fixture"
    supported_instrument = "fixture"
    extensions = dict(
        number=PropertyDefinition("A number", int),
        foo=PropertyDefinition("A string", str),
    )

    @classmethod
    def can_translate(cls, header: Mapping[str, Any], filename: str | None = None) -> bool:
        # Never used for header-based auto-detection; this translator only
        # exists so the "_translator": "fixture" reference in the archival
        # JSON can be resolved.
        return False


class ObservationInfoSerializationTestCase(unittest.TestCase):
    """Read archival ObservationInfo JSON fixtures and check the values."""

    def _check_core_properties(self, obsinfo: ObservationInfo) -> None:
        self.assertEqual(obsinfo.telescope, "Test Telescope")
        self.assertEqual(obsinfo.instrument, "TestCam")
        self.assertEqual(obsinfo.exposure_id, 987654321)
        self.assertEqual(obsinfo.visit_id, 987654000)
        self.assertEqual(obsinfo.physical_filter, "r")
        self.assertEqual(obsinfo.detector_num, 42)
        self.assertEqual(obsinfo.detector_name, "S00")
        self.assertEqual(obsinfo.detector_unique_name, "R22_S00")
        self.assertEqual(obsinfo.detector_serial, "ITL-3800C-145")
        self.assertEqual(obsinfo.detector_group, "R22")
        self.assertEqual(obsinfo.detector_exposure_id, 98765432142)
        self.assertEqual(obsinfo.object, "HD12345")
        self.assertEqual(obsinfo.science_program, "2025A-001")
        self.assertEqual(obsinfo.observation_type, "science")
        self.assertEqual(obsinfo.observation_id, "TC_O_20250521_000123")
        self.assertEqual(obsinfo.observation_reason, "science")
        self.assertEqual(obsinfo.exposure_group, "20250521_000123")
        self.assertEqual(obsinfo.observing_day, 20250520)
        self.assertEqual(obsinfo.observation_counter, 123)
        self.assertEqual(obsinfo.group_counter_start, 123)
        self.assertEqual(obsinfo.group_counter_end, 123)
        self.assertEqual(obsinfo.boresight_rotation_coord, "sky")
        self.assertEqual(obsinfo.has_simulated_content, False)
        self.assertEqual(obsinfo.can_see_sky, True)

        self.assertAlmostEqual(obsinfo.boresight_airmass, 1.1037)
        self.assertAlmostEqual(obsinfo.relative_humidity, 42.5)

        self.assertIsInstance(obsinfo.location, astropy.coordinates.EarthLocation)
        self.assertAlmostEqual(obsinfo.location.height.to_value(u.m), 2663.0, places=3)

        self.assertIsInstance(obsinfo.datetime_begin, astropy.time.Time)
        self.assertEqual(obsinfo.datetime_begin.scale, "tai")
        self.assertEqual(
            obsinfo.datetime_begin.tai.isot,
            "2025-05-21T01:23:45.000",
        )
        self.assertEqual(
            obsinfo.datetime_end.tai.isot,
            "2025-05-21T01:24:15.000",
        )

        self.assertAlmostEqual(obsinfo.exposure_time.to_value(u.s), 30.0)
        self.assertAlmostEqual(obsinfo.exposure_time_requested.to_value(u.s), 30.0)
        self.assertAlmostEqual(obsinfo.dark_time.to_value(u.s), 31.0)
        self.assertAlmostEqual(obsinfo.focus_z.to_value(u.m), 0.001)
        self.assertAlmostEqual(obsinfo.temperature.to_value(u.K), 283.15)
        self.assertAlmostEqual(obsinfo.pressure.to_value(u.hPa), 750.0)
        self.assertAlmostEqual(obsinfo.boresight_rotation_angle.to_value(u.deg), 45.0)
        self.assertEqual(round(obsinfo.observing_day_offset.to_value(u.s)), 43200)

        self.assertIsInstance(obsinfo.tracking_radec, astropy.coordinates.SkyCoord)
        self.assertAlmostEqual(obsinfo.tracking_radec.icrs.ra.to_value(u.deg), 123.456)
        self.assertAlmostEqual(obsinfo.tracking_radec.icrs.dec.to_value(u.deg), -45.678)

        self.assertIsInstance(obsinfo.altaz_begin, astropy.coordinates.AltAz)
        self.assertAlmostEqual(obsinfo.altaz_begin.az.to_value(u.deg), 110.0)
        self.assertAlmostEqual(obsinfo.altaz_begin.alt.to_value(u.deg), 65.0)
        self.assertAlmostEqual(obsinfo.altaz_end.az.to_value(u.deg), 110.5)
        self.assertAlmostEqual(obsinfo.altaz_end.alt.to_value(u.deg), 64.5)

    def test_read_full(self) -> None:
        """Read the fixture with every core property populated."""
        path = os.path.join(DATADIR, "obsinfo-full.json")
        with open(path) as fh:
            obsinfo = ObservationInfo.from_json(fh.read())

        self.assertIsInstance(obsinfo, ObservationInfo)
        self._check_core_properties(obsinfo)
        # No extensions present in this fixture.
        self.assertEqual(obsinfo.extensions, {})

        # Translator class is empty.
        self.assertIsNone(obsinfo._translator)
        self.assertEqual(obsinfo.translator_class_name, "<None>")

        # Round-tripping through JSON should be stable.
        round_tripped = ObservationInfo.from_json(obsinfo.to_json())
        self.assertEqual(round_tripped, obsinfo)

        # Should not include a spurious translator name.
        self.assertNotIn("_translator", obsinfo.model_dump())

    def test_read_full_with_extensions(self) -> None:
        """Read the fixture that also contains extension keys."""
        path = os.path.join(DATADIR, "obsinfo-full-extensions.json")
        with open(path) as fh:
            obsinfo = ObservationInfo.from_json(fh.read())

        self.assertIsInstance(obsinfo, ObservationInfo)
        self._check_core_properties(obsinfo)

        # Extension properties round-trip via the "fixture" translator name
        # stored in the JSON.
        self.assertEqual(set(obsinfo.extensions), {"number", "foo"})
        self.assertEqual(obsinfo.ext_foo, "bar")
        self.assertEqual(obsinfo.ext_number, 12345)

        self.assertIsInstance(obsinfo._translator, _FixtureTranslator)
        self.assertEqual(obsinfo.translator_name, "fixture")

        round_tripped = ObservationInfo.from_json(obsinfo.to_json())
        self.assertEqual(round_tripped, obsinfo)
        self.assertEqual(round_tripped.ext_foo, "bar")
        self.assertEqual(round_tripped.ext_number, 12345)

        # Should dump the translator name.
        self.assertIsInstance(round_tripped._translator, _FixtureTranslator)
        self.assertEqual(round_tripped.translator_name, "fixture")
        self.assertIn("_translator", obsinfo.model_dump())

        # If there are extension properties the translator must be registered
        # so that the extension definitions are available. The error must
        # identify the missing translator as the cause rather than simply
        # reporting an unknown property.
        path = os.path.join(DATADIR, "obsinfo-extensions-unknown-translator.json")
        with open(path) as fh:
            json_str = fh.read()
        with self.assertRaises(KeyError) as cm:
            ObservationInfo.from_json(json_str)
        message = str(cm.exception)
        self.assertIn("fixture_unknown", message)
        self.assertIn("not registered", message)
        self.assertIn("ext_foo", message)

    def test_read_unregistered(self) -> None:
        """Read a serialization with unregistered translator name."""
        path = os.path.join(DATADIR, "obsinfo-unregistered.json")
        with open(path) as fh:
            obsinfo = ObservationInfo.from_json(fh.read())

        self.assertIsInstance(obsinfo, ObservationInfo)
        self._check_core_properties(obsinfo)
        # No extensions present in this fixture.
        self.assertEqual(obsinfo.extensions, {})

        self.assertIsNone(obsinfo._translator, None)
        self.assertEqual(obsinfo.translator_name, "Unregistered")

        # Round-tripping through JSON should be stable.
        round_tripped = ObservationInfo.from_json(obsinfo.to_json())
        self.assertEqual(round_tripped, obsinfo)
        self.assertEqual(round_tripped.translator_name, "Unregistered")

    def test_nested_in_pydantic_model(self) -> None:
        """Nesting in another Pydantic model preserves the wire format.

        Nested serialization goes through Pydantic directly and does not
        call ``to_simple`` / ``to_json``, so the alias behavior must come
        from ``serialize_by_alias`` in the model config rather than from
        the call-site kwargs.
        """

        class Wrapper(pydantic.BaseModel):
            obs: ObservationInfo

        path = os.path.join(DATADIR, "obsinfo-full-extensions.json")
        with open(path) as fh:
            obsinfo = ObservationInfo.from_json(fh.read())

        wrapper = Wrapper(obs=obsinfo)

        nested = wrapper.model_dump()["obs"]
        self.assertIn("_translator", nested)
        self.assertNotIn("translator_name", nested)
        # Extension keys remain flat under the ext_ prefix.
        self.assertEqual(nested["ext_foo"], "bar")
        self.assertEqual(nested["ext_number"], 12345)

        nested_json = wrapper.model_dump_json()
        round_tripped = Wrapper.model_validate_json(nested_json)
        self.assertEqual(round_tripped.obs, obsinfo)


if __name__ == "__main__":
    unittest.main()
