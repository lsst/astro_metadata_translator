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

"""Tests for JSON Schema generation from the ObservationInfo model."""

import json
import os.path
import unittest
from collections.abc import Mapping
from typing import Any

import astropy.time
import astropy.units as u

try:
    import jsonschema

    HAS_JSONSCHEMA = True
except ImportError:
    HAS_JSONSCHEMA = False

from astro_metadata_translator import (
    ObservationInfo,
    PropertyDefinition,
    StubTranslator,
    makeObservationInfo,
)

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class _SchemaTranslator(StubTranslator):
    """Test translator that defines extension properties."""

    name = "schema_dummy"
    supported_instrument = "schema_dummy"
    extensions = dict(
        number=PropertyDefinition("A number", int),
        foo=PropertyDefinition("A string", str),
    )
    _const_map = {
        "ext_foo": "bar",
        "ext_number": 17,
        "detector_name": "detA",
    }
    _trivial_map = {
        "instrument": "INSTRUME",
    }

    @classmethod
    def can_translate(cls, header: Mapping[str, Any], filename: str | None = None) -> bool:
        return "INSTRUME" in header and header["INSTRUME"] == "schema_dummy"


# Mapping of internal Field-style names that must NOT appear as keys in the
# generated JSON schema (these are runtime state, not part of the wire format).
_INTERNAL_FIELDS = frozenset({"filename", "translator_class_name", "extensions", "all_properties"})


class JsonSchemaTestCase(unittest.TestCase):
    """Test JSON Schema generation."""

    def test_schema_serialization_mode_generatable(self) -> None:
        """Serialization-mode schema must be generatable and non-trivial."""
        schema = ObservationInfo.model_json_schema(mode="serialization")
        self.assertEqual(schema.get("type"), "object")
        properties = schema.get("properties", {})
        # Schema should describe each core property of the ObservationInfo,
        # not be a trivial open-ended object.
        self.assertIn("telescope", properties)
        self.assertIn("instrument", properties)
        self.assertIn("exposure_time", properties)
        self.assertIn("datetime_begin", properties)
        self.assertIn("location", properties)

    def test_schema_validation_mode_generatable(self) -> None:
        """Validation-mode schema must also be generatable."""
        schema = ObservationInfo.model_json_schema(mode="validation")
        self.assertEqual(schema.get("type"), "object")

    def test_schema_excludes_internal_fields(self) -> None:
        """Internal runtime state must not appear in the schema."""
        for mode in ("serialization", "validation"):
            with self.subTest(mode=mode):
                schema = ObservationInfo.model_json_schema(mode=mode)
                properties = set(schema.get("properties", {}).keys())
                leaked = _INTERNAL_FIELDS & properties
                self.assertEqual(leaked, set(), f"Internal fields leaked into {mode} schema: {leaked}")

    def test_schema_translator_key(self) -> None:
        """The schema should describe the _translator metadata key."""
        schema = ObservationInfo.model_json_schema(mode="serialization")
        self.assertIn("_translator", schema.get("properties", {}))

    def test_schema_allows_ext_properties(self) -> None:
        """Schema should allow ext_* keys via patternProperties."""
        schema = ObservationInfo.model_json_schema(mode="serialization")
        pattern_props = schema.get("patternProperties", {})
        # Some pattern keyed by ext_ prefix should accept arbitrary content.
        self.assertTrue(
            any("ext_" in key for key in pattern_props),
            f"Expected an ext_-prefixed patternProperties entry, got {pattern_props}",
        )

    def test_schema_astropy_field_simple_form(self) -> None:
        """Astropy-typed fields should use their simple form in the schema."""
        schema = ObservationInfo.model_json_schema(mode="serialization")
        properties = schema["properties"]

        # exposure_time: Quantity -> float (seconds)
        exposure_time = properties["exposure_time"]
        self.assertNotIn("is-instance", json.dumps(exposure_time))
        self._assert_accepts_type(exposure_time, "number")

        # datetime_begin: Time -> array of 2 numbers (TAI JD)
        datetime_begin = properties["datetime_begin"]
        self._assert_accepts_type(datetime_begin, "array")

        # location: EarthLocation -> array of 3 numbers (geocentric m)
        location = properties["location"]
        self._assert_accepts_type(location, "array")

        # boresight_rotation_angle: Angle -> float (degrees)
        rot_angle = properties["boresight_rotation_angle"]
        self._assert_accepts_type(rot_angle, "number")

        # tracking_radec: SkyCoord -> array of 2 numbers
        tracking_radec = properties["tracking_radec"]
        self._assert_accepts_type(tracking_radec, "array")

        # observing_day_offset: TimeDelta -> int (seconds)
        offset = properties["observing_day_offset"]
        self._assert_accepts_type(offset, "integer")

    @unittest.skipUnless(HAS_JSONSCHEMA, "jsonschema package is not installed")
    def test_schema_validates_real_obsinfo(self) -> None:
        """Serialized ObservationInfo must validate against the schema."""
        reference = dict(
            boresight_airmass=1.5,
            focus_z=1.0 * u.mm,
            temperature=15 * u.deg_C,
            observation_type="bias",
            exposure_time=5 * u.ks,
            detector_num=32,
            datetime_begin=astropy.time.Time("2021-02-15T12:00:00", format="isot", scale="utc"),
        )
        obsinfo = makeObservationInfo(**reference)
        data = json.loads(obsinfo.to_json())

        schema = ObservationInfo.model_json_schema(mode="serialization")
        # This will raise if invalid.
        jsonschema.validate(instance=data, schema=schema)

    @unittest.skipUnless(HAS_JSONSCHEMA, "jsonschema package is not installed")
    def test_schema_validates_with_extensions(self) -> None:
        """An ObservationInfo with ext_* properties must validate."""
        obsinfo = makeObservationInfo(
            translator_class=_SchemaTranslator,
            ext_foo="bar",
            ext_number=42,
            instrument="schema_dummy",
            detector_name="detA",
        )
        data = json.loads(obsinfo.to_json())
        self.assertIn("ext_foo", data)
        self.assertIn("ext_number", data)

        schema = ObservationInfo.model_json_schema(mode="serialization")
        jsonschema.validate(instance=data, schema=schema)

    @unittest.skipUnless(HAS_JSONSCHEMA, "jsonschema package is not installed")
    def test_schema_validates_reference_files(self) -> None:
        """Reference JSON files captured before this ticket must still
        validate against the generated schema.

        Guards against accidental wire-format regressions.
        """
        schema = ObservationInfo.model_json_schema(mode="serialization")
        for filename in ("obsinfo-full.json", "obsinfo-full-extensions.json"):
            with self.subTest(filename=filename):
                path = os.path.join(TESTDIR, "data", filename)
                with open(path) as fh:
                    data = json.load(fh)
                jsonschema.validate(instance=data, schema=schema)

    def _assert_accepts_type(self, prop_schema: dict, expected_type: str) -> None:
        """Assert that a property schema accepts the given JSON type.

        The schema may be a direct ``{"type": ...}`` or a union via
        ``anyOf``/``oneOf`` (e.g. a nullable property).

        Parameters
        ----------
        prop_schema : `dict`
            The JSON schema entry for a single property.
        expected_type : `str`
            The JSON Schema type name that the property must accept.
        """
        if prop_schema.get("type") == expected_type:
            return
        for key in ("anyOf", "oneOf"):
            for branch in prop_schema.get(key, []):
                if branch.get("type") == expected_type:
                    return
            for branch in prop_schema.get(key, []):
                # Recurse on nested unions.
                if "anyOf" in branch or "oneOf" in branch:
                    try:
                        self._assert_accepts_type(branch, expected_type)
                        return
                    except AssertionError:
                        continue
        self.fail(
            f"Property schema does not accept type {expected_type!r}: {json.dumps(prop_schema, default=str)}"
        )


if __name__ == "__main__":
    unittest.main()
