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

"""Tests for the flat dict serialization of ObservationInfo."""

import math
import os
import tempfile
import unittest
from collections.abc import Mapping
from typing import Any

import astropy.time
import astropy.units as u
from astropy.io import fits

from astro_metadata_translator import (
    FITS_FLAT_DICT_PREFIX,
    FITS_FLAT_DICT_SEP,
    ObservationInfo,
    PropertyDefinition,
    StubTranslator,
    makeObservationInfo,
)
from astro_metadata_translator.serialize.flat_dict import (
    detect_prefix_and_sep,
    find_flat_keys,
    find_translator_name,
    flatten_simple,
    make_flat_key,
    unflatten,
    validate_sep,
)
from astro_metadata_translator.tests import read_test_file

TESTDIR = os.path.abspath(os.path.dirname(__file__))
DATADIR = os.path.join(TESTDIR, "data")

HEADER_FILES = (
    "fitsheader-decam.yaml",
    "fitsheader-hsc.yaml",
    "fitsheader-megaprime.yaml",
    "fitsheader-sdss-fpC-006377-g4-0399.yaml",
    "fitsheader-suprimecam-CORR40535770.yaml",
)

PREFIX_FORMS = (("", "."), ("obs_info", "."), (FITS_FLAT_DICT_PREFIX, FITS_FLAT_DICT_SEP))


# Registers itself through StubTranslator.__init_subclass__ so that
# ``_translator: "flatdict"`` resolves and its extensions are available.
class _FlatDictTranslator(StubTranslator):
    name = "flatdict"
    supported_instrument = "flatdict"
    extensions = dict(
        number=PropertyDefinition("A number", int),
        foo=PropertyDefinition("A string", str),
        pair=PropertyDefinition("A pair of floats", tuple),
    )

    @classmethod
    def can_translate(cls, header: Mapping[str, Any], filename: str | None = None) -> bool:
        return False


def _read_obsinfo(filename: str) -> ObservationInfo:
    return ObservationInfo(read_test_file(filename, dir=DATADIR))


class FlattenTestCase(unittest.TestCase):
    """Test flattening of the simple form."""

    def test_constants(self) -> None:
        self.assertEqual(FITS_FLAT_DICT_PREFIX, "OBS_INFO")
        self.assertEqual(FITS_FLAT_DICT_SEP, " ")

    def test_validate_sep(self) -> None:
        for good in (".", " ", "-", ":"):
            validate_sep(good)
        for bad in ("", "..", "a", "Z", "1", "_"):
            with self.assertRaises(ValueError):
                validate_sep(bad)

    def test_make_flat_key(self) -> None:
        self.assertEqual(make_flat_key("instrument", prefix="", sep=".", use_upper=False), "instrument")
        self.assertEqual(make_flat_key("instrument", prefix="", sep=".", use_upper=True), "INSTRUMENT")
        self.assertEqual(
            make_flat_key("datetime_begin_0", prefix="obs_info", sep=".", use_upper=False),
            "obs_info.datetime_begin_0",
        )
        self.assertEqual(
            make_flat_key("datetime_begin_0", prefix="OBS_INFO", sep=" ", use_upper=True),
            "OBS_INFO DATETIME_BEGIN_0",
        )
        # The prefix is emitted exactly as given.
        self.assertEqual(
            make_flat_key("instrument", prefix="obs_info", sep=" ", use_upper=True),
            "obs_info INSTRUMENT",
        )

    def test_flatten_scalars_no_prefix(self) -> None:
        simple = {
            "instrument": "HSC",
            "exposure_id": 42,
            "boresight_airmass": 1.2,
            "can_see_sky": True,
            "_translator": "HSC",
        }
        self.assertEqual(
            flatten_simple(simple),
            {
                "instrument": "HSC",
                "exposure_id": 42,
                "boresight_airmass": 1.2,
                "can_see_sky": True,
                "translator": "HSC",
            },
        )

    def test_flatten_sequences(self) -> None:
        simple = {"datetime_begin": (2460000.5, 0.25), "location": [1.0, 2.0, 3.0], "ext_pair": (1, 2)}
        self.assertEqual(
            flatten_simple(simple),
            {
                "datetime_begin_0": 2460000.5,
                "datetime_begin_1": 0.25,
                "location_0": 1.0,
                "location_1": 2.0,
                "location_2": 3.0,
                "ext_pair_0": 1,
                "ext_pair_1": 2,
            },
        )

    def test_flatten_prefix_and_case(self) -> None:
        simple = {"datetime_begin": (1.0, 2.0), "instrument": "X"}
        self.assertEqual(
            flatten_simple(simple, prefix="obs_info"),
            {
                "obs_info.datetime_begin_0": 1.0,
                "obs_info.datetime_begin_1": 2.0,
                "obs_info.instrument": "X",
            },
        )
        self.assertEqual(
            flatten_simple(simple, prefix="OBS_INFO", sep=" "),
            {
                "OBS_INFO DATETIME_BEGIN_0": 1.0,
                "OBS_INFO DATETIME_BEGIN_1": 2.0,
                "OBS_INFO INSTRUMENT": "X",
            },
        )
        self.assertEqual(
            flatten_simple(simple, prefix="OBS_INFO", sep=" ", use_upper=False),
            {
                "OBS_INFO datetime_begin_0": 1.0,
                "OBS_INFO datetime_begin_1": 2.0,
                "OBS_INFO instrument": "X",
            },
        )
        self.assertEqual(
            flatten_simple(simple, use_upper=True),
            {"DATETIME_BEGIN_0": 1.0, "DATETIME_BEGIN_1": 2.0, "INSTRUMENT": "X"},
        )

    def test_flatten_rejects_bad_values(self) -> None:
        with self.assertRaises(TypeError):
            flatten_simple({"ext_bad": {"a": 1}})
        with self.assertRaises(TypeError):
            flatten_simple({"ext_nested": [[1.0, 2.0]]})
        with self.assertRaises(TypeError):
            flatten_simple({"ext_none": None})

    def test_flatten_rejects_bad_sep(self) -> None:
        with self.assertRaises(ValueError):
            flatten_simple({"instrument": "X"}, prefix="p", sep="_")
        with self.assertRaises(ValueError):
            flatten_simple({"instrument": "X"}, sep="ab")


class DetectTestCase(unittest.TestCase):
    """Test prefix detection and key collection."""

    PREFIXED = {
        "OBS_INFO DATETIME_BEGIN_0": 2460000.5,
        "OBS_INFO DATETIME_BEGIN_1": 0.25,
        "OBS_INFO INSTRUMENT": "HSC",
        "OBS_INFO OBJECT": "real",
        "OBS_INFO TRANSLATOR": "HSC",
        "OBS_INFO EXT_FOO": "bar",
        "OBS_INFO EXT_PAIR_0": 1,
        "OBJECT": "decoy",
        "TEMPERATURE": 3.0,
        "EXPTIME": 30.0,
    }

    def test_detect_nothing(self) -> None:
        self.assertIsNone(detect_prefix_and_sep({}))
        self.assertIsNone(detect_prefix_and_sep({"EXPTIME": 30.0, "FILTER": "g"}))

    def test_detect_decoys_without_anchor(self) -> None:
        # Core names appear but nothing indexed, so this is not content.
        self.assertIsNone(detect_prefix_and_sep({"OBJECT": "M31", "TEMPERATURE": 3.0}))
        self.assertIsNone(detect_prefix_and_sep({"instrument": "X", "object": "Y"}))

    def test_detect_plain(self) -> None:
        flat = {"datetime_begin_0": 1.0, "datetime_begin_1": 2.0, "instrument": "X"}
        self.assertEqual(detect_prefix_and_sep(flat), ("", ""))

    def test_detect_prefixed_with_decoys(self) -> None:
        self.assertEqual(detect_prefix_and_sep(self.PREFIXED), ("obs_info", " "))

    def test_detect_prefix_containing_sep(self) -> None:
        flat = {"LSST OBS_INFO DATETIME_BEGIN_0": 1.0, "LSST OBS_INFO EXPOSURE_ID": 5}
        self.assertEqual(detect_prefix_and_sep(flat), ("lsst obs_info", " "))

    def test_detect_dotted(self) -> None:
        flat = {"lsst.obs_info.datetime_begin_0": 1.0, "lsst.obs_info.exposure_id": 5}
        self.assertEqual(detect_prefix_and_sep(flat), ("lsst.obs_info", "."))

    def test_detect_majority(self) -> None:
        flat = {"a.datetime_begin_0": 1.0, "a.instrument": "x", "b.datetime_begin_0": 1.0}
        self.assertEqual(detect_prefix_and_sep(flat), ("a", "."))

    def test_detect_tie(self) -> None:
        flat = {"a.datetime_begin_0": 1.0, "b.datetime_begin_0": 1.0}
        with self.assertRaises(ValueError):
            detect_prefix_and_sep(flat)

    def test_find_keys_auto(self) -> None:
        found = find_flat_keys(self.PREFIXED)
        self.assertEqual(
            found,
            {
                "OBS_INFO DATETIME_BEGIN_0": "datetime_begin_0",
                "OBS_INFO DATETIME_BEGIN_1": "datetime_begin_1",
                "OBS_INFO INSTRUMENT": "instrument",
                "OBS_INFO OBJECT": "object",
                "OBS_INFO TRANSLATOR": "translator",
                "OBS_INFO EXT_FOO": "ext_foo",
                "OBS_INFO EXT_PAIR_0": "ext_pair_0",
            },
        )

    def test_find_keys_auto_nothing(self) -> None:
        self.assertEqual(find_flat_keys({"OBJECT": "M31"}), {})

    def test_find_keys_explicit(self) -> None:
        found = find_flat_keys(self.PREFIXED, prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        self.assertEqual(len(found), 7)
        self.assertNotIn("OBJECT", found)
        # Explicit empty prefix needs no indexed anchor.
        self.assertEqual(
            find_flat_keys({"instrument": "X", "EXPTIME": 1}, prefix=""), {"instrument": "instrument"}
        )
        # Explicit prefix with the wrong separator finds nothing.
        self.assertEqual(find_flat_keys(self.PREFIXED, prefix="OBS_INFO", sep="."), {})
        # Explicit prefix is matched case-insensitively.
        found = find_flat_keys(self.PREFIXED, prefix="obs_info", sep=" ")
        self.assertEqual(len(found), 7)

    def test_find_keys_ignores_non_string_keys(self) -> None:
        flat: dict[Any, Any] = {1: "x", "datetime_begin_0": 1.0, "datetime_begin_1": 2.0}
        expected = {"datetime_begin_0": "datetime_begin_0", "datetime_begin_1": "datetime_begin_1"}
        self.assertEqual(find_flat_keys(flat), expected)
        self.assertEqual(find_flat_keys(flat, prefix=""), expected)

    def test_find_keys_ignores_unknown_names(self) -> None:
        flat = {"obs_info.datetime_begin_0": 1.0, "obs_info.not_a_property": 2, "obs_info.ext_": 3}
        self.assertEqual(
            find_flat_keys(flat, prefix="obs_info"), {"obs_info.datetime_begin_0": "datetime_begin_0"}
        )

    def test_find_keys_duplicate(self) -> None:
        flat = {"OBS_INFO INSTRUMENT": "a", "obs_info instrument": "b", "OBS_INFO DATETIME_BEGIN_0": 1.0}
        with self.assertRaises(ValueError):
            find_flat_keys(flat, prefix="OBS_INFO", sep=" ")

    def test_find_keys_bad_sep(self) -> None:
        with self.assertRaises(ValueError):
            find_flat_keys({"p_instrument": "x"}, prefix="p", sep="_")


class UnflattenTestCase(unittest.TestCase):
    """Test rebuilding the simple form from flat keys."""

    def test_unflatten_core(self) -> None:
        flat = {
            "OBS_INFO DATETIME_BEGIN_1": 0.25,
            "OBS_INFO DATETIME_BEGIN_0": 2460000.5,
            "OBS_INFO LOCATION_0": 1.0,
            "OBS_INFO LOCATION_1": 2.0,
            "OBS_INFO LOCATION_2": 3.0,
            "OBS_INFO INSTRUMENT": "HSC",
            "OBS_INFO TRANSLATOR": "HSC",
            "OBJECT": "decoy",
        }
        flat_keys = find_flat_keys(flat)
        self.assertEqual(find_translator_name(flat, flat_keys), "HSC")
        self.assertEqual(
            unflatten(flat, flat_keys),
            {
                "datetime_begin": (2460000.5, 0.25),
                "location": (1.0, 2.0, 3.0),
                "instrument": "HSC",
                "_translator": "HSC",
            },
        )

    def test_unflatten_no_translator(self) -> None:
        flat = {"datetime_begin_0": 1.0, "datetime_begin_1": 2.0}
        flat_keys = find_flat_keys(flat)
        self.assertIsNone(find_translator_name(flat, flat_keys))
        self.assertEqual(unflatten(flat, flat_keys), {"datetime_begin": (1.0, 2.0)})

    def test_unflatten_gap(self) -> None:
        flat = {"datetime_begin_0": 1.0, "datetime_begin_2": 2.0}
        with self.assertRaises(ValueError):
            unflatten(flat, find_flat_keys(flat))
        flat = {"datetime_begin_1": 2.0}
        with self.assertRaises(ValueError):
            unflatten(flat, find_flat_keys(flat))

    def test_unflatten_expected_lengths(self) -> None:
        flat = {"datetime_begin_0": 1.0, "datetime_begin_1": 2.0, "location_0": 1.0, "location_1": 2.0}
        flat_keys = find_flat_keys(flat)
        lengths = {"datetime_begin": 2, "location": 3}
        # A truncated tuple is incomplete even though its indexes are
        # contiguous from 0.
        with self.assertRaises(ValueError):
            unflatten(flat, flat_keys, expected_lengths=lengths)
        # Without a declared length only contiguity is checked.
        self.assertEqual(unflatten(flat, flat_keys), {"datetime_begin": (1.0, 2.0), "location": (1.0, 2.0)})
        del flat["location_0"], flat["location_1"]
        self.assertEqual(
            unflatten(flat, find_flat_keys(flat), expected_lengths=lengths), {"datetime_begin": (1.0, 2.0)}
        )

    def test_unflatten_extensions(self) -> None:
        flat = {
            "obs_info.datetime_begin_0": 1.0,
            "obs_info.datetime_begin_1": 2.0,
            "obs_info.ext_foo": "bar",
            "obs_info.ext_pair_0": 1,
            "obs_info.ext_pair_1": 2,
            "obs_info.ext_other_1": 3,
        }
        flat_keys = find_flat_keys(flat, prefix="obs_info")
        # Known extension names restore their spelling and regroup.
        simple = unflatten(flat, flat_keys, extension_names=("Foo", "pair"))
        self.assertEqual(
            simple,
            {
                "datetime_begin": (1.0, 2.0),
                "ext_Foo": "bar",
                "ext_pair": (1, 2),
                "ext_other_1": 3,
            },
        )
        # Unknown extension names pass through as plain lower-case names.
        simple = unflatten(flat, flat_keys)
        self.assertEqual(
            simple,
            {
                "datetime_begin": (1.0, 2.0),
                "ext_foo": "bar",
                "ext_pair_0": 1,
                "ext_pair_1": 2,
                "ext_other_1": 3,
            },
        )


class ObservationInfoFlatDictTestCase(unittest.TestCase):
    """Test the flat dict methods on ObservationInfo."""

    def assert_round_trip(self, original: ObservationInfo, restored: ObservationInfo) -> None:
        # Compare the full simple forms rather than relying on ``==``,
        # which only compares the properties present on both sides.
        self.assertEqual(original.to_simple(), restored.to_simple())

    def test_round_trip_headers(self) -> None:
        for filename in HEADER_FILES:
            obsinfo = _read_obsinfo(filename)
            for prefix, sep in PREFIX_FORMS:
                with self.subTest(filename=filename, prefix=prefix):
                    flat = obsinfo.to_flat_dict(prefix=prefix, sep=sep)
                    for value in flat.values():
                        self.assertIsInstance(value, str | int | float | bool)
                    self.assertTrue(ObservationInfo.has_flat_dict(flat))
                    self.assert_round_trip(obsinfo, ObservationInfo.from_flat_dict(flat))
                    restored = ObservationInfo.from_flat_dict(flat, prefix=prefix, sep=sep)
                    self.assert_round_trip(obsinfo, restored)

    def test_fits_form_keys(self) -> None:
        obsinfo = _read_obsinfo("fitsheader-decam.yaml")
        flat = obsinfo.to_flat_dict(prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        self.assertIn("OBS_INFO DATETIME_BEGIN_0", flat)
        self.assertIn("OBS_INFO DATETIME_BEGIN_1", flat)
        self.assertIn("OBS_INFO TRACKING_RADEC_0", flat)
        self.assertIn("OBS_INFO LOCATION_2", flat)
        self.assertEqual(flat["OBS_INFO TRANSLATOR"], "DECam")
        self.assertNotIn("OBS_INFO _TRANSLATOR", flat)
        for key in flat:
            self.assertEqual(key, key.upper())

    def test_round_trip_hand_built(self) -> None:
        begin = astropy.time.Time("2025-01-01T00:00:00", scale="tai")
        obsinfo = makeObservationInfo(observation_id="A", exposure_time=30.0 * u.s, datetime_begin=begin)
        flat = obsinfo.to_flat_dict()
        self.assertNotIn("translator", flat)
        self.assert_round_trip(obsinfo, ObservationInfo.from_flat_dict(flat))

        # Without any tuple-valued property there is no indexed anchor so
        # auto-detection finds nothing and an explicit prefix is needed.
        obsinfo = makeObservationInfo(observation_id="A", exposure_time=30.0 * u.s)
        flat = obsinfo.to_flat_dict()
        self.assertFalse(ObservationInfo.has_flat_dict(flat))
        with self.assertRaises(ValueError):
            ObservationInfo.from_flat_dict(flat)
        self.assertTrue(ObservationInfo.has_flat_dict(flat, prefix=""))
        self.assert_round_trip(obsinfo, ObservationInfo.from_flat_dict(flat, prefix=""))

    def test_round_trip_extensions(self) -> None:
        begin = astropy.time.Time("2025-01-01T00:00:00", scale="tai")
        obsinfo = makeObservationInfo(
            translator_class=_FlatDictTranslator,
            instrument="flatdict",
            datetime_begin=begin,
            ext_number=42,
            ext_foo="bar",
            ext_pair=(1.5, 2.5),
        )
        for prefix, sep in PREFIX_FORMS:
            with self.subTest(prefix=prefix):
                flat = obsinfo.to_flat_dict(prefix=prefix, sep=sep)
                restored = ObservationInfo.from_flat_dict(flat)
                self.assertEqual(restored.translator_name, "flatdict")
                self.assertEqual(restored.ext_number, 42)
                self.assertEqual(restored.ext_foo, "bar")
                self.assertEqual(restored.ext_pair, (1.5, 2.5))
                self.assert_round_trip(obsinfo, restored)
        flat = obsinfo.to_flat_dict(prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        self.assertEqual(flat["OBS_INFO EXT_PAIR_0"], 1.5)
        self.assertEqual(flat["OBS_INFO EXT_PAIR_1"], 2.5)

    def test_has_and_strip(self) -> None:
        obsinfo = _read_obsinfo("fitsheader-hsc.yaml")
        flat = obsinfo.to_flat_dict(prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        n_flat = len(flat)
        header = {"SIMPLE": True, "OBJECT": "decoy", "TEMPERATURE": 3.0, **flat}
        self.assertTrue(ObservationInfo.has_flat_dict(header))
        self.assertTrue(
            ObservationInfo.has_flat_dict(header, prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        )
        ObservationInfo.strip_from_flat_dict(header)
        self.assertEqual(header, {"SIMPLE": True, "OBJECT": "decoy", "TEMPERATURE": 3.0})
        self.assertEqual(len(flat), n_flat)
        self.assertFalse(ObservationInfo.has_flat_dict(header))
        # Stripping again is a no-op.
        ObservationInfo.strip_from_flat_dict(header)
        self.assertEqual(len(header), 3)

    def test_from_flat_dict_errors(self) -> None:
        with self.assertRaises(ValueError):
            ObservationInfo.from_flat_dict({})
        with self.assertRaises(ValueError):
            ObservationInfo.from_flat_dict({"OBJECT": "M31"})

        obsinfo = _read_obsinfo("fitsheader-decam.yaml")
        flat = obsinfo.to_flat_dict(prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)

        gap = dict(flat)
        del gap["OBS_INFO DATETIME_BEGIN_1"]
        with self.assertRaises(ValueError):
            ObservationInfo.from_flat_dict(gap)

        corrupt = dict(flat)
        corrupt["OBS_INFO EXPOSURE_ID"] = "not an integer"
        # Validation errors propagate rather than being mistaken for absence.
        with self.assertRaises((TypeError, ValueError)):
            ObservationInfo.from_flat_dict(corrupt)
        self.assertTrue(ObservationInfo.has_flat_dict(corrupt))

    def test_explicit_prefix_wrong_sep(self) -> None:
        obsinfo = _read_obsinfo("fitsheader-decam.yaml")
        flat = obsinfo.to_flat_dict(prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        self.assertFalse(ObservationInfo.has_flat_dict(flat, prefix=FITS_FLAT_DICT_PREFIX, sep="."))
        with self.assertRaises(ValueError):
            ObservationInfo.from_flat_dict(flat, prefix=FITS_FLAT_DICT_PREFIX, sep=".")

    def test_to_flat_dict_bad_sep(self) -> None:
        obsinfo = _read_obsinfo("fitsheader-decam.yaml")
        with self.assertRaises(ValueError):
            obsinfo.to_flat_dict(prefix="p", sep="_")


class FitsHeaderTestCase(unittest.TestCase):
    """Test the flat dict through a real FITS header file."""

    def assert_fits_round_trip(self, original: ObservationInfo, restored: ObservationInfo) -> None:
        # FITS writers keep at least 15 significant digits of a float
        # rather than its shortest exact representation, so floats are
        # compared to that precision and everything else exactly.
        before = original.to_simple()
        after = restored.to_simple()
        self.assertEqual(set(before), set(after))
        for key, value in before.items():
            if isinstance(value, float):
                self.assertTrue(math.isclose(after[key], value, rel_tol=1e-14), key)
            elif isinstance(value, tuple) and isinstance(value[0], float):
                for element_after, element_before in zip(after[key], value, strict=True):
                    self.assertTrue(math.isclose(element_after, element_before, rel_tol=1e-14), key)
            else:
                self.assertEqual(after[key], value, key)

    def _write_and_read(self, header: fits.Header) -> fits.Header:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test.fits")
            fits.PrimaryHDU(header=header).writeto(path)
            return fits.getheader(path)

    def test_fits_header_round_trip(self) -> None:
        obsinfo = _read_obsinfo("fitsheader-decam.yaml")
        flat = obsinfo.to_flat_dict(prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)

        header = fits.Header()
        header["OBJECT"] = "decoy"
        header["HIERARCH TEMPERATURE"] = 3.0
        for key, value in flat.items():
            header[f"HIERARCH {key}"] = value
        read = self._write_and_read(header)

        self.assertTrue(ObservationInfo.has_flat_dict(read))
        self.assertTrue(
            ObservationInfo.has_flat_dict(read, prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        )

        restored = ObservationInfo.from_flat_dict(read)
        self.assert_fits_round_trip(obsinfo, restored)
        self.assertEqual(restored.object, obsinfo.object)
        self.assertNotEqual(restored.object, "decoy")
        self.assertIsInstance(restored.can_see_sky, bool)
        self.assertIsInstance(restored.has_simulated_content, bool)

        restored = ObservationInfo.from_flat_dict(read, prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        self.assert_fits_round_trip(obsinfo, restored)

        ObservationInfo.strip_from_flat_dict(read)
        self.assertFalse(ObservationInfo.has_flat_dict(read))
        self.assertFalse(any(key.startswith("OBS_INFO") for key in read))
        self.assertEqual(read["OBJECT"], "decoy")
        self.assertEqual(read["TEMPERATURE"], 3.0)

    def test_fits_header_without_content(self) -> None:
        header = fits.Header()
        header["OBJECT"] = "M31"
        header["HIERARCH TEMPERATURE"] = 3.0
        header["HIERARCH PRESSURE"] = 700.0
        read = self._write_and_read(header)
        n_cards = len(read)
        self.assertFalse(ObservationInfo.has_flat_dict(read))
        ObservationInfo.strip_from_flat_dict(read)
        self.assertEqual(len(read), n_cards)
        with self.assertRaises(ValueError):
            ObservationInfo.from_flat_dict(read)

    def test_nan_passes_through_plain_dict(self) -> None:
        obsinfo = makeObservationInfo(observation_id="A", relative_humidity=math.nan)
        flat = obsinfo.to_flat_dict()
        self.assertTrue(math.isnan(flat["relative_humidity"]))
        restored = ObservationInfo.from_flat_dict(flat, prefix="")
        self.assertIsNotNone(restored.relative_humidity)
        self.assertTrue(math.isnan(restored.relative_humidity))


if __name__ == "__main__":
    unittest.main()
