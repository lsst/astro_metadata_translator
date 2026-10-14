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

import unittest
from typing import Any

from astro_metadata_translator.serialize.flat_dict import (
    FITS_FLAT_DICT_PREFIX,
    FITS_FLAT_DICT_SEP,
    detect_prefix_and_sep,
    find_flat_keys,
    find_translator_name,
    flatten_simple,
    make_flat_key,
    unflatten,
    validate_sep,
)


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


if __name__ == "__main__":
    unittest.main()
