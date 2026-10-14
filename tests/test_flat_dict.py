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

from astro_metadata_translator.serialize.flat_dict import (
    FITS_FLAT_DICT_PREFIX,
    FITS_FLAT_DICT_SEP,
    flatten_simple,
    make_flat_key,
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


if __name__ == "__main__":
    unittest.main()
