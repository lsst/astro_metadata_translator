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

import math
import unittest

import astropy.time
import astropy.units as u
from pydantic import BaseModel, ValidationError

import astro_metadata_translator
from astro_metadata_translator import ObservationInfo, makeObservationInfo


class ModelWithObsInfo(BaseModel):
    """Pydantic model that contains an ObservationInfo."""

    obs_info: ObservationInfo
    number: int


class BasicTestCase(unittest.TestCase):
    """Test basic metadata translation functionality."""

    def test_basic(self) -> None:
        version = astro_metadata_translator.__version__
        self.assertIsNotNone(version)

    def test_obsinfo(self) -> None:
        """Test construction of ObservationInfo without header."""
        obsinfo = makeObservationInfo(boresight_airmass=1.5, tracking_radec=None)
        self.assertIsInstance(obsinfo, ObservationInfo)
        self.assertIsNone(obsinfo.tracking_radec)
        self.assertIsNotNone(obsinfo.boresight_airmass)
        self.assertAlmostEqual(obsinfo.boresight_airmass, 1.5)
        self.assertIsNone(obsinfo.observation_id)
        self.assertEqual(obsinfo.cards_used, set())
        self.assertEqual(obsinfo.stripped_header(), {})

        # Check NaN equality.
        obsinfo1 = ObservationInfo(relative_humidity=math.nan, detector_name="det1")
        self.assertEqual(obsinfo1, obsinfo1)

        with self.assertRaises(TypeError):
            _ = obsinfo > obsinfo

        with self.assertRaises(TypeError):
            _ = obsinfo < obsinfo

        with self.assertRaises(TypeError):
            _ = obsinfo > 5

        with self.assertRaises(TypeError):
            _ = obsinfo < 5

        with self.assertRaises(AttributeError):
            obsinfo.boresight_airmass = 1.0

        with self.assertRaises(ValidationError):
            ObservationInfo.makeObservationInfo(boresight_airmass=1.5, observation_id=5)

        with self.assertRaises(KeyError):
            ObservationInfo.makeObservationInfo(unrecognized=1.5, keys="unknown")

        with self.assertRaises(TypeError):
            # Translator class must be a class, not instance.
            ObservationInfo.makeObservationInfo(translator_class={})

        with self.assertRaises(TypeError):
            # Check that the translator class is checked.
            ObservationInfo(boresight_airmass=1.5, translator_class=dict)

        with self.assertRaises(TypeError):
            ObservationInfo.makeObservationInfo(boresight_airmass=1.5, extensions=[])

    def test_simple(self) -> None:
        """Test that we can simplify an ObservationInfo."""
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
        simple = obsinfo.to_simple()
        newinfo = ObservationInfo.from_simple(simple)
        self.assertEqual(obsinfo, newinfo)

        via_json = ObservationInfo.from_json(newinfo.to_json())
        self.assertEqual(via_json, obsinfo)

        # Sorting should be allowed (even if they are the same date).
        _ = sorted([obsinfo, newinfo, via_json])

    def test_pydantic(self) -> None:
        """Test pydantic serialization."""
        obsinfo = makeObservationInfo(boresight_airmass=1.5, tracking_radec=None)

        jstr = obsinfo.model_dump_json()
        from_j = ObservationInfo.model_validate_json(jstr)
        self.assertEqual(from_j, obsinfo)

        # Also from bytes.
        from_j = ObservationInfo.model_validate_json(jstr.encode())
        self.assertEqual(from_j, obsinfo)

        # Check that you get back what you gave.
        new_obsinfo = ObservationInfo.model_validate(from_j)
        self.assertEqual(id(new_obsinfo), id(from_j))

        with self.assertRaises(TypeError):
            ObservationInfo.model_validate([])

        with self.assertRaises(KeyError) as cm:
            ObservationInfo.model_validate({"_translator": "Unknown"})
        self.assertIn("Unrecognized translator", str(cm.exception))

        with self.assertRaises(KeyError) as cm:
            ObservationInfo.model_validate({"random": "Unknown"})
        self.assertIn("Unrecognized property", str(cm.exception))

        with self.assertRaises(ValidationError) as cm:
            ObservationInfo.model_validate({"detector_name": []})
        self.assertIn("should be a valid string", str(cm.exception))

    def test_embedded_pydantic(self) -> None:
        """Test that ObservationInfo can be used inside another model."""
        obsinfo = makeObservationInfo(boresight_airmass=1.5, detector_name="DETECTOR")

        test = ModelWithObsInfo(obs_info=obsinfo, number=42)

        json_str = test.model_dump_json()
        print(json_str)
        from_j = ModelWithObsInfo.model_validate_json(json_str)
        self.assertEqual(from_j.number, 42)
        self.assertEqual(from_j.obs_info.detector_name, "DETECTOR")


if __name__ == "__main__":
    unittest.main()
