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
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
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
        self.assertNotEqual(ObservationInfo(observation_id="A"), ObservationInfo(observation_id="B"))
        self.assertEqual(
            ObservationInfo(observation_id="A") == ObservationInfo(observation_id="A", instrument="HSC"),
            ObservationInfo(observation_id="A", instrument="HSC") == ObservationInfo(observation_id="A"),
        )

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

    def test_altaz_extraction(self) -> None:
        """Test that an AltAz embedded in a SkyCoord is extracted correctly."""
        # This test came originally from ip_diffim
        lsst_lat = -30.244639 * u.degree
        lsst_lon = -70.749417 * u.degree
        lsst_alt = 2663.0 * u.m
        lsst_temperature = 20.0 * u.Celsius
        lsst_humidity = 40.0  # in percent
        lsst_pressure = 73892.0 * u.pascal
        elevation = 45.0 * u.degree
        azimuth = 110.0 * u.degree
        rotangle = 30.0 * u.degree

        loc = EarthLocation(lat=lsst_lat, lon=lsst_lon, height=lsst_alt)
        airmass = 1.0 / math.sin(elevation.to_value())

        time = Time(2026.0, format="jyear", scale="tai")
        altaz = SkyCoord(
            alt=elevation,
            az=azimuth,
            obstime=time,
            frame="altaz",
            location=loc,
        )
        obsinfo = makeObservationInfo(
            location=loc,
            detector_exposure_id=12345678,
            datetime_begin=time,
            datetime_end=time,
            boresight_airmass=airmass,
            boresight_rotation_angle=rotangle,
            boresight_rotation_coord="sky",
            temperature=lsst_temperature,
            pressure=lsst_pressure,
            relative_humidity=lsst_humidity,
            tracking_radec=altaz.icrs,
            altaz_begin=altaz,
            observation_type="science",
        )
        self.assertIsInstance(obsinfo, ObservationInfo)


if __name__ == "__main__":
    unittest.main()
