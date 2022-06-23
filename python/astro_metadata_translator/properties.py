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

"""Properties calculated by this package.

Defines all properties in one place so that both `ObservationInfo` and
`MetadataTranslator` can use them.  In particular, the translator
base class can use knowledge of these properties to predefine translation
stubs with documentation attached, and `ObservationInfo` can automatically
define the getter methods.

"""
from __future__ import annotations

__all__ = (
    "PropertyDefinition",
    "PROPERTIES",
)

from dataclasses import dataclass
from typing import Any, Callable, Optional, Tuple

import astropy.coordinates
import astropy.time
import astropy.units

# Helper functions to convert complex types to simple form suitable
# for JSON serialization
# All take the complex type and return simple python form using str, float,
# int, dict, or list.
# All assume the supplied parameter is not None.


def earthlocation_to_simple(location: astropy.coordinates.EarthLocation) -> Tuple[float, ...]:
    """Convert EarthLocation to tuple.

    Parameters
    ----------
    location : `astropy.coordinates.EarthLocation`
        The location to simplify.

    Returns
    -------
    geocentric : `tuple` of (`float`, `float`, `float`)
        The geocentric location as three floats in meters.
    """
    geocentric = location.to_geocentric()
    return tuple(c.to_value(astropy.units.m) for c in geocentric)


def simple_to_earthlocation(simple: Tuple[float, ...], **kwargs: Any) -> astropy.coordinates.EarthLocation:
    """Convert simple form back to EarthLocation."""
    return astropy.coordinates.EarthLocation.from_geocentric(*simple, unit=astropy.units.m)


def datetime_to_simple(datetime: astropy.time.Time) -> Tuple[float, float]:
    """Convert Time to tuple.

    Parameters
    ----------
    datetime : `astropy.time.Time`
        The time to simplify.

    Returns
    -------
    mjds : `tuple` of (`float`, `float`)
        The two MJDs in TAI.
    """
    tai = datetime.tai
    return (tai.jd1, tai.jd2)


def simple_to_datetime(simple: Tuple[float, float], **kwargs: Any) -> astropy.time.Time:
    """Convert simple form back to astropy.time.Time"""
    return astropy.time.Time(*simple, format="jd", scale="tai")


def exptime_to_simple(exptime: astropy.units.Quantity) -> float:
    """Convert exposure time Quantity to seconds."""
    return exptime.to_value(astropy.units.s)


def simple_to_exptime(simple: float, **kwargs: Any) -> astropy.units.Quantity:
    """Convert simple form back to Quantity."""
    return simple * astropy.units.s


def angle_to_simple(angle: astropy.coordinates.Angle) -> float:
    """Convert Angle to degrees."""
    return angle.to_value(astropy.units.deg)


def simple_to_angle(simple: float, **kwargs: Any) -> astropy.coordinates.Angle:
    """Convert degrees to Angle."""
    return astropy.coordinates.Angle(simple * astropy.units.deg)


def focusz_to_simple(focusz: astropy.units.Quantity) -> float:
    """Convert focusz to meters."""
    return focusz.to_value(astropy.units.m)


def simple_to_focusz(simple: float, **kwargs: Any) -> astropy.units.Quantity:
    """Convert simple form back to Quantity."""
    return simple * astropy.units.m


def temperature_to_simple(temp: astropy.units.Quantity) -> float:
    """Convert temperature to kelvin."""
    return temp.to(astropy.units.K, equivalencies=astropy.units.temperature()).to_value()


def simple_to_temperature(simple: float, **kwargs: Any) -> astropy.units.Quantity:
    """Convert scalar kelvin value back to quantity."""
    return simple * astropy.units.K


def pressure_to_simple(press: astropy.units.Quantity) -> float:
    """Convert pressure Quantity to hPa."""
    return press.to_value(astropy.units.hPa)


def simple_to_pressure(simple: float, **kwargs: Any) -> astropy.units.Quantity:
    """Convert the pressure scalar back to Quantity."""
    return simple * astropy.units.hPa


def skycoord_to_simple(skycoord: astropy.coordinates.SkyCoord) -> Tuple[float, float]:
    """Convert SkyCoord to ICRS RA/Dec tuple"""
    icrs = skycoord.icrs
    return (icrs.ra.to_value(astropy.units.deg), icrs.dec.to_value(astropy.units.deg))


def simple_to_skycoord(simple: Tuple[float, float], **kwargs: Any) -> astropy.coordinates.SkyCoord:
    """Convert ICRS tuple to SkyCoord."""
    return astropy.coordinates.SkyCoord(*simple, unit=astropy.units.deg)


def altaz_to_simple(altaz: astropy.coordinates.AltAz) -> Tuple[float, float]:
    """Convert AltAz to Alt/Az tuple.

    Do not include obstime or location in simplification. It is assumed
    that those will be present from other properties.
    """
    return (altaz.az.to_value(astropy.units.deg), altaz.alt.to_value(astropy.units.deg))


def simple_to_altaz(simple: Tuple[float, float], **kwargs: Any) -> astropy.coordinates.AltAz:
    """Convert simple altaz tuple to AltAz.

    Will look for location and datetime_begin in kwargs.
    """
    location = kwargs.get("location")
    obstime = kwargs.get("datetime_begin")

    return astropy.coordinates.AltAz(
        simple[0] * astropy.units.deg, simple[1] * astropy.units.deg, obstime=obstime, location=location
    )


@dataclass
class PropertyDefinition:
    """Definition of an instrumental property."""

    doc: str
    """Docstring for the property."""

    str_type: str
    """Python type of property as a string (suitable for docstrings)."""

    py_type: type
    """Actual python type."""

    to_simple: Optional[Callable[[Any], Any]] = None
    """Function to convert value to simple form (can be ``None``)."""

    from_simple: Optional[Callable[[Any], Any]] = None
    """Function to convert from simple form back to required type (can be
    ``None``)."""


# Dict of properties to tuple where tuple is:
# - description of property
# - Python type of property as a string (suitable for docstrings)
# - Actual python type as a type
# - Simplification function (can be None)
# - Function to convert simple form back to required type (can be None)
PROPERTIES = {
    "telescope": PropertyDefinition("Full name of the telescope.", "str", str),
    "instrument": PropertyDefinition("The instrument used to observe the exposure.", "str", str),
    "location": PropertyDefinition(
        "Location of the observatory.",
        "astropy.coordinates.EarthLocation",
        astropy.coordinates.EarthLocation,
        earthlocation_to_simple,
        simple_to_earthlocation,
    ),
    "exposure_id": PropertyDefinition(
        "Unique (with instrument) integer identifier for this observation.", "int", int
    ),
    "visit_id": PropertyDefinition(
        """ID of the Visit this Exposure is associated with.

Science observations should essentially always be
associated with a visit, but calibration observations
may not be.""",
        "int",
        int,
    ),
    "physical_filter": PropertyDefinition("The bandpass filter used for this observation.", "str", str),
    "datetime_begin": PropertyDefinition(
        "Time of the start of the observation.",
        "astropy.time.Time",
        astropy.time.Time,
        datetime_to_simple,
        simple_to_datetime,
    ),
    "datetime_end": PropertyDefinition(
        "Time of the end of the observation.",
        "astropy.time.Time",
        astropy.time.Time,
        datetime_to_simple,
        simple_to_datetime,
    ),
    "exposure_time": PropertyDefinition(
        "Duration of the exposure with shutter open (seconds).",
        "astropy.units.Quantity",
        astropy.units.Quantity,
        exptime_to_simple,
        simple_to_exptime,
    ),
    "dark_time": PropertyDefinition(
        "Duration of the exposure with shutter closed (seconds).",
        "astropy.units.Quantity",
        astropy.units.Quantity,
        exptime_to_simple,
        simple_to_exptime,
    ),
    "boresight_airmass": PropertyDefinition("Airmass of the boresight of the telescope.", "float", float),
    "boresight_rotation_angle": PropertyDefinition(
        "Angle of the instrument in boresight_rotation_coord frame.",
        "astropy.coordinates.Angle",
        astropy.coordinates.Angle,
        angle_to_simple,
        simple_to_angle,
    ),
    "boresight_rotation_coord": PropertyDefinition(
        "Coordinate frame of the instrument rotation angle (options: sky, unknown).",
        "str",
        str,
    ),
    "detector_num": PropertyDefinition(
        "Unique (for instrument) integer identifier for the sensor.", "int", int
    ),
    "detector_name": PropertyDefinition(
        "Name of the detector within the instrument (might not be unique if there are detector groups).",
        "str",
        str,
    ),
    "detector_unique_name": PropertyDefinition(
        (
            "Unique name of the detector within the focal plane, generally combining detector_group with "
            "detector_name."
        ),
        "str",
        str,
    ),
    "detector_serial": PropertyDefinition("Serial number/string associated with this detector.", "str", str),
    "detector_group": PropertyDefinition(
        "Collection name of which this detector is a part. Can be None if there are no detector groupings.",
        "str",
        str,
    ),
    "detector_exposure_id": PropertyDefinition(
        "Unique integer identifier for this detector in this exposure.",
        "int",
        int,
    ),
    "focus_z": PropertyDefinition(
        "Defocal distance.",
        "astropy.units.Quantity",
        astropy.units.Quantity,
        focusz_to_simple,
        simple_to_focusz,
    ),
    "object": PropertyDefinition("Object of interest or field name.", "str", str),
    "temperature": PropertyDefinition(
        "Temperature outside the dome.",
        "astropy.units.Quantity",
        astropy.units.Quantity,
        temperature_to_simple,
        simple_to_temperature,
    ),
    "pressure": PropertyDefinition(
        "Atmospheric pressure outside the dome.",
        "astropy.units.Quantity",
        astropy.units.Quantity,
        pressure_to_simple,
        simple_to_pressure,
    ),
    "relative_humidity": PropertyDefinition("Relative humidity outside the dome.", "float", float),
    "tracking_radec": PropertyDefinition(
        "Requested RA/Dec to track.",
        "astropy.coordinates.SkyCoord",
        astropy.coordinates.SkyCoord,
        skycoord_to_simple,
        simple_to_skycoord,
    ),
    "altaz_begin": PropertyDefinition(
        "Telescope boresight azimuth and elevation at start of observation.",
        "astropy.coordinates.AltAz",
        astropy.coordinates.AltAz,
        altaz_to_simple,
        simple_to_altaz,
    ),
    "science_program": PropertyDefinition("Observing program (survey or proposal) identifier.", "str", str),
    "observation_type": PropertyDefinition(
        "Type of observation (currently: science, dark, flat, bias, focus).",
        "str",
        str,
    ),
    "observation_id": PropertyDefinition(
        "Label uniquely identifying this observation (can be related to 'exposure_id').",
        "str",
        str,
    ),
    "observation_reason": PropertyDefinition(
        "Reason this observation was taken, or its purpose ('science' and 'calibration' are common values)",
        "str",
        str,
    ),
    "exposure_group": PropertyDefinition(
        "Label to use to associate this exposure with others (can be related to 'exposure_id').",
        "str",
        str,
    ),
    "observing_day": PropertyDefinition(
        "Integer in YYYYMMDD format corresponding to the day of observation.", "int", int
    ),
    "observation_counter": PropertyDefinition(
        (
            "Counter of this observation. Can be counter within observing_day or a global counter. "
            "Likely to be observatory specific."
        ),
        "int",
        int,
    ),
    "has_simulated_content": PropertyDefinition(
        "Boolean indicating whether any part of this observation was simulated.", "bool", bool, None, None
    ),
    "group_counter_start": PropertyDefinition(
        "Observation counter for the start of the exposure group."
        "Depending on the instrument the relevant group may be "
        "visit_id or exposure_group.",
        "int",
        int,
        None,
        None,
    ),
    "group_counter_end": PropertyDefinition(
        "Observation counter for the end of the exposure group. "
        "Depending on the instrument the relevant group may be "
        "visit_id or exposure_group.",
        "int",
        int,
        None,
        None,
    ),
}
