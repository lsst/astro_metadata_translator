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

"""Generically useful translation helpers which translation classes
can use.

They are written as free functions.  Some of them are written
as if they are methods of `MetadataTranslator`, allowing them to be attached
to translator classes that need them.  These methods have full access to
the translator methods.

Other functions are pure helpers that can be imported and used to help
translation classes without using `MetadataTranslator` properties.
"""

from __future__ import annotations

__all__ = (
    "altitude_from_zenith_distance",
    "is_non_science",
    "to_location_via_telescope_name",
    "tracking_from_degree_headers",
)

import logging
from collections.abc import Sequence
from typing import TYPE_CHECKING

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord

if TYPE_CHECKING:
    import astropy.units

    from ..translator import MetadataTranslator

log = logging.getLogger(__name__)


def to_location_via_telescope_name(self: MetadataTranslator) -> EarthLocation:
    """Calculate the observatory location via the telescope name.

    Returns
    -------
    loc : `astropy.coordinates.EarthLocation`
        Location of the observatory.
    """
    return EarthLocation.of_site(self.to_telescope())


def is_non_science(self: MetadataTranslator) -> None:
    """Raise an exception if this is a science observation.

    Raises
    ------
    KeyError
        Is a science observation.
    """
    if self.to_observation_type() == "science":
        raise KeyError(f"{self._log_prefix}: Header represents science observation and can not default")


def altitude_from_zenith_distance(zd: astropy.units.Quantity) -> astropy.units.Quantity:
    """Convert zenith distance to altitude.

    Parameters
    ----------
    zd : `astropy.units.Quantity`
        Zenith distance as an angle.

    Returns
    -------
    alt : `astropy.units.Quantity`
        Altitude.
    """
    return 90.0 * u.deg - zd


def tracking_from_degree_headers(
    self: MetadataTranslator,
    radecsys: Sequence[str],
    radecpairs: tuple[tuple[str, str], ...],
    unit: astropy.units.Unit = u.deg,
) -> SkyCoord:
    """Calculate the tracking coordinates from lists of headers.

    Parameters
    ----------
    radecsys : `list` or `tuple`
        Header keywords to try corresponding to the tracking system.  If none
        match ICRS will be assumed.
    radecpairs : `tuple` of `tuple` of pairs of `str`
        Pairs of keywords specifying the RA/Dec in units of ``unit``.
    unit : `astropy.unit.BaseUnit` or `tuple`
        Unit definition suitable for the `~astropy.coordinate.SkyCoord`
        constructor.

    Returns
    -------
    radec = `astropy.coordinates.SkyCoord`
        The RA/Dec coordinates. None if this is a moving target or a
        non-science observation without any RA/Dec definition.

    Raises
    ------
    KeyError
        No RA/Dec keywords were found and this observation is a science
        observation.
    """
    used = []
    for k in radecsys:
        if self.is_key_ok(k):
            frame = self._header[k].strip().lower()
            used.append(k)
            if frame == "gappt":
                self._used_these_cards(*used)
                # Moving target
                return None
            break
    else:
        frame = "icrs"
    for ra_key, dec_key in radecpairs:
        if self.are_keys_ok([ra_key, dec_key]):
            radec = SkyCoord(
                self._header[ra_key],
                self._header[dec_key],
                frame=frame,
                unit=unit,
                obstime=self.to_datetime_begin(),
                location=self.to_location(),
            )
            self._used_these_cards(ra_key, dec_key, *used)
            return radec
    if self.to_observation_type() == "science":
        raise KeyError(f"{self._log_prefix}: Unable to determine tracking RA/Dec of science observation")
    return None


def altaz_from_degree_headers(
    self: MetadataTranslator,
    altazpairs: tuple[tuple[str, str], ...],
    obstime: astropy.time.Time,
    is_zd: set[str] | None = None,
) -> AltAz:
    """Calculate the altitude/azimuth coordinates from lists of headers.

    If the altitude is found but is greater than 90 deg, it will be returned
    fixed at 90 deg.
    If the altitude or azimuth are negative and this is a calibration
    observation, `None` will be returned.

    Parameters
    ----------
    altazpairs : `tuple` of `str`
        Pairs of keywords specifying Alt/Az in degrees. Each pair is tried
        in turn.
    obstime : `astropy.time.Time`
        Reference time to use for these coordinates.
    is_zd : `set`, optional
        Contains keywords that correspond to zenith distances rather than
        altitude.

    Returns
    -------
    altaz = `astropy.coordinates.AltAz`
        The AltAz coordinates associated with the telescope location
        and provided time.  Returns `None` if this observation is not
        a science observation and no AltAz keys were located.

    Raises
    ------
    KeyError
        No AltAz keywords were found and this observation is a science
        observation.
    """
    for alt_key, az_key in altazpairs:
        if self.are_keys_ok([az_key, alt_key]):
            az = self._header[az_key]
            alt = self._header[alt_key]

            # Check for zenith distance
            if is_zd and alt_key in is_zd:
                alt = altitude_from_zenith_distance(alt * u.deg).value

            if az < -360.0 or alt < 0.0:
                # Break out of loop since we have found values but
                # they are bad.
                break
            if alt > 90.0:
                log.warning("%s: Clipping altitude (%f) at 90 degrees", self._log_prefix, alt)
                alt = 90.0

            altaz = AltAz(az * u.deg, alt * u.deg, obstime=obstime, location=self.to_location())
            self._used_these_cards(az_key, alt_key)
            return altaz
    if self.to_observation_type() == "science":
        raise KeyError(f"{self._log_prefix}: Unable to determine AltAz of science observation")
    return None
