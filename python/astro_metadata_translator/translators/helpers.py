# This file is part of astro_metadata_translator.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Generically useful translation helpers which translation classes
can use.

They are written as free functions.  Some of them are written
as if they are methods of `MetadataTranslator`, allowing them to be attached
to translator classes that need them.  These methods have full access to
the translator methods.

Other functions are pure helpers that can be imported and used to help
translation classes without using `MetadataTranslator` properties.

"""

__all__ = ("to_location_via_telescope_name",
           "is_non_science",
           "tracking_from_degree_headers",
           "altitude_from_zenith_distance")

import logging
from astropy.coordinates import EarthLocation, SkyCoord, AltAz
import astropy.units as u

log = logging.getLogger(__name__)


def to_location_via_telescope_name(self):
    """Calculate the observatory location via the telescope name.

    Returns
    -------
    loc : `astropy.coordinates.EarthLocation`
        Location of the observatory.
    """
    return EarthLocation.of_site(self.to_telescope())


def is_non_science(self):
    """Raise an exception if this is a science observation.

    Raises
    ------
    KeyError
        Is a science observation.
    """
    if self.to_observation_type() == "science":
        raise KeyError("Header represents science observation and can not default")
    return


def altitude_from_zenith_distance(zd):
    """Convert zenith distance to altitude

    Parameters
    ----------
    zd : `astropy.units.Quantity`
        Zenith distance as an angle.

    Returns
    -------
    alt : `astropy.units.Quantity`
        Altitude.
    """
    return 90.*u.deg - zd


def tracking_from_degree_headers(self, radecsys, radecpairs, unit=u.deg):
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
        if k in self._header:
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
        if ra_key in self._header and dec_key in self._header:
            radec = SkyCoord(self._header[ra_key], self._header[dec_key],
                             frame=frame, unit=unit, obstime=self.to_datetime_begin(),
                             location=self.to_location())
            self._used_these_cards(ra_key, dec_key, *used)
            return radec
    if self.to_observation_type() == "science":
        raise KeyError("Unable to determine tracking RA/Dec of science observation")
    return None


def altaz_from_degree_headers(self, altazpairs, obstime, is_zd=None):
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
        if az_key in self._header and alt_key in self._header:
            az = self._header[az_key]
            alt = self._header[alt_key]

            # Check for zenith distance
            if is_zd and alt_key in is_zd:
                alt = altitude_from_zenith_distance(alt * u.deg).value

            if az < 0.0 or alt < 0.0:
                # Break out of loop since we have found values but
                # they are bad.
                break
            if alt > 90.0:
                log.warning("Clipping altitude (%f) at 90 degrees", alt)
                alt = 90.0

            altaz = AltAz(az * u.deg, alt * u.deg,
                          obstime=obstime, location=self.to_location())
            self._used_these_cards(az_key, alt_key)
            return altaz
    if self.to_observation_type() == "science":
        raise KeyError("Unable to determine AltAz of science observation")
    return None
