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

"""Metadata translation code for CFHT MegaPrime FITS headers"""

__all__ = ("MegaPrimeTranslator", )

from astropy.coordinates import EarthLocation, AltAz, Angle
import astropy.units as u

from .fits import FitsTranslator
from .helpers import tracking_from_degree_headers

filters = {'u.MP9301': 'u',
           'u.MP9302': 'u2',
           'g.MP9401': 'g',
           'g.MP9402': 'g2',
           'r.MP9601': 'r',
           'r.MP9602': 'r2',
           'i.MP9701': 'i',
           'i.MP9702': 'i2',
           'i.MP9703': 'i3',
           'z.MP9801': 'z',
           'z.MP9901': 'z2',
           }


class MegaPrimeTranslator(FitsTranslator):
    """Metadata translator for CFHT MegaPrime standard headers.
    """

    name = "MegaPrime"
    """Name of this translation class"""

    supported_instrument = "MegaPrime"
    """Supports the MegaPrime instrument."""

    _const_map = {"boresight_rotation_angle": Angle(float("nan")*u.deg),
                  "boresight_rotation_coord": "unknown"}

    _trivial_map = {"physical_filter": "FILTER",
                    "dark_time": ("DARKTIME", dict(unit=u.s)),
                    "exposure_time": ("EXPTIME", dict(unit=u.s)),
                    "observation_id": "OBSID",
                    "object": "OBJECT",
                    "science_program": "RUNID",
                    "exposure_id": "EXPNUM",
                    "visit_id": "EXPNUM",
                    "detector_name": "CCDNAME",
                    "relative_humidity": ["RELHUMID", "HUMIDITY"],
                    "temperature": (["TEMPERAT", "AIRTEMP"], dict(unit=u.deg_C)),
                    "boresight_airmass": ["AIRMASS", "BORE-AIRMASS"]}

    def to_datetime_begin(self):
        # Docstring will be inherited. Property defined in properties.py
        # We know it is UTC
        value = self._from_fits_date_string(self._header["DATE-OBS"],
                                            time_str=self._header["UTC-OBS"], scale="utc")
        self._used_these_cards("DATE-OBS", "UTC-OBS")
        return value

    def to_datetime_end(self):
        # Docstring will be inherited. Property defined in properties.py
        # Older files are missing UTCEND
        if "UTCEND" in self._header:
            # We know it is UTC
            value = self._from_fits_date_string(self._header["DATE-OBS"],
                                                time_str=self._header["UTCEND"], scale="utc")
            self._used_these_cards("DATE-OBS", "UTCEND")
        else:
            # Take a guess by adding on the exposure time
            value = self.to_datetime_begin() + self.to_exposure_time()
        return value

    def to_location(self):
        """Calculate the observatory location.

        Returns
        -------
        location : `astropy.coordinates.EarthLocation`
            An object representing the location of the telescope.
        """
        # Height is not in some MegaPrime files. Use the value from EarthLocation.of_site("CFHT")
        # Some data uses OBS-LONG, OBS-LAT, other data uses LONGITUD and LATITUDE
        for long_key, lat_key in (("LONGITUD", "LATITUDE"), ("OBS-LONG", "OBS-LAT")):
            if long_key in self._header and lat_key in self._header:
                value = EarthLocation.from_geodetic(self._header[long_key], self._header[lat_key], 4215.0)
                self._used_these_cards(long_key, lat_key)
                break
        else:
            value = EarthLocation.of_site("CFHT")
        return value

    def to_detector_num(self):
        # Docstring will be inherited. Property defined in properties.py
        try:
            extname = self._header["EXTNAME"]
            num = int(extname[3:])  # chop off "ccd"
            self._used_these_cards("EXTNAME")
            return num
        except (KeyError, ValueError):
            # Dummy value, intended for PHU (need something to get filename)
            return 99

    def to_observation_type(self):
        """Calculate the observation type.

        Returns
        -------
        typ : `str`
            Observation type. Normalized to standard set.
        """
        obstype = self._header["OBSTYPE"].strip().lower()
        self._used_these_cards("OBSTYPE")
        if obstype == "object":
            return "science"
        return obstype

    def to_tracking_radec(self):
        """Calculate the tracking RA/Dec for this observation.

        Currently will be `None` for geocentric apparent coordinates.
        Additionally, can be `None` for non-science observations.

        The method supports multiple versions of header defining tracking
        coordinates.

        Returns
        -------
        coords : `astropy.coordinates.SkyCoord`
            The tracking coordinates.
        """
        radecsys = ("RADECSYS", "OBJRADEC", "RADESYS")
        radecpairs = (("RA_DEG", "DEC_DEG"), ("BORE-RA", "BORE-DEC"))
        return tracking_from_degree_headers(self, radecsys, radecpairs)

    def to_altaz_begin(self):
        """Calculate the azimuth and elevation for the start of the
        observation.

        Can be `None` for non-science observations.
        Two different header schemes are supported.
        If headers indicate negative values it is assumed that this is
        an observation where the telescope position was not relevant.

        Returns
        -------
        altaz : `astropy.coordinates.AltAz`
            The telescope coordinates.
        """
        for az_key, alt_key in (("TELAZ", "TELALT"), ("BORE-AZ", "BORE-ALT")):
            if az_key in self._header and alt_key in self._header:
                az = self._header[az_key]
                alt = self._header[alt_key]
                if az < 1.0 or alt < 1.0:
                    # Calibrations have magic values of -9999 when telescope not
                    # involved in observation.
                    return None
                altaz = AltAz(az * u.deg, alt * u.deg,
                              obstime=self.to_datetime_begin(), location=self.to_location())
                self._used_these_cards(az_key, alt_key)
                return altaz
        if self.to_observation_type() == "science":
            raise KeyError("Unable to determine alt/az of science observation")
        return None

    def to_detector_exposure_id(self):
        # Docstring will be inherited. Property defined in properties.py
        return self.to_exposure_id() * 36 + self.to_detector_num()

    def to_pressure(self):
        # Docstring will be inherited. Property defined in properties.py
        # Can be either AIRPRESS in Pa or PRESSURE in mbar
        for key, unit in (("PRESSURE", u.hPa), ("AIRPRESS", u.Pa)):
            if key in self._header:
                return self.quantity_from_card(key, unit)
        else:
            raise KeyError("Could not find pressure keywords in header")
