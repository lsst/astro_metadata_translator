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

"""Metadata translation code for CFHT MegaPrime FITS headers"""

__all__ = ("MegaPrimeTranslator", )

import re
import posixpath
from astropy.coordinates import EarthLocation, Angle
import astropy.units as u

from ..translator import cache_translation, CORRECTIONS_RESOURCE_ROOT
from .fits import FitsTranslator
from .helpers import tracking_from_degree_headers, altaz_from_degree_headers


class MegaPrimeTranslator(FitsTranslator):
    """Metadata translator for CFHT MegaPrime standard headers.
    """

    name = "MegaPrime"
    """Name of this translation class"""

    supported_instrument = "MegaPrime"
    """Supports the MegaPrime instrument."""

    default_resource_root = posixpath.join(CORRECTIONS_RESOURCE_ROOT, "CFHT")
    """Default resource path root to use to locate header correction files."""

    # CFHT Megacam has no rotator, and the instrument angle on sky is set to
    # +Y=N, +X=W which we define as a 0 degree rotation.
    _const_map = {"boresight_rotation_angle": Angle(0*u.deg),
                  "boresight_rotation_coord": "sky",
                  "detector_group": None}

    _trivial_map = {"physical_filter": "FILTER",
                    "dark_time": ("DARKTIME", dict(unit=u.s)),
                    "exposure_time": ("EXPTIME", dict(unit=u.s)),
                    "observation_id": "OBSID",
                    "object": "OBJECT",
                    "science_program": "RUNID",
                    "exposure_id": "EXPNUM",
                    "visit_id": "EXPNUM",
                    "detector_serial": "CCDNAME",
                    "relative_humidity": ["RELHUMID", "HUMIDITY"],
                    "temperature": (["TEMPERAT", "AIRTEMP"], dict(unit=u.deg_C)),
                    "boresight_airmass": ["AIRMASS", "BORE-AIRMASS"]}

    @cache_translation
    def to_datetime_begin(self):
        # Docstring will be inherited. Property defined in properties.py
        # We know it is UTC
        value = self._from_fits_date_string(self._header["DATE-OBS"],
                                            time_str=self._header["UTC-OBS"], scale="utc")
        self._used_these_cards("DATE-OBS", "UTC-OBS")
        return value

    @cache_translation
    def to_datetime_end(self):
        # Docstring will be inherited. Property defined in properties.py
        # Older files are missing UTCEND
        if self.is_key_ok("UTCEND"):
            # We know it is UTC
            value = self._from_fits_date_string(self._header["DATE-OBS"],
                                                time_str=self._header["UTCEND"], scale="utc")
            self._used_these_cards("DATE-OBS", "UTCEND")
        else:
            # Take a guess by adding on the exposure time
            value = self.to_datetime_begin() + self.to_exposure_time()
        return value

    @cache_translation
    def to_location(self):
        """Calculate the observatory location.

        Returns
        -------
        location : `astropy.coordinates.EarthLocation`
            An object representing the location of the telescope.
        """
        # Height is not in some MegaPrime files. Use the value from
        # EarthLocation.of_site("CFHT")
        # Some data uses OBS-LONG, OBS-LAT, other data uses LONGITUD and
        # LATITUDE
        for long_key, lat_key in (("LONGITUD", "LATITUDE"), ("OBS-LONG", "OBS-LAT")):
            if self.are_keys_ok([long_key, lat_key]):
                value = EarthLocation.from_geodetic(self._header[long_key], self._header[lat_key], 4215.0)
                self._used_these_cards(long_key, lat_key)
                break
        else:
            value = EarthLocation.of_site("CFHT")
        return value

    @cache_translation
    def to_detector_name(self):
        # Docstring will be inherited. Property defined in properties.py
        if self.is_key_ok("EXTNAME"):
            name = self._header["EXTNAME"]
            # Only valid name has form "ccdNN"
            if re.match(r"ccd\d+$", name):
                self._used_these_cards("EXTNAME")
                return name

        # Dummy value, intended for PHU (need something to get filename)
        return "ccd99"

    @cache_translation
    def to_detector_num(self):
        name = self.to_detector_name()
        return int(name[3:])

    @cache_translation
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

    @cache_translation
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

    @cache_translation
    def to_altaz_begin(self):
        # Docstring will be inherited. Property defined in properties.py
        return altaz_from_degree_headers(self, (("TELALT", "TELAZ"), ("BORE-ALT", "BORE-AZ")),
                                         self.to_datetime_begin())

    @cache_translation
    def to_detector_exposure_id(self):
        # Docstring will be inherited. Property defined in properties.py
        return self.to_exposure_id() * 36 + self.to_detector_num()

    @cache_translation
    def to_pressure(self):
        # Docstring will be inherited. Property defined in properties.py
        # Can be either AIRPRESS in Pa or PRESSURE in mbar
        for key, unit in (("PRESSURE", u.hPa), ("AIRPRESS", u.Pa)):
            if self.is_key_ok(key):
                return self.quantity_from_card(key, unit)
        else:
            raise KeyError("Could not find pressure keywords in header")

    @cache_translation
    def to_observation_counter(self):
        """Return the lifetime exposure number.

        Returns
        -------
        sequence : `int`
            The observation counter.
        """
        return self.to_exposure_id()
