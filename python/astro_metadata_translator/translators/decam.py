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

"""Metadata translation code for DECam FITS headers"""

__all__ = ("DecamTranslator", )

import re

from astropy.coordinates import EarthLocation, Angle
import astropy.units as u

from ..translator import cache_translation
from .fits import FitsTranslator
from .helpers import altaz_from_degree_headers, is_non_science, \
    tracking_from_degree_headers


class DecamTranslator(FitsTranslator):
    """Metadata translator for DECam standard headers.
    """

    name = "DECam"
    """Name of this translation class"""

    supported_instrument = "DECam"
    """Supports the DECam instrument."""

    _const_map = {"boresight_rotation_angle": Angle(float("nan")*u.deg),
                  "boresight_rotation_coord": "unknown",
                  "detector_group": None}

    _trivial_map = {"exposure_time": ("EXPTIME", dict(unit=u.s)),
                    "dark_time": ("DARKTIME", dict(unit=u.s)),
                    "boresight_airmass": ("AIRMASS", dict(checker=is_non_science)),
                    "observation_id": "OBSID",
                    "object": "OBJECT",
                    "science_program": "PROPID",
                    "detector_num": "CCDNUM",
                    "detector_serial": "DETECTOR",
                    "detector_name": "DETPOS",
                    "telescope": ("TELESCOP", dict(default="CTIO 4.0-m telescope")),
                    "instrument": ("INSTRUME", dict(default="DECam")),
                    # Ensure that reasonable values are always available
                    "relative_humidity": ("HUMIDITY", dict(default=40., minimum=0, maximum=100.)),
                    "temperature": ("OUTTEMP", dict(unit=u.deg_C, default=10., minimum=-10., maximum=40.)),
                    # Header says torr but seems to be mbar. Use hPa unit
                    # which is the SI equivalent of mbar.
                    "pressure": ("PRESSURE", dict(unit=u.hPa,
                                 default=771.611, minimum=700., maximum=850.)),
                    }

    @classmethod
    def can_translate(cls, header, filename=None):
        """Indicate whether this translation class can translate the
        supplied header.

        Checks the INSTRUME and FILTER headers.

        Parameters
        ----------
        header : `dict`-like
            Header to convert to standardized form.
        filename : `str`, optional
            Name of file being translated.

        Returns
        -------
        can : `bool`
            `True` if the header is recognized by this class. `False`
            otherwise.
        """
        # Use INSTRUME. Because of defaulting behavior only do this
        # if we really have an INSTRUME header
        if "INSTRUME" in header:
            via_instrume = super().can_translate(header, filename=filename)
            if via_instrume:
                return via_instrume
        if cls.is_keyword_defined(header, "FILTER") and "DECam" in header["FILTER"]:
            return True
        return False

    @cache_translation
    def to_exposure_id(self):
        """Calculate exposure ID solely for science observations.

        Returns
        -------
        id : `int`
            ID of exposure.
        """
        if self.to_observation_type() != "science":
            return None
        value = self._header["EXPNUM"]
        self._used_these_cards("EXPNUM")
        return value

    @cache_translation
    def to_visit_id(self):
        # Docstring will be inherited. Property defined in properties.py
        return self.to_exposure_id()

    @cache_translation
    def to_datetime_end(self):
        # Docstring will be inherited. Property defined in properties.py
        return self._from_fits_date("DTUTC", scale="utc")

    def _translate_from_calib_id(self, field):
        """Fetch the ID from the CALIB_ID header.

        Calibration products made with constructCalibs have some metadata
        saved in its FITS header CALIB_ID.
        """
        data = self._header["CALIB_ID"]
        match = re.search(r".*%s=(\S+)" % field, data)
        self._used_these_cards("CALIB_ID")
        return match.groups()[0]

    @cache_translation
    def to_physical_filter(self):
        """Calculate physical filter.

        Return `None` if the keyword FILTER does not exist in the header,
        which can happen for some valid Community Pipeline products.

        Returns
        -------
        filter : `str`
            The full filter name.
        """
        if self.is_key_ok("FILTER"):
            value = self._header["FILTER"].strip()
            self._used_these_cards("FILTER")
            return value
        elif self.is_key_ok("CALIB_ID"):
            return self._translate_from_calib_id("filter")
        else:
            return None

    @cache_translation
    def to_location(self):
        """Calculate the observatory location.

        Returns
        -------
        location : `astropy.coordinates.EarthLocation`
            An object representing the location of the telescope.
        """

        if self.is_key_ok("OBS-LONG"):
            # OBS-LONG has west-positive sign so must be flipped
            lon = self._header["OBS-LONG"] * -1.0
            value = EarthLocation.from_geodetic(lon, self._header["OBS-LAT"], self._header["OBS-ELEV"])
            self._used_these_cards("OBS-LONG", "OBS-LAT", "OBS-ELEV")
        else:
            # Look up the value since some files do not have location
            value = EarthLocation.of_site("ctio")

        return value

    @cache_translation
    def to_observation_type(self):
        """Calculate the observation type.

        Returns
        -------
        typ : `str`
            Observation type. Normalized to standard set.
        """
        if not self.is_key_ok("OBSTYPE"):
            return "none"
        obstype = self._header["OBSTYPE"].strip().lower()
        self._used_these_cards("OBSTYPE")
        if obstype == "object":
            return "science"
        return obstype

    @cache_translation
    def to_tracking_radec(self):
        # Docstring will be inherited. Property defined in properties.py
        radecsys = ("RADESYS",)
        radecpairs = (("TELRA", "TELDEC"),)
        return tracking_from_degree_headers(self, radecsys, radecpairs, unit=(u.hourangle, u.deg))

    @cache_translation
    def to_altaz_begin(self):
        # Docstring will be inherited. Property defined in properties.py
        return altaz_from_degree_headers(self, (("ZD", "AZ"),),
                                         self.to_datetime_begin(), is_zd=set(["ZD"]))

    @cache_translation
    def to_detector_exposure_id(self):
        # Docstring will be inherited. Property defined in properties.py
        exposure_id = self.to_exposure_id()
        if exposure_id is None:
            return None
        return int("{:07d}{:02d}".format(exposure_id, self.to_detector_num()))
