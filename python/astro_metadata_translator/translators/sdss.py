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

"""Metadata translation code for SDSS FITS headers."""

from __future__ import annotations

__all__ = ("SdssTranslator",)

import posixpath
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

import astropy.units as u
from astropy.coordinates import AltAz, Angle, EarthLocation, UnknownSiteException

from ..translator import CORRECTIONS_RESOURCE_ROOT, cache_translation
from .fits import FitsTranslator
from .helpers import tracking_from_degree_headers

if TYPE_CHECKING:
    import astropy.coordinates
    import astropy.time


# Hard-code APO location fallback.
_APO_LOCATION = EarthLocation.from_geocentric(
    -1463969.30185172, -5166673.34223433, 3434985.71204565, unit=u.m
)


class SdssTranslator(FitsTranslator):
    """Metadata translator for SDSS standard headers.
    NB: calibration data is not handled as calibration frames were
    not available to me at time of writing.
    """

    name = "SDSS"
    """Name of this translation class"""

    supported_instrument = "Imager"
    """Supports the SDSS imager instrument."""

    default_resource_root = posixpath.join(CORRECTIONS_RESOURCE_ROOT, "SDSS")
    """Default resource path root to use to locate header correction files."""

    # SDSS has has a rotator, but in drift scan mode, the instrument
    # angle on sky is set to +X=East, +Y=North which we define as a
    # 0 degree rotation.
    _const_map = {
        "boresight_rotation_angle": Angle(0 * u.deg),
        "boresight_rotation_coord": "sky",
        "dark_time": 0.0 * u.s,  # Drift scan implies no dark time
        "instrument": "Imager on SDSS 2.5m",  # We only ever ingest data from the imager
        "telescope": "SDSS 2.5m",  # Value of TELESCOP in header is ambiguous
        "relative_humidity": None,
        "temperature": None,
        "pressure": None,
        "detector_serial": "UNKNOWN",
    }

    _trivial_map = {
        "exposure_time": ("EXPTIME", dict(unit=u.s)),
        "object": "OBJECT",
        "physical_filter": "FILTER",
        "exposure_id": "RUN",
        "visit_id": "RUN",
        "science_program": "OBJECT",  # This is the closest I can think of to a useful program
        "detector_name": "CCDLOC",  # This is a numeric incoding of the "slot", i.e. filter+camcol
    }

    #  Need a mapping from unique name to index.  The order is arbitrary.
    detector_name_id_map = {
        "g1": 0,
        "z1": 1,
        "u1": 2,
        "i1": 3,
        "r1": 4,
        "g2": 5,
        "z2": 6,
        "u2": 7,
        "i2": 8,
        "r2": 9,
        "g3": 10,
        "z3": 11,
        "u3": 12,
        "i3": 13,
        "r3": 14,
        "g4": 15,
        "z4": 16,
        "u4": 17,
        "i4": 18,
        "r4": 19,
        "g5": 20,
        "z5": 21,
        "u5": 22,
        "i5": 23,
        "r5": 24,
        "g6": 25,
        "z6": 26,
        "u6": 27,
        "i6": 28,
        "r6": 29,
    }

    @classmethod
    def can_translate(cls, header: Mapping[str, Any], filename: str | None = None) -> bool:
        """Indicate whether this translation class can translate the
        supplied header.

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
        if (
            cls.is_keyword_defined(header, "ORIGIN")
            and cls.is_keyword_defined(header, "CCDMODE")
            and cls.is_keyword_defined(header, "TELESCOP")
            and "2.5m" in header["TELESCOP"]
            and "SDSS" in header["ORIGIN"]
            and "DRIFT" in header["CCDMODE"]
        ):
            return True
        return False

    @cache_translation
    def to_detector_unique_name(self) -> str:
        # Docstring will be inherited. Property defined in properties.py
        if self.is_key_ok("CAMCOL"):
            return self.to_physical_filter() + str(self._header["CAMCOL"])
        else:
            raise ValueError(f"{self._log_prefix}: CAMCOL key is not definded")

    @cache_translation
    def to_detector_num(self) -> int:
        # Docstring will be inherited. Property defined in properties.py
        return self.detector_name_id_map[self.to_detector_unique_name()]

    @cache_translation
    def to_observation_id(self) -> str:
        """Calculate the observation ID.

        Returns
        -------
        observation_id : `str`
            A string uniquely describing the observation.
            This incorporates the run, camcol, filter and frame.
        """
        return " ".join([str(self._header[el]) for el in ["RUN", "CAMCOL", "FILTER", "FRAME"]])

    @cache_translation
    def to_datetime_begin(self) -> astropy.time.Time:
        # Docstring will be inherited. Property defined in properties.py
        # We know it is UTC
        value = self._from_fits_date_string(
            self._header["DATE-OBS"], time_str=self._header["TAIHMS"], scale="tai"
        )
        self._used_these_cards("DATE-OBS", "TAIHMS")
        return value

    @cache_translation
    def to_datetime_end(self) -> astropy.time.Time:
        # Docstring will be inherited. Property defined in properties.py
        return self.to_datetime_begin() + self.to_exposure_time()

    @cache_translation
    def to_location(self) -> EarthLocation:
        """Calculate the observatory location.

        Returns
        -------
        location : `astropy.coordinates.EarthLocation`
            An object representing the location of the telescope.
        """
        # Look up the value since files do not have location.
        # This might require a network look up if the cached database
        # is not accessible. If it fails fall back to hard-coded location.
        try:
            value = EarthLocation.of_site("apo")
        except UnknownSiteException:
            value = _APO_LOCATION

        return value

    @cache_translation
    def to_observation_type(self) -> str:
        """Calculate the observation type.

        Returns
        -------
        typ : `str`
            Observation type. Normalized to standard set.
        """
        obstype_key = "FLAVOR"
        if not self.is_key_ok(obstype_key):
            return "none"
        obstype = self._header[obstype_key].strip().lower()
        self._used_these_cards(obstype_key)
        return obstype

    @cache_translation
    def to_tracking_radec(self) -> astropy.coordinates.SkyCoord:
        # Docstring will be inherited. Property defined in properties.py
        radecsys = ("RADECSYS",)
        radecpairs = (("RA", "DEC"),)
        return tracking_from_degree_headers(self, radecsys, radecpairs, unit=u.deg)

    @cache_translation
    def to_altaz_begin(self) -> AltAz:
        # Docstring will be inherited. Property defined in properties.py
        try:
            az = self._header["AZ"]
            alt = self._header["ALT"]
            # It appears SDSS defines azimuth as increasing
            # from South through East. This translates to
            # North through East
            az = (-az + 180.0) % 360.0
            altaz = AltAz(
                az * u.deg, alt * u.deg, obstime=self.to_datetime_begin(), location=self.to_location()
            )
            self._used_these_cards("AZ", "ALT")
            return altaz
        except Exception as e:
            if self.to_observation_type() != "science":
                return None  # Allow Alt/Az not to be set for calibrations
            raise (e)

    @cache_translation
    def to_boresight_airmass(self) -> float | None:
        # Docstring will be inherited. Property defined in properties.py
        altaz = self.to_altaz_begin()
        if altaz is not None:
            return altaz.secz.value  # This is an estimate
        return None

    @cache_translation
    def to_detector_exposure_id(self) -> int | None:
        # Docstring will be inherited. Property defined in properties.py
        try:
            frame_field_map = dict(r=0, i=2, u=4, z=6, g=8)
            run = self._header["RUN"]
            filt = self._header["FILTER"]
            camcol = self._header["CAMCOL"]
            field = self._header["FRAME"] - frame_field_map[filt]
            self._used_these_cards("RUN", "FILTER", "CAMCOL", "FRAME")
        except Exception as e:
            if self.to_observation_type() != "science":
                return None
            raise (e)
        filter_id_map = dict(u=0, g=1, r=2, i=3, z=4)
        return ((int(run) * 10 + filter_id_map[filt]) * 10 + int(camcol)) * 10000 + int(field)

    @cache_translation
    def to_detector_group(self) -> str:
        # Docstring will be inherited. Property defined in properties.py
        if self.is_key_ok("CAMCOL"):
            return str(self._header["CAMCOL"])
        else:
            raise ValueError(f"{self._log_prefix}: CAMCOL key is not definded")
