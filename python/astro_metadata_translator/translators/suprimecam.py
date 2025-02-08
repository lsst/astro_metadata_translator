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

"""Metadata translation code for SuprimeCam FITS headers."""

from __future__ import annotations

__all__ = ("SuprimeCamTranslator",)

import logging
import posixpath
import re
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

import astropy.units as u
from astropy.coordinates import Angle, SkyCoord

from ..translator import CORRECTIONS_RESOURCE_ROOT, cache_translation
from .helpers import altaz_from_degree_headers
from .subaru import SubaruTranslator

if TYPE_CHECKING:
    import astropy.coordinates
    import astropy.time

log = logging.getLogger(__name__)


class SuprimeCamTranslator(SubaruTranslator):
    """Metadata translator for HSC standard headers."""

    name = "SuprimeCam"
    """Name of this translation class"""

    supported_instrument = "SuprimeCam"
    """Supports the SuprimeCam instrument."""

    default_resource_root = posixpath.join(CORRECTIONS_RESOURCE_ROOT, "SuprimeCam")
    """Default resource path root to use to locate header correction files."""

    _const_map = {"boresight_rotation_coord": "unknown", "detector_group": None}
    """Constant mappings"""

    _trivial_map: dict[str, str | list[str] | tuple[Any, ...]] = {
        "observation_id": "EXP-ID",
        "object": "OBJECT",
        "science_program": "PROP-ID",
        "detector_num": "DET-ID",
        "detector_serial": "DETECTOR",  # DETECTOR is the "call name"
        "boresight_airmass": "AIRMASS",
        "relative_humidity": "OUT-HUM",
        "temperature": ("OUT-TMP", {"unit": u.K}),
        "pressure": ("OUT-PRS", {"unit": u.hPa}),
        "exposure_time": ("EXPTIME", {"unit": u.s}),
        "dark_time": ("EXPTIME", {"unit": u.s}),  # Assume same as exposure time
    }
    """One-to-one mappings"""

    # Zero point for SuprimeCam dates: 2004-01-01
    _DAY0 = 53005

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
        if "INSTRUME" in header:
            return header["INSTRUME"] == "SuprimeCam"

        for k in ("EXP-ID", "FRAMEID"):
            if cls.is_keyword_defined(header, k):
                if header[k].startswith("SUP"):
                    return True
        return False

    def _get_adjusted_mjd(self) -> int:
        """Calculate the modified julian date offset from reference day.

        Returns
        -------
        offset : `int`
            Offset day count from reference day.
        """
        mjd = self._header["MJD"]
        self._used_these_cards("MJD")
        return int(mjd) - self._DAY0

    @cache_translation
    def to_physical_filter(self) -> str:
        # Docstring will be inherited. Property defined in properties.py
        value = self._header["FILTER01"].strip().upper()
        self._used_these_cards("FILTER01")
        # Map potential "unknown" values to standard form
        if value in {"UNRECOGNIZED", "UNRECOGNISED", "NOTSET", "UNKNOWN"}:
            value = "unknown"
        elif value == "NONE":
            value = "empty"
        return value

    @cache_translation
    def to_datetime_begin(self) -> astropy.time.Time:
        # Docstring will be inherited. Property defined in properties.py
        # We know it is UTC
        value = self._from_fits_date_string(
            self._header["DATE-OBS"], time_str=self._header["UT-STR"], scale="utc"
        )
        self._used_these_cards("DATE-OBS", "UT-STR")
        return value

    @cache_translation
    def to_datetime_end(self) -> astropy.time.Time:
        # Docstring will be inherited. Property defined in properties.py
        # We know it is UTC
        value = self._from_fits_date_string(
            self._header["DATE-OBS"], time_str=self._header["UT-END"], scale="utc"
        )
        self._used_these_cards("DATE-OBS", "UT-END")

        # Sometimes the end time is less than the begin time plus the
        # exposure time so we have to check for that.
        exposure_time = self.to_exposure_time()
        datetime_begin = self.to_datetime_begin()
        exposure_end = datetime_begin + exposure_time
        if value < exposure_end:
            value = exposure_end

        return value

    @cache_translation
    def to_exposure_id(self) -> int:
        """Calculate unique exposure integer for this observation.

        Returns
        -------
        visit : `int`
            Integer uniquely identifying this exposure.
        """
        exp_id = self._header["EXP-ID"].strip()
        m = re.search(r"^SUP[A-Z](\d{7})0$", exp_id)
        if not m:
            raise RuntimeError(f"{self._log_prefix}: Unable to interpret EXP-ID: {exp_id}")
        exposure = int(m.group(1))
        if int(exposure) == 0:
            # Don't believe it
            frame_id = self._header["FRAMEID"].strip()
            m = re.search(r"^SUP[A-Z](\d{7})\d{1}$", frame_id)
            if not m:
                raise RuntimeError(f"{self._log_prefix}: Unable to interpret FRAMEID: {frame_id}")
            exposure = int(m.group(1))
        self._used_these_cards("EXP-ID", "FRAMEID")
        return exposure

    @cache_translation
    def to_visit_id(self) -> int:
        """Calculate the unique integer ID for this visit.

        Assumed to be identical to the exposure ID in this implementation.

        Returns
        -------
        exp : `int`
            Unique visit identifier.
        """
        return self.to_exposure_id()

    @cache_translation
    def to_observation_type(self) -> str:
        """Calculate the observation type.

        Returns
        -------
        typ : `str`
            Observation type. Normalized to standard set.
        """
        obstype = self._header["DATA-TYP"].strip().lower()
        self._used_these_cards("DATA-TYP")
        if obstype == "object":
            return "science"
        return obstype

    @cache_translation
    def to_tracking_radec(self) -> SkyCoord:
        # Docstring will be inherited. Property defined in properties.py
        radec = SkyCoord(
            self._header["RA2000"],
            self._header["DEC2000"],
            frame="icrs",
            unit=(u.hourangle, u.deg),
            obstime=self.to_datetime_begin(),
            location=self.to_location(),
        )
        self._used_these_cards("RA2000", "DEC2000")
        return radec

    @cache_translation
    def to_altaz_begin(self) -> astropy.coordinates.AltAz:
        # Docstring will be inherited. Property defined in properties.py
        return altaz_from_degree_headers(self, (("ALTITUDE", "AZIMUTH"),), self.to_datetime_begin())

    @cache_translation
    def to_boresight_rotation_angle(self) -> Angle:
        # Docstring will be inherited. Property defined in properties.py
        angle = Angle(self.quantity_from_card("INR-STR", u.deg))
        angle = angle.wrap_at("360d")
        return angle

    @cache_translation
    def to_detector_exposure_id(self) -> int:
        # Docstring will be inherited. Property defined in properties.py
        return self.to_exposure_id() * 10 + self.to_detector_num()

    @cache_translation
    def to_detector_name(self) -> str:
        # Docstring will be inherited. Property defined in properties.py
        # See https://subarutelescope.org/Observing/Instruments/SCam/ccd.html
        num = self.to_detector_num()

        names = (
            "nausicaa",
            "kiki",
            "fio",
            "sophie",
            "sheeta",
            "satsuki",
            "chihiro",
            "clarisse",
            "ponyo",
            "san",
        )

        return names[num]
