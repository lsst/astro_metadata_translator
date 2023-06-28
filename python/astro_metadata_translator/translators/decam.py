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

"""Metadata translation code for DECam FITS headers."""

from __future__ import annotations

__all__ = ("DecamTranslator",)

import logging
import posixpath
import re
from collections.abc import Iterator, Mapping, MutableMapping
from typing import TYPE_CHECKING, Any

import astropy.units as u
from astropy.coordinates import Angle, EarthLocation
from astropy.io import fits

from ..translator import CORRECTIONS_RESOURCE_ROOT, cache_translation
from .fits import FitsTranslator
from .helpers import altaz_from_degree_headers, is_non_science, tracking_from_degree_headers

if TYPE_CHECKING:
    import astropy.coordinates
    import astropy.time

log = logging.getLogger(__name__)


class DecamTranslator(FitsTranslator):
    """Metadata translator for DECam standard headers."""

    name = "DECam"
    """Name of this translation class"""

    supported_instrument = "DECam"
    """Supports the DECam instrument."""

    default_resource_root = posixpath.join(CORRECTIONS_RESOURCE_ROOT, "DECam")
    """Default resource path root to use to locate header correction files."""

    # DECam has no rotator, and the instrument angle on sky is set to +Y=East,
    # +X=South which we define as a 90 degree rotation and an X-flip.
    _const_map = {
        "boresight_rotation_angle": Angle(90 * u.deg),
        "boresight_rotation_coord": "sky",
    }

    _trivial_map: dict[str, str | list[str] | tuple[Any, ...]] = {
        "exposure_time": ("EXPTIME", dict(unit=u.s)),
        "dark_time": ("DARKTIME", dict(unit=u.s)),
        "boresight_airmass": ("AIRMASS", dict(checker=is_non_science)),
        "observation_id": "OBSID",
        "object": "OBJECT",
        "science_program": "PROPID",
        "detector_num": "CCDNUM",
        "detector_serial": "DETECTOR",
        "detector_unique_name": "DETPOS",
        "telescope": ("TELESCOP", dict(default="CTIO 4.0-m telescope")),
        "instrument": ("INSTRUME", dict(default="DECam")),
        # Ensure that reasonable values are always available
        "relative_humidity": ("HUMIDITY", dict(default=40.0, minimum=0, maximum=100.0)),
        "temperature": ("OUTTEMP", dict(unit=u.deg_C, default=10.0, minimum=-10.0, maximum=40.0)),
        # Header says torr but seems to be mbar. Use hPa unit
        # which is the SI equivalent of mbar.
        "pressure": ("PRESSURE", dict(unit=u.hPa, default=771.611, minimum=700.0, maximum=850.0)),
    }

    # Unique detector names are currently not used but are read directly from
    # header.
    # The detector_group could be N or S with detector_name corresponding
    # to the number in that group.
    detector_names = {
        1: "S29",
        2: "S30",
        3: "S31",
        4: "S25",
        5: "S26",
        6: "S27",
        7: "S28",
        8: "S20",
        9: "S21",
        10: "S22",
        11: "S23",
        12: "S24",
        13: "S14",
        14: "S15",
        15: "S16",
        16: "S17",
        17: "S18",
        18: "S19",
        19: "S8",
        20: "S9",
        21: "S10",
        22: "S11",
        23: "S12",
        24: "S13",
        25: "S1",
        26: "S2",
        27: "S3",
        28: "S4",
        29: "S5",
        30: "S6",
        31: "S7",
        32: "N1",
        33: "N2",
        34: "N3",
        35: "N4",
        36: "N5",
        37: "N6",
        38: "N7",
        39: "N8",
        40: "N9",
        41: "N10",
        42: "N11",
        43: "N12",
        44: "N13",
        45: "N14",
        46: "N15",
        47: "N16",
        48: "N17",
        49: "N18",
        50: "N19",
        51: "N20",
        52: "N21",
        53: "N22",
        54: "N23",
        55: "N24",
        56: "N25",
        57: "N26",
        58: "N27",
        59: "N28",
        60: "N29",
        62: "N31",
    }

    @classmethod
    def can_translate(cls, header: Mapping[str, Any], filename: str | None = None) -> bool:
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
    def to_exposure_id(self) -> int:
        """Calculate exposure ID.

        Returns
        -------
        id : `int`
            ID of exposure.
        """
        value = self._header["EXPNUM"]
        self._used_these_cards("EXPNUM")
        return value

    @cache_translation
    def to_observation_counter(self) -> int:
        """Return the lifetime exposure number.

        Returns
        -------
        sequence : `int`
            The observation counter.
        """
        return self.to_exposure_id()

    @cache_translation
    def to_visit_id(self) -> int:
        # Docstring will be inherited. Property defined in properties.py
        return self.to_exposure_id()

    @cache_translation
    def to_datetime_end(self) -> astropy.time.Time:
        # Docstring will be inherited. Property defined in properties.py
        # Instcals have no DATE-END or DTUTC
        datetime_end = self._from_fits_date("DTUTC", scale="utc")
        if datetime_end is None:
            datetime_end = self.to_datetime_begin() + self.to_exposure_time()
        return datetime_end

    def _translate_from_calib_id(self, field: str) -> str:
        """Fetch the ID from the CALIB_ID header.

        Calibration products made with constructCalibs have some metadata
        saved in its FITS header CALIB_ID.
        """
        data = self._header["CALIB_ID"]
        match = re.search(r".*%s=(\S+)" % field, data)
        if not match:
            raise RuntimeError(f"Header CALIB_ID with value '{data}' has not field '{field}'")
        self._used_these_cards("CALIB_ID")
        return match.groups()[0]

    @cache_translation
    def to_physical_filter(self) -> str | None:
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
    def to_location(self) -> astropy.coordinates.EarthLocation:
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
    def to_observation_type(self) -> str:
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
    def to_tracking_radec(self) -> astropy.coordinates.SkyCoord:
        # Docstring will be inherited. Property defined in properties.py
        radecsys = ("RADESYS",)
        radecpairs = (("TELRA", "TELDEC"),)
        return tracking_from_degree_headers(self, radecsys, radecpairs, unit=(u.hourangle, u.deg))

    @cache_translation
    def to_altaz_begin(self) -> astropy.coordinates.AltAz:
        # Docstring will be inherited. Property defined in properties.py
        return altaz_from_degree_headers(self, (("ZD", "AZ"),), self.to_datetime_begin(), is_zd={"ZD"})

    @cache_translation
    def to_detector_exposure_id(self) -> int | None:
        # Docstring will be inherited. Property defined in properties.py
        exposure_id = self.to_exposure_id()
        if exposure_id is None:
            return None
        return int(f"{exposure_id:07d}{self.to_detector_num():02d}")

    @cache_translation
    def to_detector_group(self) -> str:
        # Docstring will be inherited. Property defined in properties.py
        name = self.to_detector_unique_name()
        return name[0]

    @cache_translation
    def to_detector_name(self) -> str:
        # Docstring will be inherited. Property defined in properties.py
        name = self.to_detector_unique_name()
        return name[1:]

    @cache_translation
    def to_focus_z(self) -> u.Quantity:
        # Docstring will be inherited. Property defined in properties.py
        # ``TELFOCUS`` is a comma-separated string with six focus offsets
        # (fx, fy, fz, tx, ty, tz) recorded in units of microns.
        tel_focus_list = self._header["TELFOCUS"].split(",")
        return float(tel_focus_list[2]) * u.um

    @classmethod
    def fix_header(
        cls, header: MutableMapping[str, Any], instrument: str, obsid: str, filename: str | None = None
    ) -> bool:
        """Fix DECam headers.

        Parameters
        ----------
        header : `dict`
            The header to update.  Updates are in place.
        instrument : `str`
            The name of the instrument.
        obsid : `str`
            Unique observation identifier associated with this header.
            Will always be provided.
        filename : `str`, optional
            Filename associated with this header. May not be set since headers
            can be fixed independently of any filename being known.

        Returns
        -------
        modified = `bool`
            Returns `True` if the header was updated.

        Notes
        -----
        Fixes the following issues:

        * If OBSTYPE contains "zero" or "bias",
          update the FILTER keyword to "solid plate 0.0 0.0".

        Corrections are reported as debug level log messages.
        """
        modified = False

        # Calculate the standard label to use for log messages
        log_label = cls._construct_log_prefix(obsid, filename)

        obstype = header.get("OBSTYPE", "unknown")

        if "bias" in obstype.lower() or "zero" in obstype.lower():
            header["FILTER"] = "solid plate 0.0 0.0"
            modified = True
            log.debug("%s: Set FILTER to %s because OBSTYPE is %s", log_label, header["FILTER"], obstype)

        return modified

    @classmethod
    def determine_translatable_headers(
        cls, filename: str, primary: MutableMapping[str, Any] | None = None
    ) -> Iterator[MutableMapping[str, Any]]:
        """Given a file return all the headers usable for metadata translation.

        DECam files are multi-extension FITS with a primary header and
        each detector stored in a subsequent extension.  DECam uses
        ``INHERIT=T`` and each detector header will be merged with the
        primary header.

        Guide headers are not returned.

        Parameters
        ----------
        filename : `str`
            Path to a file in a format understood by this translator.
        primary : `dict`-like, optional
            The primary header obtained by the caller. This is sometimes
            already known, for example if a system is trying to bootstrap
            without already knowing what data is in the file. Will be
            merged with detector headers if supplied, else will be read
            from the file.

        Yields
        ------
        headers : iterator of `dict`-like
            Each detector header in turn. The supplied header will be merged
            with the contents of each detector header.

        Notes
        -----
        This translator class is specifically tailored to raw DECam data and
        is not designed to work with general FITS files. The normal paradigm
        is for the caller to have read the first header and then called
        `determine_translator()` on the result to work out which translator
        class to then call to obtain the real headers to be used for
        translation.
        """
        # Circular dependency so must defer import.
        from ..headers import merge_headers

        # This is convoluted because we need to turn an Optional variable
        # to a Dict so that mypy is happy.
        primary_hdr = primary if primary else {}

        # Since we want to scan many HDUs we use astropy directly to keep
        # the file open rather than continually opening and closing it
        # as we go to each HDU.
        with fits.open(filename) as fits_file:
            # Astropy does not automatically handle the INHERIT=T in
            # DECam headers so the primary header must be merged.
            first_pass = True

            for hdu in fits_file:
                if first_pass:
                    if not primary_hdr:
                        primary_hdr = hdu.header
                    first_pass = False
                    continue

                header = hdu.header
                if "CCDNUM" not in header:  # Primary does not have CCDNUM
                    continue
                if header["CCDNUM"] > 62:  # ignore guide CCDs
                    continue
                yield merge_headers([primary_hdr, header], mode="overwrite")
