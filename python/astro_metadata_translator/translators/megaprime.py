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

"""Metadata translation code for CFHT MegaPrime FITS headers."""

from __future__ import annotations

__all__ = ("MegaPrimeTranslator",)

import posixpath
import re
from collections.abc import Iterator, MutableMapping
from typing import TYPE_CHECKING, Any

import astropy.time
import astropy.units as u
from astropy.coordinates import Angle, EarthLocation, UnknownSiteException
from astropy.io import fits

from ..translator import CORRECTIONS_RESOURCE_ROOT, cache_translation
from .fits import FitsTranslator
from .helpers import altaz_from_degree_headers, tracking_from_degree_headers

if TYPE_CHECKING:
    import astropy.coordinates
    import astropy.units

_MK_FALLBACK_LOCATION = EarthLocation.from_geocentric(
    -5464301.77135369, -2493489.8730419, 2151085.16950589, unit=u.m
)


class MegaPrimeTranslator(FitsTranslator):
    """Metadata translator for CFHT MegaPrime standard headers."""

    name = "MegaPrime"
    """Name of this translation class"""

    supported_instrument = "MegaPrime"
    """Supports the MegaPrime instrument."""

    default_resource_root = posixpath.join(CORRECTIONS_RESOURCE_ROOT, "CFHT")
    """Default resource path root to use to locate header correction files."""

    _observing_day_offset = astropy.time.TimeDelta(0, format="sec", scale="tai")

    # CFHT Megacam has no rotator, and the instrument angle on sky is set to
    # +Y=N, +X=W which we define as a 0 degree rotation.
    _const_map = {
        "boresight_rotation_angle": Angle(0 * u.deg),
        "boresight_rotation_coord": "sky",
        "detector_group": None,
    }

    _trivial_map: dict[str, str | list[str] | tuple[Any, ...]] = {
        "physical_filter": "FILTER",
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
        "boresight_airmass": ["AIRMASS", "BORE-AIRMASS"],
    }

    @cache_translation
    def to_datetime_begin(self) -> astropy.time.Time:
        # Docstring will be inherited. Property defined in properties.py
        # We know it is UTC
        value = self._from_fits_date_string(
            self._header["DATE-OBS"], time_str=self._header["UTC-OBS"], scale="utc"
        )
        self._used_these_cards("DATE-OBS", "UTC-OBS")
        return value

    @cache_translation
    def to_datetime_end(self) -> astropy.time.Time:
        # Docstring will be inherited. Property defined in properties.py
        # Older files are missing UTCEND
        if self.is_key_ok("UTCEND"):
            # We know it is UTC
            value = self._from_fits_date_string(
                self._header["DATE-OBS"], time_str=self._header["UTCEND"], scale="utc"
            )
            self._used_these_cards("DATE-OBS", "UTCEND")
        else:
            # Take a guess by adding on the exposure time
            value = self.to_datetime_begin() + self.to_exposure_time()
        return value

    @cache_translation
    def to_location(self) -> EarthLocation:
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
            try:
                value = EarthLocation.of_site("CFHT")
            except UnknownSiteException:
                value = _MK_FALLBACK_LOCATION
        return value

    @cache_translation
    def to_detector_name(self) -> str:
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
    def to_detector_num(self) -> int:
        name = self.to_detector_name()
        return int(name[3:])

    @cache_translation
    def to_observation_type(self) -> str:
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
    def to_tracking_radec(self) -> astropy.coordinates.SkyCoord:
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
    def to_altaz_begin(self) -> astropy.coordinates.AltAz:
        # Docstring will be inherited. Property defined in properties.py
        return altaz_from_degree_headers(
            self, (("TELALT", "TELAZ"), ("BORE-ALT", "BORE-AZ")), self.to_datetime_begin()
        )

    @cache_translation
    def to_detector_exposure_id(self) -> int:
        # Docstring will be inherited. Property defined in properties.py
        return self.to_exposure_id() * 36 + self.to_detector_num()

    @cache_translation
    def to_pressure(self) -> astropy.units.Quantity:
        # Docstring will be inherited. Property defined in properties.py
        # Can be either AIRPRESS in Pa or PRESSURE in mbar
        for key, unit in (("PRESSURE", u.hPa), ("AIRPRESS", u.Pa)):
            if self.is_key_ok(key):
                return self.quantity_from_card(key, unit)

        raise KeyError(f"{self._log_prefix}: Could not find pressure keywords in header")

    @cache_translation
    def to_observation_counter(self) -> int:
        """Return the lifetime exposure number.

        Returns
        -------
        sequence : `int`
            The observation counter.
        """
        return self.to_exposure_id()

    @classmethod
    def determine_translatable_headers(
        cls, filename: str, primary: MutableMapping[str, Any] | None = None
    ) -> Iterator[MutableMapping[str, Any]]:
        """Given a file return all the headers usable for metadata translation.

        MegaPrime files are multi-extension FITS with a primary header and
        each detector stored in a subsequent extension.  MegaPrime uses
        ``INHERIT=F`` therefore the primary header will always be ignored
        if given.

        Parameters
        ----------
        filename : `str`
            Path to a file in a format understood by this translator.
        primary : `dict`-like, optional
            The primary header obtained by the caller. This is sometimes
            already known, for example if a system is trying to bootstrap
            without already knowing what data is in the file. Will be
            ignored.

        Yields
        ------
        headers : iterator of `dict`-like
            Each detector header in turn. The supplied header will never be
            included.

        Notes
        -----
        This translator class is specifically tailored to raw MegaPrime data
        and is not designed to work with general FITS files. The normal
        paradigm is for the caller to have read the first header and then
        called `determine_translator()` on the result to work out which
        translator class to then call to obtain the real headers to be used for
        translation.
        """
        # Since we want to scan many HDUs we use astropy directly to keep
        # the file open rather than continually opening and closing it
        # as we go to each HDU.
        with fits.open(filename) as fits_file:
            for hdu in fits_file:
                # Astropy <=4.2 strips the EXTNAME header but some CFHT data
                # have two EXTNAME headers and the CCD number is in the
                # second one.
                if hdu.name == "PRIMARY":
                    continue

                if hdu.name.startswith("ccd"):
                    # It may only be some data files that are broken so
                    # handle the expected form.
                    yield hdu.header
                    continue

                # Some test data at least has the EXTNAME as
                # COMPRESSED_IMAGE but the EXTVER as the detector number.
                if hdu.name == "COMPRESSED_IMAGE":
                    header = hdu.header

                    # Astropy strips EXTNAME so put it back for the translator
                    header["EXTNAME"] = f"ccd{hdu.ver:02d}"
                    yield header

    @classmethod
    def observing_date_to_offset(cls, observing_date: astropy.time.Time) -> astropy.time.TimeDelta | None:
        """Return the offset to use when calculating the observing day.

        Parameters
        ----------
        observing_date : `astropy.time.Time`
            The date of the observation. Unused.

        Returns
        -------
        offset : `astropy.time.TimeDelta`
            The offset to apply. The offset is always 0 seconds. In Hawaii
            UTC rollover is at 2pm local time.
        """
        return cls._observing_day_offset
