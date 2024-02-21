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

"""Metadata translation code for standard FITS headers."""

from __future__ import annotations

__all__ = ("FitsTranslator",)

from collections.abc import Mapping
from typing import Any

import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time

from ..translator import MetadataTranslator, cache_translation


class FitsTranslator(MetadataTranslator):
    """Metadata translator for FITS standard headers.

    Understands:

    - DATE-OBS/MJD-OBS or DATE-BEG/MJD-BEG (-BEG is preferred).
    - INSTRUME
    - TELESCOP
    - OBSGEO-[X,Y,Z]
    """

    # Direct translation from header key to standard form
    _trivial_map: dict[str, str | list[str] | tuple[Any, ...]] = dict(
        instrument="INSTRUME",
        telescope="TELESCOP",
    )

    @classmethod
    def can_translate(cls, header: Mapping[str, Any], filename: str | None = None) -> bool:
        """Indicate whether this translation class can translate the
        supplied header.

        Checks the instrument value and compares with the supported
        instruments in the class

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
        if cls.supported_instrument is None:
            return False

        # Protect against being able to always find a standard
        # header for instrument
        try:
            translator = cls(header, filename=filename)
            instrument = translator.to_instrument()
        except KeyError:
            return False

        return instrument == cls.supported_instrument

    @classmethod
    def _from_fits_date_string(cls, date_str: str, scale: str = "utc", time_str: str | None = None) -> Time:
        """Parse standard FITS ISO-style date string and return time object.

        Parameters
        ----------
        date_str : `str`
            FITS format date string to convert to standard form. Bypasses
            lookup in the header.
        scale : `str`, optional
            Override the time scale from the TIMESYS header. Defaults to
            UTC.
        time_str : `str`, optional
            If provided, overrides any time component in the ``dateStr``,
            retaining the YYYY-MM-DD component and appending this time
            string, assumed to be of format HH:MM::SS.ss.

        Returns
        -------
        date : `astropy.time.Time`
            `~astropy.time.Time` representation of the date.
        """
        if time_str is not None:
            date_str = f"{date_str[:10]}T{time_str}"

        return Time(date_str, format="isot", scale=scale)

    def _from_fits_date(
        self, date_key: str, mjd_key: str | None = None, scale: str | None = None
    ) -> Time | None:
        """Calculate a date object from the named FITS header.

        Uses the TIMESYS header if present to determine the time scale,
        defaulting to UTC.  Can be overridden since sometimes headers
        use TIMESYS for DATE- style headers but also have headers using
        different time scales.

        Parameters
        ----------
        date_key : `str`
            The key in the header representing a standard FITS
            ISO-style date. Can be `None` to go straight to MJD key.
        mjd_key : `str`, optional
            The key in the header representing a standard FITS MJD
            style date.  This key will be tried if ``date_key`` is not
            found, is `None`, or can not be parsed.
        scale : `str`, optional
            Override value to use for the time scale in preference to
            TIMESYS or the default. Should be a form understood by
            `~astropy.time.Time`.

        Returns
        -------
        date : `astropy.time.Time`
            `~astropy.time.Time` representation of the date.
        """
        used = []
        if scale is not None:
            pass
        elif self.is_key_ok("TIMESYS"):
            scale = self._header["TIMESYS"].lower()
            used.append("TIMESYS")
        else:
            scale = "utc"
        if date_key is not None and self.is_key_ok(date_key):
            date_str = self._header[date_key]
            value = self._from_fits_date_string(date_str, scale=scale)
            used.append(date_key)
        elif self.is_key_ok(mjd_key):
            assert mjd_key is not None  # for mypy (is_key_ok checks this)
            value = Time(self._header[mjd_key], scale=scale, format="mjd")
            used.append(mjd_key)
        else:
            value = None
        self._used_these_cards(*used)
        return value

    @cache_translation
    def to_datetime_begin(self) -> Time | None:
        """Calculate start time of observation.

        Uses FITS standard ``MJD-BEG`` or ``DATE-BEG``, in conjunction
        with the ``TIMESYS`` header. Will fallback to using ``MJD-OBS``
        or ``DATE-OBS`` if the ``-BEG`` variants are not found.

        Returns
        -------
        start_time : `astropy.time.Time` or `None`
            Time corresponding to the start of the observation. Returns
            `None` if no date can be found.
        """
        # Prefer -BEG over -OBS
        begin = None
        for suffix in ("BEG", "OBS"):
            begin = self._from_fits_date(f"DATE-{suffix}", mjd_key=f"MJD-{suffix}")
            if begin is not None:
                break
        return begin

    @cache_translation
    def to_datetime_end(self) -> Time:
        """Calculate end time of observation.

        Uses FITS standard ``MJD-END`` or ``DATE-END``, in conjunction
        with the ``TIMESYS`` header.

        Returns
        -------
        start_time : `astropy.time.Time`
            Time corresponding to the end of the observation.
        """
        return self._from_fits_date("DATE-END", mjd_key="MJD-END")

    @cache_translation
    def to_location(self) -> EarthLocation:
        """Calculate the observatory location.

        Uses FITS standard ``OBSGEO-`` headers.

        Returns
        -------
        location : `astropy.coordinates.EarthLocation`
            An object representing the location of the telescope.
        """
        cards = [f"OBSGEO-{c}" for c in ("X", "Y", "Z")]
        coords = [self._header[c] for c in cards]
        value = EarthLocation.from_geocentric(*coords, unit=u.m)
        self._used_these_cards(*cards)
        return value
