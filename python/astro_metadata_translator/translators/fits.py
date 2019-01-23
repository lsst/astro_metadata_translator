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

"""Metadata translation code for standard FITS headers"""

__all__ = ("FitsTranslator", )

from astropy.time import Time
from astropy.coordinates import EarthLocation
import astropy.units as u

from ..translator import MetadataTranslator, cache_translation


class FitsTranslator(MetadataTranslator):
    """Metadata translator for FITS standard headers.

    Understands:

    - DATE-OBS
    - INSTRUME
    - TELESCOP
    - OBSGEO-[X,Y,Z]

    """

    # Direct translation from header key to standard form
    _trivial_map = dict(instrument="INSTRUME",
                        telescope="TELESCOP")

    @classmethod
    def can_translate(cls, header, filename=None):
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
    def _from_fits_date_string(cls, date_str, scale='utc', time_str=None):
        """Parse standard FITS ISO-style date string and return time object

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
            date_str = "{}T{}".format(date_str[:10], time_str)

        return Time(date_str, format="isot", scale=scale)

    def _from_fits_date(self, date_key):
        """Calculate a date object from the named FITS header

        Uses the TIMESYS header if present to determine the time scale,
        defaulting to UTC.

        Parameters
        ----------
        dateKey : `str`
            The key in the header representing a standard FITS
            ISO-style date.

        Returns
        -------
        date : `astropy.time.Time`
            `~astropy.time.Time` representation of the date.
        """
        used = [date_key, ]
        if "TIMESYS" in self._header:
            scale = self._header["TIMESYS"].lower()
            used.append("TIMESYS")
        else:
            scale = "utc"
        if date_key in self._header:
            date_str = self._header[date_key]
            value = self._from_fits_date_string(date_str, scale=scale)
            self._used_these_cards(*used)
        else:
            value = None
        return value

    @cache_translation
    def to_datetime_begin(self):
        """Calculate start time of observation.

        Uses FITS standard ``DATE-OBS`` and ``TIMESYS`` headers.

        Returns
        -------
        start_time : `astropy.time.Time`
            Time corresponding to the start of the observation.
        """
        return self._from_fits_date("DATE-OBS")

    @cache_translation
    def to_datetime_end(self):
        """Calculate end time of observation.

        Uses FITS standard ``DATE-END`` and ``TIMESYS`` headers.

        Returns
        -------
        start_time : `astropy.time.Time`
            Time corresponding to the end of the observation.
        """
        return self._from_fits_date("DATE-END")

    @cache_translation
    def to_location(self):
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
