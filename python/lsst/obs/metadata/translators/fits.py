# This file is part of obs_metadata.
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

"""Metadata translation code for standard FITS headers"""
from astropy.time import Time

from ..metadata import MetadataTranslator


class FitsTranslator(MetadataTranslator):
    """Metadata translator for FITS standard headers.

    Understands:

    - DATE-OBS
    - INSTRUME
    - TELESCOP

    """

    # Direct translation from header key to standard form
    _trivialMap = dict(instrument="INSTRUME",
                       telescope="TELESCOP")

    @classmethod
    def canTranslate(cls, header):
        """Indicate whether this translation class can translate the
        supplied header.

        Checks the instrument value and compares with the supported
        instruments in the class

        Parameters
        ----------
        header : `dict`-like
           Header to convert to standardized form.

        Returns
        -------
        can : `bool`
            `True` if the header is recognized by this class. `False`
            otherwise.
        """
        if cls.supportedInstrument is None:
            return False

        return cls.to_instrument(header) == cls.supportedInstrument

    @staticmethod
    def _from_fits_date(header, dateKey):
        """Calculate a date object from the named FITS header

        Uses the TIMESYS header if present to determine the time scale.

        Parameters
        ----------
        header : `dict`-like
            A `dict` representing a FITS-like header.
        dateKey : `str`
            The key in the header representing a standard FITS
            ISO-style date.

        Returns
        -------
        date : `astropy.time.Time`
            `~astropy.time.Time` representation of the date.
        """
        if "TIMESYS" in header:
            scale = header["TIMESYS"]
        else:
            scale = "utc"
        return Time(header[dateKey], format="isot", scale=scale)

    @staticmethod
    def to_datetime_begin(header):
        """Calculate start time of observation.

        Uses FITS standard ``DATE-OBS`` and ``TIMESYS`` headers.

        Parameters
        ----------
        header : `dict`-like
            Header values representing a FITS header.

        Returns
        -------
        start_time : `astropy.time.Time`
            Time corresponding to the start of the observation.
        """
        return FitsTranslator._from_fits_date(header, "DATE-OBS")

    @staticmethod
    def to_datetime_end(header):
        """Calculate end time of observation.

        Uses FITS standard ``DATE-END`` and ``TIMESYS`` headers.

        Parameters
        ----------
        header : `dict`-like
            Header values representing a FITS header.

        Returns
        -------
        start_time : `astropy.time.Time`
            Time corresponding to the end of the observation.
        """
        return FitsTranslator._from_fits_date(header, "DATE-END")
