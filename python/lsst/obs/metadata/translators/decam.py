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

"""Metadata translation code for DECam FITS headers"""

__all__ = ("DecamTranslator", )

import re

from .fits import FitsTranslator


class DecamTranslator(FitsTranslator):
    """Metadata translator for DECam standard headers.
    """

    name = "DECam"
    """Name of this translation class"""

    supportedInstrument = "DECam"
    """Supports the DECam instrument."""

    _trivialMap = {"exposure_time": "EXPTIME",
                   "dark_time": "DARKTIME",
                   "obsid": "OBSID",
                   "visit": "EXPNUM"}

    @classmethod
    def to_datetime_end(cls, header):
        return FitsTranslator._from_fits_date(header, "DTUTC")

    @classmethod
    def _translateFromCalibId(cls, field, header):
        """Fetch the ID from the CALIB_ID header.

        Calibration products made with constructCalibs have some metadata
        saved in its FITS header CALIB_ID.
        """
        data = header["CALIB_ID"]
        match = re.search(".*%s=(\S+)" % field, data)
        return match.groups()[0]

    @classmethod
    def to_abstract_filter(cls, header):
        """Calculate the abstract filter.

        Parameters
        ----------
        header : `dict`-like
            Representation of the FITS header.

        Returns
        -------
        filter : `str`
            The abstract filter name.

        """
        # The abstract filter can be derived from the first word in the
        # physical filter description
        physical = cls.to_physical_filter(header)
        if physical:
            return physical.split()[0]

    @classmethod
    def to_physical_filter(cls, header):
        """Calculate physical filter.

        Return `None` if the keyword FILTER does not exist in the header,
        which can happen for some valid Community Pipeline products.

        Parameters
        ----------
        header : `dict`-like
            Representation of the FITS header.

        Returns
        -------
        filter : `str`
            The full filter name.
        """
        if "FILTER" in header:
            if "OBSTYPE" in header and "zero" in header["OBSTYPE"].strip().lower():
                return "NONE"
            return header["FILTER"].strip()
        elif "CALIB_ID" in header:
            return cls._translateFromCalibId("filter", header)
        else:
            return None
