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

"""Metadata translation code for HSC FITS headers"""

__all__ = ("HscTranslator", )

from .fits import FitsTranslator


class HscTranslator(FitsTranslator):
    """Metadata translator for HSC standard headers.
    """

    name = "HSC"
    """Name of this translation class"""

    supportedInstrument = "HSC"
    """Supports the HSC instrument."""

    _constMap = {"instrument": "HSC"}
    """This translator only works for HSC."""

    _trivialMap = {"physical_filter": "FILTER01"}
    """One-to-one mappings"""

    @classmethod
    def canTranslate(cls, header):
        """Indicate whether this translation class can translate the
        supplied header.

        There is no ``INSTRUME`` header in an HSC file, so this method
        looks for HSC mentions in other headers.

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
        for k in ("EXP-ID", "FRAMEID"):
            if k in header:
                if header[k].startswith("HSC"):
                    return True
        return False

    def to_datetime_begin(self):
        # We know it is UTC
        return self._from_fits_date_string(self._header["DATE-OBS"],
                                           timeStr=self._header["UT"], scale="utc")

    def to_datetime_end(self):
        # We know it is UTC
        return self._from_fits_date_string(self._header["DATE-OBS"],
                                           timeStr=self._header["UT-END"], scale="utc")

    def to_abstract_filter(self):
        physical = self.to_physical_filter()
        if physical.startswith("HSC-"):
            return physical[4:]
        return None
