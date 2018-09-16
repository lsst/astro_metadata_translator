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

"""Metadata translation code for CFHT MegaPrime FITS headers"""

__all__ = ("MegaPrimeTranslator", )

from .fits import FitsTranslator

filters = {'u.MP9301': 'u',
           'u.MP9302': 'u2',
           'g.MP9401': 'g',
           'g.MP9402': 'g2',
           'r.MP9601': 'r',
           'r.MP9602': 'r2',
           'i.MP9701': 'i',
           'i.MP9702': 'i2',
           'i.MP9703': 'i3',
           'z.MP9801': 'z',
           'z.MP9901': 'z2',
           }


class MegaPrimeTranslator(FitsTranslator):
    """Metadata translator for CFHT MegaPrime standard headers.
    """

    name = "MegaPrime"
    """Name of this translation class"""

    supportedInstrument = "MegaPrime"
    """Supports the MegaPrime instrument."""

    @classmethod
    def to_physical_filter(cls, header):
        return header["FILTER"]

    @classmethod
    def to_abstract_filter(cls, header):
        physical = cls.to_physical_filter(header)
        if physical in filters:
            return filters[physical][0]
        return None

    @classmethod
    def to_datetime_begin(cls, header):
        # Form the FITS string and use the standard parser
        fitsStr = header["DATE-OBS"] + "T" + header["UTC-OBS"]
        # We know it is UTC
        return cls._from_fits_date_string(fitsStr, scale="utc")

    @classmethod
    def to_datetime_end(cls, header):
        # Form the FITS string and use the standard parser
        fitsStr = header["DATE-OBS"] + "T" + header["UTCEND"]
        # We know it is UTC
        return cls._from_fits_date_string(fitsStr, scale="utc")
