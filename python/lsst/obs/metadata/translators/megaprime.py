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

from astropy.coordinates import EarthLocation

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

    _trivialMap = {"physical_filter": "FILTER",
                   "dark_time": "DARKTIME",
                   "exposure_time": "EXPTIME",
                   "obsid": "OBSID",
                   "exposure": "EXPNUM",
                   "visit": "EXPNUM",
                   "detector_name": "CCDNAME",
                   "boresight_airmass": "AIRMASS"}

    def to_abstract_filter(self):
        physical = self.to_physical_filter()
        if physical in filters:
            return filters[physical][0]
        return None

    def to_datetime_begin(self):
        # We know it is UTC
        value = self._from_fits_date_string(self._header["DATE-OBS"],
                                            timeStr=self._header["UTC-OBS"], scale="utc")
        self._used_these_cards("DATE-OBS", "UTC-OBS")
        return value

    def to_datetime_end(self):
        # We know it is UTC
        value = self._from_fits_date_string(self._header["DATE-OBS"],
                                            timeStr=self._header["UTCEND"], scale="utc")
        self._used_these_cards("DATE-OBS", "UTCEND")
        return value

    def to_location(self):
        """Calculate the observatory location.

        Returns
        -------
        location : `astropy.coordinates.EarthLocation`
            An object representing the location of the telescope.
        """
        # Height is not in MegaPrime files. Use the value from EarthLocation.of_site("CFHT")
        value = EarthLocation.from_geodetic(self._header["LONGITUD"], self._header["LATITUDE"], 4215.0)
        self._used_these_cards("LONGITUD", "LATITUDE")
        return value

    def to_detector_num(self):
        try:
            extname = self._header["EXTNAME"]
            num = int(extname[3:])  # chop off "ccd"
            self._used_these_cards("EXTNAME")
            return num
        except KeyError:
            # Dummy value, intended for PHU (need something to get filename)
            return 99
