# This file is part of astro_metadata_translator.
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

"""Metadata translation code for Subaru telescope"""

__all__ = ("SubaruTranslator", )

from astropy.coordinates import EarthLocation

from ..translator import cache_translation
from .fits import FitsTranslator


class SubaruTranslator(FitsTranslator):
    """Metadata translator for Subaru telescope headers.
    """

    @cache_translation
    def to_location(self):
        """Returns the location of the Subaru telescope on Mauna Kea.

        Hardcodes the location and does not look at any headers.

        Returns
        -------
        location : `astropy.coordinates.EarthLocation`
            An object representing the location of the telescope.
        """
        return EarthLocation.from_geodetic(-155.476667, 19.825556, 4139.0)
