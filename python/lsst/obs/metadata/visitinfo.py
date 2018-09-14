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

"""Extract standard metadata from instrument FITS headers"""

__all__ = ("VisitInfo", )

import logging

from .translator import MetadataTranslator

log = logging.getLogger(__name__)


class VisitInfo:
    """Standardized representation of an instrument FITS header.

    Parameters
    ----------
    header : `dict`-like
        Representation of an instrument FITS header accessible as a `dict`.
    translator : `MetadataTranslator`-class, `optional`
        If not `None`, the class to use to translate the supplied headers
        into standard form. Otherwise each registered translator class will
        be asked in turn if it knows how to translate the supplied header.

    Raises
    ------
    ValueError
        The supplied header was not recognized by any of the registered
        translators.
    """

    def __init__(self, header, translator=None):
        if translator is None:
            translator = MetadataTranslator.determineTranslator(header)

        # Loop over each translation (not final form -- this should be
        # defined in one place and consistent with translation classes)
        translations = ("telescope", "instrument", "datetime_begin", "datetime_end")
        for t in translations:
            # prototype code
            method = f"to_{t}"
            property = f"_{t}"

            try:
                print(f"Assigning {property} via {method}")
                setattr(self, property, getattr(translator, method)(header))
            except (AttributeError, KeyError):
                # For now assign None
                setattr(self, property, None)

    @property
    def datetime_begin(self):
        return self._datetime_begin

    @property
    def datetime_end(self):
        return self._datetime_end

    @property
    def instrument(self):
        return self._instrument

    @property
    def telescope(self):
        return self._telescope
