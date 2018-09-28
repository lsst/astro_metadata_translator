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

"""Classes and support code for metadata translation"""

__all__ = ("MetadataTranslator",)

from abc import abstractmethod
import logging

import astropy.units as u

log = logging.getLogger(__name__)


class MetadataMeta(type):
    """Register all subclasses with the base class"""

    @staticmethod
    def _makeConstMapping(standardKey, constant):
        def constant_translator(self):
            return constant
        constant_translator.__doc__ = f"""Returns constant value for '{standardKey}' property

        Returns
        -------
        translation : `{type(constant).__name__}`
            Always returns {constant!r}.
        """
        return constant_translator

    @staticmethod
    def _makeTrivialMapping(standardKey, fitsKey):
        def trivial_translator(self):
            value = self._header[fitsKey]
            self._used_these_cards(fitsKey)
            return value
        trivial_translator.__doc__ = f"""Map '{fitsKey}' FITS keyword to '{standardKey}' property

        Returns
        -------
        translation : `str` or `numbers.Number`
            Translated header value derived directly from the {fitsKey}
            header value.
        """
        return trivial_translator

    def __init__(cls, name, bases, dct):
        super().__init__(name, bases, dct)

        if hasattr(cls, "name") and cls.name is not None:
            MetadataTranslator.translators[cls.name] = cls

        for standardKey, fitsKey in cls._trivialMap.items():
            translator = cls._makeTrivialMapping(standardKey, fitsKey)
            setattr(cls, f"to_{standardKey}", translator)

        for standardKey, constant in cls._constMap.items():
            translator = cls._makeConstMapping(standardKey, constant)
            setattr(cls, f"to_{standardKey}", translator)


class MetadataTranslator(metaclass=MetadataMeta):
    """Per-instrument metadata translation support

    Parameters
    ----------
    header : `dict`-like
        Representation of an instrument FITS header that can be manipulated
        as if it was a `dict`.
    """

    _trivialMap = {}
    """Dict of one-to-one mappings for header translation from standard
    property to corresponding FITS keyword."""

    _constMap = {}
    """Dict defining a constant for specified standard properties."""

    translators = dict()
    """All registered metadata translation classes."""

    supportedInstrument = None
    """Name of instrument understood by this translation class."""

    def __init__(self, header):
        self._header = header
        self._used_cards = set()

    @classmethod
    @abstractmethod
    def canTranslate(cls, header):
        """Indicate whether this translation class can translate the
        supplied header.

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
        raise NotImplementedError()

    @classmethod
    def determineTranslator(cls, header):
        """Determine a translation class by examining the header

        Parameters
        ----------
        header : `dict`-like
            Representation of a FITS header.

        Returns
        -------
        translator : `MetadataTranslator`
            Translation class that knows how to extract metadata from
            the supplied header.
        """
        for name, trans in cls.translators.items():
            if trans.canTranslate(header):
                log.debug(f"Using translation class {name}")
                return trans
        else:
            raise ValueError("None of the registered translation classes understood this header")

    def _used_these_cards(self, *args):
        """Indicate that the supplied cards have been used for translation.

        Parameters
        ----------
        args : sequence of `str`
            Keywords used to process a translation.
        """
        self._used_cards.update(set(args))

    def cards_used(self):
        """Cards used during metadata extraction.

        Returns
        -------
        used : `frozenset`
            Cards used when extracting metadata.
        """
        return frozenset(self._used_cards)

    def quantity_from_card(self, keyword, unit):
        """Calculate a Astropy Quantity from a header card and a unit.

        Parameters
        ----------
        keyword : `str`
            Keyword to use from header.
        unit : `astropy.units.UnitBase`
            Unit of the item in the header.

        Returns
        -------
        q : `astropy.units.Quantity`
            Quantity representing the header value.
        """
        value = self._header[keyword]
        self._used_these_cards(keyword)
        return u.Quantity(value, unit=unit)
