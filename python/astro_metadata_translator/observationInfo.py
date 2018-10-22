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

"""Represent standard metadata from instrument headers"""

__all__ = ("ObservationInfo", )

import itertools
import logging
import copy

from .translator import MetadataTranslator
from .properties import PROPERTIES

log = logging.getLogger(__name__)


class ObservationInfo:
    """Standardized representation of an instrument header for a single
    exposure observation.

    Parameters
    ----------
    header : `dict`-like
        Representation of an instrument header accessible as a `dict`.
    translator_class : `MetadataTranslator`-class, `optional`
        If not `None`, the class to use to translate the supplied headers
        into standard form. Otherwise each registered translator class will
        be asked in turn if it knows how to translate the supplied header.

    Raises
    ------
    ValueError
        The supplied header was not recognized by any of the registered
        translators.
    TypeError
        The supplied translator class was not a MetadataTranslator.
    """

    _PROPERTIES = PROPERTIES
    """All the properties supported by this class with associated
    documentation."""

    def __init__(self, header, translator_class=None):

        # Store the supplied header for later stripping
        self._header = header

        # PropertyList is not dict-like so force to a dict here to simplify
        # the translation code.
        if hasattr(header, "toOrderedDict"):
            header = header.toOrderedDict()

        if translator_class is None:
            translator_class = MetadataTranslator.determine_translator(header)
        elif not issubclass(translator_class, MetadataTranslator):
            raise TypeError(f"Translator class must be a MetadataTranslator, not {translator_class}")

        # Create an instance for this header
        translator = translator_class(header)

        # Store the translator
        self._translator = translator

        # Loop over each translation (not final form -- this should be
        # defined in one place and consistent with translation classes)

        for t in self._PROPERTIES:
            # prototype code
            method = f"to_{t}"
            property = f"_{t}"

            try:
                setattr(self, property, getattr(translator, method)())
            except NotImplementedError as e:
                raise NotImplementedError(f"No translation exists for property '{t}'"
                                          f" using translator {translator.__class__}") from e
            except KeyError as e:
                raise KeyError(f"Error calculating property '{t}'"
                               f" using translator {translator.__class__}") from e

    @property
    def cards_used(self):
        """Header cards used for the translation.

        Returns
        -------
        used : `frozenset` of `str`
            Set of card used.
        """
        return self._translator.cards_used()

    def stripped_header(self):
        """Return a copy of the supplied header with used keywords removed.

        Returns
        -------
        stripped : `dict`-like
            Same class as header supplied to constructor, but with the
            headers used to calculate the generic information removed.
        """
        hdr = copy.copy(self._header)
        used = self._translator.cards_used()
        for c in used:
            del hdr[c]
        return hdr

    def __str__(self):
        # Put more interesting answers at front of list
        # and then do remainder
        priority = ("instrument", "telescope", "datetime_begin")
        properties = sorted(set(self._PROPERTIES.keys()) - set(priority))

        result = ""
        for p in itertools.chain(priority, properties):
            result += f"{p}: {getattr(self, p)}\n"

        return result

    def __eq__(self, other):
        """Compares equal if standard properties are equal
        """
        if type(self) != type(other):
            return False

        for p in self._PROPERTIES:
            # Use string comparison since SkyCoord.__eq__ seems unreliable
            # otherwise.  Should have per-type code so that floats and
            # quantities can be compared properly.
            v1 = f"{getattr(self, p)}"
            v2 = f"{getattr(other, p)}"
            if v1 != v2:
                return False

        return True

    def __getstate__(self):
        """Get pickleable state

        Returns the properties, the name of the translator, and the
        cards that were used.  Does not return the full header.

        Returns
        -------
        state : `dict`
            Dict containing items that can be persisted.
        """
        state = dict()
        for p in self._PROPERTIES:
            property = f"_{p}"
            state[p] = getattr(self, property)

        return state

    def __setstate__(self, state):
        for p in self._PROPERTIES:
            property = f"_{p}"
            setattr(self, property, state[p])


# Method to add the standard properties
def _make_property(property, doc, return_type):
    """Create a getter method with associated docstring.

    Parameters
    ----------
    property : `str`
        Name of the property getter to be created.
    doc : `str`
        Description of this property.
    return_type : `str`
        Type of this property (used in the doc string).

    Returns
    -------
    p : `function`
        Getter method for this property.
    """
    def getter(self):
        return getattr(self, f"_{property}")

    getter.__doc__ = f"""{doc}

    Returns
    -------
    {property} : `{return_type}`
        Access the property.
    """
    return getter


# Initialize the internal properties (underscored) and add the associated
# getter methods.
for name, description in ObservationInfo._PROPERTIES.items():
    setattr(ObservationInfo, f"_{name}", None)
    setattr(ObservationInfo, name, property(_make_property(name, *description)))
