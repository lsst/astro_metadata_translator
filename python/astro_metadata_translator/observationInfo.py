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

"""Represent standard metadata from instrument headers"""

__all__ = ("ObservationInfo", "makeObservationInfo")

import itertools
import logging
import copy

import astropy.time

from .translator import MetadataTranslator
from .properties import PROPERTIES
from .headers import fix_header

log = logging.getLogger(__name__)


class ObservationInfo:
    """Standardized representation of an instrument header for a single
    exposure observation.

    Parameters
    ----------
    header : `dict`-like
        Representation of an instrument header accessible as a `dict`.
        May be updated with header corrections if corrections are found.
    filename : `str`, optional
        Name of the file whose header is being translated.  For some
        datasets with missing header information this can sometimes
        allow for some fixups in translations.
    translator_class : `MetadataTranslator`-class, optional
        If not `None`, the class to use to translate the supplied headers
        into standard form. Otherwise each registered translator class will
        be asked in turn if it knows how to translate the supplied header.
    pedantic : `bool`, optional
        If True the translation must succeed for all properties.  If False
        individual property translations must all be implemented but can fail
        and a warning will be issued.
    search_path : iterable, optional
        Override search paths to use during header fix up.

    Raises
    ------
    ValueError
        Raised if the supplied header was not recognized by any of the
        registered translators.
    TypeError
        Raised if the supplied translator class was not a MetadataTranslator.
    KeyError
        Raised if a translation fails and pedantic mode is enabled.
    NotImplementedError
        Raised if the selected translator does not support a required
        property.

    Notes
    -----
    Headers will be corrected if correction files are located and this will
    modify the header provided to the constructor.
    """

    _PROPERTIES = PROPERTIES
    """All the properties supported by this class with associated
    documentation."""

    def __init__(self, header, filename=None, translator_class=None, pedantic=False,
                 search_path=None):

        # Initialize the empty object
        self._header = {}
        self.filename = filename
        self._translator = None
        self.translator_class_name = "<None>"

        # To allow makeObservationInfo to work, we special case a None
        # header
        if header is None:
            return

        # Fix up the header (if required)
        fix_header(header, translator_class=translator_class, filename=filename,
                   search_path=search_path)

        # Store the supplied header for later stripping
        self._header = header

        # PropertyList is not dict-like so force to a dict here to simplify
        # the translation code.
        if hasattr(header, "toOrderedDict"):
            header = header.toOrderedDict()

        if translator_class is None:
            translator_class = MetadataTranslator.determine_translator(header, filename=filename)
        elif not issubclass(translator_class, MetadataTranslator):
            raise TypeError(f"Translator class must be a MetadataTranslator, not {translator_class}")

        # Create an instance for this header
        translator = translator_class(header, filename=filename)

        # Store the translator
        self._translator = translator
        self.translator_class_name = translator_class.__name__

        # Form file information string in case we need an error message
        if filename:
            file_info = f" and file {filename}"
        else:
            file_info = ""

        # Loop over each property and request the translated form
        for t in self._PROPERTIES:
            # prototype code
            method = f"to_{t}"
            property = f"_{t}"

            try:
                value = getattr(translator, method)()
            except NotImplementedError as e:
                raise NotImplementedError(f"No translation exists for property '{t}'"
                                          f" using translator {translator.__class__}") from e
            except KeyError as e:
                err_msg = f"Error calculating property '{t}' using translator {translator.__class__}" \
                    f"{file_info}"
                if pedantic:
                    raise KeyError(err_msg) from e
                else:
                    log.debug("Calculation of property '%s' failed with header: %s", t, header)
                    log.warning(f"Ignoring {err_msg}: {e}")
                    continue

            if value is not None and not isinstance(value, self._PROPERTIES[t][2]):
                err_msg = f"Value calculated for property '{t}' is wrong type " \
                    f"({type(value)} != {self._PROPERTIES[t][1]}) using translator {translator.__class__}" \
                    f"{file_info}"
                if pedantic:
                    raise TypeError(err_msg)
                else:
                    log.debug("Calcuation of property '%s' had unexpected type with header: %s", t, header)
                    log.warning(f"Ignoring {err_msg}")

            setattr(self, property, value)

    @property
    def cards_used(self):
        """Header cards used for the translation.

        Returns
        -------
        used : `frozenset` of `str`
            Set of card used.
        """
        if not self._translator:
            return frozenset()
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
        used = self.cards_used
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
            value = getattr(self, p)
            if isinstance(value, astropy.time.Time):
                value.format = "isot"
                value = str(value.value)
            result += f"{p}: {value}\n"

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

    def __lt__(self, other):
        return self.datetime_begin < other.datetime_begin

    def __gt__(self, other):
        return self.datetime_begin > other.datetime_begin

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

    @classmethod
    def makeObservationInfo(cls, **kwargs):  # noqa: N802
        """Construct an `ObservationInfo` from the supplied parameters.

        Notes
        -----
        The supplied parameters should use names matching the property.
        The type of the supplied value will be checked against the property.
        Any properties not supplied will be assigned a value of `None`.

        Raises
        ------
        KeyError
            Raised if a supplied parameter key is not a known property.
        TypeError
            Raised if a supplied value does not match the expected type
            of the property.
        """

        obsinfo = cls(None)

        unused = set(kwargs)

        for p in cls._PROPERTIES:
            if p in kwargs:
                property = f"_{p}"
                value = kwargs[p]
                if not isinstance(value, cls._PROPERTIES[p][2]):
                    raise TypeError(f"Supplied value {value} for property {p} "
                                    f"should be of class {cls._PROPERTIES[p][1]} not {value.__class__}")
                setattr(obsinfo, property, value)
                unused.remove(p)

        if unused:
            raise KeyError(f"Unrecognized properties provided: {', '.join(unused)}")

        return obsinfo


# Method to add the standard properties
def _make_property(property, doc, return_typedoc, return_type):
    """Create a getter method with associated docstring.

    Parameters
    ----------
    property : `str`
        Name of the property getter to be created.
    doc : `str`
        Description of this property.
    return_typedoc : `str`
        Type string of this property (used in the doc string).
    return_type : `class`
        Type of this property.

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
    {property} : `{return_typedoc}`
        Access the property.
    """
    return getter


# Initialize the internal properties (underscored) and add the associated
# getter methods.
for name, description in ObservationInfo._PROPERTIES.items():
    setattr(ObservationInfo, f"_{name}", None)
    setattr(ObservationInfo, name, property(_make_property(name, *description)))


def makeObservationInfo(**kwargs):  # noqa: N802
    """Construct an `ObservationInfo` from the supplied parameters.

    Notes
    -----
    The supplied parameters should use names matching the property.
    The type of the supplied value will be checked against the property.
    Any properties not supplied will be assigned a value of `None`.

    Raises
    ------
    KeyError
        Raised if a supplied parameter key is not a known property.
    TypeError
        Raised if a supplied value does not match the expected type
        of the property.
    """
    return ObservationInfo.makeObservationInfo(**kwargs)
