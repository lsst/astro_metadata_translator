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
import json
import math

import astropy.time
from astropy.coordinates import SkyCoord, AltAz

from .translator import MetadataTranslator
from .properties import PROPERTIES as _PROPERTIES  # renamed because we'll use the original name later
from .headers import fix_header

log = logging.getLogger(__name__)


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


class ObservationInfoMeta(type):
    """Metaclass for `ObservationInfo`

    We use the list of properties (in the ``PROPERTIES`` class variable) to set
    up internal properties (each of the properties with a leading underscore)
    and associated getter methods.

    This is done via a metaclass instead of a class decorator or other
    post-processing of the class in order to make subclassing of
    `ObservationInfo` convenient: the user writing the subclass doesn't need to
    worry about this detail, as the metaclass is inherited and fires
    automatically.
    """
    def __new__(cls, name, bases, attrs):
        assert "PROPERTIES" in attrs, "The PROPERTIES class variable must be defined"
        # Initialize the internal properties (underscored) and add the
        # associated getter methods.
        for name, description in attrs["PROPERTIES"].items():
            attrs[f"_{name}"] = None
            attrs[name] = property(_make_property(name, *description[:3]))
        return super().__new__(cls, name, bases, attrs)


class ObservationInfo(metaclass=ObservationInfoMeta):
    """Standardized representation of an instrument header for a single
    exposure observation.

    When subclassing in order to extend the list of properties, the only
    required update is the ``PROPERTIES`` class variable dict, e.g.:

        class MyObsInfo(ObservationInfo):
            PROPERTIES = dict(
                myInfo=("My special information", "str", str, None, None),
                **ObservationInfo.PROPERTIES
            )

    The tuple that must be provided for each property includes:
    - description of property
    - Python type of property as a string (suitable for docstrings)
    - Actual python type as a type
    - Simplification function (can be None)
    - Function to convert simple form back to required type (can be None)

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
    required : `set`, optional
        This parameter can be used to confirm that all properties contained
        in the set must translate correctly and also be non-None.  For the case
        where ``pedantic`` is `True` this will still check that the resulting
        value is not `None`.
    subset : `set`, optional
        If not `None`, controls the translations that will be performed
        during construction. This can be useful if the caller is only
        interested in a subset of the properties and knows that some of
        the others might be slow to compute (for example the airmass if it
        has to be derived).

    Raises
    ------
    ValueError
        Raised if the supplied header was not recognized by any of the
        registered translators. Also raised if the request property subset
        is not a subset of the known properties.
    TypeError
        Raised if the supplied translator class was not a MetadataTranslator.
    KeyError
        Raised if a required property cannot be calculated, or if pedantic
        mode is enabled and any translations fails.
    NotImplementedError
        Raised if the selected translator does not support a required
        property.

    Notes
    -----
    Headers will be corrected if correction files are located and this will
    modify the header provided to the constructor.
    """

    PROPERTIES = _PROPERTIES
    """All the properties supported by this class with associated
    documentation."""

    def __init__(self, header, filename=None, translator_class=None, pedantic=False,
                 search_path=None, required=None, subset=None):

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

        # Determine the properties of interest
        all_properties = set(self.PROPERTIES)
        if subset is not None:
            if not subset:
                raise ValueError("Cannot request no properties be calculated.")
            if not subset.issubset(all_properties):
                raise ValueError("Requested subset is not a subset of known properties. "
                                 f"Got extra: {subset - all_properties}")
            properties = subset
        else:
            properties = all_properties

        if required is None:
            required = set()
        else:
            if not required.issubset(all_properties):
                raise ValueError("Requested required properties include unknowns: "
                                 f"{required - all_properties}")

        # Loop over each property and request the translated form
        for t in properties:
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
                if pedantic or t in required:
                    raise KeyError(err_msg) from e
                else:
                    log.debug("Calculation of property '%s' failed with header: %s", t, header)
                    log.warning(f"Ignoring {err_msg}: {e}")
                    continue

            if not self._is_property_ok(t, value):
                err_msg = f"Value calculated for property '{t}' is wrong type " \
                    f"({type(value)} != {self.PROPERTIES[t][1]}) using translator {translator.__class__}" \
                    f"{file_info}"
                if pedantic or t in required:
                    raise TypeError(err_msg)
                else:
                    log.debug("Calcuation of property '%s' had unexpected type with header: %s", t, header)
                    log.warning(f"Ignoring {err_msg}")

            if value is None and t in required:
                raise KeyError(f"Calculation of required property {t} resulted in a value of None")

            setattr(self, property, value)

    @classmethod
    def _is_property_ok(cls, property, value):
        """Compare the supplied value against the expected type as defined
        for the corresponding property.

        Parameters
        ----------
        property : `str`
            Name of property.
        value : `object`
            Value of the property to validate.

        Returns
        -------
        is_ok : `bool`
            `True` if the value is of an appropriate type.

        Notes
        -----
        Currently only the type of the property is validated. There is no
        attempt to check bounds or determine that a Quantity is compatible
        with the property.
        """
        if value is None:
            return True

        property_type = cls.PROPERTIES[property][2]

        # For AltAz coordinates, they can either arrive as AltAz or
        # as SkyCoord(frame=AltAz) so try to find the frame inside
        # the SkyCoord.
        if issubclass(property_type, AltAz) and isinstance(value, SkyCoord):
            value = value.frame

        if not isinstance(value, property_type):
            return False

        return True

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
        properties = sorted(set(self.PROPERTIES.keys()) - set(priority))

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
        if not isinstance(other, ObservationInfo):
            return NotImplemented

        # Compare simplified forms.
        # Cannot compare directly because nan will not equate as equal
        # whereas they should be equal for our purposes
        self_simple = self.to_simple()
        other_simple = other.to_simple()

        for k, self_value in self_simple.items():
            other_value = other_simple[k]
            if self_value != other_value:
                if math.isnan(self_value) and math.isnan(other_value):
                    # If both are nan this is fine
                    continue
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
        for p in self.PROPERTIES:
            property = f"_{p}"
            state[p] = getattr(self, property)

        return state

    def __setstate__(self, state):
        for p in self.PROPERTIES:
            property = f"_{p}"
            setattr(self, property, state[p])

    def to_simple(self):
        """Convert the contents of this object to simple dict form.

        The keys of the dict are the standard properties but the values
        can be simplified to support JSON serialization. For example a
        SkyCoord might be represented as an ICRS RA/Dec tuple rather than
        a full SkyCoord representation.

        Any properties with `None` value will be skipped.

        Can be converted back to an `ObservationInfo` using `from_simple()`.

        Returns
        -------
        simple : `dict` of [`str`, `Any`]
            Simple dict of all properties.
        """
        simple = {}

        for p in self.PROPERTIES:
            property = f"_{p}"
            value = getattr(self, property)
            if value is None:
                continue

            # Access the function to simplify the property
            simplifier = self.PROPERTIES[p][3]

            if simplifier is None:
                simple[p] = value
                continue

            simple[p] = simplifier(value)

        return simple

    def to_json(self):
        """Serialize the object to JSON string.

        Returns
        -------
        j : `str`
            The properties of the ObservationInfo in JSON string form.
        """
        return json.dumps(self.to_simple())

    @classmethod
    def from_simple(cls, simple):
        """Convert the entity returned by `to_simple` back into an
        `ObservationInfo`.

        Parameters
        ----------
        simple : `dict` [`str`, `Any`]
            The dict returned by `to_simple()`

        Returns
        -------
        obsinfo : `ObservationInfo`
            New object constructed from the dict.
        """
        processed = {}
        for k, v in simple.items():

            if v is None:
                continue

            # Access the function to convert from simple form
            complexifier = cls.PROPERTIES[k][4]

            if complexifier is not None:
                v = complexifier(v, **processed)

            processed[k] = v

        return cls.makeObservationInfo(**processed)

    @classmethod
    def from_json(cls, json_str):
        """Create `ObservationInfo` from JSON string.

        Parameters
        ----------
        json_str : `str`
            The JSON representation.

        Returns
        -------
        obsinfo : `ObservationInfo`
            Reconstructed object.
        """
        simple = json.loads(json_str)
        return cls.from_simple(simple)

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

        for p in cls.PROPERTIES:
            if p in kwargs:
                property = f"_{p}"
                value = kwargs[p]
                if not cls._is_property_ok(p, value):
                    raise TypeError(f"Supplied value {value} for property {p} "
                                    f"should be of class {cls.PROPERTIES[p][1]} not {value.__class__}")
                setattr(obsinfo, property, value)
                unused.remove(p)

        if unused:
            n = len(unused)
            raise KeyError(f"Unrecognized propert{'y' if n == 1 else 'ies'} provided: {', '.join(unused)}")

        return obsinfo


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
