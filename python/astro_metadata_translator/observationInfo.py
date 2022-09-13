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

from __future__ import annotations

__all__ = ("ObservationInfo", "makeObservationInfo")

import copy
import itertools
import json
import logging
import math
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    FrozenSet,
    MutableMapping,
    Optional,
    Sequence,
    Set,
    Tuple,
    Type,
)

import astropy.time
from astropy.coordinates import AltAz, SkyCoord

from .headers import fix_header
from .properties import PROPERTIES, PropertyDefinition
from .translator import MetadataTranslator

if TYPE_CHECKING:
    import astropy.coordinates
    import astropy.units

log = logging.getLogger(__name__)


class ObservationInfo:
    """Standardized representation of an instrument header for a single
    exposure observation.

    There is a core set of instrumental properties that are pre-defined.
    Additional properties may be defined, either through the
    ``makeObservationInfo`` factory function by providing the ``extensions``
    definitions, or through the regular ``ObservationInfo`` constructor when
    the extensions have been defined in the ``MetadataTranslator`` for the
    instrument of interest (or in the provided ``translator_class``).

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

    Values of the properties are read-only.
    """

    # Static typing requires that we define the standard dynamic properties
    # statically.
    if TYPE_CHECKING:
        telescope: int
        instrument: str
        location: astropy.coordinates.EarthLocation
        exposure_id: int
        visit_id: int
        physical_filter: str
        datetime_begin: astropy.time.Time
        datetime_end: astropy.time.Time
        exposure_group: str
        exposure_time: astropy.units.Quantity
        dark_time: astropy.units.Quantity
        boresight_airmass: float
        boresight_rotation_angle: astropy.units.Quantity
        boresight_rotation_coord: str
        detector_num: int
        detector_name: str
        detector_serial: str
        detector_group: str
        detector_exposure_id: int
        focus_z: astropy.units.Quantity
        object: str
        temperature: astropy.units.Quantity
        pressure: astropy.units.Quantity
        relative_humidity: float
        tracking_radec: astropy.coordinates.SkyCoord
        altaz_begin: astropy.coordinates.AltAz
        science_program: str
        observation_counter: int
        observation_reason: str
        observation_type: str
        observation_id: str
        observing_day: int
        group_counter_start: int
        group_counter_end: int
        has_simulated_content: bool

    def __init__(
        self,
        header: Optional[MutableMapping[str, Any]],
        filename: Optional[str] = None,
        translator_class: Optional[Type[MetadataTranslator]] = None,
        pedantic: bool = False,
        search_path: Optional[Sequence[str]] = None,
        required: Optional[Set[str]] = None,
        subset: Optional[Set[str]] = None,
    ) -> None:

        # Initialize the empty object
        self._header: MutableMapping[str, Any] = {}
        self.filename = filename
        self._translator = None
        self.translator_class_name = "<None>"

        # To allow makeObservationInfo to work, we special case a None
        # header
        if header is None:
            return

        # Fix up the header (if required)
        fix_header(header, translator_class=translator_class, filename=filename, search_path=search_path)

        # Store the supplied header for later stripping
        self._header = header

        if translator_class is None:
            translator_class = MetadataTranslator.determine_translator(header, filename=filename)
        elif not issubclass(translator_class, MetadataTranslator):
            raise TypeError(f"Translator class must be a MetadataTranslator, not {translator_class}")

        self._declare_extensions(translator_class.extensions)

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
        full_set = set(self.all_properties)
        if subset is not None:
            if not subset:
                raise ValueError("Cannot request no properties be calculated.")
            if not subset.issubset(full_set):
                raise ValueError(
                    "Requested subset is not a subset of known properties. " f"Got extra: {subset - full_set}"
                )
            properties = subset
        else:
            properties = full_set

        if required is None:
            required = set()
        else:
            if not required.issubset(full_set):
                raise ValueError("Requested required properties include unknowns: " f"{required - full_set}")

        # Loop over each property and request the translated form
        for t in properties:
            # prototype code
            method = f"to_{t}"
            property = f"_{t}" if not t.startswith("ext_") else t

            try:
                value = getattr(translator, method)()
            except NotImplementedError as e:
                raise NotImplementedError(
                    f"No translation exists for property '{t}'" f" using translator {translator.__class__}"
                ) from e
            except Exception as e:
                err_msg = (
                    f"Error calculating property '{t}' using translator {translator.__class__}" f"{file_info}"
                )
                if pedantic or t in required:
                    raise KeyError(err_msg) from e
                else:
                    log.debug("Calculation of property '%s' failed with header: %s", t, header)
                    log.warning(f"Ignoring {err_msg}: {e}")
                    continue

            definition = self.all_properties[t]
            if not self._is_property_ok(definition, value):
                err_msg = (
                    f"Value calculated for property '{t}' is wrong type "
                    f"({type(value)} != {definition.str_type}) using translator {translator.__class__}"
                    f"{file_info}"
                )
                if pedantic or t in required:
                    raise TypeError(err_msg)
                else:
                    log.debug("Calcuation of property '%s' had unexpected type with header: %s", t, header)
                    log.warning(f"Ignoring {err_msg}")

            if value is None and t in required:
                raise KeyError(f"Calculation of required property {t} resulted in a value of None")

            super().__setattr__(property, value)  # allows setting even write-protected extensions

    @staticmethod
    def _get_all_properties(
        extensions: Optional[Dict[str, PropertyDefinition]] = None
    ) -> Dict[str, PropertyDefinition]:
        """Return the definitions of all properties

        Parameters
        ----------
        extensions : `dict` [`str`: `PropertyDefinition`]
            List of extension property definitions, indexed by name (with no
            "ext_" prefix).

        Returns
        -------
        properties : `dict` [`str`: `PropertyDefinition`]
            Merged list of all property definitions, indexed by name. Extension
            properties will be listed with an ``ext_`` prefix.
        """
        properties = dict(PROPERTIES)
        if extensions:
            properties.update({"ext_" + pp: dd for pp, dd in extensions.items()})
        return properties

    def _declare_extensions(self, extensions: Optional[Dict[str, PropertyDefinition]]) -> None:
        """Declare and set up extension properties

        This should always be called internally as part of the creation of a
        new `ObservationInfo`.

        The core set of properties each have a python ``property`` that makes
        them read-only, and serves as a useful place to hang the docstring.
        However, the core set are set up at compile time, whereas the extension
        properties have to be configured at run time (because we don't know
        what they will be until we look at the header and figure out what
        instrument we're dealing with) when we have an instance rather than a
        class (and python ``property`` doesn't work on instances; only on
        classes). We therefore use a separate scheme for the extension
        properties: we write them directly to their associated instance
        variable, and we use ``__setattr__`` to protect them as read-only.
        Unfortunately, with this scheme, we can't give extension properties a
        docstring; but we're setting them up at runtime, so maybe that's not
        terribly important.

        Parameters
        ----------
        extensions : `dict` [`str`: `PropertyDefinition`]
            List of extension property definitions, indexed by name (with no
            "ext_" prefix).
        """
        if not extensions:
            extensions = {}
        for name in extensions:
            super().__setattr__("ext_" + name, None)
        self.extensions = extensions
        self.all_properties = self._get_all_properties(extensions)

    def __setattr__(self, name: str, value: Any) -> Any:
        """Set attribute

        This provides read-only protection for the extension properties. The
        core set of properties have read-only protection via the use of the
        python ``property``.
        """
        if hasattr(self, "extensions") and name.startswith("ext_") and name[4:] in self.extensions:
            raise AttributeError(f"Attribute {name} is read-only")
        return super().__setattr__(name, value)

    @classmethod
    def _is_property_ok(cls, definition: PropertyDefinition, value: Any) -> bool:
        """Compare the supplied value against the expected type as defined
        for the corresponding property.

        Parameters
        ----------
        definition : `PropertyDefinition`
            Property definition.
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

        # For AltAz coordinates, they can either arrive as AltAz or
        # as SkyCoord(frame=AltAz) so try to find the frame inside
        # the SkyCoord.
        if issubclass(definition.py_type, AltAz) and isinstance(value, SkyCoord):
            value = value.frame

        if not isinstance(value, definition.py_type):
            return False

        return True

    @property
    def cards_used(self) -> FrozenSet[str]:
        """Header cards used for the translation.

        Returns
        -------
        used : `frozenset` of `str`
            Set of card used.
        """
        if not self._translator:
            return frozenset()
        return self._translator.cards_used()

    def stripped_header(self) -> MutableMapping[str, Any]:
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

    def __str__(self) -> str:
        # Put more interesting answers at front of list
        # and then do remainder
        priority = ("instrument", "telescope", "datetime_begin")
        properties = sorted(set(self.all_properties) - set(priority))

        result = ""
        for p in itertools.chain(priority, properties):
            value = getattr(self, p)
            if isinstance(value, astropy.time.Time):
                value.format = "isot"
                value = str(value.value)
            result += f"{p}: {value}\n"

        return result

    def __eq__(self, other: Any) -> bool:
        """Compares equal if standard properties are equal"""
        if not isinstance(other, ObservationInfo):
            return NotImplemented

        # Compare simplified forms.
        # Cannot compare directly because nan will not equate as equal
        # whereas they should be equal for our purposes
        self_simple = self.to_simple()
        other_simple = other.to_simple()

        # We don't care about the translator internal detail
        self_simple.pop("_translator", None)
        other_simple.pop("_translator", None)

        for k, self_value in self_simple.items():
            other_value = other_simple[k]
            if self_value != other_value:
                if math.isnan(self_value) and math.isnan(other_value):
                    # If both are nan this is fine
                    continue
                return False
        return True

    def __lt__(self, other: Any) -> bool:
        if not isinstance(other, ObservationInfo):
            return NotImplemented
        return self.datetime_begin < other.datetime_begin

    def __gt__(self, other: Any) -> bool:
        if not isinstance(other, ObservationInfo):
            return NotImplemented
        return self.datetime_begin > other.datetime_begin

    def __getstate__(self) -> Tuple[Any, ...]:
        """Get pickleable state

        Returns the properties.  Deliberately does not preserve the full
        current state; in particular, does not return the full header or
        translator.

        Returns
        -------
        state : `tuple`
            Pickled state.
        """
        state = dict()
        for p in self.all_properties:
            state[p] = getattr(self, p)

        return state, self.extensions

    def __setstate__(self, state: Tuple[Any, ...]) -> None:
        """Set object state from pickle

        Parameters
        ----------
        state : `tuple`
            Pickled state.
        """
        try:
            state, extensions = state
        except ValueError:
            # Backwards compatibility for pickles generated before DM-34175
            extensions = {}
        self._declare_extensions(extensions)
        for p in self.all_properties:
            if p.startswith("ext_"):
                # allows setting even write-protected extensions
                super().__setattr__(p, state[p])  # type: ignore
            else:
                property = f"_{p}"
                setattr(self, property, state[p])  # type: ignore

    def to_simple(self) -> MutableMapping[str, Any]:
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

        Notes
        -----
        Round-tripping of extension properties requires that the
        `ObservationInfo` was created with the help of a registered
        `MetadataTranslator` (which contains the extension property
        definitions).
        """
        simple = {}
        if hasattr(self, "_translator") and self._translator and self._translator.name:
            simple["_translator"] = self._translator.name

        for p in self.all_properties:
            property = f"_{p}" if not p.startswith("ext_") else p
            value = getattr(self, property)
            if value is None:
                continue

            # Access the function to simplify the property
            simplifier = self.all_properties[p].to_simple

            if simplifier is None:
                simple[p] = value
                continue

            simple[p] = simplifier(value)

        return simple

    def to_json(self) -> str:
        """Serialize the object to JSON string.

        Returns
        -------
        j : `str`
            The properties of the ObservationInfo in JSON string form.

        Notes
        -----
        Round-tripping of extension properties requires that the
        `ObservationInfo` was created with the help of a registered
        `MetadataTranslator` (which contains the extension property
        definitions).
        """
        return json.dumps(self.to_simple())

    @classmethod
    def from_simple(cls, simple: MutableMapping[str, Any]) -> ObservationInfo:
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

        Notes
        -----
        Round-tripping of extension properties requires that the
        `ObservationInfo` was created with the help of a registered
        `MetadataTranslator` (which contains the extension property
        definitions).
        """
        extensions = {}
        translator = simple.pop("_translator", None)
        if translator:
            if translator not in MetadataTranslator.translators:
                raise KeyError(f"Unrecognised translator: {translator}")
            extensions = MetadataTranslator.translators[translator].extensions

        properties = cls._get_all_properties(extensions)

        processed: Dict[str, Any] = {}
        for k, v in simple.items():

            if v is None:
                continue

            # Access the function to convert from simple form
            complexifier = properties[k].from_simple

            if complexifier is not None:
                v = complexifier(v, **processed)

            processed[k] = v

        return cls.makeObservationInfo(extensions=extensions, **processed)

    @classmethod
    def from_json(cls, json_str: str) -> ObservationInfo:
        """Create `ObservationInfo` from JSON string.

        Parameters
        ----------
        json_str : `str`
            The JSON representation.

        Returns
        -------
        obsinfo : `ObservationInfo`
            Reconstructed object.

        Notes
        -----
        Round-tripping of extension properties requires that the
        `ObservationInfo` was created with the help of a registered
        `MetadataTranslator` (which contains the extension property
        definitions).
        """
        simple = json.loads(json_str)
        return cls.from_simple(simple)

    @classmethod
    def makeObservationInfo(  # noqa: N802
        cls, *, extensions: Optional[Dict[str, PropertyDefinition]] = None, **kwargs: Any
    ) -> ObservationInfo:
        """Construct an `ObservationInfo` from the supplied parameters.

        Parameters
        ----------
        extensions : `dict` [`str`: `PropertyDefinition`], optional
            Optional extension definitions, indexed by extension name (without
            the ``ext_`` prefix, which will be added by `ObservationInfo`).
        **kwargs
            Name-value pairs for any properties to be set. In the case of
            extension properties, the names should include the ``ext_`` prefix.

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
        obsinfo._declare_extensions(extensions)

        unused = set(kwargs)

        for p in obsinfo.all_properties:
            if p in kwargs:
                property = f"_{p}" if not p.startswith("ext_") else p
                value = kwargs[p]
                definition = obsinfo.all_properties[p]
                if not cls._is_property_ok(definition, value):
                    raise TypeError(
                        f"Supplied value {value} for property {p} "
                        f"should be of class {definition.str_type} not {value.__class__}"
                    )
                super(cls, obsinfo).__setattr__(property, value)  # allows setting write-protected extensions
                unused.remove(p)

        # Recent additions to ObservationInfo may not be present in
        # serializations. In theory they can be derived from other
        # values in the default case. This might not be the right thing
        # to do.
        for k in ("group_counter_start", "group_counter_end"):
            if k not in kwargs and "observation_counter" in kwargs:
                super(cls, obsinfo).__setattr__(f"_{k}", obsinfo.observation_counter)
        if (k := "has_simulated_content") not in kwargs:
            super(cls, obsinfo).__setattr__(f"_{k}", False)

        if unused:
            n = len(unused)
            raise KeyError(f"Unrecognized propert{'y' if n == 1 else 'ies'} provided: {', '.join(unused)}")

        return obsinfo


# Method to add the standard properties
def _make_property(property: str, doc: str, return_typedoc: str, return_type: Type) -> Callable:
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

    def getter(self: ObservationInfo) -> Any:
        return getattr(self, f"_{property}")

    getter.__doc__ = f"""{doc}

    Returns
    -------
    {property} : `{return_typedoc}`
        Access the property.
    """
    return getter


# Set up the core set of properties
# In order to provide read-only protection, each attribute is hidden behind a
# python "property" wrapper.
for name, definition in PROPERTIES.items():
    setattr(ObservationInfo, f"_{name}", None)
    setattr(
        ObservationInfo,
        name,
        property(_make_property(name, definition.doc, definition.str_type, definition.py_type)),
    )


def makeObservationInfo(  # noqa: N802
    *, extensions: Optional[Dict[str, PropertyDefinition]] = None, **kwargs: Any
) -> ObservationInfo:
    """Construct an `ObservationInfo` from the supplied parameters.

    Parameters
    ----------
    extensions : `dict` [`str`: `PropertyDefinition`], optional
        Optional extension definitions, indexed by extension name (without
        the ``ext_`` prefix, which will be added by `ObservationInfo`).
    **kwargs
        Name-value pairs for any properties to be set. In the case of
        extension properties, the names should include the ``ext_`` prefix.

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
    return ObservationInfo.makeObservationInfo(extensions=extensions, **kwargs)
