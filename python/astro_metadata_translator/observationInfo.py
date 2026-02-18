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

"""Represent standard metadata from instrument headers."""

from __future__ import annotations

__all__ = ("ObservationInfo", "makeObservationInfo")

import copy
import itertools
import logging
import math
from collections.abc import MutableMapping, Sequence
from typing import TYPE_CHECKING, Any, cast, overload

import astropy.time
from lsst.resources import ResourcePath
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    PrivateAttr,
    ValidationInfo,
    field_validator,
    model_serializer,
)

from .headers import fix_header
from .properties import (
    PROPERTIES,
    PropertyDefinition,
)
from .translator import MetadataTranslator

if TYPE_CHECKING:
    import astropy.coordinates
    import astropy.units

log = logging.getLogger(__name__)
_CORE_FROM_SIMPLE_FIELDS = tuple(name for name, definition in PROPERTIES.items() if definition.from_simple)


class ObservationInfo(BaseModel):
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
        and a warning will be issued. Only used if a ``header`` is specified.
    search_path : `~collections.abc.Iterable`, optional
        Override search paths to use during header fix up. Only used if a
        ``header`` is specified.
    required : `set`, optional
        This parameter can be used to confirm that all properties contained
        in the set must translate correctly and also be non-None.  For the case
        where ``pedantic`` is `True` this will still check that the resulting
        value is not `None`. Only used if a ``header`` is specified.
    subset : `set`, optional
        If not `None`, controls the translations that will be performed
        during construction. This can be useful if the caller is only
        interested in a subset of the properties and knows that some of
        the others might be slow to compute (for example the airmass if it
        has to be derived). Only used if a ``header`` is specified.
    **kwargs : `typing.Any`
        Property name/value pairs for kwargs-based construction mode. This
        mode creates an `ObservationInfo` directly from supplied properties
        rather than by translating a header. If ``header`` is provided it is
        an error to also provide ``kwargs``.

    Raises
    ------
    ValueError
        Raised if the supplied header was not recognized by any of the
        registered translators. Also raised if the request property subset
        is not a subset of the known properties or if a header is given along
        with kwargs.
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
    There is a core set of instrumental properties that are pre-defined.
    Additional properties may be defined, either through the
    `makeObservationInfo` factory function by providing the ``extensions``
    definitions, or through the regular `ObservationInfo` constructor when
    the extensions have been defined in the `MetadataTranslator` for the
    instrument of interest (or in the provided ``translator_class``).

    There are two forms of the constructor. If the ``header`` is given
    then a translator will be determined and the properties will be populated
    accordingly. No generic keyword arguments will be expected and the
    remaining parameters control the behavior of the translator.

    If the header is not given it is assumed that the keyword arguments
    are direct specifications of observation properties. In this mode only
    the ``filename`` and ``translator_class`` parameters will be used. The
    latter is used to determine any extensions that are being provided,
    although when using standard serializations the special ``_translator``
    key will be used instead to specify the name of the registered translator
    from which to extract extension definitions.

    Headers will be corrected if correction files are located and this will
    modify the header provided to the constructor. Modifying the supplied
    header after construction will modify the internal cached header.

    Values of the properties are read-only.
    """

    model_config = ConfigDict(
        extra="forbid",
        arbitrary_types_allowed=True,
        validate_assignment=False,
        ser_json_inf_nan="constants",  # Allow for inf and nan to round trip.
    )

    filename: str | None = Field(default=None, exclude=True)
    translator_class_name: str = Field(default="<None>", exclude=True)
    extensions: dict[str, PropertyDefinition] = Field(default_factory=dict, exclude=True)
    all_properties: dict[str, PropertyDefinition] = Field(default_factory=dict, exclude=True)
    telescope: str | None = Field(default=None, description=PROPERTIES["telescope"].doc)
    instrument: str | None = Field(default=None, description=PROPERTIES["instrument"].doc)
    location: astropy.coordinates.EarthLocation | None = Field(
        default=None, description=PROPERTIES["location"].doc
    )
    exposure_id: int | None = Field(default=None, description=PROPERTIES["exposure_id"].doc)
    visit_id: int | None = Field(default=None, description=PROPERTIES["visit_id"].doc)
    physical_filter: str | None = Field(default=None, description=PROPERTIES["physical_filter"].doc)
    datetime_begin: astropy.time.Time | None = Field(
        default=None, description=PROPERTIES["datetime_begin"].doc
    )
    datetime_end: astropy.time.Time | None = Field(default=None, description=PROPERTIES["datetime_end"].doc)
    exposure_time: astropy.units.Quantity | None = Field(
        default=None, description=PROPERTIES["exposure_time"].doc
    )
    exposure_time_requested: astropy.units.Quantity | None = Field(
        default=None, description=PROPERTIES["exposure_time_requested"].doc
    )
    dark_time: astropy.units.Quantity | None = Field(default=None, description=PROPERTIES["dark_time"].doc)
    boresight_airmass: float | None = Field(default=None, description=PROPERTIES["boresight_airmass"].doc)
    boresight_rotation_angle: astropy.coordinates.Angle | None = Field(
        default=None, description=PROPERTIES["boresight_rotation_angle"].doc
    )
    boresight_rotation_coord: str | None = Field(
        default=None, description=PROPERTIES["boresight_rotation_coord"].doc
    )
    detector_num: int | None = Field(default=None, description=PROPERTIES["detector_num"].doc)
    detector_name: str | None = Field(default=None, description=PROPERTIES["detector_name"].doc)
    detector_unique_name: str | None = Field(default=None, description=PROPERTIES["detector_unique_name"].doc)
    detector_serial: str | None = Field(default=None, description=PROPERTIES["detector_serial"].doc)
    detector_group: str | None = Field(default=None, description=PROPERTIES["detector_group"].doc)
    detector_exposure_id: int | None = Field(default=None, description=PROPERTIES["detector_exposure_id"].doc)
    focus_z: astropy.units.Quantity | None = Field(default=None, description=PROPERTIES["focus_z"].doc)
    object: str | None = Field(default=None, description=PROPERTIES["object"].doc)
    temperature: astropy.units.Quantity | None = Field(
        default=None, description=PROPERTIES["temperature"].doc
    )
    pressure: astropy.units.Quantity | None = Field(default=None, description=PROPERTIES["pressure"].doc)
    relative_humidity: float | None = Field(default=None, description=PROPERTIES["relative_humidity"].doc)
    tracking_radec: astropy.coordinates.SkyCoord | None = Field(
        default=None, description=PROPERTIES["tracking_radec"].doc
    )
    altaz_begin: astropy.coordinates.AltAz | None = Field(
        default=None, description=PROPERTIES["altaz_begin"].doc
    )
    altaz_end: astropy.coordinates.AltAz | None = Field(default=None, description=PROPERTIES["altaz_end"].doc)
    science_program: str | None = Field(default=None, description=PROPERTIES["science_program"].doc)
    observation_type: str | None = Field(default=None, description=PROPERTIES["observation_type"].doc)
    observation_id: str | None = Field(default=None, description=PROPERTIES["observation_id"].doc)
    observation_reason: str | None = Field(default=None, description=PROPERTIES["observation_reason"].doc)
    exposure_group: str | None = Field(default=None, description=PROPERTIES["exposure_group"].doc)
    observing_day: int | None = Field(default=None, description=PROPERTIES["observing_day"].doc)
    observing_day_offset: astropy.time.TimeDelta | None = Field(
        default=None, description=PROPERTIES["observing_day_offset"].doc
    )
    observation_counter: int | None = Field(default=None, description=PROPERTIES["observation_counter"].doc)
    has_simulated_content: bool | None = Field(
        default=None, description=PROPERTIES["has_simulated_content"].doc
    )
    group_counter_start: int | None = Field(default=None, description=PROPERTIES["group_counter_start"].doc)
    group_counter_end: int | None = Field(default=None, description=PROPERTIES["group_counter_end"].doc)
    can_see_sky: bool | None = Field(default=None, description=PROPERTIES["can_see_sky"].doc)

    _header: MutableMapping[str, Any] = PrivateAttr(default_factory=dict)
    _translator: MetadataTranslator | None = PrivateAttr(default=None)
    _sealed: bool = PrivateAttr(default=False)

    @field_validator(*_CORE_FROM_SIMPLE_FIELDS, mode="before")
    @classmethod
    def _before_core_from_simple(cls, value: Any, info: ValidationInfo) -> Any:
        assert info.field_name is not None
        definition = PROPERTIES[info.field_name]
        context = info.data if isinstance(info.data, dict) else {}
        return cls._coerce_from_simple(definition, value, context)

    @overload
    def __init__(
        self,
        header: MutableMapping[str, Any],
        filename: str | ResourcePath | None = None,
        translator_class: type[MetadataTranslator] | None = None,
        pedantic: bool = False,
        search_path: Sequence[str] | None = None,
        required: set[str] | None = None,
        subset: set[str] | None = None,
    ) -> None: ...

    @overload
    def __init__(
        self,
        header: None = None,
        filename: str | ResourcePath | None = None,
        translator_class: type[MetadataTranslator] | None = None,
        **kwargs: Any,
    ) -> None: ...

    def __init__(
        self,
        header: MutableMapping[str, Any] | None = None,
        filename: str | ResourcePath | None = None,
        translator_class: type[MetadataTranslator] | None = None,
        pedantic: bool = False,
        search_path: Sequence[str] | None = None,
        required: set[str] | None = None,
        subset: set[str] | None = None,
        **kwargs: Any,
    ) -> None:
        if filename is not None:
            filename = str(ResourcePath(filename, forceAbsolute=True))
        if header is not None:
            if kwargs:
                raise ValueError(
                    "kwargs not allowed if constructor given a header to translate. "
                    f"Unrecognized keys: {[k for k in kwargs]}"
                )
            self._init_from_header(
                header,
                filename=filename,
                translator_class=translator_class,
                pedantic=pedantic,
                search_path=search_path,
                required=required,
                subset=subset,
            )
            return

        self._init_from_kwargs(filename=filename, translator_class=translator_class, **kwargs)

    def _init_from_header(
        self,
        header: MutableMapping[str, Any],
        *,
        filename: str | None,
        translator_class: type[MetadataTranslator] | None,
        pedantic: bool,
        search_path: Sequence[str] | None,
        required: set[str] | None,
        subset: set[str] | None,
    ) -> None:
        super().__init__(filename=filename)
        self._sealed = False
        # Initialize the empty object
        self._header = {}
        self._translator = None

        # Look for translator class before header fixup. fix_header calls
        # determine_translator immediately on the basis that you need to know
        # enough of the header to work out the translator before you can fix
        # it up. There is no gain in asking fix_header to determine the
        # translator and then trying to work it out again here.
        if translator_class is None:
            translator_class = MetadataTranslator.determine_translator(header, filename=filename)
        elif not issubclass(translator_class, MetadataTranslator):
            raise TypeError(f"Translator class must be a MetadataTranslator, not {translator_class}")

        # Fix up the header (if required)
        fix_header(header, translator_class=translator_class, filename=filename, search_path=search_path)

        # Store the supplied header for later stripping
        self._header = header

        # This configures both self.extensions and self.all_properties.
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
                    f"Requested subset is not a subset of known properties. Got extra: {subset - full_set}"
                )
            properties = subset
        else:
            properties = full_set

        if required is None:
            required = set()
        else:
            if not required.issubset(full_set):
                raise ValueError(f"Requested required properties include unknowns: {required - full_set}")

        # Loop over each property and request the translated form
        for property in properties:
            # prototype code
            method = f"to_{property}"

            try:
                value = getattr(translator, method)()
            except NotImplementedError as e:
                raise NotImplementedError(
                    f"No translation exists for property '{property}' using translator {translator.__class__}"
                ) from e
            except Exception as e:
                err_msg = (
                    f"Error calculating property '{property}' using "
                    f"translator {translator.__class__}{file_info}"
                )
                if pedantic or property in required:
                    raise KeyError(err_msg) from e
                else:
                    log.debug("Calculation of property '%s' failed with header: %s", property, header)
                    log.warning(f"Ignoring {err_msg}: {e}")
                    continue

            definition = self.all_properties[property]
            # Some translators can return a compatible form that needs to
            # be coerced to the correct type (e.g., returning SkyCoord when you
            # need AltAz). In theory we could patch the translators to return
            # AltAz but code has historically not been as picky about this
            # until pydantic turned up.
            value = self._coerce_from_simple(definition, value, {})
            if not definition.is_value_conformant(value):
                err_msg = (
                    f"Value calculated for property '{property}' is wrong type "
                    f"({type(value)} != {definition.str_type}) using translator {translator.__class__}"
                    f"{file_info}"
                )
                if pedantic or property in required:
                    raise TypeError(err_msg)
                else:
                    log.debug(
                        "Calculation of property '%s' had unexpected type with header: %s", property, header
                    )
                    log.warning(f"Ignoring {err_msg}")

            if value is None and property in required:
                raise KeyError(f"Calculation of required property {property} resulted in a value of None")

            object.__setattr__(self, property, value)  # allows setting even write-protected extensions

        self._sealed = True

    def _init_from_kwargs(
        self,
        *,
        filename: str | None,
        translator_class: type[MetadataTranslator] | None,
        **kwargs: Any,
    ) -> None:
        supplied_keys = set(kwargs)
        translator_name = kwargs.pop("_translator", None)
        supplied_extensions = kwargs.pop("_extensions", None)
        if translator_name is not None:
            if translator_name not in MetadataTranslator.translators:
                raise KeyError(f"Unrecognized translator: {translator_name}")
            translator_class = MetadataTranslator.translators[translator_name]

        if translator_class is not None and not issubclass(translator_class, MetadataTranslator):
            raise TypeError(f"Translator class must be a MetadataTranslator, not {translator_class}")

        if supplied_extensions is not None:
            if translator_class is not None:
                raise ValueError("Provide either translator_class or _extensions, not both.")
            if not isinstance(supplied_extensions, dict):
                raise TypeError("_extensions must be a dictionary of PropertyDefinition entries.")
            extensions = supplied_extensions
        else:
            extensions = translator_class.extensions if translator_class is not None else {}

        all_properties = self._get_all_properties(extensions)
        for key in kwargs:
            if key not in all_properties:
                raise KeyError(f"Unrecognized property '{key}' provided")

        processed = {k: v for k, v in kwargs.items() if k in PROPERTIES and v is not None}
        processed = self._apply_constructor_defaults(processed, supplied_keys)

        super().__init__(filename=filename, **processed)
        self._sealed = False

        # This configures both self.extensions and self.all_properties.
        self._declare_extensions(extensions)

        # Handle extensions.
        ext_input = {k: v for k, v in kwargs.items() if k.startswith("ext_")}
        processed_ext = self._validate_property_mapping(ext_input, extensions)
        for key, value in processed_ext.items():
            object.__setattr__(self, key, value)

        if translator_class is not None:
            self._translator = translator_class({})
            self.translator_class_name = translator_class.__name__

        self._sealed = True

    @staticmethod
    def _apply_constructor_defaults(processed: dict[str, Any], supplied_keys: set[str]) -> dict[str, Any]:
        """Apply derived/default values for kwargs-style construction.

        Parameters
        ----------
        processed : `dict` [`str`, `typing.Any`]
            Properties validated from kwargs input.
        supplied_keys : `set` [`str`]
            Property names explicitly supplied by the caller.

        Returns
        -------
        updated : `dict` [`str`, `typing.Any`]
            Updated property mapping with defaults/backfilled values applied.
        """
        updated = dict(processed)
        for key in ("group_counter_start", "group_counter_end"):
            if (
                key not in supplied_keys
                and "observation_counter" in supplied_keys
                and "observation_counter" in updated
            ):
                updated[key] = updated["observation_counter"]
        if "has_simulated_content" not in supplied_keys:
            updated["has_simulated_content"] = False
        return updated

    @classmethod
    def from_header(
        cls,
        header: MutableMapping[str, Any],
        *,
        filename: str | None = None,
        translator_class: type[MetadataTranslator] | None = None,
        pedantic: bool = False,
        search_path: Sequence[str] | None = None,
        required: set[str] | None = None,
        subset: set[str] | None = None,
    ) -> ObservationInfo:
        """Create an `ObservationInfo` by translating a metadata header.

        Parameters
        ----------
        header : `dict`-like
            Header mapping to translate.
        filename : `str`, optional
            Name of file associated with this header.
        translator_class : `MetadataTranslator`-class, optional
            Translator class to use. If `None`, translator will be
            auto-determined.
        pedantic : `bool`, optional
            If `True`, translation failures are fatal.
        search_path : `~collections.abc.Sequence` [`str`], optional
            Search paths for header corrections.
        required : `set` [`str`], optional
            Properties that must be translated and non-`None`.
        subset : `set` [`str`], optional
            Restrict translation to this subset of properties.

        Returns
        -------
        obsinfo : `ObservationInfo`
            Translated observation metadata.
        """
        return cls(
            header=header,
            filename=filename,
            translator_class=translator_class,
            pedantic=pedantic,
            search_path=search_path,
            required=required,
            subset=subset,
        )

    @staticmethod
    def _get_all_properties(
        extensions: dict[str, PropertyDefinition] | None = None,
    ) -> dict[str, PropertyDefinition]:
        """Return the definitions of all properties.

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

    def _declare_extensions(self, extensions: dict[str, PropertyDefinition] | None) -> None:
        """Declare and set up extension properties.

        This should always be called internally as part of the creation of a
        new `ObservationInfo`.

        The core set of properties are declared as model fields at import
        time. Extension properties have to be configured at runtime (because
        we don't know what they will be until we look at the header and figure
        out what instrument we're dealing with), so we add them to the model
        and then use ``__setattr__`` to protect them as read-only. All
        extension properties are set to `None`.

        Parameters
        ----------
        extensions : `dict` [`str`: `PropertyDefinition`]
            List of extension property definitions, indexed by name (with no
            "ext_" prefix).
        """
        if extensions:
            for name in extensions:
                field_name = "ext_" + name
                if not hasattr(self, field_name):
                    object.__setattr__(self, field_name, None)
            self.extensions = extensions
        self.all_properties = self._get_all_properties(self.extensions)

    def __setattr__(self, name: str, value: Any) -> Any:
        """Set attribute.

        This provides read-only protection for all properties once the
        instance has been sealed.

        Parameters
        ----------
        name : `str`
            Name of attribute to set.
        value : `typing.Any`
            Value to set it to.
        """
        if (
            getattr(self, "_sealed", False)
            and hasattr(self, "all_properties")
            and name in self.all_properties
        ):
            raise AttributeError(f"Attribute {name} is read-only")
        return super().__setattr__(name, value)

    @model_serializer(mode="plain")
    def _serialize(self) -> dict[str, Any]:
        simple: dict[str, Any] = {}
        if self._translator and self._translator.name:
            simple["_translator"] = self._translator.name

        for p, definition in self.all_properties.items():
            value = getattr(self, p)
            if value is None:
                continue
            simplifier = definition.to_simple
            if simplifier is None:
                simple[p] = value
            else:
                simple[p] = simplifier(value)

        return simple

    @classmethod
    def _validate_property_mapping(
        cls,
        data: MutableMapping[str, Any],
        extensions: dict[str, PropertyDefinition] | None,
    ) -> dict[str, Any]:
        # Validate extension properties.
        properties = {f"ext_{name}": definition for name, definition in (extensions or {}).items()}
        processed: dict[str, Any] = {}

        for key, value in data.items():
            if key not in properties:
                raise KeyError(f"Unrecognized property '{key}' provided")
            if value is None:
                continue
            processed[key] = cls._coerce_property_value(key, value, properties, processed)
        return processed

    @classmethod
    def _coerce_property_value(
        cls,
        key: str,
        value: Any,
        properties: dict[str, PropertyDefinition],
        processed: dict[str, Any],
    ) -> Any:
        definition = properties[key]
        converted = cls._coerce_from_simple(definition, value, processed)
        if not definition.is_value_conformant(converted):
            raise TypeError(
                f"Supplied value {value} for property {key} "
                f"should be of class {definition.str_type} not {converted.__class__}"
            )
        return converted

    @classmethod
    def _coerce_from_simple(
        cls,
        definition: PropertyDefinition,
        value: Any,
        processed: dict[str, Any],
    ) -> Any:
        if definition.is_value_conformant(value):
            return value
        complexifier = definition.from_simple
        if complexifier is None:
            # Not the correct type, assumes caller will check.
            return value
        return complexifier(value, **processed)

    @property
    def cards_used(self) -> frozenset[str]:
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
            if c in hdr:
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
        """Check equality with another object.

        Compares equal if standard properties are equal.

        Parameters
        ----------
        other : `typing.Any`
            Thing to compare with.
        """
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

        for k in self_simple.keys() & other_simple.keys():
            self_value = self_simple[k]
            other_value = other_simple[k]
            if self_value != other_value:
                try:
                    both_nan = math.isnan(self_value) and math.isnan(other_value)
                except TypeError:
                    both_nan = False
                if both_nan:
                    # If both are nan this is fine
                    continue
                return False
        return True

    def __lt__(self, other: Any) -> bool:
        if not isinstance(other, ObservationInfo):
            return NotImplemented
        if self.datetime_begin is None or other.datetime_begin is None:
            raise TypeError("Cannot compare ObservationInfo without datetime_begin values")
        return self.datetime_begin < other.datetime_begin

    def __gt__(self, other: Any) -> bool:
        if not isinstance(other, ObservationInfo):
            return NotImplemented
        if self.datetime_begin is None or other.datetime_begin is None:
            raise TypeError("Cannot compare ObservationInfo without datetime_begin values")
        return self.datetime_begin > other.datetime_begin

    def __getstate__(self) -> dict[str, Any]:
        """Get pickleable state.

        Returns the properties.  Deliberately does not preserve the full
        current state; in particular, does not return the full header or
        translator.

        Returns
        -------
        state : `dict`
            Pickled state.
        """
        state: dict[str, Any] = {}
        for p in self.all_properties:
            state[p] = getattr(self, p)

        return {"state": state, "extensions": self.extensions}

    def __setstate__(self, state: dict[Any, Any]) -> None:
        """Set object state from pickle.

        Parameters
        ----------
        state : `dict`
            Pickled state.
        """
        state_any = cast(Any, state)
        if isinstance(state_any, dict) and "state" in state_any:
            state = state_any["state"]
            extensions = state_any.get("extensions", {})
        else:
            try:
                state, extensions = state_any
            except ValueError:
                # Backwards compatibility for pickles generated before DM-34175
                extensions = {}
        super().__init__()
        self._sealed = False
        self._declare_extensions(extensions)
        for p in self.all_properties:
            # allows setting even write-protected extensions
            object.__setattr__(self, p, state[p])
        self._sealed = True

    def to_simple(self) -> MutableMapping[str, Any]:
        """Convert the contents of this object to simple dict form.

        The keys of the dict are the standard properties but the values
        can be simplified to support JSON serialization. For example a
        SkyCoord might be represented as an ICRS RA/Dec tuple rather than
        a full SkyCoord representation.

        Any properties with `None` value will be skipped.

        Can be converted back to an `ObservationInfo` using `from_simple`.

        Returns
        -------
        simple : `dict` of [`str`, `~typing.Any`]
            Simple dict of all properties.

        Notes
        -----
        Round-tripping of extension properties requires that the
        `ObservationInfo` was created with the help of a registered
        `MetadataTranslator` (which contains the extension property
        definitions).
        """
        return self.model_dump(mode="python")

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
        return self.model_dump_json()

    @classmethod
    def from_simple(cls, simple: MutableMapping[str, Any]) -> ObservationInfo:
        """Convert the entity returned by `to_simple` back into an
        `ObservationInfo`.

        Parameters
        ----------
        simple : `dict` [`str`, `~typing.Any`]
            The dict returned by `to_simple`.

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
        return cls.model_validate(simple)

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
        return cls.model_validate_json(json_str)

    @classmethod
    def makeObservationInfo(  # noqa: N802
        cls,
        *,
        extensions: dict[str, PropertyDefinition] | None = None,
        translator_class: type[MetadataTranslator] | None = None,
        **kwargs: Any,
    ) -> ObservationInfo:
        """Construct an `ObservationInfo` from the supplied parameters.

        Parameters
        ----------
        extensions : `dict` [`str`: `PropertyDefinition`], optional
            Optional extension definitions, indexed by extension name (without
            the ``ext_`` prefix, which will be added by `ObservationInfo`).
        translator_class : `MetadataTranslator`-class, optional
            Optional translator class defining the extension properties. If
            provided, this can be used instead of ``extensions`` and will be
            stored in the instance for JSON round-tripping.
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
        return cls(filename=None, translator_class=translator_class, _extensions=extensions, **kwargs)


def makeObservationInfo(  # noqa: N802
    *,
    extensions: dict[str, PropertyDefinition] | None = None,
    translator_class: type[MetadataTranslator] | None = None,
    **kwargs: Any,
) -> ObservationInfo:
    """Construct an `ObservationInfo` from the supplied parameters.

    Parameters
    ----------
    extensions : `dict` [`str`: `PropertyDefinition`], optional
        Optional extension definitions, indexed by extension name (without
        the ``ext_`` prefix, which will be added by `ObservationInfo`).
    translator_class : `MetadataTranslator`-class, optional
        Optional translator class defining the extension properties. If
        provided, this can be used instead of ``extensions`` and will be
        stored in the instance for JSON round-tripping.
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
    return ObservationInfo.makeObservationInfo(
        extensions=extensions, translator_class=translator_class, **kwargs
    )
