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

"""Classes and support code for metadata translation"""

__all__ = ("MetadataTranslator", "StubTranslator", "cache_translation")

from abc import abstractmethod
import logging
import warnings
import math

import astropy.units as u
from astropy.coordinates import Angle

from .properties import PROPERTIES


log = logging.getLogger(__name__)


def cache_translation(func, method=None):
    """Decorator to cache the result of a translation method.

    Especially useful when a translation uses many other translation
    methods.  Should be used only on ``to_x()`` methods.

    Parameters
    ----------
    func : `function`
        Translation method to cache.
    method : `str`, optional
        Name of the translation method to cache.  Not needed if the decorator
        is used around a normal method, but necessary when the decorator is
        being used in a metaclass.

    Returns
    -------
    wrapped : `function`
        Method wrapped by the caching function.
    """
    name = func.__name__ if method is None else method

    def func_wrapper(self):
        if name not in self._translation_cache:
            self._translation_cache[name] = func(self)
        return self._translation_cache[name]
    func_wrapper.__doc__ = func.__doc__
    func_wrapper.__name__ = f"{name}_cached"
    return func_wrapper


class MetadataTranslator:
    """Per-instrument metadata translation support

    Parameters
    ----------
    header : `dict`-like
        Representation of an instrument header that can be manipulated
        as if it was a `dict`.
    filename : `str`, optional
        Name of the file whose header is being translated.  For some
        datasets with missing header information this can sometimes
        allow for some fixups in translations.
    """

    # These are all deliberately empty in the base class.

    _trivial_map = {}
    """Dict of one-to-one mappings for header translation from standard
    property to corresponding keyword."""

    _const_map = {}
    """Dict defining a constant for specified standard properties."""

    translators = dict()
    """All registered metadata translation classes."""

    supported_instrument = None
    """Name of instrument understood by this translation class."""

    @staticmethod
    def _make_const_mapping(property_key, constant):
        """Make a translator method that returns a constant value.

        Parameters
        ----------
        property_key : `str`
            Name of the property to be calculated (for the docstring).
        constant : `str` or `numbers.Number`
            Value to return for this translator.

        Returns
        -------
        f : `function`
            Function returning the constant.
        """
        def constant_translator(self):
            return constant

        if property_key in PROPERTIES:
            property_doc, return_type = PROPERTIES[property_key]
        else:
            return_type = type(constant).__name__
            property_doc = f"Returns constant value for '{property_key}' property"

        constant_translator.__doc__ = f"""{property_doc}

        Returns
        -------
        translation : `{return_type}`
            Translated property.
        """
        return constant_translator

    @staticmethod
    def _make_trivial_mapping(property_key, header_key, default=None, minimum=None, maximum=None,
                              unit=None, checker=None):
        """Make a translator method returning a header value.

        The header value can be converted to a `~astropy.units.Quantity`
        if desired, and can also have its value validated.

        See `MetadataTranslator.validate_value()` for details on the use
        of default parameters.

        Parameters
        ----------
        property_key : `str`
            Name of the translator to be constructed (for the docstring).
        header_key : `str` or `list` of `str`
            Name of the key to look up in the header. If a `list` each
            key will be tested in turn until one matches.  This can deal with
            header styles that evolve over time.
        default : `numbers.Number` or `astropy.units.Quantity`, `str`, optional
            If not `None`, default value to be used if the parameter read from
            the header is not defined or if the header is missing.
        minimum : `numbers.Number` or `astropy.units.Quantity`, optional
            If not `None`, and if ``default`` is not `None`, minimum value
            acceptable for this parameter.
        maximum : `numbers.Number` or `astropy.units.Quantity`, optional
            If not `None`, and if ``default`` is not `None`, maximum value
            acceptable for this parameter.
        unit : `astropy.units.Unit`, optional
            If not `None`, the value read from the header will be converted
            to a `~astropy.units.Quantity`.  Only supported for numeric values.
        checker : `function`, optional
            Callback function to be used by the translator method in case the
            keyword is not present.  Function will be executed as if it is
            a method of the translator class.  Running without raising an
            exception will allow the default to be used. Should usually raise
            `KeyError`.

        Returns
        -------
        t : `function`
            Function implementing a translator with the specified
            parameters.
        """
        if property_key in PROPERTIES:
            property_doc, return_type = PROPERTIES[property_key]
        else:
            return_type = "str` or `numbers.Number"
            property_doc = f"Map '{header_key}' header keyword to '{property_key}' property"

        def trivial_translator(self):
            if unit is not None:
                q = self.quantity_from_card(header_key, unit,
                                            default=default, minimum=minimum, maximum=maximum,
                                            checker=checker)
                # Convert to Angle if this quantity is an angle
                if return_type == "astropy.coordinates.Angle":
                    q = Angle(q)
                return q

            keywords = header_key if isinstance(header_key, list) else [header_key]
            for key in keywords:
                if key in self._header:
                    value = self._header[key]
                    if default is not None and not isinstance(value, str):
                        value = self.validate_value(value, default, minimum=minimum, maximum=maximum)
                    self._used_these_cards(key)
                    break
            else:
                # No keywords found, use default, checking first, or raise
                # A None default is only allowed if a checker is provided.
                if checker is not None:
                    try:
                        checker(self)
                        return default
                    except Exception:
                        raise KeyError(f"Could not find {keywords} in header")
                    value = default
                elif default is not None:
                    value = default
                else:
                    raise KeyError(f"Could not find {keywords} in header")

            # If we know this is meant to be a string, force to a string.
            # Sometimes headers represent items as integers which generically
            # we want as strings (eg OBSID).  Sometimes also floats are
            # written as "NaN" strings.
            if return_type == "str" and not isinstance(value, str):
                value = str(value)
            elif return_type == "float" and not isinstance(value, float):
                value = float(value)
            return value

        # Docstring inheritance means it is confusing to specify here
        # exactly which header value is being used.
        trivial_translator.__doc__ = f"""{property_doc}

        Returns
        -------
        translation : `{return_type}`
            Translated value derived from the header.
        """
        return trivial_translator

    @classmethod
    def __init_subclass__(cls, **kwargs):
        """Register all subclasses with the base class and create dynamic
        translator methods.

        The method provides two facilities.  Firstly, every subclass
        of `MetadataTranslator` that includes a ``name`` class property is
        registered as a translator class that could be selected when automatic
        header translation is attempted.  Only name translator subclasses that
        correspond to a complete instrument.  Translation classes providing
        generic translation support for multiple instrument translators should
        not be named.

        The second feature of this method is to convert simple translations
        to full translator methods.  Sometimes a translation is fixed (for
        example a specific instrument name should be used) and rather than provide
        a full ``to_property()`` translation method the mapping can be defined
        in a class variable named ``_constMap``.  Similarly, for one-to-one
        trivial mappings from a header to a property, ``_trivialMap`` can be
        defined.  Trivial mappings are a dict mapping a generic property
        to either a header keyword, or a tuple consisting of the header keyword
        and a dict containing key value pairs suitable for the
        `MetadataTranslator.quantity_from_card()` method.
        """
        super().__init_subclass__(**kwargs)

        # Only register classes with declared names
        if hasattr(cls, "name") and cls.name is not None:
            MetadataTranslator.translators[cls.name] = cls

        # Go through the trival mappings for this class and create
        # corresponding translator methods
        for property_key, header_key in cls._trivial_map.items():
            kwargs = {}
            if type(header_key) == tuple:
                kwargs = header_key[1]
                header_key = header_key[0]
            translator = cls._make_trivial_mapping(property_key, header_key, **kwargs)
            method = f"to_{property_key}"
            translator.__name__ = f"{method}_trivial"
            setattr(cls, method, cache_translation(translator, method=method))
            if property_key not in PROPERTIES:
                log.warning(f"Unexpected trivial translator for '{property_key}' defined in {cls}")

        # Go through the constant mappings for this class and create
        # corresponding translator methods
        for property_key, constant in cls._const_map.items():
            translator = cls._make_const_mapping(property_key, constant)
            method = f"to_{property_key}"
            translator.__name__ = f"{method}_constant"
            setattr(cls, method, translator)
            if property_key not in PROPERTIES:
                log.warning(f"Unexpected constant translator for '{property_key}' defined in {cls}")

    def __init__(self, header, filename=None):
        self._header = header
        self.filename = filename
        self._used_cards = set()

        # Cache assumes header is read-only once stored in object
        self._translation_cache = {}

    @classmethod
    @abstractmethod
    def can_translate(cls, header, filename=None):
        """Indicate whether this translation class can translate the
        supplied header.

        Parameters
        ----------
        header : `dict`-like
            Header to convert to standardized form.
        filename : `str`, optional
            Name of file being translated.

        Returns
        -------
        can : `bool`
            `True` if the header is recognized by this class. `False`
            otherwise.
        """
        raise NotImplementedError()

    @classmethod
    def can_translate_with_options(cls, header, options, filename=None):
        """Helper method for `can_translate` allowing options.

        Parameters
        ----------
        header : `dict`-like
            Header to convert to standardized form.
        options : `dict`
            Headers to try to determine whether this header can
            be translated by this class.  If a card is found it will
            be compared with the expected value and will return that
            comparison.  Each card will be tried in turn until one is
            found.
        filename : `str`, optional
            Name of file being translated.

        Returns
        -------
        can : `bool`
            `True` if the header is recognized by this class. `False`
            otherwise.

        Notes
        -----
        Intended to be used from within `can_translate` implementations
        for specific translators.  Is not intended to be called directly
        from `determine_translator`.
        """
        for card, value in options.items():
            if card in header:
                return header[card] == value
        return False

    @classmethod
    def determine_translator(cls, header, filename=None):
        """Determine a translation class by examining the header

        Parameters
        ----------
        header : `dict`-like
            Representation of a header.
        filename : `str`, optional
            Name of file being translated.

        Returns
        -------
        translator : `MetadataTranslator`
            Translation class that knows how to extract metadata from
            the supplied header.

        Raises
        ------
        ValueError
            None of the registered translation classes understood the supplied
            header.
        """
        for name, trans in cls.translators.items():
            if trans.can_translate(header, filename=filename):
                log.debug(f"Using translation class {name}")
                return trans
        else:
            raise ValueError(f"None of the registered translation classes {list(cls.translators.keys())}"
                             " understood this header")

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
        used : `frozenset` of `str`
            Cards used when extracting metadata.
        """
        return frozenset(self._used_cards)

    @staticmethod
    def validate_value(value, default, minimum=None, maximum=None):
        """Validate the supplied value, returning a new value if out of range

        Parameters
        ----------
        value : `float`
            Value to be validated.
        default : `float`
            Default value to use if supplied value is invalid or out of range.
            Assumed to be in the same units as the value expected in the
            header.
        minimum : `float`
            Minimum possible valid value, optional.  If the calculated value
            is below this value, the default value will be used.
        maximum : `float`
            Maximum possible valid value, optional.  If the calculated value
            is above this value, the default value will be used.

        Returns
        -------
        value : `float`
            Either the supplied value, or a default value.
        """
        if value is None or math.isnan(value):
            value = default
        else:
            if minimum is not None and value < minimum:
                value = default
            elif maximum is not None and value > maximum:
                value = default
        return value

    def quantity_from_card(self, keywords, unit, default=None, minimum=None, maximum=None, checker=None):
        """Calculate a Astropy Quantity from a header card and a unit.

        Parameters
        ----------
        keywords : `str` or `list` of `str`
            Keyword to use from header.  If a list each keyword will be tried
            in turn until one matches.
        unit : `astropy.units.UnitBase`
            Unit of the item in the header.
        default : `float`, optional
            Default value to use if the header value is invalid.  Assumed
            to be in the same units as the value expected in the header.  If
            None, no default value is used.
        minimum : `float`, optional
            Minimum possible valid value, optional.  If the calculated value
            is below this value, the default value will be used.
        maximum : `float`, optional
            Maximum possible valid value, optional.  If the calculated value
            is above this value, the default value will be used.
        checker : `function`, optional
            Callback function to be used by the translator method in case the
            keyword is not present.  Function will be executed as if it is
            a method of the translator class.  Running without raising an
            exception will allow the default to be used. Should usually raise
            `KeyError`.

        Returns
        -------
        q : `astropy.units.Quantity`
            Quantity representing the header value.

        Raises
        ------
        KeyError
            The supplied header key is not present.
        """
        keywords = keywords if isinstance(keywords, list) else [keywords]
        for k in keywords:
            if k in self._header:
                value = self._header[k]
                keyword = k
                break
        else:
            if checker is not None:
                try:
                    checker(self)
                    value = default
                    if value is not None:
                        value = u.Quantity(value, unit=unit)
                    return value
                except Exception:
                    pass
            raise KeyError(f"Could not find {keywords} in header")
        if isinstance(value, str):
            # Sometimes the header has the wrong type in it but this must
            # be a number if we are creating a quantity.
            value = float(value)
        self._used_these_cards(keyword)
        if default is not None:
            value = self.validate_value(value, default, maximum=maximum, minimum=minimum)
        return u.Quantity(value, unit=unit)


def _make_abstract_translator_method(property, doc, return_type):
    """Create a an abstract translation method for this property.

    Parameters
    ----------
    property : `str`
        Name of the translator for property to be created.
    doc : `str`
        Description of the property.
    return_type : `str`
        Type of this property (used in the doc string).

    Returns
    -------
    m : `function`
        Translator method for this property.
    """
    def to_property(self):
        raise NotImplementedError(f"Translator for '{property}' undefined.")

    to_property.__doc__ = f"""Return value of {property} from headers.

    {doc}

    Returns
    -------
    {property} : `{return_type}`
        The translated property.
    """
    return to_property


# Make abstract methods for all the translators methods.
# Unfortunately registering them as abstractmethods does not work
# as these assignments come after the class has been created.
# Assigning to __abstractmethods__ directly does work but interacts
# poorly with the metaclass automatically generating methods from
# _trivialMap and _constMap.

for name, description in PROPERTIES.items():
    setattr(MetadataTranslator, f"to_{name}",
            abstractmethod(_make_abstract_translator_method(name, *description)))


class StubTranslator(MetadataTranslator):
    """Translator where all the translations are stubbed out and issue
    warnings.

    This translator can be used as a base class whilst developing a new
    translator.  It allows testing to proceed without being required to fully
    define all translation methods.  Once complete the class should be
    removed from the inheritance tree.

    """
    pass


def _make_stub_translator_method(property, doc, return_type):
    """Create a an stub translation method for this property.

    Parameters
    ----------
    property : `str`
        Name of the translator for property to be created.
    doc : `str`
        Description of the property.
    return_type : `str`
        Type of this property (used in the doc string).

    Returns
    -------
    m : `function`
        Stub translator method for this property.
    """
    def to_stub(self):
        warnings.warn(f"Please implement translator for property '{property}' for translator {self}",
                      stacklevel=3)
        return None

    to_stub.__doc__ = f"""Unimplemented translator for {property}.

    {doc}

    Issues a warning reminding the implementer to override this method.

    Returns
    -------
    {property} : `None`
        Always returns `None`.
    """
    return to_stub


# Create stub translation methods
for name, description in PROPERTIES.items():
    setattr(StubTranslator, f"to_{name}", _make_stub_translator_method(name, *description))
