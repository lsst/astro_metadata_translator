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

"""Flatten the simple form of an ObservationInfo into scalar key/value pairs.

The simple form produced by `ObservationInfo.to_simple` contains tuples
for some properties. A FITS header can only hold scalars, so each
element of a tuple is written to its own key with an ``_N`` index
suffix. Keys can carry an opaque prefix joined by a separator and can
be upper-cased to suit FITS conventions. The reverse direction finds
the keys again case-insensitively, with or without knowing the prefix.
"""

from __future__ import annotations

__all__ = (
    "FITS_FLAT_DICT_PREFIX",
    "FITS_FLAT_DICT_SEP",
    "FlatValue",
    "flatten_simple",
    "make_flat_key",
    "validate_sep",
)

import re
from collections.abc import Mapping
from typing import Any

from ..properties import PROPERTIES

FITS_FLAT_DICT_PREFIX = "OBS_INFO"
"""Recommended key prefix when writing to a FITS header."""

FITS_FLAT_DICT_SEP = " "
"""Recommended key separator when writing to a FITS header."""

FlatValue = str | int | float | bool
"""Scalar value types permitted in a flat dict."""

_FLAT_TYPES = (str, int, float, bool)

_TRANSLATOR_FLAT = "translator"
"""Flat name used for the ``_translator`` entry of the simple form."""

_TRANSLATOR_SIMPLE = "_translator"
_EXT_PREFIX = "ext_"
_CORE_NAMES = frozenset(PROPERTIES) | {_TRANSLATOR_FLAT}


def validate_sep(sep: str) -> None:
    """Check that a separator is usable in flat keys.

    Parameters
    ----------
    sep : `str`
        Separator to check.

    Raises
    ------
    ValueError
        Raised if the separator is not exactly one character or is a
        letter, digit, or underscore. Those characters appear inside
        property names and would make keys ambiguous.
    """
    if len(sep) != 1:
        raise ValueError(f"Separator for flat keys must be a single character. Got {sep!r}.")
    if re.match(r"\w$", sep):
        raise ValueError(f"Separator for flat keys can not be a letter, digit, or underscore. Got {sep!r}.")


def make_flat_key(flat_name: str, *, prefix: str, sep: str, use_upper: bool) -> str:
    """Build the full key for a flat name.

    Parameters
    ----------
    flat_name : `str`
        Flat property name, for example ``datetime_begin_0``.
    prefix : `str`
        Prefix emitted exactly as given. Empty for no prefix.
    sep : `str`
        Separator placed between the prefix and the flat name. Unused
        when the prefix is empty.
    use_upper : `bool`
        If `True` the flat name is upper-cased.

    Returns
    -------
    key : `str`
        The full key.
    """
    name = flat_name.upper() if use_upper else flat_name
    if prefix:
        return f"{prefix}{sep}{name}"
    return name


def _resolve_use_upper(prefix: str, use_upper: bool | None) -> bool:
    """Decide the flat name case, following the prefix when unspecified."""
    if use_upper is not None:
        return use_upper
    return prefix[0].isupper() if prefix else False


def flatten_simple(
    simple: Mapping[str, Any],
    *,
    prefix: str = "",
    sep: str = ".",
    use_upper: bool | None = None,
) -> dict[str, FlatValue]:
    """Flatten a simple-form mapping into scalar key/value pairs.

    Parameters
    ----------
    simple : `~collections.abc.Mapping` [`str`, `~typing.Any`]
        Mapping as returned by `ObservationInfo.to_simple`.
    prefix : `str`, optional
        Prefix for every key. Emitted exactly as given.
    sep : `str`, optional
        Separator between prefix and flat name. Must be one character
        that is not a letter, digit, or underscore.
    use_upper : `bool` or `None`, optional
        If `None` the flat names are upper-cased when the first
        character of ``prefix`` is upper case, otherwise left lower
        case. `True` or `False` force the case.

    Returns
    -------
    flat : `dict` [`str`, `str` | `int` | `float` | `bool`]
        Flat mapping. Each element of a list or tuple value becomes its
        own entry named ``{name}_{index}`` counting from 0. The
        ``_translator`` entry is renamed ``translator``.

    Raises
    ------
    ValueError
        Raised if ``sep`` is not valid.
    TypeError
        Raised if a value is not a `str`, `int`, `float`, `bool`, or a
        one-level list or tuple of those.
    """
    validate_sep(sep)
    upper = _resolve_use_upper(prefix, use_upper)
    flat: dict[str, FlatValue] = {}

    def _add(flat_name: str, value: Any, property_name: str) -> None:
        if not isinstance(value, _FLAT_TYPES):
            raise TypeError(
                f"Value for property {property_name!r} must be str, int, float, or bool "
                f"(or a flat sequence of those). Got {type(value).__name__}."
            )
        flat[make_flat_key(flat_name, prefix=prefix, sep=sep, use_upper=upper)] = value

    for name, value in simple.items():
        flat_name = _TRANSLATOR_FLAT if name == _TRANSLATOR_SIMPLE else name
        if isinstance(value, list | tuple):
            for i, element in enumerate(value):
                _add(f"{flat_name}_{i}", element, name)
        else:
            _add(flat_name, value, name)
    return flat
