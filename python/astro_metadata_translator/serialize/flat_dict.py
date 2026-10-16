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
    "detect_prefix_and_sep",
    "find_flat_keys",
    "find_translator_name",
    "flatten_simple",
    "make_flat_key",
    "unflatten",
    "validate_sep",
)

import math
import re
from collections.abc import Collection, Mapping
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

# Longest names first so the regex alternation prefers the full name.
_NAME_ALTERNATION = "|".join(re.escape(n) for n in sorted(_CORE_NAMES, key=len, reverse=True))

# A lower-cased key is an optional prefix and one-character non-word
# separator, then a core name, then an optional ``_N`` index. The greedy
# prefix group and the end anchor let a prefix contain the separator.
_DETECT_RE = re.compile(
    rf"^(?:(?P<prefix>.+)(?P<sep>[^\w]))?(?P<name>{_NAME_ALTERNATION})(?:_(?P<idx>[0-9]+))?$"
)

_INDEX_RE = re.compile(r"^(?P<stem>.+)_(?P<idx>[0-9]+)$")


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


def flatten_simple(
    simple: Mapping[str, Any],
    *,
    prefix: str = "",
    sep: str = ".",
    use_upper: bool | None = None,
    encode_non_finite: bool = False,
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
    encode_non_finite : `bool`, optional
        If `True`, ``NaN`` and infinite floats are written as the strings
        ``NAN``, ``+INFINITY``, and ``-INFINITY``, the spellings cfitsio
        and afw use for such values in FITS headers. By default they are
        passed through as floats.

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
    if use_upper is None:
        # Follow the case of the prefix when the caller has no preference.
        use_upper = prefix[0].isupper() if prefix else False
    flat: dict[str, FlatValue] = {}

    def _add(flat_name: str, value: Any, property_name: str) -> None:
        if not isinstance(value, _FLAT_TYPES):
            raise TypeError(
                f"Value for property {property_name!r} must be str, int, float, or bool "
                f"(or a flat sequence of those). Got {type(value).__name__}."
            )
        if encode_non_finite and isinstance(value, float) and not math.isfinite(value):
            value = _encode_non_finite(value)
        flat[make_flat_key(flat_name, prefix=prefix, sep=sep, use_upper=use_upper)] = value

    for name, value in simple.items():
        flat_name = _TRANSLATOR_FLAT if name == _TRANSLATOR_SIMPLE else name
        if isinstance(value, list | tuple):
            for i, element in enumerate(value):
                _add(f"{flat_name}_{i}", element, name)
        else:
            _add(flat_name, value, name)
    return flat


def _split_index(flat_name: str) -> tuple[str, int | None]:
    """Split a trailing ``_N`` index from a flat name, if present.

    Parameters
    ----------
    flat_name : `str`
        Lower-case flat name such as ``datetime_begin_0``.

    Returns
    -------
    stem : `str`
        The name without the index, or the whole name if there is none.
    index : `int` or `None`
        The index, or `None` if the name has no ``_N`` suffix.
    """
    match = _INDEX_RE.match(flat_name)
    if match is None:
        return flat_name, None
    return match.group("stem"), int(match.group("idx"))


def _is_flat_name(flat_name: str) -> bool:
    """Report whether a lower-case name could have come from flattening.

    Parameters
    ----------
    flat_name : `str`
        Lower-case candidate with any prefix already removed.

    Returns
    -------
    is_flat : `bool`
        `True` for a core property name, ``translator``, an indexed core
        property name, or any name starting with ``ext_``.
    """
    if flat_name.startswith(_EXT_PREFIX):
        return len(flat_name) > len(_EXT_PREFIX)
    if flat_name in _CORE_NAMES:
        return True
    stem, idx = _split_index(flat_name)
    return idx is not None and stem in PROPERTIES


def detect_prefix_and_sep(
    flat: Mapping[str, Any], *, anchor_names: Collection[str] | None = None
) -> tuple[str, str] | None:
    """Work out the prefix and separator used for flat keys in a mapping.

    Parameters
    ----------
    flat : `~collections.abc.Mapping` [`str`, `~typing.Any`]
        Mapping that may contain keys written by `flatten_simple`
        alongside unrelated keys, for example a FITS header.
    anchor_names : `~collections.abc.Collection` [`str`] or `None`, optional
        Property names whose indexed form, such as ``datetime_begin_0``,
        may anchor detection. `None` allows any core property name.
        Callers that know which properties are tuple-valued should pass
        them so that an unrelated key such as ``OBJECT_1`` can not
        anchor.

    Returns
    -------
    detected : `tuple` [`str`, `str`] or `None`
        Lower-cased ``(prefix, sep)``. Both are empty strings when the
        keys have no prefix. `None` if no flat keys were found.

    Raises
    ------
    ValueError
        Raised if two prefixes match equally well.

    Notes
    -----
    This is a heuristic. Core property names such as ``object`` and
    ``temperature`` also occur as ordinary FITS cards, so every key that
    looks like a core name votes for its ``(prefix, sep)`` and the group
    with the most distinct names wins. Only groups that include an
    indexed name such as ``datetime_begin_0`` are considered, because
    every flattened observation has one and ordinary headers do not.
    """
    votes: dict[tuple[str, str], set[str]] = {}
    anchored: set[tuple[str, str]] = set()
    anchors = None if anchor_names is None else {name.lower() for name in anchor_names}
    for key in flat:
        # Non-string keys can never match a flat name.
        match = _DETECT_RE.match(str(key).lower())
        if match is None:
            continue
        group = (match.group("prefix") or "", match.group("sep") or "")
        votes.setdefault(group, set()).add(match.group("name"))
        if match.group("idx") is not None and (anchors is None or match.group("name") in anchors):
            anchored.add(group)

    candidates = {group: names for group, names in votes.items() if group in anchored}
    if not candidates:
        return None
    ranked = sorted(candidates.items(), key=lambda item: len(item[1]), reverse=True)
    if len(ranked) > 1 and len(ranked[0][1]) == len(ranked[1][1]):
        raise ValueError(
            "Can not determine the ObservationInfo key prefix: "
            f"{ranked[0][0]!r} and {ranked[1][0]!r} match equally well."
        )
    return ranked[0][0]


def find_flat_keys(
    flat: Mapping[str, Any],
    *,
    prefix: str | None = None,
    sep: str = ".",
    anchor_names: Collection[str] | None = None,
    unique: bool = True,
) -> dict[str, str]:
    """Find the keys in a mapping that were written by `flatten_simple`.

    Parameters
    ----------
    flat : `~collections.abc.Mapping` [`str`, `~typing.Any`]
        Mapping to scan. Keys are matched case-insensitively.
    prefix : `str` or `None`, optional
        Prefix used when the keys were written. `None` detects the
        prefix and separator with `detect_prefix_and_sep`, in which case
        ``sep`` is ignored. Any other value, including the empty string,
        is used as given.
    sep : `str`, optional
        Separator used when the keys were written. Only used with an
        explicit non-empty ``prefix``.
    anchor_names : `~collections.abc.Collection` [`str`] or `None`, optional
        Passed to `detect_prefix_and_sep` when detecting the prefix.
    unique : `bool`, optional
        If `True` two keys that collapse to the same flat name, for
        example differing only in case, raise `ValueError`. If `False`
        both are returned, which suits callers that remove keys.

    Returns
    -------
    found : `dict` [`str`, `str`]
        Each matching key exactly as it appears in ``flat``, mapped to
        its lower-case flat name with the prefix removed. Empty if
        nothing was found.

    Raises
    ------
    ValueError
        Raised if ``sep`` is invalid, if detection finds a tie, or if
        ``unique`` is `True` and two keys collapse to the same flat name.
    """
    if prefix is None:
        detected = detect_prefix_and_sep(flat, anchor_names=anchor_names)
        if detected is None:
            return {}
        prefix, sep = detected
    elif prefix:
        validate_sep(sep)

    lead = f"{prefix}{sep}".lower() if prefix else ""
    found: dict[str, str] = {}
    seen: dict[str, str] = {}
    for key in flat:
        # Non-string keys can never match a flat name.
        lower = str(key).lower()
        if lead:
            if not lower.startswith(lead):
                continue
            flat_name = lower[len(lead) :]
        else:
            flat_name = lower
        if not _is_flat_name(flat_name):
            continue
        if unique and flat_name in seen:
            raise ValueError(f"Keys {seen[flat_name]!r} and {key!r} both map to flat name {flat_name!r}.")
        seen[flat_name] = key
        found[key] = flat_name
    return found


def find_translator_name(flat: Mapping[str, Any], flat_keys: Mapping[str, str]) -> str | None:
    """Return the translator name recorded in flat keys, if any.

    Parameters
    ----------
    flat : `~collections.abc.Mapping` [`str`, `~typing.Any`]
        The mapping holding the values.
    flat_keys : `~collections.abc.Mapping` [`str`, `str`]
        Keys found by `find_flat_keys`.

    Returns
    -------
    name : `str` or `None`
        The translator name, or `None` if the key is absent.
    """
    for key, flat_name in flat_keys.items():
        if flat_name == _TRANSLATOR_FLAT:
            return str(flat[key])
    return None


def _encode_non_finite(value: float) -> str:
    """Return the FITS string form of a non-finite float.

    Parameters
    ----------
    value : `float`
        A ``NaN`` or infinite value.

    Returns
    -------
    encoded : `str`
        ``NAN``, ``+INFINITY``, or ``-INFINITY``, as written by cfitsio.
    """
    if math.isnan(value):
        return "NAN"
    return "+INFINITY" if value > 0 else "-INFINITY"


def _decode_float(value: Any, name: str) -> Any:
    """Convert a string holding a float back to a float.

    Parameters
    ----------
    value : `~typing.Any`
        Value read for a float property. Anything other than a `str` is
        returned unchanged.
    name : `str`
        Property name, used in the error message.

    Returns
    -------
    decoded : `~typing.Any`
        The float parsed from the string, or the original value.

    Raises
    ------
    ValueError
        Raised if the string does not parse as a float.
    """
    if not isinstance(value, str):
        return value
    try:
        return float(value)
    except ValueError:
        raise ValueError(f"Value {value!r} for float property {name!r} is not a number.") from None


def unflatten(
    flat: Mapping[str, Any],
    flat_keys: Mapping[str, str],
    *,
    extension_names: Collection[str] = (),
    expected_lengths: Mapping[str, int] | None = None,
    float_names: Collection[str] = (),
) -> dict[str, Any]:
    """Rebuild the simple form from flat keys.

    Parameters
    ----------
    flat : `~collections.abc.Mapping` [`str`, `~typing.Any`]
        The mapping holding the values.
    flat_keys : `~collections.abc.Mapping` [`str`, `str`]
        Keys found by `find_flat_keys`, mapping each original key to
        its lower-case flat name.
    extension_names : `~collections.abc.Collection` [`str`], optional
        Extension property names (without the ``ext_`` prefix) defined
        by the relevant translator. A trailing ``_N`` on an ``ext_`` key
        is treated as an index only when the stem is one of these, and
        the original spelling is restored in the result.
    expected_lengths : `~collections.abc.Mapping` [`str`, `int`], optional
        Required tuple length for each indexed name. A name listed here
        must have exactly that many elements. Names not listed only need
        a contiguous run of indexes starting at 0.
    float_names : `~collections.abc.Collection` [`str`], optional
        Simple-form names of properties whose values, or tuple elements,
        are floats. A `str` found for one of these is converted with
        `float`, which accepts the ``NAN``, ``+INFINITY``, and
        ``-INFINITY`` strings written by `flatten_simple`, cfitsio, and
        afw, as well as the ``NaN`` and ``Infinity`` spellings used in
        JSON.

    Returns
    -------
    simple : `dict` [`str`, `~typing.Any`]
        Simple form suitable for `ObservationInfo.from_simple`. Indexed
        entries are gathered into tuples and ``translator`` is renamed
        back to ``_translator``.

    Raises
    ------
    ValueError
        Raised if the indexes for a name do not form a contiguous run
        starting at 0, or do not match the expected length. Also raised
        if a name appears both as a single value and as indexed values,
        or if a string for a name in ``float_names`` is not a number.
    """
    ext_by_lower = {f"{_EXT_PREFIX}{name.lower()}": f"{_EXT_PREFIX}{name}" for name in extension_names}

    scalars: dict[str, Any] = {}
    indexed: dict[str, dict[int, Any]] = {}
    for key, flat_name in flat_keys.items():
        value = flat[key]
        stem, idx = _split_index(flat_name)
        if idx is not None and (stem in PROPERTIES or stem in ext_by_lower):
            indexed.setdefault(stem, {})[idx] = value
        else:
            scalars[flat_name] = value

    if conflicts := sorted(scalars.keys() & indexed.keys()):
        raise ValueError(f"Properties {conflicts} appear both as a single value and as indexed values.")

    for stem, parts in indexed.items():
        n_expected = len(parts) if expected_lengths is None else expected_lengths.get(stem, len(parts))
        expected = list(range(n_expected))
        if sorted(parts) != expected:
            raise ValueError(
                f"Indexed values for {stem!r} must have indexes {expected}. Got indexes {sorted(parts)}."
            )
        scalars[stem] = tuple(parts[i] for i in expected)

    simple: dict[str, Any] = {}
    for flat_name, value in scalars.items():
        if flat_name == _TRANSLATOR_FLAT:
            simple[_TRANSLATOR_SIMPLE] = value
        else:
            simple[ext_by_lower.get(flat_name, flat_name)] = value

    for name in float_names:
        if name not in simple:
            continue
        value = simple[name]
        if isinstance(value, tuple):
            simple[name] = tuple(_decode_float(element, name) for element in value)
        else:
            simple[name] = _decode_float(value, name)
    return simple
