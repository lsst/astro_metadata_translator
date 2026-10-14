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


def _split_index(flat_name: str) -> tuple[str, int | None]:
    """Split a trailing ``_N`` index from a flat name, if present."""
    match = _INDEX_RE.match(flat_name)
    if match is None:
        return flat_name, None
    return match.group("stem"), int(match.group("idx"))


def _is_flat_name(flat_name: str) -> bool:
    """Report whether a lower-case name could have come from flattening."""
    if flat_name.startswith(_EXT_PREFIX):
        return len(flat_name) > len(_EXT_PREFIX)
    if flat_name in _CORE_NAMES:
        return True
    stem, idx = _split_index(flat_name)
    return idx is not None and stem in PROPERTIES


def detect_prefix_and_sep(flat: Mapping[str, Any]) -> tuple[str, str] | None:
    """Work out the prefix and separator used for flat keys in a mapping.

    Parameters
    ----------
    flat : `~collections.abc.Mapping` [`str`, `~typing.Any`]
        Mapping that may contain keys written by `flatten_simple`
        alongside unrelated keys, for example a FITS header.

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
    for key in flat:
        # Non-string keys can never match a flat name.
        match = _DETECT_RE.match(str(key).lower())
        if match is None:
            continue
        group = (match.group("prefix") or "", match.group("sep") or "")
        votes.setdefault(group, set()).add(match.group("name"))
        if match.group("idx") is not None:
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


def find_flat_keys(flat: Mapping[str, Any], *, prefix: str | None = None, sep: str = ".") -> dict[str, str]:
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

    Returns
    -------
    found : `dict` [`str`, `str`]
        Each matching key exactly as it appears in ``flat``, mapped to
        its lower-case flat name with the prefix removed. Empty if
        nothing was found.

    Raises
    ------
    ValueError
        Raised if ``sep`` is invalid, if detection finds a tie, or if two
        keys collapse to the same flat name.
    """
    if prefix is None:
        detected = detect_prefix_and_sep(flat)
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
        if flat_name in seen:
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


def unflatten(
    flat: Mapping[str, Any],
    flat_keys: Mapping[str, str],
    *,
    extension_names: Collection[str] = (),
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
        starting at 0.
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

    for stem, parts in indexed.items():
        expected = list(range(len(parts)))
        if sorted(parts) != expected:
            raise ValueError(
                f"Indexed values for {stem!r} are not contiguous from 0. Got indexes {sorted(parts)}."
            )
        scalars[stem] = tuple(parts[i] for i in expected)

    simple: dict[str, Any] = {}
    for flat_name, value in scalars.items():
        if flat_name == _TRANSLATOR_FLAT:
            simple[_TRANSLATOR_SIMPLE] = value
        else:
            simple[ext_by_lower.get(flat_name, flat_name)] = value
    return simple
