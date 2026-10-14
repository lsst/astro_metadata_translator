# ObservationInfo Flat Dict Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `to_flat_dict`, `from_flat_dict`, `has_flat_dict`, and `strip_from_flat_dict` to `ObservationInfo` so an observation can be stored as scalar key/value pairs in a FITS header and rebuilt from it losslessly.

**Architecture:** A new pure-function module `serialize/flat_dict.py` flattens the existing `to_simple()` dict (tuples become `name_0`, `name_1`, ...; `_translator` becomes `translator`; an optional opaque prefix joined by a separator; optional upper-casing) and reverses it with case-insensitive, majority-vote prefix detection anchored on an indexed key. `ObservationInfo` gains four thin methods that delegate to the helper and to the existing `to_simple`/`from_simple`, so all property type knowledge stays in `properties.py`.

**Tech Stack:** Python 3.11+, pydantic (already used by `ObservationInfo`), astropy (`astropy.io.fits` for the integration test), unittest, ruff, mypy.

**Spec:** `docs/superpowers/specs/2026-10-14-observation-info-flat-dict-design.md`

## Global Constraints

- Only astro_metadata_translator changes on this ticket. obs_base and lsst.images are follow-up tickets.
- Flat value types are exactly `str`, `int`, `float`, `bool`. Anything else raises `TypeError`.
- `sep` must be exactly one character and must not be a letter, digit, or underscore. Violations raise `ValueError`.
- Index suffix is always `_<digits>` joined with an underscore, in every prefix form.
- `_translator` in simple form is `translator` in flat form.
- Constants: `FITS_FLAT_DICT_PREFIX = "OBS_INFO"`, `FITS_FLAT_DICT_SEP = " "`, re-exported from the top-level `astro_metadata_translator` package.
- Auto-detection (`prefix=None`) only accepts a candidate prefix group that contains at least one indexed match such as `datetime_begin_0`.
- Style: ruff clean (`line-length = 110`, `max-doc-length = 79`, numpydoc docstrings on every public function), mypy clean under the repo's strict settings for `astro_metadata_translator.*` (`disallow_untyped_defs`).
- Docstrings and comments describe the code as it is; never mention plans, history, or alternatives.
- Commit after each task. Never push.
- Run tests with `python -m pytest tests/test_flat_dict.py -q` from the repo root; lint with `ruff check python tests && ruff format --check python tests && mypy python`.

## Deviations during execution

The code as merged differs from the task text below in these ways; the spec reflects all of them.

- Tests run through the EUPS runner (`lsst-run -l . -- pytest ...`) because the stack's installed copy of the package shadows the checkout on PYTHONPATH.
- `detect_prefix_and_sep` and `find_flat_keys` normalize keys with `str(key).lower()` instead of an `isinstance` guard, which mypy rejected as unreachable.
- `unflatten` takes `expected_lengths` and `from_flat_dict` passes lengths derived from the model JSON schema, so a truncated tuple is rejected.
- `detect_prefix_and_sep` and `find_flat_keys` take `anchor_names`; the ObservationInfo methods pass the tuple-valued property names so an unrelated `OBJECT_1` card can not anchor detection.
- `unflatten` rejects a property present both as a single value and as indexed values.
- `find_flat_keys` takes `unique`; `has_flat_dict` and `strip_from_flat_dict` pass `False` so keys differing only in case are counted and removed rather than raising.
- `serialize/__init__.py` exports only the two constants from `flat_dict`.
- The real-FITS test compares floats to 15 significant digits because astropy truncates long floats to a 20-character field.
- Two private helpers in `flat_dict.py` carry full numpydoc Parameters and Returns sections to satisfy the pre-commit hook.
- `flatten_simple` resolves `use_upper` inline instead of through a `_resolve_use_upper` helper.

## Review Focus

1. A plain FITS header with `OBJECT` and `TEMPERATURE` cards but no ObservationInfo content: `has_flat_dict(hdr)` must be `False` and `strip_from_flat_dict(hdr)` must delete nothing. Pinned in Task 2 (`test_detect_decoys_without_anchor`) and Task 5 (`test_fits_header_without_content`).
2. A header containing both the prefixed set and an unprefixed decoy `OBJECT`: the rebuilt `object` must come from `OBS_INFO OBJECT`, not from the decoy. Pinned in Task 5 (`test_fits_header_round_trip`).
3. A mapping with non-string keys (an `int` key in a plain dict): detection and collection must skip it rather than crash on `.lower()`. Pinned in Task 2 (`test_find_keys_ignores_non_string_keys`).
4. Explicit prefix with the wrong separator (`prefix="OBS_INFO", sep="."` against space-separated keys): `has_flat_dict` is `False`, `from_flat_dict` raises `ValueError`. Pinned in Task 4 (`test_explicit_prefix_wrong_sep`).
5. `NaN` float values such as `relative_humidity`: `to_flat_dict` passes them through unchanged and a plain dict round trip preserves them, but astropy refuses `NaN` header cards, so FITS writers must drop them first. Pinned in Task 5 (`test_nan_passes_through_plain_dict`) and documented in the `to_flat_dict` docstring.

---

### Task 1: Helper module, key building and flattening

**Files:**
- Create: `python/astro_metadata_translator/serialize/flat_dict.py`
- Modify: `python/astro_metadata_translator/serialize/__init__.py` (add `from .flat_dict import *`)
- Create: `tests/test_flat_dict.py`

**Interfaces:**
- Consumes: `astro_metadata_translator.properties.PROPERTIES` (dict of core property names).
- Produces:
  - `FITS_FLAT_DICT_PREFIX: str = "OBS_INFO"`, `FITS_FLAT_DICT_SEP: str = " "`
  - `FlatValue = str | int | float | bool`
  - `validate_sep(sep: str) -> None`
  - `make_flat_key(flat_name: str, *, prefix: str, sep: str, use_upper: bool) -> str`
  - `flatten_simple(simple: Mapping[str, Any], *, prefix: str = "", sep: str = ".", use_upper: bool | None = None) -> dict[str, FlatValue]`

- [ ] **Step 1: Create the test file with failing tests for flattening**

Create `tests/test_flat_dict.py`. Copy the ten-line license comment block from the top of `tests/test_serialization.py` verbatim, then:

```python
"""Tests for the flat dict serialization of ObservationInfo."""

import unittest

from astro_metadata_translator.serialize.flat_dict import (
    FITS_FLAT_DICT_PREFIX,
    FITS_FLAT_DICT_SEP,
    flatten_simple,
    make_flat_key,
    validate_sep,
)


class FlattenTestCase(unittest.TestCase):
    """Test flattening of the simple form."""

    def test_constants(self) -> None:
        self.assertEqual(FITS_FLAT_DICT_PREFIX, "OBS_INFO")
        self.assertEqual(FITS_FLAT_DICT_SEP, " ")

    def test_validate_sep(self) -> None:
        for good in (".", " ", "-", ":"):
            validate_sep(good)
        for bad in ("", "..", "a", "Z", "1", "_"):
            with self.assertRaises(ValueError):
                validate_sep(bad)

    def test_make_flat_key(self) -> None:
        self.assertEqual(make_flat_key("instrument", prefix="", sep=".", use_upper=False), "instrument")
        self.assertEqual(make_flat_key("instrument", prefix="", sep=".", use_upper=True), "INSTRUMENT")
        self.assertEqual(
            make_flat_key("datetime_begin_0", prefix="obs_info", sep=".", use_upper=False),
            "obs_info.datetime_begin_0",
        )
        self.assertEqual(
            make_flat_key("datetime_begin_0", prefix="OBS_INFO", sep=" ", use_upper=True),
            "OBS_INFO DATETIME_BEGIN_0",
        )
        # The prefix is emitted exactly as given.
        self.assertEqual(
            make_flat_key("instrument", prefix="obs_info", sep=" ", use_upper=True),
            "obs_info INSTRUMENT",
        )

    def test_flatten_scalars_no_prefix(self) -> None:
        simple = {
            "instrument": "HSC",
            "exposure_id": 42,
            "boresight_airmass": 1.2,
            "can_see_sky": True,
            "_translator": "HSC",
        }
        self.assertEqual(
            flatten_simple(simple),
            {
                "instrument": "HSC",
                "exposure_id": 42,
                "boresight_airmass": 1.2,
                "can_see_sky": True,
                "translator": "HSC",
            },
        )

    def test_flatten_sequences(self) -> None:
        simple = {"datetime_begin": (2460000.5, 0.25), "location": [1.0, 2.0, 3.0], "ext_pair": (1, 2)}
        self.assertEqual(
            flatten_simple(simple),
            {
                "datetime_begin_0": 2460000.5,
                "datetime_begin_1": 0.25,
                "location_0": 1.0,
                "location_1": 2.0,
                "location_2": 3.0,
                "ext_pair_0": 1,
                "ext_pair_1": 2,
            },
        )

    def test_flatten_prefix_and_case(self) -> None:
        simple = {"datetime_begin": (1.0, 2.0), "instrument": "X"}
        self.assertEqual(
            flatten_simple(simple, prefix="obs_info"),
            {
                "obs_info.datetime_begin_0": 1.0,
                "obs_info.datetime_begin_1": 2.0,
                "obs_info.instrument": "X",
            },
        )
        self.assertEqual(
            flatten_simple(simple, prefix="OBS_INFO", sep=" "),
            {
                "OBS_INFO DATETIME_BEGIN_0": 1.0,
                "OBS_INFO DATETIME_BEGIN_1": 2.0,
                "OBS_INFO INSTRUMENT": "X",
            },
        )
        self.assertEqual(
            flatten_simple(simple, prefix="OBS_INFO", sep=" ", use_upper=False),
            {
                "OBS_INFO datetime_begin_0": 1.0,
                "OBS_INFO datetime_begin_1": 2.0,
                "OBS_INFO instrument": "X",
            },
        )
        self.assertEqual(
            flatten_simple(simple, use_upper=True),
            {"DATETIME_BEGIN_0": 1.0, "DATETIME_BEGIN_1": 2.0, "INSTRUMENT": "X"},
        )

    def test_flatten_rejects_bad_values(self) -> None:
        with self.assertRaises(TypeError):
            flatten_simple({"ext_bad": {"a": 1}})
        with self.assertRaises(TypeError):
            flatten_simple({"ext_nested": [[1.0, 2.0]]})
        with self.assertRaises(TypeError):
            flatten_simple({"ext_none": None})

    def test_flatten_rejects_bad_sep(self) -> None:
        with self.assertRaises(ValueError):
            flatten_simple({"instrument": "X"}, prefix="p", sep="_")
        with self.assertRaises(ValueError):
            flatten_simple({"instrument": "X"}, sep="ab")


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_flat_dict.py -q`
Expected: FAIL at import with `ModuleNotFoundError: No module named 'astro_metadata_translator.serialize.flat_dict'`.

- [ ] **Step 3: Create the helper module**

Create `python/astro_metadata_translator/serialize/flat_dict.py`. Copy the ten-line license comment block from the top of `python/astro_metadata_translator/serialize/fits.py` verbatim, then:

```python
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
        Opaque prefix emitted exactly as given. Empty for no prefix.
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
        Opaque prefix for every key. Emitted exactly as given.
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
```

Then edit `python/astro_metadata_translator/serialize/__init__.py` so the body reads:

```python
from .fits import *
from .flat_dict import *
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_flat_dict.py -q`
Expected: 8 passed.

- [ ] **Step 5: Lint and type-check**

Run: `ruff check python tests && ruff format --check python tests && mypy python`
Expected: no errors. If `ruff format --check` complains, run `ruff format python tests` and re-run the tests.

- [ ] **Step 6: Commit**

```bash
git add python/astro_metadata_translator/serialize/flat_dict.py python/astro_metadata_translator/serialize/__init__.py tests/test_flat_dict.py
git commit -m "Add flat dict key building and flattening helpers.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 2: Helper module, prefix detection and key collection

**Files:**
- Modify: `python/astro_metadata_translator/serialize/flat_dict.py`
- Modify: `tests/test_flat_dict.py`

**Interfaces:**
- Consumes: `validate_sep`, `_CORE_NAMES`, `_EXT_PREFIX` from Task 1.
- Produces:
  - `detect_prefix_and_sep(flat: Mapping[str, Any]) -> tuple[str, str] | None` returning lower-cased `(prefix, sep)`, `("", "")` for no prefix, or `None` when nothing is found. Raises `ValueError` on a tie.
  - `find_flat_keys(flat: Mapping[str, Any], *, prefix: str | None = None, sep: str = ".") -> dict[str, str]` mapping each original key to its lower-case flat name with the prefix removed. Empty when nothing is found.
  - `_split_index(flat_name: str) -> tuple[str, int | None]` (module private, reused by Task 3).

- [ ] **Step 1: Add failing detection and collection tests**

Add to the import block in `tests/test_flat_dict.py`:

```python
from astro_metadata_translator.serialize.flat_dict import (
    FITS_FLAT_DICT_PREFIX,
    FITS_FLAT_DICT_SEP,
    detect_prefix_and_sep,
    find_flat_keys,
    flatten_simple,
    make_flat_key,
    validate_sep,
)
```

Append this class before the `if __name__ == "__main__":` line:

```python
class DetectTestCase(unittest.TestCase):
    """Test prefix detection and key collection."""

    PREFIXED = {
        "OBS_INFO DATETIME_BEGIN_0": 2460000.5,
        "OBS_INFO DATETIME_BEGIN_1": 0.25,
        "OBS_INFO INSTRUMENT": "HSC",
        "OBS_INFO OBJECT": "real",
        "OBS_INFO TRANSLATOR": "HSC",
        "OBS_INFO EXT_FOO": "bar",
        "OBS_INFO EXT_PAIR_0": 1,
        "OBJECT": "decoy",
        "TEMPERATURE": 3.0,
        "EXPTIME": 30.0,
    }

    def test_detect_nothing(self) -> None:
        self.assertIsNone(detect_prefix_and_sep({}))
        self.assertIsNone(detect_prefix_and_sep({"EXPTIME": 30.0, "FILTER": "g"}))

    def test_detect_decoys_without_anchor(self) -> None:
        # Core names appear but nothing indexed, so this is not content.
        self.assertIsNone(detect_prefix_and_sep({"OBJECT": "M31", "TEMPERATURE": 3.0}))
        self.assertIsNone(detect_prefix_and_sep({"instrument": "X", "object": "Y"}))

    def test_detect_plain(self) -> None:
        flat = {"datetime_begin_0": 1.0, "datetime_begin_1": 2.0, "instrument": "X"}
        self.assertEqual(detect_prefix_and_sep(flat), ("", ""))

    def test_detect_prefixed_with_decoys(self) -> None:
        self.assertEqual(detect_prefix_and_sep(self.PREFIXED), ("obs_info", " "))

    def test_detect_prefix_containing_sep(self) -> None:
        flat = {"LSST OBS_INFO DATETIME_BEGIN_0": 1.0, "LSST OBS_INFO EXPOSURE_ID": 5}
        self.assertEqual(detect_prefix_and_sep(flat), ("lsst obs_info", " "))

    def test_detect_dotted(self) -> None:
        flat = {"lsst.obs_info.datetime_begin_0": 1.0, "lsst.obs_info.exposure_id": 5}
        self.assertEqual(detect_prefix_and_sep(flat), ("lsst.obs_info", "."))

    def test_detect_majority(self) -> None:
        flat = {"a.datetime_begin_0": 1.0, "a.instrument": "x", "b.datetime_begin_0": 1.0}
        self.assertEqual(detect_prefix_and_sep(flat), ("a", "."))

    def test_detect_tie(self) -> None:
        flat = {"a.datetime_begin_0": 1.0, "b.datetime_begin_0": 1.0}
        with self.assertRaises(ValueError):
            detect_prefix_and_sep(flat)

    def test_find_keys_auto(self) -> None:
        found = find_flat_keys(self.PREFIXED)
        self.assertEqual(
            found,
            {
                "OBS_INFO DATETIME_BEGIN_0": "datetime_begin_0",
                "OBS_INFO DATETIME_BEGIN_1": "datetime_begin_1",
                "OBS_INFO INSTRUMENT": "instrument",
                "OBS_INFO OBJECT": "object",
                "OBS_INFO TRANSLATOR": "translator",
                "OBS_INFO EXT_FOO": "ext_foo",
                "OBS_INFO EXT_PAIR_0": "ext_pair_0",
            },
        )

    def test_find_keys_auto_nothing(self) -> None:
        self.assertEqual(find_flat_keys({"OBJECT": "M31"}), {})

    def test_find_keys_explicit(self) -> None:
        found = find_flat_keys(self.PREFIXED, prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        self.assertEqual(len(found), 7)
        self.assertNotIn("OBJECT", found)
        # Explicit empty prefix needs no indexed anchor.
        self.assertEqual(
            find_flat_keys({"instrument": "X", "EXPTIME": 1}, prefix=""), {"instrument": "instrument"}
        )
        # Explicit prefix with the wrong separator finds nothing.
        self.assertEqual(find_flat_keys(self.PREFIXED, prefix="OBS_INFO", sep="."), {})
        # Explicit prefix is matched case-insensitively.
        found = find_flat_keys(self.PREFIXED, prefix="obs_info", sep=" ")
        self.assertEqual(len(found), 7)

    def test_find_keys_ignores_non_string_keys(self) -> None:
        flat: dict[Any, Any] = {1: "x", "datetime_begin_0": 1.0, "datetime_begin_1": 2.0}
        expected = {"datetime_begin_0": "datetime_begin_0", "datetime_begin_1": "datetime_begin_1"}
        self.assertEqual(find_flat_keys(flat), expected)
        self.assertEqual(find_flat_keys(flat, prefix=""), expected)

    def test_find_keys_ignores_unknown_names(self) -> None:
        flat = {"obs_info.datetime_begin_0": 1.0, "obs_info.not_a_property": 2, "obs_info.ext_": 3}
        self.assertEqual(
            find_flat_keys(flat, prefix="obs_info"), {"obs_info.datetime_begin_0": "datetime_begin_0"}
        )

    def test_find_keys_duplicate(self) -> None:
        flat = {"OBS_INFO INSTRUMENT": "a", "obs_info instrument": "b", "OBS_INFO DATETIME_BEGIN_0": 1.0}
        with self.assertRaises(ValueError):
            find_flat_keys(flat, prefix="OBS_INFO", sep=" ")

    def test_find_keys_bad_sep(self) -> None:
        with self.assertRaises(ValueError):
            find_flat_keys({"p_instrument": "x"}, prefix="p", sep="_")
```

Also add `from typing import Any` to the imports at the top of the test file.

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_flat_dict.py -q`
Expected: FAIL at import with `ImportError: cannot import name 'detect_prefix_and_sep'`.

- [ ] **Step 3: Implement detection and collection**

In `python/astro_metadata_translator/serialize/flat_dict.py`, extend `__all__`:

```python
__all__ = (
    "FITS_FLAT_DICT_PREFIX",
    "FITS_FLAT_DICT_SEP",
    "FlatValue",
    "detect_prefix_and_sep",
    "find_flat_keys",
    "flatten_simple",
    "make_flat_key",
    "validate_sep",
)
```

After the `_CORE_NAMES` definition add:

```python
# Longest names first so the regex alternation prefers the full name.
_NAME_ALTERNATION = "|".join(re.escape(n) for n in sorted(_CORE_NAMES, key=len, reverse=True))

# A lower-cased key is an optional prefix and one-character non-word
# separator, then a core name, then an optional ``_N`` index. The greedy
# prefix group and the end anchor let a prefix contain the separator.
_DETECT_RE = re.compile(
    rf"^(?:(?P<prefix>.+)(?P<sep>[^\w]))?(?P<name>{_NAME_ALTERNATION})(?:_(?P<idx>[0-9]+))?$"
)

_INDEX_RE = re.compile(r"^(?P<stem>.+)_(?P<idx>[0-9]+)$")
```

Append these functions at the end of the module:

```python
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
        if not isinstance(key, str):
            continue
        match = _DETECT_RE.match(key.lower())
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


def find_flat_keys(
    flat: Mapping[str, Any], *, prefix: str | None = None, sep: str = "."
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
        if not isinstance(key, str):
            continue
        lower = key.lower()
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
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_flat_dict.py -q`
Expected: 23 passed.

- [ ] **Step 5: Lint and type-check**

Run: `ruff check python tests && ruff format --check python tests && mypy python`
Expected: no errors. If `ruff format --check` complains, run `ruff format python tests` and re-run the tests.

- [ ] **Step 6: Commit**

```bash
git add python/astro_metadata_translator/serialize/flat_dict.py tests/test_flat_dict.py
git commit -m "Add flat dict prefix detection and key collection.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 3: Helper module, unflattening

**Files:**
- Modify: `python/astro_metadata_translator/serialize/flat_dict.py`
- Modify: `tests/test_flat_dict.py`

**Interfaces:**
- Consumes: `find_flat_keys`, `_split_index`, `_TRANSLATOR_FLAT`, `_TRANSLATOR_SIMPLE`, `_EXT_PREFIX` from Tasks 1 and 2.
- Produces:
  - `find_translator_name(flat: Mapping[str, Any], flat_keys: Mapping[str, str]) -> str | None`
  - `unflatten(flat: Mapping[str, Any], flat_keys: Mapping[str, str], *, extension_names: Collection[str] = ()) -> dict[str, Any]` returning a simple-form dict ready for `ObservationInfo.from_simple`.

- [ ] **Step 1: Add failing unflatten tests**

Extend the import block in `tests/test_flat_dict.py` with `find_translator_name` and `unflatten` (keep the list alphabetical). Append this class before the `if __name__ == "__main__":` line:

```python
class UnflattenTestCase(unittest.TestCase):
    """Test rebuilding the simple form from flat keys."""

    def test_unflatten_core(self) -> None:
        flat = {
            "OBS_INFO DATETIME_BEGIN_1": 0.25,
            "OBS_INFO DATETIME_BEGIN_0": 2460000.5,
            "OBS_INFO LOCATION_0": 1.0,
            "OBS_INFO LOCATION_1": 2.0,
            "OBS_INFO LOCATION_2": 3.0,
            "OBS_INFO INSTRUMENT": "HSC",
            "OBS_INFO TRANSLATOR": "HSC",
            "OBJECT": "decoy",
        }
        flat_keys = find_flat_keys(flat)
        self.assertEqual(find_translator_name(flat, flat_keys), "HSC")
        self.assertEqual(
            unflatten(flat, flat_keys),
            {
                "datetime_begin": (2460000.5, 0.25),
                "location": (1.0, 2.0, 3.0),
                "instrument": "HSC",
                "_translator": "HSC",
            },
        )

    def test_unflatten_no_translator(self) -> None:
        flat = {"datetime_begin_0": 1.0, "datetime_begin_1": 2.0}
        flat_keys = find_flat_keys(flat)
        self.assertIsNone(find_translator_name(flat, flat_keys))
        self.assertEqual(unflatten(flat, flat_keys), {"datetime_begin": (1.0, 2.0)})

    def test_unflatten_gap(self) -> None:
        flat = {"datetime_begin_0": 1.0, "datetime_begin_2": 2.0}
        with self.assertRaises(ValueError):
            unflatten(flat, find_flat_keys(flat))
        flat = {"datetime_begin_1": 2.0}
        with self.assertRaises(ValueError):
            unflatten(flat, find_flat_keys(flat))

    def test_unflatten_extensions(self) -> None:
        flat = {
            "obs_info.datetime_begin_0": 1.0,
            "obs_info.datetime_begin_1": 2.0,
            "obs_info.ext_foo": "bar",
            "obs_info.ext_pair_0": 1,
            "obs_info.ext_pair_1": 2,
            "obs_info.ext_other_1": 3,
        }
        flat_keys = find_flat_keys(flat, prefix="obs_info")
        # Known extension names restore their spelling and regroup.
        simple = unflatten(flat, flat_keys, extension_names=("Foo", "pair"))
        self.assertEqual(
            simple,
            {
                "datetime_begin": (1.0, 2.0),
                "ext_Foo": "bar",
                "ext_pair": (1, 2),
                "ext_other_1": 3,
            },
        )
        # Unknown extension names pass through as plain lower-case names.
        simple = unflatten(flat, flat_keys)
        self.assertEqual(
            simple,
            {
                "datetime_begin": (1.0, 2.0),
                "ext_foo": "bar",
                "ext_pair_0": 1,
                "ext_pair_1": 2,
                "ext_other_1": 3,
            },
        )
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_flat_dict.py -q`
Expected: FAIL at import with `ImportError: cannot import name 'find_translator_name'`.

- [ ] **Step 3: Implement unflattening**

In `python/astro_metadata_translator/serialize/flat_dict.py`, extend `__all__`:

```python
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
```

Change the collections import to `from collections.abc import Collection, Mapping`. Append at the end of the module:

```python
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
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_flat_dict.py -q`
Expected: 27 passed.

- [ ] **Step 5: Lint and type-check**

Run: `ruff check python tests && ruff format --check python tests && mypy python`
Expected: no errors.

- [ ] **Step 6: Commit**

```bash
git add python/astro_metadata_translator/serialize/flat_dict.py tests/test_flat_dict.py
git commit -m "Add flat dict unflattening helper.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 4: ObservationInfo methods and round trips

**Files:**
- Modify: `python/astro_metadata_translator/observationInfo.py` (imports near line 10 and 29; new methods after `from_json`, which ends around line 1035)
- Modify: `python/astro_metadata_translator/__init__.py`
- Modify: `tests/test_flat_dict.py`

**Interfaces:**
- Consumes: `FlatValue`, `find_flat_keys`, `find_translator_name`, `flatten_simple`, `unflatten` from Tasks 1 to 3; `ObservationInfo.to_simple`, `ObservationInfo.from_simple`, `MetadataTranslator.translators`.
- Produces on `ObservationInfo`:
  - `to_flat_dict(self, *, prefix: str = "", sep: str = ".", use_upper: bool | None = None) -> dict[str, FlatValue]`
  - `has_flat_dict(cls, flat: Mapping[str, Any], *, prefix: str | None = None, sep: str = ".") -> bool`
  - `strip_from_flat_dict(cls, flat: MutableMapping[str, Any], *, prefix: str | None = None, sep: str = ".") -> None`
  - `from_flat_dict(cls, flat: Mapping[str, Any], *, prefix: str | None = None, sep: str = ".") -> ObservationInfo`
  - `astro_metadata_translator.FITS_FLAT_DICT_PREFIX` and `astro_metadata_translator.FITS_FLAT_DICT_SEP`

- [ ] **Step 1: Add failing round-trip tests**

Add these imports to `tests/test_flat_dict.py` (keep groups sorted; ruff will tell you if not):

```python
import os
from collections.abc import Mapping

import astropy.time
import astropy.units as u

from astro_metadata_translator import (
    FITS_FLAT_DICT_PREFIX,
    FITS_FLAT_DICT_SEP,
    ObservationInfo,
    PropertyDefinition,
    StubTranslator,
    makeObservationInfo,
)
from astro_metadata_translator.tests import read_test_file
```

Remove `FITS_FLAT_DICT_PREFIX` and `FITS_FLAT_DICT_SEP` from the `serialize.flat_dict` import since they now come from the top-level package. After the imports add:

```python
TESTDIR = os.path.abspath(os.path.dirname(__file__))
DATADIR = os.path.join(TESTDIR, "data")

HEADER_FILES = (
    "fitsheader-decam.yaml",
    "fitsheader-hsc.yaml",
    "fitsheader-megaprime.yaml",
    "fitsheader-sdss-fpC-006377-g4-0399.yaml",
    "fitsheader-suprimecam-CORR40535770.yaml",
)

PREFIX_FORMS = (("", "."), ("obs_info", "."), (FITS_FLAT_DICT_PREFIX, FITS_FLAT_DICT_SEP))


# Registers itself through StubTranslator.__init_subclass__ so that
# ``_translator: "flatdict"`` resolves and its extensions are available.
class _FlatDictTranslator(StubTranslator):
    name = "flatdict"
    supported_instrument = "flatdict"
    extensions = dict(
        number=PropertyDefinition("A number", int),
        foo=PropertyDefinition("A string", str),
        pair=PropertyDefinition("A pair of floats", tuple),
    )

    @classmethod
    def can_translate(cls, header: Mapping[str, Any], filename: str | None = None) -> bool:
        return False


def _read_obsinfo(filename: str) -> ObservationInfo:
    return ObservationInfo(read_test_file(filename, dir=DATADIR))
```

Append this class before the `if __name__ == "__main__":` line:

```python
class ObservationInfoFlatDictTestCase(unittest.TestCase):
    """Test the flat dict methods on ObservationInfo."""

    def assert_round_trip(self, original: ObservationInfo, restored: ObservationInfo) -> None:
        # Compare the full simple forms rather than relying on ``==``,
        # which only compares the properties present on both sides.
        self.assertEqual(original.to_simple(), restored.to_simple())

    def test_round_trip_headers(self) -> None:
        for filename in HEADER_FILES:
            obsinfo = _read_obsinfo(filename)
            for prefix, sep in PREFIX_FORMS:
                with self.subTest(filename=filename, prefix=prefix):
                    flat = obsinfo.to_flat_dict(prefix=prefix, sep=sep)
                    for value in flat.values():
                        self.assertIsInstance(value, str | int | float | bool)
                    self.assertTrue(ObservationInfo.has_flat_dict(flat))
                    self.assert_round_trip(obsinfo, ObservationInfo.from_flat_dict(flat))
                    restored = ObservationInfo.from_flat_dict(flat, prefix=prefix, sep=sep)
                    self.assert_round_trip(obsinfo, restored)

    def test_fits_form_keys(self) -> None:
        obsinfo = _read_obsinfo("fitsheader-decam.yaml")
        flat = obsinfo.to_flat_dict(prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        self.assertIn("OBS_INFO DATETIME_BEGIN_0", flat)
        self.assertIn("OBS_INFO DATETIME_BEGIN_1", flat)
        self.assertIn("OBS_INFO TRACKING_RADEC_0", flat)
        self.assertIn("OBS_INFO LOCATION_2", flat)
        self.assertEqual(flat["OBS_INFO TRANSLATOR"], "DECam")
        self.assertNotIn("OBS_INFO _TRANSLATOR", flat)
        for key in flat:
            self.assertEqual(key, key.upper())

    def test_round_trip_hand_built(self) -> None:
        begin = astropy.time.Time("2025-01-01T00:00:00", scale="tai")
        obsinfo = makeObservationInfo(observation_id="A", exposure_time=30.0 * u.s, datetime_begin=begin)
        flat = obsinfo.to_flat_dict()
        self.assertNotIn("translator", flat)
        self.assert_round_trip(obsinfo, ObservationInfo.from_flat_dict(flat))

        # Without any tuple-valued property there is no indexed anchor so
        # auto-detection finds nothing and an explicit prefix is needed.
        obsinfo = makeObservationInfo(observation_id="A", exposure_time=30.0 * u.s)
        flat = obsinfo.to_flat_dict()
        self.assertFalse(ObservationInfo.has_flat_dict(flat))
        with self.assertRaises(ValueError):
            ObservationInfo.from_flat_dict(flat)
        self.assertTrue(ObservationInfo.has_flat_dict(flat, prefix=""))
        self.assert_round_trip(obsinfo, ObservationInfo.from_flat_dict(flat, prefix=""))

    def test_round_trip_extensions(self) -> None:
        begin = astropy.time.Time("2025-01-01T00:00:00", scale="tai")
        obsinfo = makeObservationInfo(
            translator_class=_FlatDictTranslator,
            instrument="flatdict",
            datetime_begin=begin,
            ext_number=42,
            ext_foo="bar",
            ext_pair=(1.5, 2.5),
        )
        for prefix, sep in PREFIX_FORMS:
            with self.subTest(prefix=prefix):
                flat = obsinfo.to_flat_dict(prefix=prefix, sep=sep)
                restored = ObservationInfo.from_flat_dict(flat)
                self.assertEqual(restored.translator_name, "flatdict")
                self.assertEqual(restored.ext_number, 42)
                self.assertEqual(restored.ext_foo, "bar")
                self.assertEqual(restored.ext_pair, (1.5, 2.5))
                self.assert_round_trip(obsinfo, restored)
        flat = obsinfo.to_flat_dict(prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        self.assertEqual(flat["OBS_INFO EXT_PAIR_0"], 1.5)
        self.assertEqual(flat["OBS_INFO EXT_PAIR_1"], 2.5)

    def test_has_and_strip(self) -> None:
        obsinfo = _read_obsinfo("fitsheader-hsc.yaml")
        flat = obsinfo.to_flat_dict(prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        n_flat = len(flat)
        header = {"SIMPLE": True, "OBJECT": "decoy", "TEMPERATURE": 3.0, **flat}
        self.assertTrue(ObservationInfo.has_flat_dict(header))
        self.assertTrue(
            ObservationInfo.has_flat_dict(header, prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        )
        ObservationInfo.strip_from_flat_dict(header)
        self.assertEqual(header, {"SIMPLE": True, "OBJECT": "decoy", "TEMPERATURE": 3.0})
        self.assertEqual(len(flat), n_flat)
        self.assertFalse(ObservationInfo.has_flat_dict(header))
        # Stripping again is a no-op.
        ObservationInfo.strip_from_flat_dict(header)
        self.assertEqual(len(header), 3)

    def test_from_flat_dict_errors(self) -> None:
        with self.assertRaises(ValueError):
            ObservationInfo.from_flat_dict({})
        with self.assertRaises(ValueError):
            ObservationInfo.from_flat_dict({"OBJECT": "M31"})

        obsinfo = _read_obsinfo("fitsheader-decam.yaml")
        flat = obsinfo.to_flat_dict(prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)

        gap = dict(flat)
        del gap["OBS_INFO DATETIME_BEGIN_1"]
        with self.assertRaises(ValueError):
            ObservationInfo.from_flat_dict(gap)

        corrupt = dict(flat)
        corrupt["OBS_INFO EXPOSURE_ID"] = "not an integer"
        # Validation errors propagate rather than being mistaken for absence.
        with self.assertRaises((TypeError, ValueError)):
            ObservationInfo.from_flat_dict(corrupt)
        self.assertTrue(ObservationInfo.has_flat_dict(corrupt))

    def test_explicit_prefix_wrong_sep(self) -> None:
        obsinfo = _read_obsinfo("fitsheader-decam.yaml")
        flat = obsinfo.to_flat_dict(prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        self.assertFalse(ObservationInfo.has_flat_dict(flat, prefix=FITS_FLAT_DICT_PREFIX, sep="."))
        with self.assertRaises(ValueError):
            ObservationInfo.from_flat_dict(flat, prefix=FITS_FLAT_DICT_PREFIX, sep=".")

    def test_to_flat_dict_bad_sep(self) -> None:
        obsinfo = _read_obsinfo("fitsheader-decam.yaml")
        with self.assertRaises(ValueError):
            obsinfo.to_flat_dict(prefix="p", sep="_")
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_flat_dict.py -q`
Expected: FAIL at import with `ImportError: cannot import name 'FITS_FLAT_DICT_PREFIX' from 'astro_metadata_translator'`.

- [ ] **Step 3: Export the constants from the top-level package**

In `python/astro_metadata_translator/__init__.py`, after the line `from .properties import *`, add:

```python
from .serialize.flat_dict import FITS_FLAT_DICT_PREFIX, FITS_FLAT_DICT_SEP
```

- [ ] **Step 4: Add the methods to ObservationInfo**

In `python/astro_metadata_translator/observationInfo.py`:

Change the `collections.abc` import to:

```python
from collections.abc import Generator, Mapping, MutableMapping, Sequence
```

After `from .properties import (...)` and before `from .translator import MetadataTranslator`, add:

```python
from .serialize.flat_dict import FlatValue, find_flat_keys, find_translator_name, flatten_simple, unflatten
```

Immediately after the `from_json` classmethod (the one ending with `return cls.model_validate_json(json_str)`), add:

```python
    def to_flat_dict(
        self,
        *,
        prefix: str = "",
        sep: str = ".",
        use_upper: bool | None = None,
    ) -> dict[str, FlatValue]:
        """Convert the contents of this object to a flat dict of scalars.

        Parameters
        ----------
        prefix : `str`, optional
            Opaque prefix placed before every key. Emitted exactly as
            given. Use `FITS_FLAT_DICT_PREFIX` for FITS headers.
        sep : `str`, optional
            Separator between the prefix and the property name. Must be
            one character that is not a letter, digit, or underscore.
            Unused when ``prefix`` is empty. Use `FITS_FLAT_DICT_SEP` for
            FITS headers.
        use_upper : `bool` or `None`, optional
            If `None` the property names are upper-cased when the first
            character of ``prefix`` is upper case and left lower case
            otherwise. `True` or `False` force the case.

        Returns
        -------
        flat : `dict` [`str`, `str` | `int` | `float` | `bool`]
            Flat mapping of all properties that are not `None`.

        Raises
        ------
        ValueError
            Raised if ``sep`` is not valid.
        TypeError
            Raised if a property value can not be represented as a
            scalar or a flat sequence of scalars.

        Notes
        -----
        Each key is the property name as used by `to_simple`, optionally
        preceded by the prefix and separator. A property whose simple
        form is a tuple, such as ``datetime_begin``, is written as one
        key per element with an index suffix: ``datetime_begin_0`` and
        ``datetime_begin_1``. The ``_translator`` entry is written as
        ``translator``. With ``prefix="OBS_INFO"`` and ``sep=" "`` the
        keys are suitable for FITS ``HIERARCH`` cards, for example
        ``OBS_INFO DATETIME_BEGIN_0``.

        Values are passed through unchanged. FITS can not represent
        ``NaN`` so a writer must remove such values before storing the
        result in a header.

        The result can be converted back with `from_flat_dict`.
        Round-tripping of extension properties requires that the
        translator is registered, and an extension name must not end in
        an underscore followed by digits if a tuple-valued extension
        with the shorter name also exists.
        """
        return flatten_simple(self.to_simple(), prefix=prefix, sep=sep, use_upper=use_upper)

    @classmethod
    def has_flat_dict(cls, flat: Mapping[str, Any], *, prefix: str | None = None, sep: str = ".") -> bool:
        """Report whether a mapping contains keys written by `to_flat_dict`.

        Parameters
        ----------
        flat : `~collections.abc.Mapping` [`str`, `~typing.Any`]
            Mapping to inspect, for example a FITS header. Keys are
            matched case-insensitively.
        prefix : `str` or `None`, optional
            Prefix given to `to_flat_dict`. `None` detects the prefix
            and separator automatically, in which case ``sep`` is
            ignored.
        sep : `str`, optional
            Separator given to `to_flat_dict`.

        Returns
        -------
        present : `bool`
            `True` if at least one key was found.

        Raises
        ------
        ValueError
            Raised if automatic detection finds two prefixes that match
            equally well.

        Notes
        -----
        Use this before `from_flat_dict` to tell the absence of content
        from content that fails validation. Automatic detection requires
        an indexed property such as ``datetime_begin_0`` to be present.
        """
        return bool(find_flat_keys(flat, prefix=prefix, sep=sep))

    @classmethod
    def strip_from_flat_dict(
        cls, flat: MutableMapping[str, Any], *, prefix: str | None = None, sep: str = "."
    ) -> None:
        """Remove keys written by `to_flat_dict` from a mapping in place.

        Parameters
        ----------
        flat : `~collections.abc.MutableMapping` [`str`, `~typing.Any`]
            Mapping to modify, for example a FITS header. Keys are
            matched case-insensitively.
        prefix : `str` or `None`, optional
            Prefix given to `to_flat_dict`. `None` detects the prefix
            and separator automatically, in which case ``sep`` is
            ignored.
        sep : `str`, optional
            Separator given to `to_flat_dict`.

        Notes
        -----
        Nothing is removed if no keys are found. Other keys in the
        mapping are left untouched.
        """
        for key in find_flat_keys(flat, prefix=prefix, sep=sep):
            del flat[key]

    @classmethod
    def from_flat_dict(
        cls, flat: Mapping[str, Any], *, prefix: str | None = None, sep: str = "."
    ) -> ObservationInfo:
        """Create an `ObservationInfo` from the output of `to_flat_dict`.

        Parameters
        ----------
        flat : `~collections.abc.Mapping` [`str`, `~typing.Any`]
            Mapping containing the flat keys, for example a FITS header.
            Keys are matched case-insensitively and unrelated keys are
            ignored.
        prefix : `str` or `None`, optional
            Prefix given to `to_flat_dict`. `None` detects the prefix
            and separator automatically, in which case ``sep`` is
            ignored.
        sep : `str`, optional
            Separator given to `to_flat_dict`.

        Returns
        -------
        obsinfo : `ObservationInfo`
            New object constructed from the mapping.

        Raises
        ------
        ValueError
            Raised if no keys are found, if automatic detection is
            ambiguous, or if the indexed values of a property are
            incomplete. Also raised, as `pydantic.ValidationError`, if a
            value fails validation.

        Notes
        -----
        Round-tripping of extension properties requires that the
        translator named in the mapping is registered, because its
        extension definitions are needed to interpret ``ext_`` keys.
        """
        flat_keys = find_flat_keys(flat, prefix=prefix, sep=sep)
        if not flat_keys:
            raise ValueError("No ObservationInfo content found in the given mapping.")
        extension_names: tuple[str, ...] = ()
        translator_name = find_translator_name(flat, flat_keys)
        if translator_name is not None:
            translator_class = MetadataTranslator.translators.get(translator_name)
            if translator_class is not None:
                extension_names = tuple(translator_class.extensions)
        return cls.from_simple(unflatten(flat, flat_keys, extension_names=extension_names))
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python -m pytest tests/test_flat_dict.py -q`
Expected: 35 passed. Then run the whole suite to check nothing else broke: `python -m pytest -q`. Expected: all pass.

- [ ] **Step 6: Lint and type-check**

Run: `ruff check python tests && ruff format --check python tests && mypy python`
Expected: no errors. If the `serialize.flat_dict` import in `observationInfo.py` causes a circular import at runtime, the test run in Step 5 fails at import; in that case confirm that `serialize/fits.py` still imports `ObservationInfo` only under `TYPE_CHECKING`.

- [ ] **Step 7: Commit**

```bash
git add python/astro_metadata_translator/observationInfo.py python/astro_metadata_translator/__init__.py tests/test_flat_dict.py
git commit -m "Add ObservationInfo flat dict methods.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 5: FITS integration tests, NaN pin, and news fragment

**Files:**
- Modify: `tests/test_flat_dict.py`
- Create: `doc/changes/DM-56697.feature.rst`

**Interfaces:**
- Consumes: everything from Task 4.
- Produces: nothing new; verifies the real FITS path and documents the feature.

- [ ] **Step 1: Add the FITS and NaN tests**

Add to the imports in `tests/test_flat_dict.py`:

```python
import math
import tempfile

from astropy.io import fits
```

Append this class before the `if __name__ == "__main__":` line:

```python
class FitsHeaderTestCase(unittest.TestCase):
    """Test the flat dict through a real FITS header file."""

    def _write_and_read(self, header: fits.Header) -> fits.Header:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test.fits")
            fits.PrimaryHDU(header=header).writeto(path)
            return fits.getheader(path)

    def test_fits_header_round_trip(self) -> None:
        obsinfo = _read_obsinfo("fitsheader-decam.yaml")
        flat = obsinfo.to_flat_dict(prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)

        header = fits.Header()
        header["OBJECT"] = "decoy"
        header["HIERARCH TEMPERATURE"] = 3.0
        for key, value in flat.items():
            header[f"HIERARCH {key}"] = value
        read = self._write_and_read(header)

        self.assertTrue(ObservationInfo.has_flat_dict(read))
        self.assertTrue(
            ObservationInfo.has_flat_dict(read, prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        )

        restored = ObservationInfo.from_flat_dict(read)
        self.assertEqual(obsinfo.to_simple(), restored.to_simple())
        self.assertEqual(restored.object, obsinfo.object)
        self.assertNotEqual(restored.object, "decoy")
        self.assertIsInstance(restored.can_see_sky, bool)
        self.assertIsInstance(restored.has_simulated_content, bool)

        restored = ObservationInfo.from_flat_dict(read, prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
        self.assertEqual(obsinfo.to_simple(), restored.to_simple())

        ObservationInfo.strip_from_flat_dict(read)
        self.assertFalse(ObservationInfo.has_flat_dict(read))
        self.assertFalse(any(key.startswith("OBS_INFO") for key in read))
        self.assertEqual(read["OBJECT"], "decoy")
        self.assertEqual(read["TEMPERATURE"], 3.0)

    def test_fits_header_without_content(self) -> None:
        header = fits.Header()
        header["OBJECT"] = "M31"
        header["HIERARCH TEMPERATURE"] = 3.0
        header["HIERARCH PRESSURE"] = 700.0
        read = self._write_and_read(header)
        n_cards = len(read)
        self.assertFalse(ObservationInfo.has_flat_dict(read))
        ObservationInfo.strip_from_flat_dict(read)
        self.assertEqual(len(read), n_cards)
        with self.assertRaises(ValueError):
            ObservationInfo.from_flat_dict(read)

    def test_nan_passes_through_plain_dict(self) -> None:
        obsinfo = makeObservationInfo(observation_id="A", relative_humidity=math.nan)
        flat = obsinfo.to_flat_dict()
        self.assertTrue(math.isnan(flat["relative_humidity"]))
        restored = ObservationInfo.from_flat_dict(flat, prefix="")
        self.assertIsNotNone(restored.relative_humidity)
        self.assertTrue(math.isnan(restored.relative_humidity))
```

Note: `relative_humidity` is a `float | None` property so the `assertIsNotNone` keeps mypy happy before `math.isnan`; if mypy still complains, use `assert restored.relative_humidity is not None`.

- [ ] **Step 2: Run the tests to verify the new ones pass**

Run: `python -m pytest tests/test_flat_dict.py -q`
Expected: 38 passed. These tests exercise code written in Task 4, so they pass immediately; if any fails, the failure is a real defect in Task 4 code and must be fixed there, not by weakening the test.

- [ ] **Step 3: Write the news fragment**

Create `doc/changes/DM-56697.feature.rst`:

```rst
Added ``ObservationInfo.to_flat_dict``, ``from_flat_dict``, ``has_flat_dict``, and ``strip_from_flat_dict`` for storing an ``ObservationInfo`` as scalar key/value pairs, for example in a FITS header, and rebuilding it losslessly.
Tuple-valued properties are written as indexed keys such as ``datetime_begin_0`` and keys can carry an optional prefix; the recommended FITS form uses the new ``FITS_FLAT_DICT_PREFIX`` and ``FITS_FLAT_DICT_SEP`` constants.
```

- [ ] **Step 4: Full verification**

Run: `python -m pytest -q && ruff check python tests && ruff format --check python tests && mypy python`
Expected: all tests pass, no lint or type errors.

- [ ] **Step 5: Commit**

```bash
git add tests/test_flat_dict.py doc/changes/DM-56697.feature.rst
git commit -m "Add FITS header tests and news fragment for flat dict.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```
