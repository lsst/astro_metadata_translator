# ObservationInfo flat dict serialization

Ticket: DM-56697.
Date: 2026-10-14.

## Goal

Provide a lossless way to store an `ObservationInfo` as scalar key/value pairs that can live in a FITS header, and to rebuild the `ObservationInfo` from that header later.

The motivating flow is Rubin raw ingest and image conversion.
obs_base builds an `ObservationInfo` when it reads a raw file and can write the flat dict into the exposure metadata.
Those cards ride along in the primary header through the pipeline.
`VisitImage.from_legacy` in lsst.images can then rebuild the exact `ObservationInfo` from the header instead of reconstructing an incomplete one from `VisitInfo` and butler dimension records.

The API mirrors `lsst.daf.butler.DatasetProvenance.to_flat_dict` and `from_flat_dict` so the two mechanisms feel the same to callers.

## Scope

In scope for this ticket, all in astro_metadata_translator:

- `ObservationInfo.to_flat_dict`, `ObservationInfo.from_flat_dict`, `ObservationInfo.has_flat_dict`, and `ObservationInfo.strip_from_flat_dict`.
- A new helper module `python/astro_metadata_translator/serialize/flat_dict.py` holding the key logic as pure functions.
- Two shared constants naming the recommended FITS prefix and separator.
- Tests and a news fragment.

Out of scope, to be done on follow-up tickets:

- obs_base writing the flat dict into raw exposure metadata.
- lsst.images reading it in `VisitImage.from_legacy`.

## Public API

All four methods live on `ObservationInfo` in `python/astro_metadata_translator/observationInfo.py` and delegate to the helper module.

### `to_flat_dict(self, *, prefix: str = "", sep: str = ".", use_upper: bool | None = None) -> dict[str, str | int | float | bool]`

Calls `to_simple()` and flattens the result.

- `prefix` is an opaque string placed before every key.
  It is never split into components.
- `sep` joins the prefix to the flat name.
  It must be exactly one character and must not be a word character or underscore, the same rule as provenance.
  It is unused when `prefix` is empty.
- `use_upper=None` follows the case of the first character of `prefix`, so an upper-case prefix yields upper-case flat names and an empty prefix yields lower-case flat names.
  `True` and `False` force the case.
  The prefix itself is always emitted exactly as given, as provenance does.
- Properties whose value is `None` are omitted, exactly as `to_simple` omits them.

### `from_flat_dict(cls, flat: Mapping[str, Any], *, prefix: str | None = None, sep: str = ".") -> ObservationInfo`

Classmethod.

- `prefix=None` auto-detects the prefix and separator using the majority vote described below.
  The `sep` argument is ignored in that mode.
- Any other `prefix`, including `""`, skips detection and is used together with `sep`.
- Keys are matched case-insensitively so an `astropy.io.fits.Header` or an `lsst.daf.base.PropertyList` can be passed directly.
- Raises `ValueError` when no keys are found.

### `has_flat_dict(cls, flat: Mapping[str, Any], *, prefix: str | None = None, sep: str = ".") -> bool`

Classmethod.
Returns `True` when the detection described below finds at least one key.
It exists so a reader can distinguish "absent" from "present but corrupt": check first, then let any error from `from_flat_dict` propagate.

### `strip_from_flat_dict(cls, flat: MutableMapping[str, Any], *, prefix: str | None = None, sep: str = ".") -> None`

Classmethod.
Deletes every matched key from the mapping in place.
Uses the same detection as `from_flat_dict` but is a silent no-op when nothing is found, so a writer can call it unconditionally before adding fresh cards.

### Constants

`FITS_FLAT_DICT_PREFIX = "OBS_INFO"` and `FITS_FLAT_DICT_SEP = " "` are defined in the helper module and re-exported from the top-level `astro_metadata_translator` package.
obs_base and lsst.images should import these rather than spelling the prefix themselves.

## Key format

### Flat names

The flat name for a scalar property is the property name, for example `exposure_time`.

A property whose simple form is a list or tuple is expanded one level into `{name}_{i}` for each element, counting from 0.
Today that is `datetime_begin`, `datetime_end`, `location`, `tracking_radec`, `altaz_begin`, and `altaz_end`, giving names such as `datetime_begin_0` and `location_2`.
Extension values are expanded by the same rule, so an extension named `foo` with a two-element simple form becomes `ext_foo_0` and `ext_foo_1`.
The rule is applied to whatever `to_simple` returns, so a new tuple-valued property needs no change in this code.

The `_translator` key from `to_simple` is renamed to `translator` so the FITS card has no leading underscore.
It is restored to `_translator` on read.

### Full keys

With a prefix the key is `f"{prefix}{sep}{flat_name}"`.
Without a prefix the key is the flat name alone.
The flat name is upper-cased when `use_upper` resolves to true; the prefix is never changed.

Examples for `datetime_begin`:

| Call | Key |
|---|---|
| `to_flat_dict()` | `datetime_begin_0` |
| `to_flat_dict(prefix="obs_info")` | `obs_info.datetime_begin_0` |
| `to_flat_dict(prefix="OBS_INFO", sep=" ")` | `OBS_INFO DATETIME_BEGIN_0` |

The last form is what goes into a FITS header, where the writer adds the `HIERARCH` keyword.

### Naming constraint

Because the index is joined with an underscore, the reader must decide whether a trailing `_<digits>` is an index or part of the name.
The rule is that it is an index only when the stem is a known property or extension name.
No core property ends in `_<digits>`, so core names are never ambiguous.
An extension is ambiguous only when both `foo` and `foo_1` are defined on the same translator, and that combination is unsupported.
This constraint is documented in the `to_flat_dict` docstring.

## Detection

Detection runs when `from_flat_dict` or `strip_from_flat_dict` is called with `prefix=None`.

Each key is lower-cased and matched against this pattern, where `{core}` is the alternation of all `PROPERTIES` names plus `translator`:

```
^(?:(?P<prefix>.+)(?P<sep>[^\w]))?(?P<name>{core})(?:_(?P<idx>\d+))?$
```

The greedy `prefix` group and the end anchor mean a prefix that contains the separator, such as `LSST OBS_INFO` with `sep=" "`, is recovered intact.
An empty prefix group means no prefix and no separator.

Matches are grouped by `(prefix, sep)` and each group counts the distinct names it matched.
Only groups that contain an indexed match for a tuple-valued core property, such as `datetime_begin_0`, are candidates; the indexed form of a scalar property such as `object_1` does not count.
Every real `ObservationInfo` has `datetime_begin`, and ordinary FITS headers never carry `name_<digits>` cards for these names, so this anchor stops a header with only an `OBJECT` card from being mistaken for ObservationInfo content and stops `strip_from_flat_dict` from deleting it.
A hand-built `ObservationInfo` with no tuple-valued property can still be read by passing `prefix` explicitly.
Among the candidates, the group with the highest count wins.
This is a heuristic chosen because core property names collide with ordinary FITS cards such as `OBJECT`, `TEMPERATURE`, and `PRESSURE`.
A real header always contains `datetime_begin_0` and a couple of dozen other names under one prefix, so a stray `OBJECT` card loses the vote.
A tie between the top two groups raises `ValueError`.
Extension keys never vote because their names are open-ended.

## Reading

Once the prefix and separator are known, either by detection or from the arguments, every key is lower-cased.
With a non-empty prefix, only keys whose lower-cased form starts with the lower-cased prefix followed by the separator are considered, and that leading part is removed.
With an empty prefix every key is considered as is.
A remainder is kept when it is a core name, `translator`, or starts with `ext_`, with or without an `_<digits>` suffix.
Every other key in the mapping is ignored, so other header cards are untouched.

Indexed entries are gathered per stem, sorted by index, checked to be a contiguous run from 0, and turned into a tuple.
For core properties the run must also have the length the `ObservationInfo` JSON schema declares for that property, so a truncated set such as `datetime_begin_0` without `datetime_begin_1` is rejected rather than passed on.
For `ext_` keys the known names come from the extension definitions of the translator named by the `translator` value, looked up in `MetadataTranslator.translators`.
Extension names are matched case-insensitively against those definitions and the original spelling from the definition is used in the result.
If the translator is absent or unregistered, `ext_` keys are passed through as plain lower-case names and `from_simple` reports them as it does today.

`translator` is renamed back to `_translator` and the regrouped dict is passed to `from_simple`.
All type coercion and validation stays there.

## Error handling

Writing:

- A `sep` that is not exactly one character, or that is a letter, digit, or underscore, raises `ValueError`.
- A value that is not `str`, `int`, `float`, `bool`, or a one-level list or tuple of those raises `TypeError` naming the property.

Reading, all `ValueError`:

- No keys found.
- Vote tie between prefixes.
- An index gap, for example `datetime_begin_0` present but `datetime_begin_1` missing.
- Two keys that collapse to the same lower-case name.
- An explicit prefix that matches nothing, which is what a `sep` mismatch looks like.

Errors raised by `from_simple` propagate unchanged so corrupt cards are never mistaken for absent ones.

Stripping never raises for an empty result.

## FITS considerations

These do not change the astro_metadata_translator code but shaped the choices above and matter to the downstream tickets.

Card budget for the longest flat name, `BORESIGHT_ROTATION_ANGLE`, counting `HIERARCH `, the keyword, and ` = ` on an 80-character card:

| Prefix | Characters left for value |
|---|---|
| `OBSERVATION_INFO` | 27 |
| `LSST OBS_INFO` | 30 |
| `OBS_INFO` | 35 |

Floats are at most 20 characters from astropy and 23 from cfitsio, so no float or integer card overflows with any of these prefixes.
`OBS_INFO` was chosen because it buys extra room for strings, keeps cards readable, and matches the `VisitImage.obs_info` attribute name.

Long strings are written with `CONTINUE` cards by cfitsio through afw's `fits_write_key_longstr` path and by astropy from version 7.1.
obs_base already elides provenance strings on older astropy and can apply the same `allow_long_headers` treatment here.

Floats do not always round-trip bit for bit.
Astropy writes a float exactly when its shortest representation fits in 20 characters and otherwise truncates the formatted value to that 20-character field, which can leave 15 significant digits; cfitsio always writes 16 significant digits.
A value such as a declination of `-0.0015811524500870524` therefore comes back from astropy as `-0.00158115245008705`.
Round-trip tests through a real FITS file, here and in lsst.images, compare floats to 15 significant digits and everything else exactly.

## Intended downstream use

obs_base, when building raw exposure metadata:

```python
ObservationInfo.strip_from_flat_dict(md, prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
md.update(obs_info.to_flat_dict(prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP))
```

lsst.images, in `VisitImage.from_legacy`, before the current `obs_info_from_legacy` fallback:

```python
if ObservationInfo.has_flat_dict(primary_header, prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP):
    obs_info = ObservationInfo.from_flat_dict(primary_header, prefix=FITS_FLAT_DICT_PREFIX, sep=FITS_FLAT_DICT_SEP)
else:
    obs_info = obs_info_from_legacy(...)
```

Catching `ValueError` around `from_flat_dict` would be wrong because pydantic validation errors are also `ValueError`, and corrupt cards must not fall back silently.

## Testing

New file `tests/test_flat_dict.py`.

Helper module unit tests:

- Key building for empty, lower-case, and upper-case prefixes and for `use_upper` overrides.
- Separator validation.
- Detection with `OBJECT` and `TEMPERATURE` decoy cards, a prefix containing the separator, an empty prefix, a tie, and a header with decoys but no indexed anchor.
- Regrouping of contiguous indexes, a gap, and an `ext_` stem that is and is not a known extension.

Round trips on `ObservationInfo`:

- The DECam, HSC, MegaPrime, SDSS, and SuprimeCam headers in `tests/data`, each through no prefix, `obs_info` with `.`, and `OBS_INFO` with a space, compared with the existing `__eq__`.
- The same headers read back with an explicit prefix.
- One real FITS pass: put the `OBS_INFO` dict into an `astropy.io.fits.Header`, write it to a temporary file, read it back, and rebuild from the `Header` object.
- The dummy translator from `tests/test_extensions.py` with `ext_number` and `ext_foo`, plus a tuple-valued extension to exercise `ext_*_N` regrouping.
- A hand-built `ObservationInfo` from `makeObservationInfo` with no translator.

Presence and strip:

- `has_flat_dict` is true for each prefix form and false for a header with only decoy cards.
- A mixed header where only the `OBS_INFO` cards are removed.
- An empty result is a no-op.

Every error case in the error handling section.

## Housekeeping

- News fragment `doc/changes/DM-56697.feature.rst`.
- ruff and mypy clean.
- Docstrings follow the numpydoc style used elsewhere in the package.
