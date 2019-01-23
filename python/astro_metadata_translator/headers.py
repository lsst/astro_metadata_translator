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

"""Code to support header manipulation operations."""

__all__ = ("merge_headers",)

import itertools
import copy


def merge_headers(headers, mode="overwrite"):
    """Merge multiple headers into a single dict.

    Given a list of dict-like data headers, combine them following the
    specified mode.

    Parameters
    ----------
    headers : `list` of `dict` (or `dict`-like)
        Collection of headers to combine. `~lsst.daf.base.PropertyList`
        is supported.
    mode : `str`
        Scheme to use when a header has the same key as another header
        but different value. Options are:

        - ``'overwrite'`` : Value in later header overwrites earlier value.
        - ``'drop'`` : Entire key is dropped.
        - ``'first'`` : Retain first value encountered.
        - ``'append'`` : Convert value to list with a value for each header
          (`None` if the key was not present). If the value is
          identical in multiple headers but key is missing in
          some, then the single identical header is stored.

    Returns
    -------
    merged : `dict`
        Single `dict` combining all the headers using the specified
        combination mode.
    """
    if not headers:
        raise ValueError("No headers supplied.")

    # Force PropertyList to OrderedDict
    # In python 3.7 dicts are guaranteed to retain order
    headers = [h.toOrderedDict() if hasattr(h, "toOrderedDict") else h for h in headers]

    if mode == "overwrite":
        merged = copy.deepcopy(headers.pop(0))
        for h in headers:
            merged.update(h)

    elif mode == "first":
        # Reversing the headers and using overwrite mode would result in the
        # header order being inconsistent dependent on mode.
        merged = copy.deepcopy(headers.pop(0))
        for hdr in headers:
            for key in hdr:
                if key not in merged:
                    merged[key] = hdr[key]

    elif mode == "drop":
        merged = copy.deepcopy(headers.pop(0))
        drop = set()
        for hdr in headers:
            for key in hdr:
                if key not in merged:
                    merged[key] = hdr[key]
                elif merged[key] != hdr[key]:
                    # Key should be dropped later (not in loop since removing
                    # the key now might add it back for the next header).
                    drop.add(key)

        for key in drop:
            del merged[key]

    elif mode == "append":
        first = headers.pop(0)
        merged = copy.deepcopy(first)
        fill = set()
        for hdr in headers:
            for key in hdr:
                if key not in merged:
                    merged[key] = hdr[key]
                elif not isinstance(merged[key], list) and merged[key] != hdr[key]:
                    # If we detect different values, store an empty list
                    # in the slot and fill it later.  Do it at end so
                    # we can pick up earlier values and fill empty with None.
                    merged[key] = []
                    fill.add(key)

        # Fill the entries that have multiple differing values
        for key in fill:
            merged[key] = [h[key] if key in h else None
                           for h in itertools.chain([first], headers)]

    else:
        raise ValueError(f"Unsupported value of '{mode}' for mode parameter.")

    return merged
