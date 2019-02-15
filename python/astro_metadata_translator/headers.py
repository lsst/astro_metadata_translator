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

import logging
import itertools
import copy

from .translator import MetadataTranslator
from .translators import FitsTranslator

log = logging.getLogger(__name__)


def merge_headers(headers, mode="overwrite", sort=False, first=None, last=None):
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
    sort : `bool`
        If `True`, sort the supplied headers into date order if possible.
        This affects how the resulting merged output depending on the requested
        merge mode.  An attempt will be made to extract a date from the
        headers.
    first : `list` or `tuple`
        Keys to retain even if they differ.  The value in the merged header
        will always be the value first encountered, independently of ``mode``.
        This is usually to allow time-dependent headers such as ``DATE-OBS``
        and ``AZSTART`` to be retained to allow the header to indicate the
        range of values.  Not used if ``mode`` is ``append``.    No exception
        is raised if a key can not be found in a header since this allows a
        range of expected headers to be listed covering multiple instruments.
    last : `list` or `tuple`
        Keys to retain even if they differ.  The value in the merged header
        will always be the final value encountered, independently of ``mode``.
        This is usually to allow time-dependent headers such as ``DATE-END``
        and ``AZEND`` to be retained to allow the header to indicate the
        range of values.  Not used if ``mode`` is ``append``.  No exception
        is raised if a key can not be found in a header since this allows a
        range of expected headers to be listed covering multiple instruments.

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

    # With a single header provided return a copy immediately
    if len(headers) == 1:
        return copy.deepcopy(headers[0])

    if sort:
        def key_func(hdr):
            translator_class = None
            try:
                translator_class = MetadataTranslator.determine_translator(hdr)
            except ValueError:
                # Try the FITS translator
                translator_class = FitsTranslator
            translator = translator_class(hdr)
            return translator.to_datetime_begin()

        headers = sorted(headers, key=key_func)

    log.debug("Received %d headers for merging", len(headers))

    # Take copy of first header
    first_hdr = headers.pop(0)

    # Seed the merged header
    merged = copy.deepcopy(first_hdr)

    if mode == "overwrite":
        for h in headers:
            merged.update(h)

    elif mode == "first":
        # Reversing the headers and using overwrite mode would result in the
        # header order being inconsistent dependent on mode.
        for hdr in headers:
            for key in hdr:
                if key not in merged:
                    merged[key] = hdr[key]

    elif mode == "drop":
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
                           for h in itertools.chain([first_hdr], headers)]

    else:
        raise ValueError(f"Unsupported value of '{mode}' for mode parameter.")

    # Force the first and last values to be inserted
    #
    if mode != "append":
        def retain_value(to_receive, to_retain, sources):
            if to_retain:
                for k in to_retain:
                    # Look for values until we find one
                    for h in sources:
                        if k in h:
                            to_receive[k] = h[k]
                            break

        all_headers = (first_hdr, *headers)
        retain_value(merged, first, all_headers)
        retain_value(merged, last, tuple(reversed(all_headers)))

    return merged
