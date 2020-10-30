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

__all__ = ("merge_headers", "fix_header")

import pkg_resources
import posixpath
import logging
import itertools
import copy
import os
import yaml
from collections.abc import Mapping

from .translator import MetadataTranslator
from .translators import FitsTranslator

log = logging.getLogger(__name__)

ENV_VAR_NAME = "METADATA_CORRECTIONS_PATH"
"""Name of environment variable containing search path for header fix up."""


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
    sort : `bool`, optional
        If `True`, sort the supplied headers into date order if possible.
        This affects the resulting merged output depending on the requested
        merge mode.  An attempt will be made to extract a date from the
        headers.
    first : `list` or `tuple`, optional
        Keys to retain even if they differ.  For all modes excepting ``append``
        (where it is ignored) the value in the merged header will always be
        the value first encountered.  This is usually to allow time-dependent
        headers such as ``DATE-OBS`` and ``AZSTART`` to be retained to allow
        the header to indicate the range of values.  No exception is raised if
        a key can not be found in a header since this allows a range of
        expected headers to be listed covering multiple instruments.
    last : `list` or `tuple`, optional
        Keys to retain even if they differ.  For all modes excepting ``append``
        (where it is ignored) the value in the merged header will always be
        the final value encountered.  This is usually to allow time-dependent
        headers such as ``DATE-END`` and ``AZEND`` to be retained to allow
        the header to indicate the range of values.  No exception is raised if
        a key can not be found in a header since this allows a range of
        expected headers to be listed covering multiple instruments.

    Returns
    -------
    merged : `dict`
        Single `dict` combining all the headers using the specified
        combination mode.

    Notes
    -----
    If ``first`` and ``last`` are supplied, the keys from ``first`` are
    handled first, followed by the keys from ``last``.  No check is made to
    ensure that the keys do not overlap.
    """
    if not headers:
        raise ValueError("No headers supplied.")

    # Copy the input list because we will be reorganizing it
    headers = list(headers)

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

    # Pull out first header
    first_hdr = headers.pop(0)

    # Seed the merged header with a copy
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


def _read_yaml(fh, msg):
    """Read YAML from file descriptor.

    Parameters
    ----------
    fh : `io.IOBase`
        Open file handle containing the YAML stream
    msg : `str`
        Text to include in log file when referring to this stream. Examples
        could be "file something.yaml" or "resource module:resource".

    Returns
    -------
    parsed : `dict` or `None`
        The contents of the YAML file if it was a `dict`, else `None` if
        the contents could not be parsed or the contents were YAML but
        not a mapping.
    """
    try:
        content = yaml.safe_load(fh)
    except Exception as e:
        log.warning("Error parsing YAML header corrections from %s: %s", msg, str(e))
        return None

    if not isinstance(content, Mapping):
        log.warning("YAML Mapping not found in %s. Ignoring contents.", msg)
        return None

    return content


def _find_from_file(header, paths, target_file):
    """Search file system for matching correction files.

    Parameters
    ----------
    header : `dict`
        Header to update.
    paths : `list`
        Paths to search.
    target_file : `str`
        File to locate in the path.

    Returns
    -------
    modified : `bool`
        `True` if a correction was found. Only the first correction located
        in a path is used.
    """
    for p in paths:
        correction_file = os.path.join(p, target_file)
        if os.path.exists(correction_file):
            with open(correction_file) as fh:
                log.debug("Applying header corrections from file %s", correction_file)
                corrections = _read_yaml(fh, f"file {correction_file}")

            if corrections is None:
                continue

            # Apply corrections
            header.update(corrections)

            return True
    return False


def _find_from_resource(header, package, resource_root, target_file):
    """Search package resource for correction information.

    Parameters
    ----------
    header : `dict`
        Header to update.
    package : `str`
        Package resource to search.
    resource_root : `str`
        Resource root.
    target_file : `str`
        Resource to locate.

    Returns
    -------
    modified : `bool`
        `True` if a correction was found.
    """
    if package is not None and resource_root is not None:
        resource_name = posixpath.join(resource_root, target_file)
        if pkg_resources.resource_exists(package, resource_name):
            log.debug("Applying header corrections from package resource %s:%s", package, resource_name)
            with pkg_resources.resource_stream(package, resource_name) as fh:
                corrections = _read_yaml(fh, f"package resource {package}:{resource_name}")

            if corrections is None:
                return False

            header.update(corrections)

            return True
    return False


def fix_header(header, search_path=None, translator_class=None, filename=None):
    """Update, in place, the supplied header with known corrections.

    Parameters
    ----------
    header : `dict`-like
        Header to correct.
    search_path : `list` or `str`, optional
        Explicit directory paths to search for correction files.
        A single directory path can be given as a string.
    translator_class : `MetadataTranslator`-class, optional
        If not `None`, the class to use to translate the supplied headers
        into standard form. Otherwise each registered translator class will
        be asked in turn if it knows how to translate the supplied header.
    filename : `str`, optional
        Name of the file whose header is being translated.  For some
        datasets with missing header information this can sometimes
        allow for some fixups in translations.

    Returns
    -------
    fixed : `bool`
        `True` if the header was updated.

    Raises
    ------
    TypeError
        Raised if the supplied translation class is not a `MetadataTranslator`.

    Notes
    -----
    In order to determine that a header update is required it is
    necessary for the header to be handled by the supplied translator
    class or else support automatic translation class determination.
    It is also required that the ``observation_id`` and ``instrument``
    be calculable prior to header fix up.  If a translator class can not
    be found or if there is a problem determining the instrument or
    observation ID, the function will return without action.

    Correction files use names of the form ``instrument-obsid.yaml`` (for
    example ``LATISS-AT_O_20190329_000022.yaml``).
    The YAML file should have the format of:

    .. code-block:: yaml

        EXPTIME: 30.0
        IMGTYPE: bias

    where each key/value pair is copied directly into the supplied header,
    overwriting any previous values.

    This function searches a number of locations for such a correction file.
    The search order is:

    - Any paths explicitly supplied through ``search_path``.
    - The contents of the PATH-like environment variable
      ``$METADATA_CORRECTIONS_PATH``.
    - Any search paths supplied by the matching translator class.

    The first file located in the search path is used for the correction.
    """

    if translator_class is None:
        try:
            translator_class = MetadataTranslator.determine_translator(header,
                                                                       filename=filename)
        except ValueError:
            # if the header is not recognized, we should not complain
            # and should not proceed further.
            return False
    elif not issubclass(translator_class, MetadataTranslator):
        raise TypeError(f"Translator class must be a MetadataTranslator, not {translator_class}")

    # Create an instance for this header
    translator = translator_class(header, filename=filename)

    # To determine the file look up we need the observation_id and instrument
    try:
        obsid = translator.to_observation_id()
        instrument = translator.to_instrument()
    except Exception:
        # Return without comment if these translations failed
        return False

    target_file = f"{instrument}-{obsid}.yaml"
    log.debug("Checking for header correction file named %s", target_file)

    # Work out the search path
    paths = []
    if search_path is not None:
        if isinstance(search_path, str):
            # Allow a single path to be given as a string
            search_path = [search_path]
        paths.extend(search_path)
    if ENV_VAR_NAME in os.environ and os.environ[ENV_VAR_NAME]:
        paths.extend(os.environ[ENV_VAR_NAME].split(os.path.pathsep))

    paths.extend(translator.search_paths())

    # Prioritize file system overrides
    modified = _find_from_file(header, paths, target_file)

    # Apply updates from resources only if none found in files
    if not modified:
        package, resource_root = translator.resource_root()
        modified = _find_from_resource(header, package, resource_root, target_file)

    # Allow a translation class to do local fixups
    # Allow it to fail but log the failure
    try:
        translator_modified = translator_class.fix_header(header, instrument, obsid, filename=filename)
    except Exception as e:
        log.fatal("Ignoring translator header fixup of %s %s: %s",
                  instrument, obsid, e)
        translator_modified = False

    return modified or translator_modified
