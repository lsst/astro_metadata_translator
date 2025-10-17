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

"""Support functions for script implementations.

These functions should not be treated as part of the public API.
"""

from __future__ import annotations

__all__ = ("find_files", "read_basic_metadata_from_file", "read_file_info")

import json
import logging
import re
import traceback
from collections.abc import Iterable, MutableMapping
from typing import IO, TYPE_CHECKING, Any

from astropy.io import fits
from lsst.resources import ResourcePath

from .headers import merge_headers
from .observationInfo import ObservationInfo
from .tests import read_test_file

if TYPE_CHECKING:
    from lsst.resources import ResourcePathExpression

log = logging.getLogger(__name__)

# Prefer afw over Astropy
try:
    import lsst.daf.base  # noqa: F401 need PropertyBase for readMetadata
    from lsst.afw.fits import FitsError, readMetadata

    have_afw = True

    def _read_fits_metadata_afw(
        file: str, hdu: int, can_raise: bool = False
    ) -> MutableMapping[str, Any] | None:
        # Only works with local files.
        # Tries to catch a FitsError 104 and convert to `FileNotFoundError`.
        # For detailed docstrings see the _read_fits_metadata implementation
        # below.
        try:
            return readMetadata(file, hdu=hdu)
        except FitsError as e:
            if can_raise:
                # Try to convert a basic fits error code
                if "(104)" in str(e):
                    raise FileNotFoundError(f"No such file or directory: {file}") from e
                raise e
            return None

except ImportError:
    have_afw = False


def _read_fits_metadata_astropy(
    file: ResourcePathExpression, hdu: int, can_raise: bool = False
) -> MutableMapping[str, Any] | None:
    """Read a FITS header using astropy.

    Parameters
    ----------
    file : `str` or `lsst.resources.ResourcePathExpression`
        The file to read.
    hdu : `int`
        The header number to read.
    can_raise : `bool`, optional
        Indicate whether the function can raise and exception (default)
        or should return `None` on error. Can still raise if an unexpected
        error is encountered.

    Returns
    -------
    md : `dict`
        The requested header. `None` if it could not be read and
        ``can_raise`` is `False`.
    """
    header = None
    uri = ResourcePath(file, forceDirectory=False)
    try:
        fs, fspath = uri.to_fsspec()
        with fs.open(fspath) as f, fits.open(f) as fits_file:
            try:
                # Copy forces a download of the remote resource.
                header = fits_file[hdu].header.copy()
            except IndexError as e:
                if can_raise:
                    raise e
    except Exception as e:
        if can_raise:
            raise e
    return header


def _read_fits_metadata(
    file: ResourcePathExpression, hdu: int, can_raise: bool = False
) -> MutableMapping[str, Any] | None:
    """Read a FITS header using afw or astropy.

    Prefer afw for local reads if available.

    Parameters
    ----------
    file : `str` or `lsst.resources.ResourcePathExpression`
        The file to read.
    hdu : `int`
        The header number to read.
    can_raise : `bool`, optional
        Indicate whether the function can raise and exception (default)
        or should return `None` on error. Can still raise if an unexpected
        error is encountered.

    Returns
    -------
    md : `dict`
        The requested header. `None` if it could not be read and
        ``can_raise`` is `False`.
    """
    uri = ResourcePath(file, forceAbsolute=False)
    if have_afw and uri.isLocal:
        return _read_fits_metadata_afw(uri.ospath, hdu, can_raise=can_raise)
    return _read_fits_metadata_astropy(uri, hdu, can_raise=can_raise)


def find_files(files: Iterable[ResourcePathExpression], regex: str) -> list[ResourcePath]:
    """Find files for processing.

    Parameters
    ----------
    files : iterable of `lsst.resources.ResourcePathExpression`
        The files or directories from which the headers are to be read.
    regex : `str`
        Regular expression string used to filter files when a directory is
        scanned.

    Returns
    -------
    found_files : `list` of `lsst.resources.ResourcePath`
        The files that were found.
    """
    file_regex = re.compile(regex)
    found_files: list[ResourcePath] = []

    # Find all the files of interest
    for candidate in files:
        uri = ResourcePath(candidate, forceAbsolute=False)
        if uri.isdir():
            found_files.extend(ResourcePath.findFileResources([uri], file_filter=file_regex, grouped=False))
        else:
            found_files.append(uri)

    return found_files


def read_basic_metadata_from_file(
    file: ResourcePathExpression, hdrnum: int, can_raise: bool = True
) -> MutableMapping[str, Any] | None:
    """Read a raw header from a file, merging if necessary.

    Parameters
    ----------
    file : `str` or `lsst.resources.ResourcePathExpression`
        Name of file to read. Can be FITS, YAML or JSON. YAML or JSON must be
        a simple top-level dict.
    hdrnum : `int`
        Header number to read. Only relevant for FITS. If greater than 1
        it will be merged with the primary header. If a negative number is
        given the second header, if present, will be merged with the primary
        header. If there is only a primary header a negative number behaves
        identically to specifying 0 for the HDU number.
    can_raise : `bool`, optional
        Indicate whether the function can raise an exception (default)
        or should return `None` on error. Can still raise if an unexpected
        error is encountered.

    Returns
    -------
    header : `dict`
        The header as a dict. Can be `None` if there was a problem reading
        the file.
    """
    uri = ResourcePath(file, forceAbsolute=False)
    if uri.getExtension() in (".yaml", ".json"):
        try:
            md = read_test_file(
                uri,
            )
        except Exception as e:
            if not can_raise:
                md = None
            else:
                raise e
        if hdrnum != 0:
            # YAML can't have HDUs so skip merging below
            hdrnum = 0
    else:
        md = _read_fits_metadata(uri, 0, can_raise=can_raise)
    if md is None:
        log.warning("Unable to open file %s", file)
        return None
    if hdrnum < 0:
        if "EXTEND" in md and md["EXTEND"]:
            hdrnum = 1
    if hdrnum > 0:
        # Allow this to fail
        mdn = _read_fits_metadata(uri, int(hdrnum), can_raise=False)
        # Astropy does not allow append mode since it does not
        # convert lists to multiple cards. Overwrite for now
        if mdn is not None:
            md = merge_headers([md, mdn], mode="overwrite")
        else:
            log.warning("HDU %d was not found in file %s. Ignoring request.", hdrnum, uri)

    return md


def read_file_info(
    file: ResourcePathExpression,
    hdrnum: int,
    print_trace: bool | None = None,
    content_mode: str = "translated",
    content_type: str = "simple",
    outstream: IO | None = None,
) -> str | MutableMapping[str, Any] | ObservationInfo | None:
    """Read information from file.

    Parameters
    ----------
    file : `str`
        The file from which the header is to be read.
    hdrnum : `int`
        The HDU number to read. The primary header is always read and
        merged with the header from this HDU.
    print_trace : `bool` or `None`
        If there is an error reading the file and this parameter is `True`,
        a full traceback of the exception will be reported. If `False` prints
        a one line summary of the error condition. If `None` the exception
        will be allowed to propagate.
    content_mode : `str`
        Content returned. This can be: ``metadata`` to return the unfixed
        metadata headers; ``translated`` to return the output from metadata
        translation.
    content_type : `str`, optional
        Form of content to be returned. Can be either ``json`` to return a
        JSON string, ``simple`` to always return a `dict`, or ``native`` to
        return either a `dict` (for ``metadata``) or `.ObservationInfo` for
        ``translated``.
    outstream : `io.StringIO` or `None`, optional
        Output stream to use for standard messages. Defaults to `None` which
        uses the default output stream.

    Returns
    -------
    simple : `dict` of `str` or `.ObservationInfo`
        The return value of `.ObservationInfo.to_simple`. Returns `None`
        if there was a problem and ``print_trace`` is not `None`.
    """
    if content_mode not in ("metadata", "translated"):
        raise ValueError(f"Unrecognized content mode request: {content_mode}")

    if content_type not in ("native", "simple", "json"):
        raise ValueError(f"Unrecognized content type request {content_type}")

    uri = ResourcePath(file, forceAbsolute=False)
    try:
        # Calculate the JSON from the file
        md = read_basic_metadata_from_file(uri, hdrnum, can_raise=True if print_trace is None else False)
        if md is None:
            return None
        if content_mode == "metadata":
            # Do not fix the header
            if content_type == "json":
                # Add a key to tell the reader whether this is md or translated
                md["__CONTENT__"] = content_mode
                try:
                    json_str = json.dumps(md)
                except TypeError:
                    # Cast to dict and try again -- PropertyList is a problem
                    json_str = json.dumps(dict(md))
                return json_str
            return md
        obs_info = ObservationInfo(md, pedantic=True, filename=str(uri))
        if content_type == "native":
            return obs_info
        simple = obs_info.to_simple()
        if content_type == "simple":
            return simple
        if content_type == "json":
            # Add a key to tell the reader if this is metadata or translated
            simple["__CONTENT__"] = content_mode
            return json.dumps(simple)
        raise RuntimeError(f"Logic error. Unrecognized mode for reading file: {content_mode}/{content_type}")
    except Exception as e:
        if print_trace is None:
            raise e
        if print_trace:
            traceback.print_exc(file=outstream)
        else:
            print(repr(e), file=outstream)
    return None
