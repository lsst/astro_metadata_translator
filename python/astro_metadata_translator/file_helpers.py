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

"""Support functions for script implementations."""

__all__ = ("read_fits_metadata", "find_files", "read_basic_metadata_from_file", "read_file_info")

import re
import os
import sys
import traceback

from astro_metadata_translator import merge_headers, ObservationInfo
from .tests import read_test_file


# Prefer afw over Astropy
try:
    from lsst.afw.fits import readMetadata
    import lsst.daf.base  # noqa: F401 need PropertyBase for readMetadata

    def read_fits_metadata(file, hdu, can_raise=False):
        """Read a FITS header using afw.

        Parameters
        ----------
        file : `str`
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
        try:
            return readMetadata(file, hdu=hdu)
        except lsst.afw.fits.FitsError as e:
            if can_raise:
                # Try to convert a basic fits error code
                if "(104)" in str(e):
                    raise FileNotFoundError(f"No such file or directory: {file}") from e
                raise e
            return None

except ImportError:
    from astropy.io import fits

    def read_fits_metadata(file, hdu, can_raise=False):
        """Read a FITS header using astropy.

        Parameters
        ----------
        file : `str`
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
        fits_file = fits.open(file)
        try:
            header = fits_file[hdu].header
        except IndexError as e:
            if can_raise:
                raise e
            header = None
        return header


def find_files(files, regex):
    """Find files for processing.

    Parameters
    ----------
    files : iterable of `str`
        The files or directories from which the headers are to be read.
    regex : `str`
        Regular expression string used to filter files when a directory is
        scanned.
    """
    file_regex = re.compile(regex)
    found_files = []

    # Find all the files of interest
    for file in files:
        if os.path.isdir(file):
            for root, dirs, files in os.walk(file):
                for name in files:
                    path = os.path.join(root, name)
                    if os.path.isfile(path) and file_regex.search(name):
                        found_files.append(path)
        else:
            found_files.append(file)

    return found_files


def read_basic_metadata_from_file(file, hdrnum, errstream=sys.stderr, can_raise=False):
    """Read a raw header from a file, merging if necessary

    Parameters
    ----------
    file : `str`
        Name of file to read. Can be FITS or YAML. YAML must be a simple
        top-level dict.
    hdrnum : `int`
        Header number to read. Only relevant for FITS. If greater than 1
        it will be merged with the primary header.
    errstream : `io.StringIO`, optional
        Stream to send messages that would normally be sent to standard
        error. Defaults to `sys.stderr`. Only used if exceptions are disabled.
    can_raise : `bool`, optional
        Indicate whether the function can raise and exception (default)
        or should return `None` on error. Can still raise if an unexpected
        error is encountered.

    Returns
    -------
    header : `dict`
        The header as a dict. Can be `None` if there was a problem reading
        the file.
    """
    if file.endswith(".yaml"):
        md = read_test_file(file,)
        if hdrnum != 0:
            # YAML can't have HDUs
            hdrnum = 0
    else:
        md = read_fits_metadata(file, 0, can_raise=can_raise)
    if md is None:
        print(f"Unable to open file {file}", file=errstream)
        return None
    if hdrnum != 0:
        mdn = read_fits_metadata(file, int(hdrnum), can_raise=can_raise)
        # Astropy does not allow append mode since it does not
        # convert lists to multiple cards. Overwrite for now
        if mdn is not None:
            md = merge_headers([md, mdn], mode="overwrite")
        else:
            print(f"HDU {hdrnum} was not found in file {file}. Ignoring request.", file=errstream)

    return md


def read_file_info(file, hdrnum, print_trace=None, mode="simple",
                   outstream=sys.stdout, errstream=sys.stderr):
    """Read information from file

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
    mode : `str`
        Type of content to return. Options are:
        ``simple`` returns the simplified form of an ObservationInfo (default);
        ``obsInfo`` returns an ObsrvationInfo;
        ``json`` returns a JSON string of the ObservationInfo;
        ``metadata`` returns the metadata (will be unfixed form);.
    outstream : `io.StringIO`, optional
        Output stream to use for standard messages. Defaults to `sys.stdout`.
    errstream : `io.StringIO`, optional
        Stream to send messages that would normally be sent to standard
        error. Defaults to `sys.stderr`.

    Returns
    -------
    simple : `dict` of `str`
        The return value of `ObservationInfo.to_simple()`. Returns `None`
        if there was a problem and `print_trace` is not `None`.
    """

    if mode not in ("metadata", "obsInfo", "simple", "json"):
        raise ValueError(f"Unrecognized mode {mode}")

    try:
        # Calculate the JSON from the file
        md = read_basic_metadata_from_file(file, hdrnum, errstream=errstream,
                                           can_raise=True if print_trace is None else False)
        if md is None:
            return None
        if mode == "metadata":
            # Do not fix the header
            return md
        obs_info = ObservationInfo(md, pedantic=True, filename=file)
        if mode == "obsInfo":
            return obs_info
        if mode == "simple":
            return obs_info.to_simple()
        if mode == "json":
            return obs_info.to_json()
        raise RuntimeError(f"Logic error. Unrecognized mode for reading file: {mode}")
    except Exception as e:
        if print_trace is None:
            raise e
        if print_trace:
            traceback.print_exc(file=outstream)
        else:
            print(repr(e), file=outstream)
    return None
