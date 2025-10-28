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

from __future__ import annotations

__all__ = ("write_sidecar_file", "write_sidecar_files")

import logging
import traceback
from collections.abc import Sequence
from typing import IO, TYPE_CHECKING

from lsst.resources import ResourcePath

from ..file_helpers import find_files, read_file_info

if TYPE_CHECKING:
    from lsst.resources import ResourcePathExpression

log = logging.getLogger(__name__)


def write_sidecar_file(
    file: ResourcePathExpression,
    hdrnum: int,
    content_mode: str,
    print_trace: bool,
    outstream: IO | None = None,
) -> bool:
    """Write JSON summary to sidecar file.

    Parameters
    ----------
    file : `str` or `lsst.resources.ResourcePathExpression`
        The file from which the header is to be read.
    hdrnum : `int`
        The HDU number to read. The primary header is always read and
        merged with the header from this HDU.
    content_mode : `str`
        Content mode to use when writing JSON to sidecar. Options are:
        ``metadata`` to write the unmodified header;
        ``translated`` to write the translated ObservationInfo.
    print_trace : `bool`
        If there is an error reading the file and this parameter is `True`,
        a full traceback of the exception will be reported. If `False` prints
        a one line summary of the error condition. If `None` the exception
        will be allowed.
    outstream : `io.StringIO`, optional
        Output stream to use for standard messages. Defaults to `None` which
        uses the default output stream.

    Returns
    -------
    success : `bool`
        `True` if the file was handled successfully, `False` if the file
        could not be processed.
    """
    if content_mode not in ("metadata", "translated"):
        raise ValueError(f"Specified content mode '{content_mode}' is not understood.")

    try:
        # Calculate the JSON from the file
        json_str = read_file_info(
            file,
            hdrnum,
            content_mode=content_mode,
            content_type="json",
            print_trace=print_trace,
            outstream=outstream,
        )
        if json_str is None:
            return False

        # Calculate sidecar file name derived from this file.
        # .fits.gz should be replaced
        # with .json, and not resulting in .fits.json.
        uri = ResourcePath(file)
        newfile = uri.updatedExtension(".json")
        newfile.write(str(json_str).encode())
        log.debug("Writing sidecar file %s", newfile)

    except Exception as e:
        if print_trace is None:
            raise e
        if print_trace:
            traceback.print_exc(file=outstream)
        else:
            print(repr(e), file=outstream)
        return False
    return True


def write_sidecar_files(
    files: Sequence[str],
    regex: str,
    hdrnum: int,
    content_mode: str,
    print_trace: bool,
    outstream: IO | None = None,
) -> tuple[list[str], list[str]]:
    """Process each file and create sidecar file.

    Parameters
    ----------
    files : iterable of `str`
        The files or directories from which the headers are to be read.
    regex : `str`
        Regular expression string used to filter files when a directory is
        scanned.
    hdrnum : `int`
        The HDU number to read. The primary header is always read and
        merged with the header from this HDU.
    content_mode : `str`
        Content mode to use when writing JSON to sidecar. Options are:
        ``metadata`` to write the unmodified header;
        ``translated`` to write the translated ObservationInfo.
    print_trace : `bool`
        If there is an error reading the file and this parameter is `True`,
        a full traceback of the exception will be reported. If `False` prints
        a one line summary of the error condition.
    outstream : `io.StringIO`, optional
        Output stream to use for standard messages. Defaults to `None` which
        uses the default output stream.

    Returns
    -------
    okay : `list` of `str`
        All the files that were processed successfully.
    failed : `list` of `str`
        All the files that could not be processed.
    """
    found_files = find_files(files, regex)

    # Process each file
    failed = []
    okay = []
    for path in sorted(found_files):
        isok = write_sidecar_file(path, hdrnum, content_mode, print_trace, outstream)
        if isok:
            okay.append(str(path))
        else:
            failed.append(str(path))

    return okay, failed
