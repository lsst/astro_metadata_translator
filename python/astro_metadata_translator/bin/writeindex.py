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

__all__ = ("write_index_files")

import logging
import json
import os
import sys

from ..file_helpers import find_files
from ..indexing import index_files

log = logging.getLogger(__name__)


def write_index_files(files, regex, hdrnum, print_trace, mode="obsInfo",
                      outstream=sys.stdout, errstream=sys.stderr):
    """Process each file and create JSON index file.

    The index file will have common information in the toplevel.
    There is then a ``__DIFF__`` key that is a dictionary with file
    names as keys and per-file differences as the values in a dict.

    Parameters
    ----------
    files : iterable of `str`
        The files or directories from which the headers are to be read.
    regex : `str`
        Regular expression string used to filter files when a directory is
        scanned.
    hdrnum : `int`
        The HDU number to read. The primary header is always read and
    print_trace : `bool`
        If there is an error reading the file and this parameter is `True`,
        a full traceback of the exception will be reported. If `False` prints
        a one line summary of the error condition.
    mode : `str`
        Form of data to write in index file. Options are:
        ``obsInfo`` (default) to write ObservationInfo to the index;
        ``metadata`` to write native metadata headers to the index.
        The index file is called ``{mode}_index.json``
    outstream : `io.StringIO`, optional
        Output stream to use for standard messages. Defaults to `sys.stdout`.
    errstream : `io.StringIO`, optional
        Stream to send messages that would normally be sent to standard
        error. Defaults to `sys.stderr`.

    Returns
    -------
    okay : `list` of `str`
        All the files that were processed successfully.
    failed : `list` of `str`
        All the files that could not be processed.
    """
    if mode not in ("obsInfo", "metadata"):
        raise ValueError("Unrecognized mode {mode}")

    found_files = find_files(files, regex)

    failed = []
    okay = []
    by_directory = {}

    # Group each file by directory
    for path in found_files:
        head, tail = os.path.split(path)
        by_directory.setdefault(head, []).append(tail)

    # Extract translated metadata for each file in each directory
    for directory, files_in_dir in by_directory.items():
        output, this_okay, this_failed = index_files(files_in_dir, directory, hdrnum, print_trace, mode,
                                                     outstream, errstream)

        failed.extend(this_failed)
        okay.extend(this_okay)

        # Write the index file
        with open(outfile := os.path.join(directory, f"{mode}_index.json"), "w") as fd:
            print(json.dumps(output), file=fd)
            log.info("Wrote index file to %s", outfile)

    return okay, failed
