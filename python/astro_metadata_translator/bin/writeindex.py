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

__all__ = "write_index_files"

import json
import logging
import os
import sys
from typing import IO, List, MutableMapping, Optional, Sequence, Tuple

from ..file_helpers import find_files
from ..indexing import index_files

log = logging.getLogger(__name__)


def write_index_files(
    files: Sequence[str],
    regex: str,
    hdrnum: int,
    print_trace: bool,
    content_mode: str = "translated",
    outpath: Optional[str] = None,
    outstream: IO = sys.stdout,
    errstream: IO = sys.stderr,
) -> Tuple[List[str], List[str]]:
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
        The HDU number to read. The primary header is always read and merged
        with the specified header.
    print_trace : `bool`
        If there is an error reading the file and this parameter is `True`,
        a full traceback of the exception will be reported. If `False` prints
        a one line summary of the error condition.
    content_mode : `str`
        Form of data to write in index file. Options are:
        ``translated`` (default) to write ObservationInfo to the index;
        ``metadata`` to write native metadata headers to the index.
        The index file is called ``_index.json``
    outpath : `str`, optional
        If specified a single index file will be written to this location
        combining all the information from all files. If `None`, the default,
        and index file will be written to each directory in which files
        are found.
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
    if content_mode not in ("translated", "metadata"):
        raise ValueError(f"Unrecognized content mode {content_mode}")

    if outpath is not None:
        _, ext = os.path.splitext(outpath)
        if ext != ".json":
            raise ValueError(f"Override output file must end in .json but given {outpath}")

    found_files = find_files(files, regex)

    failed = []
    okay = []
    files_per_directory: MutableMapping[str, List[str]] = {}

    # Group each file by directory if no explicit output path
    if outpath is None:
        for path in found_files:
            head, tail = os.path.split(path)
            files_per_directory.setdefault(head, []).append(tail)
    else:
        files_per_directory["."] = list(found_files)

    # Extract translated metadata for each file in each directory
    for directory, files_in_dir in files_per_directory.items():
        output, this_okay, this_failed = index_files(
            files_in_dir, directory, hdrnum, print_trace, content_mode, outstream, errstream
        )

        failed.extend(this_failed)
        okay.extend(this_okay)

        # Write the index file
        if outpath is None:
            index_file = os.path.join(directory, "_index.json")
        else:
            index_file = outpath
        with open(index_file, "w") as fd:
            print(json.dumps(output), file=fd)
            log.info("Wrote index file to %s", index_file)

    return okay, failed
