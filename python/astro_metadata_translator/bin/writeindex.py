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
import traceback

from astro_metadata_translator import ObservationInfo, merge_headers, fix_header

from .helper import find_files, read_metadata_from_file

log = logging.getLogger(__name__)


def read_file(file, hdrnum, print_trace, mode="simple", outstream=sys.stdout, errstream=sys.stderr):
    """Read information from file

    Parameters
    ----------
    file : `str`
        The file from which the header is to be read.
    hdrnum : `int`
        The HDU number to read. The primary header is always read and
        merged with the header from this HDU.
    print_trace : `bool`
        If there is an error reading the file and this parameter is `True`,
        a full traceback of the exception will be reported. If `False` prints
        a one line summary of the error condition.
    mode : `str`
        Type of content to return. Options are:
        ``simple`` returns the simplified form of an ObservationInfo (default);
        ``obsInfo`` returns an ObsrvationInfo;
        ``metadata`` returns the metadata (will be fixed form).
    outstream : `io.StringIO`, optional
        Output stream to use for standard messages. Defaults to `sys.stdout`.
    errstream : `io.StringIO`, optional
        Stream to send messages that would normally be sent to standard
        error. Defaults to `sys.stderr`.

    Returns
    -------
    simple : `dict` of `str`
        The return value of `ObservationInfo.to_simple()`.
    """

    if mode not in ("metadata", "obsInfo", "simple"):
        raise ValueError(f"Unrecognized mode {mode}")

    try:
        # Calculate the JSON from the file
        md = read_metadata_from_file(file, hdrnum, errstream=errstream)
        if md is None:
            return None
        if mode == "metadata":
            fix_header(md)
            return md
        obs_info = ObservationInfo(md, pedantic=True, filename=file)
        if mode == "obsInfo":
            return obs_info
        return obs_info.to_simple()
    except Exception as e:
        if print_trace:
            traceback.print_exc(file=outstream)
        else:
            print(repr(e), file=outstream)
    return None


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

    # We want the reader to return simple dict to us here
    read_mode = mode
    if read_mode == "obsInfo":
        read_mode = "simple"

    # Extract translated metadata for each file in each directory
    for directory, files_in_dir in by_directory.items():
        by_file = {}
        for file in sorted(files_in_dir):
            path = os.path.join(directory, file)
            simple = read_file(path, hdrnum, print_trace, read_mode, outstream, errstream)
            if simple is None:
                failed.append(path)
                continue
            else:
                okay.append(path)

            # Store the information indexed by the filename within dir
            by_file[file] = simple

        # Merge all the information into a primary plus diff
        merged = merge_headers(by_file.values(), mode="diff")

        # The structure to write to file is intended to look like (in YAML):
        # __COMMON__:
        #   KEY1: value1
        #   KEY2: value2
        # FILE1:
        #   KEY3: value3a
        # FILE2:
        #   KEY3: value3b

        # if there was only one file there will not be a diff but we
        # want it to look like there was.
        diff_dict = merged.pop("__DIFF__", [dict()])

        # Put the common headers first in the output.
        output = {"__COMMON__": merged}
        for file, diff in zip(by_file, diff_dict):
            output[file] = diff

        # Write the index file
        with open(outfile := os.path.join(directory, f"{mode}_index.json"), "w") as fd:
            print(json.dumps(output), file=fd)
            log.info("Wrote index file to %s", outfile)

    return okay, failed
