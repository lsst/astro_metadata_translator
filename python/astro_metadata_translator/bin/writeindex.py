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

from astro_metadata_translator import ObservationInfo, merge_headers

from .helper import find_files, read_metadata_from_file

log = logging.getLogger(__name__)


def read_simple(file, hdrnum, print_trace, outstream=sys.stdout, errstream=sys.stderr):
    """Read simplified translated information from file

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

    try:
        # Calculate the JSON from the file
        md = read_metadata_from_file(file, hdrnum, errstream=errstream)
        if md is None:
            return None
        obs_info = ObservationInfo(md, pedantic=True, filename=file)
        return obs_info.to_simple()
    except Exception as e:
        if print_trace:
            traceback.print_exc(file=outstream)
        else:
            print(repr(e), file=outstream)
    return None


def write_index_files(files, regex, hdrnum, print_trace,
                      outstream=sys.stdout, errstream=sys.stderr):
    """Process each file and create JSON index file.

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
        by_file = {}
        for file in files_in_dir:
            path = os.path.join(directory, file)
            simple = read_simple(path, hdrnum, print_trace, outstream, errstream)
            if simple is None:
                failed.append(path)
                continue
            else:
                okay.append(path)

            # Store the information indexed by the filename within dir
            by_file[file] = simple

        # Merge all the information into a primary plus diff
        merged = merge_headers(by_file.values(), mode="diff")

        # Convert the diff into a dict indexed by file name
        diff_dict = {}
        for file, diff in zip(by_file, merged["__DIFF__"]):
            diff_dict[file] = diff
        merged["__DIFF__"] = diff_dict

        # Write the index file
        with open(outfile := os.path.join(directory, "obsinfo_index.json"), "w") as fd:
            print(json.dumps(merged), file=fd)
            log.info("Wrote index file to %s", outfile)

    return okay, failed
