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

__all__ = ("write_sidecar_files", "write_sidecar_file")

import os
import sys
import traceback

from ..file_helpers import find_files, read_file_info


def write_sidecar_file(file, hdrnum, content_mode, print_trace, outstream=sys.stdout, errstream=sys.stderr):
    """Write JSON summary to sidecar file.

    Parameters
    ----------
    file : `str`
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
        Output stream to use for standard messages. Defaults to `sys.stdout`.
    errstream : `io.StringIO`, optional
        Stream to send messages that would normally be sent to standard
        error. Defaults to `sys.stderr`.

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
        json_str = read_file_info(file, hdrnum, content_mode=content_mode, content_type="json",
                                  print_trace=print_trace, outstream=outstream, errstream=errstream)
        if json_str is None:
            return False

        # Calculate file name derived on this file
        root, ext = os.path.splitext(file)
        newfile = root + ".json"
        with open(newfile, "w") as fd:
            print(json_str, file=fd)

    except Exception as e:
        if print_trace is None:
            raise e
        if print_trace:
            traceback.print_exc(file=outstream)
        else:
            print(repr(e), file=outstream)
        return False
    return True


def write_sidecar_files(files, regex, hdrnum, content_mode, print_trace,
                        outstream=sys.stdout, errstream=sys.stderr):
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
    content_mode : `str`
        Content mode to use when writing JSON to sidecar. Options are:
        ``metadata`` to write the unmodified header;
        ``translated`` to write the translated ObservationInfo.
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

    # Process each file
    failed = []
    okay = []
    for path in sorted(found_files):
        isok = write_sidecar_file(path, hdrnum, content_mode, print_trace, outstream, errstream)
        if isok:
            okay.append(path)
        else:
            failed.append(path)

    return okay, failed
