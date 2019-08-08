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

"""Implementation of the ``translate_header.py`` script.

Read file metadata from the specified files and report the translated content.
"""

__all__ = ("main", "process_files")

import argparse

import os
import re
import sys
import traceback
import importlib
import yaml
from astro_metadata_translator import ObservationInfo, merge_headers
from astro_metadata_translator.tests import read_test_file

# Prefer afw over Astropy
try:
    from lsst.afw.fits import readMetadata
    import lsst.daf.base  # noqa: F401 need PropertyBase for readMetadata

    def read_metadata(file, hdu):
        try:
            return readMetadata(file, hdu=hdu)
        except lsst.afw.fits.FitsError:
            return None

except ImportError:
    from astropy.io import fits

    def read_metadata(file, hdu):
        fits_file = fits.open(file)
        try:
            header = fits_file[hdu].header
        except IndexError:
            header = None
        return header


def build_argparser():
    """Construct an argument parser for the ``translate_header.py`` script.

    Returns
    -------
    argparser : `argparse.ArgumentParser`
        The argument parser that defines the ``translate_header.py``
        command-line interface.
    """

    parser = argparse.ArgumentParser(description="Summarize headers from astronomical data files")
    parser.add_argument("files", metavar="file", type=str, nargs="+",
                        help="File(s) from which headers will be parsed."
                        " If a directory is given it will be scanned for files matching the regular"
                        " expression defined in --regex.")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Do not report the translation content from each header")
    parser.add_argument("-d", "--dumphdr", action="store_true",
                        help="Dump the header in YAML format to standard output rather than translating it")
    parser.add_argument("--traceback", action="store_true",
                        help="Give detailed trace back when any errors encountered")
    parser.add_argument("-n", "--hdrnum", default=1,
                        help="HDU number to read.  If the HDU can not be found, a warning is issued but "
                             "translation is attempted using the primary header.  "
                             "The primary header is always read and merged with this header.")

    re_default = r"\.fit[s]?\b"
    parser.add_argument("-r", "--regex", default=re_default,
                        help="When looking in a directory, regular expression to use to determine whether"
                        f" a file should be examined. Default: '{re_default}'")

    parser.add_argument("-p", "--packages", action="append", type=str,
                        help="Python packages to import to register additional translators")

    return parser


def read_file(file, hdrnum, dumphdr, quiet, print_trace,
              outstream=sys.stdout, errstream=sys.stderr):
    """Read the specified file and process it.

    Parameters
    ----------
    file : `str`
        The file from which the header is to be read.
    hdrnum : `int`
        The HDU number to read. The primary header is always read and
        merged with the header from this HDU.
    dumphdr : `bool`
        If `True` dump the merged header to standard output rather than
        translating it.
    quiet : `bool`
        If `True` do not report the translated values for the file.
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
    success : `bool`
        `True` if the file was handled successfully, `False` if the file
        could not be processed.
    """
    print(f"Analyzing {file}...", file=errstream)
    try:
        if file.endswith(".yaml"):
            md = read_test_file(file,)
            if hdrnum != 0:
                # YAML can't have HDUs
                hdrnum = 0
        else:
            md = read_metadata(file, 0)
        if md is None:
            print(f"Unable to open file {file}", file=errstream)
            return False
        if hdrnum != 0:
            mdn = read_metadata(file, int(hdrnum))
            # Astropy does not allow append mode since it does not
            # convert lists to multiple cards. Overwrite for now
            if mdn is not None:
                md = merge_headers([md, mdn], mode="overwrite")
            else:
                print(f"HDU {hdrnum} was not found. Ignoring request.", file=errstream)

        if dumphdr:
            # The header should be written out in the insertion order
            print(yaml.dump(md, sort_keys=False), file=outstream)
            return True
        obs_info = ObservationInfo(md, pedantic=True, filename=file)
        if not quiet:
            print(f"{obs_info}", file=outstream)
    except Exception as e:
        if print_trace:
            traceback.print_exc(file=outstream)
        else:
            print(repr(e), file=outstream)
        return False
    return True


def process_files(files, regex, hdrnum, dumphdr, quiet, print_trace,
                  outstream=sys.stdout, errstream=sys.stderr):
    """Read and translate metadata from the specified files.

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
    dumphdr : `bool`
        If `True` dump the merged header to standard output rather than
        translating it.
    quiet : `bool`
        If `True` do not report the translated values for the file.
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
    file_regex = re.compile(regex)
    failed = []
    okay = []
    for file in files:
        if os.path.isdir(file):
            for root, dirs, files in os.walk(file):
                for name in files:
                    path = os.path.join(root, name)
                    if os.path.isfile(path) and file_regex.search(name):
                        isok = read_file(path, hdrnum, dumphdr, quiet, print_trace, outstream, errstream)
                        if isok:
                            okay.append(path)
                        else:
                            failed.append(path)
        else:
            isok = read_file(file, hdrnum, dumphdr, quiet, print_trace, outstream, errstream)
            if isok:
                okay.append(file)
            else:
                failed.append(file)

    return okay, failed


def main():
    """Read metadata from the supplied files and translate the content to
    standard form.
    """
    args = build_argparser().parse_args()

    # Process import requests
    if args.packages:
        for m in args.packages:
            importlib.import_module(m)

    # Main loop over files
    _, failed = process_files(args.files, args.regex, args.hdrnum,
                              args.dumphdr, args.quiet, args.traceback)

    if failed:
        print("Files with failed translations:", file=sys.stderr)
        for f in failed:
            print(f"\t{f}", file=sys.stderr)
