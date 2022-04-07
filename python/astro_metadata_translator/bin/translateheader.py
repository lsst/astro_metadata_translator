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

from __future__ import annotations

__all__ = ("main", "process_files")

import argparse
import importlib
import logging
import sys
import traceback
from typing import IO, List, Sequence, Tuple

import yaml

from astro_metadata_translator import MetadataTranslator, ObservationInfo, fix_header

from ..file_helpers import find_files, read_basic_metadata_from_file

# Output mode choices
OUTPUT_MODES = ("auto", "verbose", "table", "yaml", "fixed", "yamlnative", "fixednative", "none")

# Definitions for table columns
TABLE_COLUMNS = (
    {"format": "32.32s", "attr": "observation_id", "label": "ObsId"},
    {
        "format": "8.8s",
        "attr": "observation_type",
        "label": "ImgType",
    },
    {
        "format": "16.16s",
        "attr": "object",
        "label": "Object",
    },
    {
        "format": "16.16s",
        "attr": "physical_filter",
        "label": "Filter",
    },
    {"format": ">8.8s", "attr": "detector_unique_name", "label": "Detector"},
    {
        "format": "5.1f",
        "attr": "exposure_time",
        "label": "ExpTime",
    },
)


def build_argparser() -> argparse.ArgumentParser:
    """Construct an argument parser for the ``translate_header.py`` script.

    Returns
    -------
    argparser : `argparse.ArgumentParser`
        The argument parser that defines the ``translate_header.py``
        command-line interface.
    """

    parser = argparse.ArgumentParser(description="Summarize headers from astronomical data files")
    parser.add_argument(
        "files",
        metavar="file",
        type=str,
        nargs="+",
        help="File(s) from which headers will be parsed."
        " If a directory is given it will be scanned for files matching the regular"
        " expression defined in --regex.",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Do not report the translation content from each header. This forces output mode 'none'.",
    )
    parser.add_argument(
        "-d",
        "--dumphdr",
        action="store_true",
        help="Dump the header in YAML format to standard output rather than translating it."
        " This is the same as using mode=yaml",
    )
    parser.add_argument(
        "--traceback", action="store_true", help="Give detailed trace back when any errors encountered"
    )
    parser.add_argument(
        "-n",
        "--hdrnum",
        default=1,
        help="HDU number to read.  If the HDU can not be found, a warning is issued but "
        "translation is attempted using the primary header.  "
        "The primary header is always read and merged with this header.",
    )
    parser.add_argument(
        "-m",
        "--mode",
        default="auto",
        choices=OUTPUT_MODES,
        help="Display mode for translated parameters. 'verbose' displays all the information"
        " available. 'table' displays important information in tabular form."
        " 'yaml' dumps the header in YAML format (this is equivalent to -d option)."
        " 'fixed' dumps the header in YAML after it has had corrections applied."
        " Add 'native' suffix to dump YAML in PropertyList or Astropy native form."
        " 'none' displays no translated header information and is an alias for the "
        " '--quiet' option."
        " 'auto' mode is 'verbose' for a single file and 'table' for multiple files.",
    )
    parser.add_argument("-l", "--log", default="warn", help="Python logging level to use.")

    re_default = r"\.fit[s]?\b"
    parser.add_argument(
        "-r",
        "--regex",
        default=re_default,
        help="When looking in a directory, regular expression to use to determine whether"
        f" a file should be examined. Default: '{re_default}'",
    )

    parser.add_argument(
        "-p",
        "--packages",
        action="append",
        type=str,
        help="Python packages to import to register additional translators",
    )

    return parser


def read_file(
    file: str,
    hdrnum: int,
    print_trace: bool,
    outstream: IO = sys.stdout,
    errstream: IO = sys.stderr,
    output_mode: str = "verbose",
    write_heading: bool = False,
) -> bool:
    """Read the specified file and process it.

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
    output_mode : `str`, optional
        Output mode to use. Must be one of "verbose", "none", "table",
        "yaml", or "fixed".  "yaml" and "fixed" can be modified with a
        "native" suffix to indicate that the output should be a representation
        of the native object type representing the header (which can be
        PropertyList or an Astropy header).  Without this modify headers
        will be dumped as simple `dict` form.
        "auto" is used to indicate that a single file has been specified
        but the output will depend on whether the file is a multi-extension
        FITS file or not.
    write_heading: `bool`, optional
        If `True` and in table mode, write a table heading out before writing
        the content.

    Returns
    -------
    success : `bool`
        `True` if the file was handled successfully, `False` if the file
        could not be processed.
    """
    if output_mode not in OUTPUT_MODES:
        raise ValueError(f"Output mode of '{output_mode}' is not understood.")

    # This gets in the way in tabular mode
    if output_mode != "table":
        print(f"Analyzing {file}...", file=errstream)

    try:
        md = read_basic_metadata_from_file(file, hdrnum, errstream=errstream, can_raise=True)
        if md is None:
            raise RuntimeError(f"Failed to read file {file} HDU={hdrnum}")

        if output_mode.endswith("native"):
            # Strip native and don't change type of md
            output_mode = output_mode[: -len("native")]
        else:
            # Rewrite md as simple dict for output
            md = {k: v for k, v in md.items()}

        if output_mode in ("yaml", "fixed"):

            if output_mode == "fixed":
                fix_header(md, filename=file)

            # The header should be written out in the insertion order
            print(yaml.dump(md, sort_keys=False), file=outstream)
            return True

        # Try to work out a translator class.
        translator_class = MetadataTranslator.determine_translator(md, filename=file)

        # Work out which headers to translate, assuming the default if
        # we have a YAML test file.
        if file.endswith(".yaml"):
            headers = [md]
        else:
            headers = list(translator_class.determine_translatable_headers(file, md))
        if output_mode == "auto":
            output_mode = "table" if len(headers) > 1 else "verbose"

        wrote_heading = False
        for md in headers:
            obs_info = ObservationInfo(md, pedantic=True, filename=file)
            if output_mode == "table":
                columns = [
                    "{:{fmt}}".format(getattr(obs_info, c["attr"]), fmt=c["format"]) for c in TABLE_COLUMNS
                ]

                if write_heading and not wrote_heading:
                    # Construct headings of the same width as the items
                    # we have calculated.  Doing this means we don't have to
                    # work out for ourselves how many characters will be used
                    # for non-strings (especially Quantity)
                    headings = []
                    separators = []
                    for thiscol, defn in zip(columns, TABLE_COLUMNS):
                        width = len(thiscol)
                        headings.append("{:{w}.{w}}".format(defn["label"], w=width))
                        separators.append("-" * width)
                    print(" ".join(headings), file=outstream)
                    print(" ".join(separators), file=outstream)
                    wrote_heading = True

                row = " ".join(columns)
                print(row, file=outstream)
            elif output_mode == "verbose":
                print(f"{obs_info}", file=outstream)
            elif output_mode == "none":
                pass
            else:
                raise RuntimeError(f"Output mode of '{output_mode}' not recognized but should be known.")
    except Exception as e:
        if print_trace:
            traceback.print_exc(file=outstream)
        else:
            print(f"Failure processing {file}: {e}", file=outstream)
        return False
    return True


def process_files(
    files: Sequence[str],
    regex: str,
    hdrnum: int,
    print_trace: bool,
    outstream: IO = sys.stdout,
    errstream: IO = sys.stderr,
    output_mode: str = "auto",
) -> Tuple[List[str], List[str]]:
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
    print_trace : `bool`
        If there is an error reading the file and this parameter is `True`,
        a full traceback of the exception will be reported. If `False` prints
        a one line summary of the error condition.
    outstream : `io.StringIO`, optional
        Output stream to use for standard messages. Defaults to `sys.stdout`.
    errstream : `io.StringIO`, optional
        Stream to send messages that would normally be sent to standard
        error. Defaults to `sys.stderr`.
    output_mode : `str`, optional
        Output mode to use for the translated information.
        "auto" switches based on how many files are found.

    Returns
    -------
    okay : `list` of `str`
        All the files that were processed successfully.
    failed : `list` of `str`
        All the files that could not be processed.
    """
    found_files = find_files(files, regex)

    # Convert "auto" to correct mode but for a single file keep it
    # auto in case that file has multiple headers
    if output_mode == "auto":
        if len(found_files) > 1:
            output_mode = "table"

    # Process each file
    failed = []
    okay = []
    heading = True
    for path in sorted(found_files):
        isok = read_file(path, hdrnum, print_trace, outstream, errstream, output_mode, heading)
        heading = False
        if isok:
            okay.append(path)
        else:
            failed.append(path)

    return okay, failed


def main() -> int:
    """Read metadata from the supplied files and translate the content to
    standard form.

    Returns
    -------
    status : `int`
        Exit status to be passed to `sys.exit()`. 0 if any of the files
        could be translated. 1 otherwise.
    """

    logging.warn(
        "This command is deprecated. Please use 'astrometadata translate' "
        " or 'astrometadata dump' instead. See 'astrometadata -h' for more details."
    )

    args = build_argparser().parse_args()

    # Process import requests
    if args.packages:
        for m in args.packages:
            importlib.import_module(m)

    output_mode = args.mode
    if args.quiet:
        output_mode = "none"
    elif args.dumphdr:
        output_mode = "yaml"

    # Set the log level. Convert to upper case to allow the user to
    # specify --log=DEBUG or --log=debug
    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {args.log}")
    logging.basicConfig(level=numeric_level)

    # Main loop over files
    okay, failed = process_files(args.files, args.regex, args.hdrnum, args.traceback, output_mode=output_mode)

    if failed:
        print("Files with failed translations:", file=sys.stderr)
        for f in failed:
            print(f"\t{f}", file=sys.stderr)

    if okay:
        # Good status if anything was returned in okay
        return 0
    else:
        return 1
