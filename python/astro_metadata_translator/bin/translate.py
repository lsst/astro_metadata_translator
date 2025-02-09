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

"""Implementation of the ``astrometadata translate`` command-line.

Read file metadata from the specified files and report the translated content.
"""

from __future__ import annotations

__all__ = ("translate_or_dump_headers",)

import logging
import math
import traceback
from collections import defaultdict
from collections.abc import Sequence
from typing import IO, Any

import astropy.time
import astropy.units as u
import yaml
from astropy.table import Column, MaskedColumn, QTable

from astro_metadata_translator import MetadataTranslator, ObservationInfo, fix_header

from ..file_helpers import find_files, read_basic_metadata_from_file
from ..properties import PROPERTIES

log = logging.getLogger(__name__)

# Output mode choices
OUTPUT_MODES = ("auto", "verbose", "table", "yaml", "fixed", "yamlnative", "fixednative", "none")

# Number of rows per table page.
# This is a minimum given that DECam data files can include multiple headers.
MAX_TABLE_PAGE_SIZE = 500

# Definitions for table columns
TABLE_COLUMNS = (
    {"format": "<", "attr": "observation_id", "label": "ObsId"},
    {
        "format": "<",
        "attr": "observation_type",
        "label": "ImgType",
    },
    {
        "format": "<",
        "attr": "object",
        "label": "Object",
    },
    {
        "format": "<",
        "attr": "physical_filter",
        "label": "Filter",
    },
    {"format": ">8.8s", "attr": "detector_unique_name", "label": "Detector"},
    {
        "format": None,
        "attr": "exposure_time",
        "label": "ExpTime",
        "bad": math.nan * u.s,
    },
    {"attr": "datetime_begin", "label": "Observation Date", "bad": astropy.time.Time(0.0, format="jd")},
)


def read_file(
    file: str,
    hdrnum: int,
    print_trace: bool,
    output_columns: defaultdict[str, list],
    outstream: IO | None = None,
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
    output_columns : `collections.defaultdict` [ `str`, `list` ]
        For table mode, a place to store the columns.
    outstream : `io.StringIO`, optional
        Output stream to use for standard messages. Defaults to `None` which
        uses the default output stream.
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
    write_heading : `bool`, optional
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
        log.info("Analyzing %s...", file)

    try:
        md = read_basic_metadata_from_file(file, hdrnum, can_raise=True)
        if md is None:
            raise RuntimeError(f"Failed to read file {file} HDU={hdrnum}")

        if output_mode.endswith("native"):
            # Strip native and don't change type of md
            output_mode = output_mode[: -len("native")]
        else:
            # Rewrite md as simple dict for output
            md = dict(md.items())

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

        if output_mode == "table" and output_columns is None:
            raise ValueError("Table mode requires output columns")

        for md in headers:
            obs_info = ObservationInfo(md, pedantic=False, filename=file)
            if output_mode == "table":
                for c in TABLE_COLUMNS:
                    output_columns[c["label"]].append(getattr(obs_info, c["attr"]))
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


def _fill_bad_values(value: Any, fillvalue: Any) -> Any:
    """Convert None values to the fill value.

    Parameters
    ----------
    value : `Any`
        Value to check.
    fillvalue : `Any`
        Value to use if `None`.

    Returns
    -------
    filled : `Any`
        Original value or the fill value.
    """
    return value if value is not None else fillvalue


def _dump_columns(output_columns: dict[str, list], outstream: IO | None = None) -> None:
    """Write columns to output stream as a table.

    Parameters
    ----------
    output_columns : `dict` [`str`, `list`]
        The columns to be written, indexed by column name.
    outstream : `io.StringIO` or `None`, optional
        Output stream to use for standard messages. Defaults to `None` which
        uses the default output stream.
    """
    if not output_columns:
        return

    qt = QTable()
    for c in TABLE_COLUMNS:
        data = output_columns[c["label"]]
        # If the column has no content, assume no columns have content and
        # abandon table output.
        if not len(data):
            return

        need_mask = False
        mask = []
        for v in data:
            is_none = v is None
            if is_none:
                need_mask = True
            mask.append(is_none)
        col_format = c.get("format")

        if need_mask:
            data = [_fill_bad_values(v, c.get("bad", "-")) for v in output_columns[c["label"]]]

        # Quantity have to be handled in special way since they need
        # to be merged into a single entity before they can be stored
        # in a column.
        if issubclass(PROPERTIES[c["attr"]].py_type, u.Quantity):
            data = u.Quantity(data)
        elif issubclass(PROPERTIES[c["attr"]].py_type, astropy.time.Time):
            # Force to ISO string.
            data = astropy.time.Time(data).isot

        if need_mask:
            data = MaskedColumn(name=c["label"], data=data, mask=mask, format=col_format)
        else:
            data = Column(data=data, name=c["label"], format=col_format)
        qt[c["label"]] = data

    print("\n".join(qt.pformat(max_lines=-1, max_width=-1)), file=outstream)


def translate_or_dump_headers(
    files: Sequence[str],
    regex: str,
    hdrnum: int,
    print_trace: bool,
    outstream: IO | None = None,
    output_mode: str = "auto",
) -> tuple[list[str], list[str]]:
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
    outstream : `io.StringIO` or `None`, optional
        Output stream to use for standard messages. Defaults to `None` which
        uses the default output stream.
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
    output_columns: defaultdict[str, list] = defaultdict(list)
    for path in sorted(found_files):
        isok = read_file(path, hdrnum, print_trace, output_columns, outstream, output_mode, heading)
        heading = False
        if isok:
            okay.append(path)
        else:
            failed.append(path)

        # Check if we have exceeded the page size and should dump the table.
        if output_mode == "table" and len(output_columns[TABLE_COLUMNS[0]["label"]]) >= MAX_TABLE_PAGE_SIZE:
            _dump_columns(output_columns, outstream)
            output_columns = defaultdict(list)

    if output_columns:
        _dump_columns(output_columns, outstream)

    return okay, failed
