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

__all__ = ("read_index", "calculate_index", "index_files", "process_index_data")

"""Functions to support file indexing."""

import json
import logging
import os
import sys
from copy import deepcopy

from .observationInfo import ObservationInfo
from .observationGroup import ObservationGroup
from .headers import merge_headers
from .file_helpers import read_file_info

log = logging.getLogger(__name__)

COMMON_KEY = "__COMMON__"
CONTENT_KEY = "__CONTENT__"


def index_files(files, root, hdrnum, print_trace, content, outstream=sys.stdout, errstream=sys.stderr):
    """Create an index from the supplied files.

    No file is written. The Python structure returned is suitable
    for writing.

    Parameters
    ----------
    files : iterable of `str`
        Paths to the files to be indexed. They do not have to all be
        in a single directory but all content will be indexed into a single
        index.
    root : `str`
        Directory root that can be combined with each file (if the supplied)
        file is relative. Will be ignored if `None`.
    hdrnum : `int`
        The HDU number to read. The primary header is always read and
    print_trace : `bool`
        If there is an error reading the file and this parameter is `True`,
        a full traceback of the exception will be reported. If `False` prints
        a one line summary of the error condition. If `None` the exception
        will be allowed.
    content : `str`
        Form of data to write in index file. Options are:
        ``translated`` (default) to write ObservationInfo to the index;
        ``metadata`` to write native metadata headers to the index.
        The index file is called ``{mode}_index.json``
    outstream : `io.StringIO`, optional
        Output stream to use for standard messages. Defaults to `sys.stdout`.
    errstream : `io.StringIO`, optional
        Stream to send messages that would normally be sent to standard
        error. Defaults to `sys.stderr`.

    Returns
    -------
    file_index : `dict` of [`str`, `dict`]
        The headers in form suitable for writing to an index. The keys will
        be ``__COMMON__`` for shared content, ``__CONTENT_`` to record the
        content mode used to construct the index, and paths to the files. The
        paths will be the supplied paths and will not include any supplied
        ``root``.
    okay : `list` of `str`
        All the files that were processed successfully.
    failed : `list` of `str`
        All the files that could not be processed. Will be empty if
        ``print_trace`` is not `None`.
    """
    if content not in ("translated", "metadata"):
        raise ValueError("Unrecognized mode {mode}")

    failed = []
    okay = []

    content_by_file = {}  # Mapping of path to file content
    for file in sorted(files):
        if root is not None:
            path = os.path.join(root, file)
        else:
            path = file
        simple = read_file_info(path, hdrnum, print_trace, content, "simple", outstream, errstream)
        if simple is None:
            failed.append(path)
            continue
        else:
            okay.append(path)

        # Store the information indexed by the filename within dir
        content_by_file[file] = simple

    output = calculate_index(content_by_file, content)

    return output, okay, failed


def calculate_index(headers, content):
    """Calculate an index data structure from the supplied headers.

    Parameters
    ----------
    headers : `dict` of [`str`, `dict`]
        The headers indexed by filename.
    content : `str`
        The mode associated with these headers. Not used other than to
        store the information in the data structure for later use on
        deserialization.

    Returns
    -------
    index_ : `dict` of [`str`, `dict`]
        The headers in form suitable for writing to an index.
    """
    if content not in ("metadata", "translated"):
        raise ValueError(f"Unrecognized mode for index creation: {content}")

    # Merge all the information into a primary plus diff
    merged = merge_headers(headers.values(), mode="diff")

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
    # Store the mode so that we can work out how to read the file in
    output = {CONTENT_KEY: content, COMMON_KEY: merged}
    for file, diff in zip(headers, diff_dict):
        output[file] = diff

    return output


def read_index(path, force_dict=False):
    """Read an index file.

    Parameters
    ----------
    path : `str`
        Path to the index file.
    force_dict : `bool`, optional
        If `True` the structure returned will always be a dict keyed
        by filename.

    Returns
    -------
    index_ : `ObservationGroup` or `dict[str, Union[dict, ObservaitonInfo]]`
       If the index file referred to `ObservationInfo` this will return
       an `ObservationGroup`, otherwise a `dict` will be returned with the
       keys being paths to files and the values being the keys and values
       stored in the index (with common information merged in). If
       ``force_dict`` is `True` a `dict` will be returned with filename
       keys even if the index file refers to `ObservationInfo` (the values
       will be `ObservationInfo`).
    """
    if not path.endswith(".json"):
        raise ValueError(f"Index files must be in .json format; got {path}")

    with open(path, "r") as fd:
        content = json.loads(fd.read())

    return process_index_data(content, force_dict=force_dict)


def process_index_data(content, force_metadata=False, force_dict=False):
    """Process the content read from a JSON index file.

    Parameters
    ----------
    content : `dict`
        Data structure stored in JSON index file converted to simple python
        form.
    force_metadata : `bool`, optional
        By default the content returned will match the original form that
        was used for the index. If this parameter is `True` an index of
        `ObservationInfo` will be returned as if it was simple dict content.
    force_dict : `bool`, optional
        If `True` the structure returned will always be a dict keyed
        by filename.

    Returns
    -------
    index : `ObservationGroup` or `dict` of [`str`, `dict`]
       If the index file referred to `ObservationInfo` this will return
       an `ObservationGroup`, otherwise a `dict` will be returned with the
       keys being paths to files and the values being the keys and values
       stored in the index (with common information merged in). This
       can be overridden using the ``force_metadata`` parameter. If
       ``force_dict`` is `True` a `dict` will be returned with filename
       keys even if the index file refers to `ObservationInfo` (the values
       will be `ObservationInfo` unless ``force_metadata`` is `True`).

    Notes
    -----
    File keys will be relative to the location of the index file.
    """

    if COMMON_KEY not in content:
        raise ValueError(f"No '{COMMON_KEY}' key found in dict. Does not look like an index data structure.")

    # Copy the input structure so we can update in place
    unpacked = deepcopy(content)

    content = unpacked.pop(CONTENT_KEY, None)
    if force_metadata:
        content = "metadata"
    elif content is None:
        log.warning("No '%s' key in data structure, assuming 'metadata'", CONTENT_KEY)
        content = "metadata"

    # The common headers will be copied into each header
    common = unpacked.pop(COMMON_KEY)

    for file in unpacked:
        unpacked[file].update(common)

    if content == "metadata":
        # nothing more to be done
        return unpacked

    obs_infos = []
    by_file = {}
    for file, hdr in unpacked.items():
        info = ObservationInfo.from_simple(hdr)
        info.filename = file
        obs_infos.append(info)
        by_file[file] = info

    if force_dict:
        return by_file
    return ObservationGroup(obs_infos)


def read_sidecar(path):
    """Read a metadata sidecar file.

    Parameters
    ----------
    path : `str`
        Path to the sidecar file.

    Returns
    -------
    info : `ObservationInfo` or `dict` of [`str`, `dict`]
        If the sidecar file referred to `ObservationInfo` this will return
        an `ObservationInfo`, otherwise a `dict` will be returned.
    """
    if not path.endswith(".json"):
        raise ValueError(f"Sidecar files must be in .json format; got {path}")

    with open(path, "r") as fd:
        content = json.loads(fd.read())

    return process_sidecar_data(content)


def process_sidecar_data(content, force_metadata=False):
    """Process the content read from a JSON sidecar file.

    Parameters
    ----------
    content : `dict`
        Data structure stored in JSON sidecar file converted to simple python
        form.
    force_metadata : `bool`, optional
        By default the content returned will match the original form that
        was used for the sidecar. If this parameter is `True` a sidecar of
        `ObservationInfo` will be returned as if it was simple dict content.

    Returns
    -------
    info : `ObservationInfo` or `dict` of [`str`, `dict`]
        If the sidecar file referred to `ObservationInfo` this will return
        an `ObservationGroup`, otherwise a `dict` will be returned. This
        can be overridden using the ``force_metadata`` parameter.
    """

    if not isinstance(content, dict):
        raise TypeError(f"Content of sidecar must be a dict, not {type(content)}")

    # Copy the input structure so we can update in place
    content = deepcopy(content)

    guessing = False
    content_mode = content.pop(CONTENT_KEY, None)
    if force_metadata:
        content_mode = "metadata"
    elif content is None:
        # All ObservationInfo objects will have observation_id and instrument
        # so if they are there we can guess
        guessing = True
        if "observation_id" in content and "instrument" in content_mode:
            content_mode = "translated"
        else:
            content_mode = "metadata"
        log.warning("No '%s' key in data structure, assuming '%s'", CONTENT_KEY, content_mode)

    if content_mode == "metadata":
        # nothing more to be done
        return content

    try:
        info = ObservationInfo.from_simple(content)
    except Exception as e:
        if guessing:
            # We were guessing so seems like this is not ObservationInfo
            return content
        raise e

    return info
