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

"""Support functions for script implementations."""

__all__ = ("read_metadata", "find_files", "read_metadata_from_file")

import re
import os
import sys

from astro_metadata_translator import merge_headers
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


def find_files(files, regex):
    """Find files for processing.

    Parameters
    ----------
    files : iterable of `str`
        The files or directories from which the headers are to be read.
    regex : `str`
        Regular expression string used to filter files when a directory is
        scanned.
    """
    file_regex = re.compile(regex)
    found_files = []

    # Find all the files of interest
    for file in files:
        if os.path.isdir(file):
            for root, dirs, files in os.walk(file):
                for name in files:
                    path = os.path.join(root, name)
                    if os.path.isfile(path) and file_regex.search(name):
                        found_files.append(path)
        else:
            found_files.append(file)

    return found_files


def read_metadata_from_file(file, hdrnum, errstream=sys.stderr):
    if file.endswith(".yaml"):
        md = read_test_file(file,)
        if hdrnum != 0:
            # YAML can't have HDUs
            hdrnum = 0
    else:
        md = read_metadata(file, 0)
    if md is None:
        print(f"Unable to open file {file}", file=errstream)
        return None
    if hdrnum != 0:
        mdn = read_metadata(file, int(hdrnum))
        # Astropy does not allow append mode since it does not
        # convert lists to multiple cards. Overwrite for now
        if mdn is not None:
            md = merge_headers([md, mdn], mode="overwrite")
        else:
            print(f"HDU {hdrnum} was not found in file {file}. Ignoring request.", file=errstream)

    return md
