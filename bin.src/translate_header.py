#!/usr/bin/env python3

import argparse
import os
import re
import sys
import traceback
import importlib
import yaml
from astro_metadata_translator import ObservationInfo, merge_headers

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

args = parser.parse_args()

# Process import requests
if args.packages:
    for m in args.packages:
        importlib.import_module(m)


def read_file(file, failed):
    print(f"Analyzing {file}...", file=sys.stderr)
    try:
        md = read_metadata(file, 0)
        if args.hdrnum != 0:
            mdn = read_metadata(file, int(args.hdrnum))
            # Astropy does not allow append mode since it does not
            # convert lists to multiple cards. Overwrite for now
            if mdn is not None:
                md = merge_headers([md, mdn], mode="overwrite")
            else:
                print(f"HDU {args.hdrnum} was not found. Ignoring request.", file=sys.stderr)

        if args.dumphdr:
            print(yaml.dump(md))
            return
        obs_info = ObservationInfo(md, pedantic=True, filename=file)
        if not args.quiet:
            print(f"{obs_info}")
    except Exception as e:
        if args.traceback:
            traceback.print_exc(file=sys.stdout)
        else:
            print(repr(e))
        failed.append(file)


file_regex = re.compile(args.regex)
failed = []
for file in args.files:
    if os.path.isdir(file):
        for root, dirs, files in os.walk(file):
            for name in files:
                path = os.path.join(root, name)
                if os.path.isfile(path) and file_regex.search(name):
                    read_file(path, failed)
    else:
        read_file(file, failed)

if failed:
    print("Files with failed translations:", file=sys.stderr)
    for f in failed:
        print(f"\t{f}", file=sys.stderr)
