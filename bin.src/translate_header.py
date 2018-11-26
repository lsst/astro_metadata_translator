#!/usr/bin/env python3

import argparse
import os
import re
import sys
import traceback
import importlib
import yaml
from astro_metadata_translator import ObservationInfo

# Prefer afw over Astropy
try:
    from lsst.afw.fits import readMetadata as read_metadata  # noqa: N813
    import lsst.daf.base  # noqa: F401 need PropertyBase for readMetadata
except ImportError:
    from astropy.io import fits

    def read_metadata(file, hdu=1):
        fits_file = fits.open(file)
        return fits_file[hdu].header

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
        md = read_metadata(file)
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
        for name in os.listdir(file):
            path = os.path.join(file, name)
            if os.path.isfile(path) and file_regex.search(name):
                read_file(path, failed)
    else:
        read_file(file, failed)

if failed:
    print("Files with failed translations:", file=sys.stderr)
    for f in failed:
        print(f"\t{f}", file=sys.stderr)
