#!/usr/bin/env python3

import sys
import traceback
import argparse
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
                    help="File(s) from which headers will be parsed")
parser.add_argument("-q", "--quiet", const=True, default=False, action="store_const",
                    help="Do not report the translation content from each header")
parser.add_argument("-d", "--dumphdr", const=True, default=False, action="store_const",
                    help="Dump the header in YAML format to standard output rather than translating it")
parser.add_argument("--traceback", const=True, default=False, action="store_const",
                    help="Give detailed trace back when any errors encountered")

args = parser.parse_args()

failed = []
for file in args.files:
    print(f"Analyzing {file}...", file=sys.stderr)
    try:
        md = read_metadata(file)
        if args.dumphdr:
            print(yaml.dump(md))
            continue
        obs_info = ObservationInfo(md, pedantic=True)
        if not args.quiet:
            print(f"{obs_info}")
    except Exception as e:
        if args.traceback:
            traceback.print_exc(file=sys.stdout)
        else:
            print(repr(e))
        failed.append(file)

if failed:
    print("Files with failed translations:", file=sys.stderr)
    for f in failed:
        print(f"\t{f}", file=sys.stderr)
