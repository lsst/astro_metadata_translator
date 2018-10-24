#!/usr/bin/env python3

import argparse
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
                    help="Do not report the summary for each file")

args = parser.parse_args()
print(repr(args))

failed = []
for file in args.files:
    print(f"Analyzing {file}...")
    try:
        md = read_metadata(file)
        obs_info = ObservationInfo(md, pedantic=True)
        if not args.quiet:
            print(f"{obs_info}")
    except Exception as e:
        print(e)
        failed.append(file)

if failed:
    print("Files with failed translations:")
    for f in failed:
        print(f"\t{f}")
