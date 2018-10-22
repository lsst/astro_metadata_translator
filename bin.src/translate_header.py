#!/usr/bin/env python3

import sys
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


failed = []
for file in sys.argv[1:]:
    print(f"{file}:")
    md = read_metadata(file)
    try:
        obs_info = ObservationInfo(md, pedantic=True)
        print(f"{obs_info}")
    except Exception:
        failed.append(file)

if failed:
    print("Files with failed translations:")
    for f in failed:
        print(f"\t{f}")
