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

__all__ = ("read_test_file", "MetadataAssertHelper")

import astropy.units as u
import os
import pickle
import yaml
from collections import OrderedDict

from astro_metadata_translator import ObservationInfo

# PropertyList is optional
try:
    import lsst.daf.base as daf_base
except ImportError:
    daf_base = None


# For YAML >= 5.1 need a different Loader for the constructor
try:
    Loader = yaml.FullLoader
except AttributeError:
    Loader = yaml.Loader


# Define a YAML loader for lsst.daf.base.PropertySet serializations that
# we can use if daf_base is not available.
def pl_constructor(loader, node):
    """Construct an OrderedDict from a YAML file containing a PropertyList."""
    pl = OrderedDict()
    yield pl
    state = loader.construct_sequence(node, deep=True)
    for key, dtype, value, comment in state:
        if dtype == "Double":
            pl[key] = float(value)
        elif dtype == "Int":
            pl[key] = int(value)
        elif dtype == "Bool":
            pl[key] = True if value == "true" else False
        else:
            pl[key] = value


if daf_base is None:
    yaml.add_constructor("lsst.daf.base.PropertyList", pl_constructor, Loader=Loader)


def read_test_file(filename, dir=None):
    """Read the named test file relative to the location of this helper

    Parameters
    ----------
    filename : `str`
        Name of file in the data directory.
    dir : `str`, optional.
        Directory from which to read file. Current directory used if none
        specified.

    Returns
    -------
    header : `dict`-like
        Header read from file.
    """
    if dir is not None and not os.path.isabs(filename):
        filename = os.path.join(dir, filename)
    with open(filename) as fd:
        header = yaml.load(fd, Loader=Loader)
    return header


class MetadataAssertHelper:
    """Class with helpful asserts that can be used for testing metadata
    translations.
    """

    def assertCoordinatesConsistent(self, obsinfo, max_sep=1.0, amdelta=0.01):  # noqa: N802
        """Check that SkyCoord, AltAz, and airmass are self consistent.

        Parameters
        ----------
        obsinfo : `ObservationInfo`
            Object to check.
        max_sep : `float`, optional
            Maximum separation between AltAz derived from RA/Dec headers
            and that found in the AltAz headers.
        amdelta : `float`, optional
            Max difference between header airmass and derived airmass.

        Raises
        ------
        AssertionError
            Inconsistencies found.
        """
        self.assertIsNotNone(obsinfo.tracking_radec)
        self.assertIsNotNone(obsinfo.altaz_begin)

        # Is airmass from header close to airmass from AltAz headers?
        # In most cases there is uncertainty over whether the elevation
        # and airmass in the header are for the start, end, or middle
        # of the observation.  Sometimes the AltAz is from the start
        # but the airmass is from the middle so accuracy is not certain.
        self.assertAlmostEqual(obsinfo.altaz_begin.secz.to_value(), obsinfo.boresight_airmass, delta=amdelta)

        # Is AltAz from headers close to AltAz from RA/Dec headers?
        sep = obsinfo.altaz_begin.separation(obsinfo.tracking_radec.altaz)
        self.assertLess(sep.to_value(unit="arcmin"), max_sep, msg="AltAz inconsistent with RA/Dec")

    def assertObservationInfoFromYaml(self, file, dir=None, check_wcs=True,  # noqa: N802
                                      wcs_params=None, **kwargs):
        """Check contents of an ObservationInfo.

        Parameters
        ----------
        file : `str`
            Path to YAML file representing the header.
        dir : `str`, optional
            Optional directory from which to read ``file``.
        check_wcs : `bool`, optional
            Check the consistency of the RA/Dec and AltAz values.
        wcs_params : `dict`, optional
            Parameters to pass to `assertCoordinatesConsistent`.
        kwargs : `dict`
            Keys matching `ObservationInfo` properties with values
            to be tested.

        Raises
        ------
        AssertionError
            A value in the ObservationInfo derived from the file is
            inconsistent.
        """
        header = read_test_file(file, dir=dir)
        self.assertObservationInfo(header, check_wcs=check_wcs, wcs_params=wcs_params, **kwargs)

    def assertObservationInfo(self, header, check_wcs=True, wcs_params=None, **kwargs):  # noqa: N802
        """Check contents of an ObservationInfo.

        Parameters
        ----------
        header : `dict`-like
            Header to be checked.
        check_wcs : `bool`, optional
            Check the consistency of the RA/Dec and AltAz values.  Checks
            are automatically disabled if the translated header does
            not appear to be "science".
        wcs_params : `dict`, optional
            Parameters to pass to `assertCoordinatesConsistent`.
        kwargs : `dict`
            Keys matching `ObservationInfo` properties with values
            to be tested.

        Raises
        ------
        AssertionError
            A value in the ObservationInfo derived from the file is
            inconsistent.
        """
        # For testing we force pedantic mode since we are in charge
        # of all the translations
        obsinfo = ObservationInfo(header, pedantic=True)
        translator = obsinfo.translator_class_name

        # Check that we can pickle and get back the same properties
        newinfo = pickle.loads(pickle.dumps(obsinfo))
        self.assertEqual(obsinfo, newinfo)

        # Check the properties
        for property, expected in kwargs.items():
            calculated = getattr(obsinfo, property)
            msg = f"Comparing property {property} using translator {translator}"
            if isinstance(expected, u.Quantity) and calculated is not None:
                calculated = calculated.to_value(unit=expected.unit)
                expected = expected.to_value()
                self.assertAlmostEqual(calculated, expected, msg=msg)
            elif isinstance(calculated, u.Quantity):
                # Only happens if the test is not a quantity when it should be
                self.fail(f"Expected {expected!r} but got Quantity '{calculated}': {msg}")
            elif isinstance(expected, float) and calculated is not None:
                self.assertAlmostEqual(calculated, expected, msg=msg)
            else:
                self.assertEqual(calculated, expected, msg=msg)

        # Check the WCS consistency
        if check_wcs and obsinfo.observation_type == "science":
            if wcs_params is None:
                wcs_params = {}
            self.assertCoordinatesConsistent(obsinfo, **wcs_params)
