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

from __future__ import annotations

__all__ = ("read_test_file", "MetadataAssertHelper")

import json
import os
import pickle
import warnings
from collections.abc import MutableMapping
from typing import TYPE_CHECKING, Any, NoReturn

import astropy.units as u
import astropy.utils.exceptions
import yaml
from astropy.io.fits import Header
from astropy.time import Time

from astro_metadata_translator import ObservationInfo

# PropertyList is optional
try:
    import lsst.daf.base as daf_base
except ImportError:
    daf_base = None


# For YAML >= 5.1 need a different Loader for the constructor
Loader: type[yaml.Loader] | type[yaml.FullLoader]
try:
    Loader = yaml.FullLoader
except AttributeError:
    Loader = yaml.Loader


# Define a YAML loader for lsst.daf.base.PropertySet serializations that
# we can use if daf_base is not available.
def pl_constructor(loader: yaml.Loader, node: yaml.SequenceNode) -> Any:
    """Construct an OrderedDict from a YAML file containing a PropertyList."""
    pl: dict[str, Any] = {}
    yield pl
    state = loader.construct_sequence(node, deep=True)
    for key, dtype, value, comment in state:
        if dtype == "Double":
            pl[key] = float(value)
        elif dtype == "Int":
            pl[key] = int(value)
        elif dtype == "Bool":
            pl[key] = value
        else:
            pl[key] = value


if daf_base is None:
    yaml.add_constructor("lsst.daf.base.PropertyList", pl_constructor, Loader=Loader)  # type: ignore


def read_test_file(filename: str, dir: str | None = None) -> MutableMapping[str, Any]:
    """Read the named test file relative to the location of this helper.

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
        if filename.endswith(".yaml"):
            header = yaml.load(fd, Loader=Loader)
        elif filename.endswith(".json"):
            header = json.load(fd)
        else:
            raise RuntimeError(f"Unrecognized file format: {filename}")

    # Cannot directly check for Mapping because PropertyList is not one
    if not hasattr(header, "items"):
        raise ValueError(f"Contents of YAML file {filename} are not a mapping, they are {type(header)}")
    return header


class MetadataAssertHelper:
    """Class with helpful asserts that can be used for testing metadata
    translations.
    """

    # This class is assumed to be combined with unittest.TestCase but mypy
    # does not know this. We need to teach mypy about the APIs we are using.
    if TYPE_CHECKING:

        def assertAlmostEqual(  # noqa: N802
            self,
            a: float,
            b: float,
            places: int | None = None,
            msg: str | None = None,
            delta: float | None = None,
        ) -> None:
            pass

        def assertIsNotNone(self, a: Any, msg: str | None = None) -> None:  # noqa: N802
            pass

        def assertEqual(self, a: Any, b: Any, msg: str | None = None) -> None:  # noqa: N802
            pass

        def assertLess(self, a: Any, b: Any, msg: str | None = None) -> None:  # noqa: N802
            pass

        def assertLessEqual(self, a: Any, b: Any, msg: str | None = None) -> None:  # noqa: N802
            pass

        def fail(self, msg: str) -> NoReturn:
            pass

    def assertCoordinatesConsistent(  # noqa: N802
        self, obsinfo: ObservationInfo, max_sep: float = 1.0, amdelta: float = 0.01
    ) -> None:
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
        # Can trigger warnings from Astropy if date is in future
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=astropy.utils.exceptions.AstropyWarning)
            sep = obsinfo.altaz_begin.separation(obsinfo.tracking_radec.altaz)
        self.assertLess(sep.to_value(unit="arcmin"), max_sep, msg="AltAz inconsistent with RA/Dec")

    def assertObservationInfoFromYaml(  # noqa: N802
        self,
        file: str,
        dir: str | None = None,
        check_wcs: bool = True,
        wcs_params: dict[str, Any] | None = None,
        check_altaz: bool = False,
        **kwargs: Any,
    ) -> None:
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
        check_altaz : `bool`, optional
            Check that an alt/az value has been calculated.
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

        # DM-30093: Astropy Header does not always behave dict-like
        astropy_header = Header()
        for key, val in header.items():
            values = val if isinstance(val, list) else [val]
            for v in values:
                # appending ensures *all* duplicated keys are also preserved
                astropy_header.append((key, v))

        for hdr in (header, astropy_header):
            try:
                self.assertObservationInfo(
                    header,
                    filename=file,
                    check_wcs=check_wcs,
                    wcs_params=wcs_params,
                    check_altaz=check_altaz,
                    **kwargs,
                )
            except AssertionError as e:
                raise AssertionError(
                    f"ObservationInfo derived from {type(hdr)} type is inconsistent: {e}"
                ) from e

    def assertObservationInfo(  # noqa: N802
        self,
        header: MutableMapping[str, Any],
        filename: str | None = None,
        check_wcs: bool = True,
        wcs_params: dict[str, Any] | None = None,
        check_altaz: bool = False,
        **kwargs: Any,
    ) -> None:
        """Check contents of an ObservationInfo.

        Parameters
        ----------
        header : `dict`-like
            Header to be checked.
        filename : `str`, optional
            Name of the filename associated with this header if known.
        check_wcs : `bool`, optional
            Check the consistency of the RA/Dec and AltAz values.  Checks
            are automatically disabled if the translated header does
            not appear to be "science".
        wcs_params : `dict`, optional
            Parameters to pass to `assertCoordinatesConsistent`.
        check_altaz : `bool`, optional
            Check that an alt/az value has been calculated.
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
        obsinfo = ObservationInfo(header, pedantic=True, filename=filename)
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

        # Date comparison error reports will benefit by specifying ISO
        # format.  Generate a new Time object at fixed precision
        # to work around the fact that (as of astropy 3.1) adding 0.0 seconds
        # to a Time results in a new Time object that is a few picoseconds in
        # the past.
        def _format_date_for_testing(date: Time | None) -> Time | None:
            if date is not None:
                date.format = "isot"
                date.precision = 9
                date = Time(str(date), scale=date.scale, format="isot")
            return date

        datetime_begin = _format_date_for_testing(obsinfo.datetime_begin)
        datetime_end = _format_date_for_testing(obsinfo.datetime_end)

        # Check that dates are defined and end is the same or after beginning
        self.assertLessEqual(datetime_begin, datetime_end)

        # Check that exposure time is not outside datetime_end
        self.assertLessEqual(obsinfo.datetime_begin + obsinfo.exposure_time, obsinfo.datetime_end)

        # Do we expect an AltAz value or not.
        if check_altaz:
            self.assertIsNotNone(obsinfo.altaz_begin, "altaz_begin is None but should have a value")

        # Check the WCS consistency
        if check_wcs and obsinfo.observation_type == "science":
            if wcs_params is None:
                wcs_params = {}
            self.assertCoordinatesConsistent(obsinfo, **wcs_params)
