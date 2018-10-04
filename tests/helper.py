# This file is part of obs_metadata.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import yaml
from collections import OrderedDict

# PropertyList is optional
try:
    import lsst.daf.base as dafBase
except ImportError:
    dafBase = None


TESTDIR = os.path.abspath(os.path.dirname(__file__))


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


if dafBase is None:
    yaml.add_constructor("lsst.daf.base.PropertyList", pl_constructor)


def readTestFile(filename):
    """Read the named test file relative to the location of this helper

    Parameters
    ----------
    filename : `str`
        Name of file in the data directory.

    Returns
    -------
    header : `dict`-like
        Header read from file.
    """
    with open(os.path.join(TESTDIR, "data", filename)) as fd:
        header = yaml.load(fd)
    return header


class UsefulAsserts:
    """Class with helpful asserts"""

    def assertCoordinatesConsistent(self, obsinfo, max_sep=1.0, amdelta=0.01):
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
        self.assertLess(sep.to_value(unit="arcmin"), max_sep)
