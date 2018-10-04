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

import unittest

from helper import readTestFile, UsefulAsserts
from lsst.obs.metadata import ObservationInfo


class HscTestCase(unittest.TestCase, UsefulAsserts):

    def testHscTranslator(self):
        for file, filter in (("fitsheader-hsc.yaml", "i"),
                             ("fitsheader-hsc-HSCA04090107.yaml", "r")):
            header = readTestFile(file)
            v1 = ObservationInfo(header)
            self.assertEqual(v1.instrument, "HSC")
            self.assertEqual(v1.telescope, "Subaru")
            self.assertEqual(v1.abstract_filter, filter)

            # Sanity check WCS
            self.assertCoordinatesConsistent(v1)


if __name__ == "__main__":
    unittest.main()
