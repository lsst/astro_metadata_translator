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

from helper import readTestFile
from lsst.obs.metadata import ObservationInfo


class DecamTestCase(unittest.TestCase):

    def testDecamTranslator(self):
        header = readTestFile("fitsheader-decam.yaml")
        v1 = ObservationInfo(header)
        self.assertEqual(v1.instrument, "DECam")
        self.assertEqual(v1.telescope, "CTIO 4.0-m telescope")
        self.assertEqual(v1.abstract_filter, "z")


if __name__ == "__main__":
    unittest.main()
