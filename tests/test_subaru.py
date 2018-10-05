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

from helper import UsefulAsserts


class HscTestCase(unittest.TestCase, UsefulAsserts):

    def testHscTranslator(self):
        test_data = (("fitsheader-hsc.yaml", dict(telescope="Subaru", instrument="HSC",
                                                  abstract_filter="i")),
                     ("fitsheader-hsc-HSCA04090107.yaml", dict(telescope="Subaru", instrument="HSC",
                                                               abstract_filter="r")),
                     )
        for file, expected in test_data:
            self.assertObservationInfo(file, **expected)

    def testSuprimeCamTranslator(self):
        # In this case the airmass is average during observation
        # but it looks like ALTITUDE is from a different time so loosen amdelta
        test_data = (("fitsheader-suprimecam-CORR40535770.yaml",
                      dict(telescope="Subaru", instrument="SuprimeCam",
                           abstract_filter="r", wcsParams=dict(amdelta=0.015))),
                     )
        for file, expected in test_data:
            self.assertObservationInfo(file, **expected)


if __name__ == "__main__":
    unittest.main()
