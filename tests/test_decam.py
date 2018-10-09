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

from helper import MetadataAssertHelper


class DecamTestCase(unittest.TestCase, MetadataAssertHelper):

    def testDecamTranslator(self):
        test_data = (("fitsheader-decam.yaml",
                      dict(telescope="CTIO 4.0-m telescope", instrument="DECam",
                           abstract_filter="z", wcsParams=dict(max_sep=1.5))),
                     )
        for file, expected in test_data:
            self.assertObservationInfoFromYaml(file, **expected)


if __name__ == "__main__":
    unittest.main()
